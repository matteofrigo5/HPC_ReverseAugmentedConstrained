#include "RevAugm_params.h"
#include "RevAugm_default_params.h"
#include "RevAugm.h"

#include "Chronos.h"
#include "SadPointMat.h"

#include "cpt_D_coef_loc.h"

#define PRINT_D false

#if PRINT_D
namespace
{
std::string toString( int const i, int const w = 2 )
{
   char fmt[200];
   std::snprintf( fmt, 200, "%%%i.%ii", w, w );
   char buf[200];
   std::snprintf( buf, 200, fmt, i );
   std::string str( buf );
   return str;
}
}
#endif

//----------------------------------------------------------------------------------------

// Creates an empty object.
RevAugm::RevAugm(){
   setDefaultParameters( params );
}

//----------------------------------------------------------------------------------------

// Deletes the object.
RevAugm::~RevAugm(){

   if ( prec_K != nullptr ) {
      delete prec_K;
      prec_K = nullptr;
   }

   D_coef.~VEC();

   MatC=nullptr;
   MatCt=nullptr;

}

//----------------------------------------------------------------------------------------

// Prepares the data structure necessary for MxV.
void RevAugm::Prepare_MxV(const iReg nRHS){

   // Retrieve current communicator
   MPI_Comm const currComm = Chronos.Get_currComm();

   // Check if aAMG is prepared
   if ( Prepared_flag == true && nRHS == get_nRHS_prep() ) return;

   prec_K->Prepare_MxV(nRHS);

   // set flag
   Prepared_flag = true;
   nRHS_prep = nRHS;

   // Syncronize MPI processes before exit
   MPI_Barrier(currComm);

}

//----------------------------------------------------------------------------------------

// Compute MxV
void RevAugm::MxV(DDMat& __restrict__ x, DDMat& __restrict__ b,
                           bool Barrier_flag){

   iReg n1,n2;
   iReg nRHS;
   DDMat b1,b2,x1,x2;

   // Retrieve number of threads
   type_OMP_iReg nthreads = Chronos.Get_inpNthreads();

   // Retrive dimensions
   n1   = MatC->get_nrows();
   n2   = MatCt->get_nrows();
   nRHS = b.get_ncols();

   rExt* ptr_b = b.get_ptr_coef_data();
   rExt* ptr_x = x.get_ptr_coef_data();
   b1.overlap(n1,nRHS,ptr_b);
   b2.overlap(n2,nRHS,ptr_b+n1);
   x1.overlap(n1,nRHS,ptr_x);
   x2.overlap(n2,nRHS,ptr_x+n1);

   DDMat tmp1=DDMat(n1,nRHS);
   DDMat tmp2=DDMat(n2,nRHS);

   nthreads = static_cast<type_OMP_iReg>(std::min(static_cast<iReg>(nthreads),n1*nRHS));
   rExt* ptr_x1  = x1.get_ptr_coef_data();
   rExt* ptr_x2  = x2.get_ptr_coef_data();
   rExt* ptr_b2  = b2.get_ptr_coef_data();
   rExt* ptr_tmp1 = tmp1.get_ptr_coef_data();
   rExt* ptr_tmp2 = tmp2.get_ptr_coef_data();
   rExt* ptr_D_coef = D_coef.data();

   //  Apply D_coef to x2 to get tmp2
   #pragma omp parallel for num_threads(nthreads)
   for (iReg irow = 0; irow < n2; irow++){
      iReg ind = irow*nRHS;
      for (iReg j = 0; j < nRHS; j++) ptr_tmp2[ind+j] = ptr_D_coef[irow]*ptr_x2[ind+j];
   }

   // tmp1 = C x tmp2
   MatC->MxV(tmp2,tmp1,Barrier_flag);

   // tmp1 <- tmp1 + x1
   #pragma omp parallel for num_threads(nthreads)
   for(iReg k=0; k<n1*nRHS; ++k){
      ptr_tmp1[k]+=ptr_x1[k];
   }

   //  Apply prec_K to tmp1 to get b1
   prec_K->MxV(tmp1,b1,Barrier_flag);

   // tmp2 = Ct x b1
   MatCt->MxV(b1,tmp2,Barrier_flag);

   // tmp2 <- tmp2 - x2
   #pragma omp parallel for num_threads(nthreads)
   for(iReg k=0; k<n2*nRHS; ++k){
      ptr_tmp2[k]-=ptr_x2[k];
   }

   //  Apply D_coef to tmp2 to get b2
   #pragma omp parallel for num_threads(nthreads)
   for (iReg irow = 0; irow < n2; irow++){
      iReg ind = irow*nRHS;
      for (iReg j = 0; j < nRHS; j++) ptr_b2[ind+j] = ptr_D_coef[irow]*ptr_tmp2[ind+j];
   }

   x1.assign_null();
   x2.assign_null();
   b1.assign_null();
   b2.assign_null();

}

//----------------------------------------------------------------------------------------

// Undo the data structure necessary for MxV.
void RevAugm::UndoPrepare_MxV(){

   prec_K->UndoPrepare_MxV();
   Prepared_flag=false;
}

//----------------------------------------------------------------------------------------

void RevAugm::Compute ( DSMat& __restrict__ Mat_in , DDMat *const __restrict__ V0 /*= nullptr*/ ){

   // Retrieve current communicator
   MPI_Comm const currComm = Chronos.Get_currComm();

   type_MPI_iReg rank;
   MPI_Comm_rank(currComm,&rank);

   SadPointMat *Mat;
   try
   {
      Mat = dynamic_cast< SadPointMat* >( &Mat_in );
   }
   catch (bad_cast const &)
   {
      throw linsol_error( "RevAugm::Compute", "dynamic_cast error" );
   }

   // Set pointer to matrices MatC and MatCt
   iReg n2 = Mat->MatCt.get_nrows();
   D_coef.assign(n2,params.gamma);

   DSMat & MatK = Mat->MatK;
   MatC = &(Mat->MatC);
   MatCt = &(Mat->MatCt);

   // Set time to compute the augmentation
   rExt time = MPI_Wtime();

   // Compute the diagonal matrix D
   cptAugMat( *MatCt, MatK, D_coef );

   // Get time to compute the augmentation
   time = MPI_Wtime() - time;

#if PRINT_D
   std::ofstream file( "D_" + toString( rank ) );
   for( int i = 0; i < n2; ++i )
   {
      file << D_coef[i] << std::endl;
   }
   file.close();
#endif

   if( rank == 0 )
   {
     std::cout << "Time to compute augmentation matrix [s] " << time << std::endl;
   }

   DSMat MatDCt;
   MatDCt = *MatCt;

   // Compute the product MatDCt = D*Ct
   // Retrieve number of threads
   type_OMP_iReg nthreads = Chronos.Get_inpNthreads();
   if( n2 == 0 )
   {
      nthreads = 1;
   }
   else
   {
      nthreads = static_cast< type_OMP_iReg >( std::min( static_cast< iReg >( nthreads ), n2 ) );
   }

   #pragma omp parallel num_threads(nthreads)
   {
      // Retrieve processor ID
      iReg const myid = static_cast< iReg >( omp_get_thread_num() );

      // Get the first and last row of the current thread
      iReg const bsize = n2/nthreads;
      iReg const resto = n2%nthreads;
      iReg nrows_th, firstrow;
      if (myid <= resto)
      {
         nrows_th = bsize+1;
         firstrow = myid*nrows_th;
         if (myid == resto)
         {
            nrows_th--;
         }
      }
      else
      {
         nrows_th = bsize;
         firstrow = myid*bsize + resto;
      }
      iReg const lastrow = firstrow + nrows_th;

      iReg const nBlk = MatDCt.MatGraph->get_ncon_L() + 1 + MatDCt.MatGraph->get_ncon_R();
      rExt const * const pt_D = D_coef.data();
      for (iReg ib = 0; ib < nBlk; ib++)
      {
         iExt const *iat = MatDCt.csr[ib].get_ptr_iat_data();
         rExt *coef_out = MatDCt.csr[ib].get_ptr_coef_data();
         for (iReg i = firstrow; i < lastrow; i++)
         {
            for (iReg j = iat[i]; j < iat[i+1]; j++)
            {
               coef_out[j] *= pt_D[i];
            }
         }
      }

   } // end parallel

   // Compute the augmented block
   DSMat CDCt;

   try
   {
      Mat->MatC.MxM( MatDCt, CDCt );
   }
   catch (linsol_error )
   {
      throw linsol_error( "RevAugm::Compute", "computing MxM" );
   }

   CDCt.set_sym_flag( true );
   delete CDCt.OutptDom;
   delete CDCt.OutMatGraph;
   CDCt.mk_OutptDom();
   CDCt.OutMatGraph = CDCt.MatGraph;

   try
   {
      Mat->MatK.sum( CDCt, Schur,params.gamma, 1.0 );
   }
   catch (linsol_error )
   {
      throw linsol_error( "RevAugm::Compute", "computing M+M" );
   }
   Schur.set_sym_flag(true);

   // Compute preconditioner for block [0,0]
   try
   {
      prec_K->Compute(Schur,V0);
   }
   catch (linsol_error)
   {
      throw linsol_error( "RevAugm::Compute", "in smoother computation" );
   }

}

//----------------------------------------------------------------------------------------

void  RevAugm::PrecK_set_type(const iReg type){

   params.type_precK=type;

   if (type==0){
      // Allocate and set default parameters
      prec_K = new(nothrow) aFSAI;
      if ( prec_K == nullptr ) throw linsol_error ("RevAugm::PrecK_set_type","allocating prec_K");

   } else if (type==1) {

      // Allocate and set default parameters
      prec_K = new(nothrow) aAMG;
      if ( prec_K == nullptr ) throw linsol_error ("RevAugm::PrecK_set_type","allocating prec_K");

   } else if (type==2) {

      // Allocate and set default parameters
      prec_K = new(nothrow) aAMG;
      if ( prec_K == nullptr ) throw linsol_error ("RevAugm::PrecK_set_type","allocating prec_K");

   }

   setParametersBlock11( params );
}

//----------------------------------------------------------------------------------------

void RevAugm::set_gamma(const rExt gamma_in){
   params.gamma=gamma_in;
}

//----------------------------------------------------------------------------------------

void RevAugm::cpt_D_coef( DSMat const & MatK,
                                   DSMat const * const MatC,
                                   DSMat const * const MatCt )
{
   // Retrive dimensions
   iReg const n1   = MatC->get_nrows();
   iReg const n2   = MatCt->get_nrows();

   // Retrieve number of threads
   type_OMP_iReg nthreads = Chronos.Get_inpNthreads();
   nthreads = static_cast<type_OMP_iReg>(min(static_cast<iReg>(nthreads),n1));

   // vecstart for MatCt [n2 x n1]
   iReg const bsize = n2 / nthreads;
   iReg const rest = n2 % nthreads;
   VEC<iReg> vec_vecstart;
   try
   {
      vec_vecstart.resize(nthreads+1);
   }
   catch( linsol_error )
   {
      throw linsol_error ("EXTI","allocationg vecstart");
   }
   iReg* vecstart = vec_vecstart.data();

   vecstart[0] = 0;
   vecstart[1] = bsize + rest;
   for (iReg i=2; i<=nthreads; i++)
   {
      vecstart[i] = vecstart[i-1] + bsize;
   }

   iExt* iat_Ct  = MatCt->matC->get_ptr_iat_data();
   iReg* ja_Ct   = MatCt->matC->get_ptr_ja_data();
   rExt* coef_Ct = MatCt->matC->get_ptr_coef_data();

   iExt* iat_K   = MatK.matC->get_ptr_iat_data();
   iReg* ja_K    = MatK.matC->get_ptr_ja_data();
   rExt* coef_K  = MatK.matC->get_ptr_coef_data();

   rExt* ptr_D  = D_coef.data();

   #pragma omp parallel num_threads(nthreads)
   {
      // Retrive processor ID
      iReg myid = omp_get_thread_num();

      // Get the first row of the current processor
      iReg firstrow = vecstart[myid];

      // Get the number of rows of the current processor
      iReg nrows_L = vecstart[myid+1] - firstrow;
      try
      {
         cpt_D_coef_loc(firstrow,nrows_L,iat_K,ja_K,coef_K,iat_Ct,ja_Ct,coef_Ct,ptr_D+firstrow);
      }
      catch( linsol_error )
      {
         throw linsol_error("RevAugm::cpt_D_coef","Error in cpt_D_coef_loc");
      }
   }
}

//----------------------------------------------------------------------------------------

void RevAugm::print_OperStats( FILE *ofile )
{
   if( params.type_precK == 0 )
   {
      aFSAI* prec = dynamic_cast<aFSAI*>( prec_K );
      prec->print_OperStats( ofile );
   }
   else if( params.type_precK == 1 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->print_OperStats( ofile );
   }
   else if( params.type_precK == 2 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->print_OperStats( ofile );
   }
}

//----------------------------------------------------------------------------------------

void RevAugm::print_Times( FILE *ofile )
{
   if( params.type_precK == 0 )
   {
      aFSAI* prec = dynamic_cast<aFSAI*>( prec_K );
      prec->print_Times( ofile );
   }
   else if( params.type_precK == 1 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->print_Times( ofile );
   }
   else if( params.type_precK == 2 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->print_Times( ofile );
   }
}

//----------------------------------------------------------------------------------------

void RevAugm::setDefaultParametersBlock11( RevAugm_params & paramsOut )
{
   // FSAI case
   {
     aFSAI* prec = dynamic_cast<aFSAI*>( prec_K );
     prec->setDefaultParameters( paramsOut.precK_aFSAI_params );
   }

   // AMG case
   {
     aAMG* prec = dynamic_cast<aAMG*>( prec_K );
     prec->setDefaultParameters( paramsOut.precK_aAMG_params );
   }
}

//----------------------------------------------------------------------------------------

void RevAugm::setParametersBlock11( RevAugm_params const & paramsIn )
{
   if( params.type_precK == 0 )
   {
      aFSAI* prec = dynamic_cast<aFSAI*>( prec_K );
      prec->Set_Parameters( &( paramsIn.precK_aFSAI_params ) );
   }
   else if( params.type_precK == 1 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->Set_Parameters( &( paramsIn.precK_aAMG_params ) );
   }
   else if( params.type_precK == 2 )
   {
      aAMG* prec = dynamic_cast<aAMG*>( prec_K );
      prec->Set_Parameters( &( paramsIn.precK_aAMG_params ) );
   }
}

//----------------------------------------------------------------------------------------

void RevAugm::setDefaultParameters( RevAugm_params & paramsOut )
{
   paramsOut.type_precK = 1;
   paramsOut.gamma = 1.0;
   setDefaultParametersBlock11( paramsOut );
}

//----------------------------------------------------------------------------------------

void RevAugm::Set_Parameters( const void *input_params )
{
   // cast the pointer
   RevAugm_params* input_RA_params = (RevAugm_params*) input_params;

   params.type_precK = input_RA_params->type_precK;
   params.gamma = input_RA_params->gamma;
   setParametersBlock11( *input_RA_params );
}

//----------------------------------------------------------------------------------------

iExt RevAugm::Get_FlopEst() const
{
  return prec_K->Get_FlopEst() +
         MatC->Get_Globalnterm() +
         MatCt->Get_Globalnterm() +
         2*MatCt->get_nrows();
}

//----------------------------------------------------------------------------------------
