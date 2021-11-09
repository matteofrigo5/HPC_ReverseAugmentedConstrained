#include "SadPointMat.h"

//----------------------------------------------------------------------------------------

// Creates an empty object.
SadPointMat::SadPointMat(){
}

//----------------------------------------------------------------------------------------

// Deletes the object.
SadPointMat::~SadPointMat(){

   MatK.~DSMat();
   MatC.~DSMat();
   MatCt.~DSMat();

}

//----------------------------------------------------------------------------------------

// Prepares the data structure necessary for MxV.
void SadPointMat::Prepare_MxV(const iReg nRHS){

   // Check if the matrix is prepared
   if ( Prepared_flag == true ){
      if (nRHS <= get_nRHS_prep()){
         return;
      } else {
         UndoPrepare_MxV();
      }
   }

   MatK.Prepare_MxV(nRHS);
   MatC.Prepare_MxV(nRHS);
   MatCt.Prepare_MxV(nRHS);

   // set flag
   Prepared_flag = true;
   nRHS_prep = nRHS;

}

//----------------------------------------------------------------------------------------

// Compute MxV
void SadPointMat::MxV(DDMat& __restrict__ x, DDMat& __restrict__ b,
                bool Barrier_flag){

   iReg n1,n2;//,nn;
   iReg nRHS;
   DDMat b1,b2,x1,x2;

   // retrieve number of threads
   type_OMP_iReg nthreads = Chronos.Get_inpNthreads();

   // Retrive dimensions
   n1   = MatK.get_nrows();
   n2   = MatCt.get_nrows();
   nRHS = b.get_ncols();
   //nn   = b.get_nrows();
   DDMat tmp=DDMat(n1,nRHS);

   rExt* ptr_b = b.get_ptr_coef_data();
   b1.overlap(n1,nRHS,ptr_b);
   b2.overlap(n2,nRHS,ptr_b+n1);

   rExt* ptr_x = x.get_ptr_coef_data();
   x1.overlap(n1,nRHS,ptr_x);
   x2.overlap(n2,nRHS,ptr_x+n1);


   MatK.MxV(x1,b1,Barrier_flag);
   MatC.MxV(x2,tmp,Barrier_flag);
   MatCt.MxV(x1,b2,Barrier_flag);

   nthreads = static_cast<type_OMP_iReg>(min(static_cast<iReg>(nthreads),n1*nRHS));
   rExt* ptr_b1  = b1.get_ptr_coef_data();
   rExt* ptr_tmp = tmp.get_ptr_coef_data();

   // b1 <- b1 + tmp
   #pragma omp parallel for num_threads(nthreads)
   for(iReg k=0; k<n1*nRHS; ++k){
      ptr_b1[k]+=ptr_tmp[k];
   }

   b1.assign_null();
   b2.assign_null();
   x1.assign_null();
   x2.assign_null();

}

//----------------------------------------------------------------------------------------

// Undo the data structure necessary for MxV.
void SadPointMat::UndoPrepare_MxV(){

   MatK.UndoPrepare_MxV();
   MatC.UndoPrepare_MxV();
   MatCt.UndoPrepare_MxV();
   Prepared_flag = false;
}

//----------------------------------------------------------------------------------------
