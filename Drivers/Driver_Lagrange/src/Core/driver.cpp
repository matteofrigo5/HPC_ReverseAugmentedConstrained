#include <iostream>
#include <mpi.h>
#include <omp.h>

// include Chronos methods
#include "precision.h"        // use
#include "linsol_error.h"     // to throw general errors
#include "mpi_error.h"        // to throw mpi errors
#include "Chronos.h"          // to set the library environment
#include "DSMat.h"            // to use the DSMat object
#include "SadPointMat.h"      // to use the SadPointMat object
#include "DDMat.h"            // to use the DDMat object
#include "IDMat.h"            // to use the IDMat object
#include "RevAugm.h"          // to use the Reverse Augmented preconditioner
#include "BiCGstab.h"         // to use BiCGstab solver
#include "GMRES.h"            // to use GMRES solver
#include "dnrm2_par.h"        // to use euclidean norm

// include Chronos methods to read parameters and objects
//#include "next_line.h"
#include "check_file_exist.h"
#include "read_DSMat.h"
#include "read_DSMat_ascii.h"
#include "read_DDMat.h"
#include "read_DDMat_ascii.h"

// include Chronos methods to partitioning matrices
#include "Create_InvPart.h"

// include pugixml for XML parsing
#include "pugixml.hpp"

// include files for reading functions
#include "check_XML.h"
#include "read_fnames.h"
#include "read_general_parms.h"
#include "set_solver_parms.h"
#include "get_solver_kind.h"
#include "set_prec_parms.h"

//----------------------------------------------------------------------------------------

int main(int argc,char **argv){

   // --- MPI initialization -------------------------------------------------------------
   type_MPI_iReg mpisupport,rank,size;

   type_MPI_iReg ierr = MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&mpisupport);
   if(ierr != 0) {
      mpi_error("driver","MPI_Init_thread","Unknown");
      return 0;
   }
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

   // --- Init Chronos -------------------------------------------------------------------
   if (!Chronos.init()){
      cout << "Error in Chronos initialization" << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
   }

   // --- Open Log file ------------------------------------------------------------------
   FILE *logfile = fopen("CHRONOS.log","w");

   if(rank==0){
      cout << endl << "MPI Initialization done (" << size << " processes)" << endl;
      fprintf(logfile,"\nMPI Initialization done\n%i processes\n", size);
   }

   // --- Check Multi-Thread Communications ------------------------------------------------
   #if defined (MultThrCom)
   if(rank==0){cout << endl << "MULTIPLE THREAD COMMUNICATIONS" << endl;}
   if(rank==0){cout << endl << "NOT YET IMPLEMENTED" << endl;}
   MPI_Finalize();
   return 0;
   #else
   if(rank==0){cout << endl << "FUNNELED THREAD COMMUNICATIONS" << endl;}
   #endif

   // --- Read elaboration parameters ----------------------------------------------------
   iReg ierr_xml;

   if(rank==0){
      cout<< endl << "Read elaboration parameters"<<endl;
      fprintf(logfile,"\nRead elaboration parameters\n");
   }
   MPI_Barrier(MPI_COMM_WORLD);

   // Parse XML configuration file
   pugi::xml_node node;
   pugi::xml_attribute attr;
   pugi::xml_document input_xml;
   pugi::xml_parse_result result = input_xml.load_file(argv[1]);
   if (result){
      if (rank==0)
         std::cout << "XML [" << argv[1] << "] parsed without errors\n" << endl;
   } else {
      if (rank==0){
         std::cout << "XML [" << "config.xml" << "] parsed with errors, attr value: ["
                   << input_xml.child("node").attribute("attr").value() << "]\n";
         std::cout << "Error description: " << result.description() << "\n";
         std::cout << "Error offset: " << result.offset << " (error at [..."
                   << (argv[1] + result.offset) << "]\n\n";
      }
      linsol_error("Driver","wrong XML format");
      MPI_Finalize();
      return 0;
   }

   // Check correctness of XML inputs
   ierr_xml = check_XML(input_xml);
   if (ierr_xml != 0) return 0;

   // Read names of the other input files
   string K0_fname, C0_fname, Ct0_fname, RHS0_fname, TS0_fname;
   ierr_xml = read_fnames(rank,input_xml,K0_fname,C0_fname,Ct0_fname, RHS0_fname, TS0_fname);
   if (ierr_xml != 0) return 0;

   // Read general parameters
   type_OMP_iReg nthreads;
   iReg verbosity;
   bool PART;
   ierr_xml = read_general_parms(rank,input_xml,nthreads,verbosity,PART);
   if (ierr_xml != 0) return 0;

   // Get solver kind
   iReg const solverID = get_solver_kind( input_xml );

   // Set solver parameters
   iReg RHS_Flag;

   LinSolver * solver;
   switch( solverID )
   {
     case( 1 ):
     {
       solver = new BiCGstab();
       ierr_xml = set_solver_parms(rank,input_xml,*solver,RHS_Flag);
     }
     case( 2 ):
     {
       solver = new GMRES();
       ierr_xml = set_solver_parms(rank,input_xml,*solver,RHS_Flag);
     }
   }
   if (ierr_xml != 0) return 0;

   // Set preconditioner parameters
   RevAugm PREC_RA;
   try {
      set_prec_parms(rank,input_xml,PREC_RA);
   } catch (linsol_error){
      linsol_error("Driver","set prec parms");
      MPI_Finalize();
      return 0;
   }

   // --- Set-up of the Chronos Environment ----------------------------------------------
   Chronos.Set_inpNthreads(nthreads);
   Chronos.Set_Verbosity(verbosity);
   bool const printTimes = true;
   Chronos.Set_GetTime(printTimes);

   // +++ TRUE execution SCOPE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   {

      rExt T_start;

      // --- Read Matrix -----------------------------------------------------------------
      // NOTES: this function reads a sparse matrix stored in coordinate format (row,col,coef)
      // and transform it into a DSMat object
      if(rank==0){
         cout << endl << "Read Matrices" << endl;
         fprintf(logfile,"\nRead Matrices\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      iExt globalNterm = 0;

      SadPointMat A0_SP;
      try {

         if(rank==0){
            cout  << "   Read MatK" << endl;
            fprintf(logfile,"   Read MatK\n");
         }
         //DSMat Kscr = read_DSMat_ascii(K0_fname);
         DSMat Kscr = read_DSMat(K0_fname);
         A0_SP.MatK = Kscr;
         globalNterm += Kscr.Get_Globalnterm();
         //A0_SP.MatK.fPrintAll_PARDDMAT_ASCII("K_prova.ascii");

         if(rank==0){
            cout << "   Read MatC" << endl;
            fprintf(logfile,"   Read MatC\n");
         }
         //DSMat Cscr = read_DSMat_ascii(C0_fname,true);
         DSMat Cscr = read_DSMat(C0_fname,true);
         A0_SP.MatC = Cscr;
         globalNterm += Cscr.Get_Globalnterm();
         //A0_SP.MatC.fPrintAll_PARDDMAT_ASCII("C_prova.ascii");

         if(rank==0){
            cout << "   Read MatCt" << endl << endl;
            fprintf(logfile,"   Read MatCt\n\n");
         }
         //DSMat Ctscr = read_DSMat_ascii(Ct0_fname,true);
         DSMat Ctscr = read_DSMat(Ct0_fname,true);
         A0_SP.MatCt = Ctscr;
         globalNterm += Ctscr.Get_Globalnterm();
         //A0_SP.MatCt.fPrintAll_PARDDMAT_ASCII("Ct_prova.ascii");
      } catch (linsol_error) {
         linsol_error("Driver","reading A0_SP");
         MPI_Finalize();
         return 0;
      }

      iGlo nrGlo = A0_SP.Get_Globalnrows();
      iGlo ntGlo = A0_SP.Get_Globalnterm();
      if (rank == 0 ) {
         cout << "+ nrGlo    = " << nrGlo << endl;
         cout << "+ ntGlo    = " << ntGlo << endl;
         fprintf(logfile,"+ nrGlo    = %12ld\n",nrGlo);
         fprintf(logfile,"+ ntGlo    = %12ld\n",ntGlo);
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "reading time [s] = " << T_start << endl;
         fprintf(logfile,"reading time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // retrieve number of rows belonging to the MPI rank
      iReg my_nrows_11 = A0_SP.MatK.get_nrows();
      iReg my_nrows = A0_SP.get_nrows();

      // --- Treat Boundary conditions ---------------------------------------------------
      // NOTES: this function modifies a DSMat object by detecting rows/columns associated to
      // Dirichlet boundary conditions and setting to those rows/columns a diagonal entry equal
      // to the maximum diagonal entry of the whole matrix. In such a way the minimum eigenvalue
      // of the matrix is not erroneously associated to the arbitrary value given to the
      // diagonal entries of boundary dofs.
      if(rank==0){
         cout << endl << "Treat BC" << endl;
         fprintf(logfile,"\nTreat BC\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      A0_SP.MatK.TreatBC();

      //-------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "BC time [s] = " << T_start << endl;
         fprintf(logfile,"BC time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // --- Read the Test Space ---------------------------------------------------------
      // NOTES: this function reads the test space as a dense matrix and transform it
      // into a DDMat object
      if(rank==0){
         cout << endl << "Read Test Space" << endl;
         fprintf(logfile,"\nRead Test Space\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      DDMat TS0;
      if ( check_file_exist(TS0_fname) ){
         // If the file exists read it
         try {
            //DDMat TSscr = read_DDMat_ascii(my_nrows_11,TS0_fname);
            DDMat TSscr = read_DDMat(my_nrows_11,TS0_fname);
            TS0 = TSscr;
         } catch (linsol_error) {
            linsol_error("Driver","reading TS0");
            MPI_Finalize();
            return 0;
         }
      } else {
         // If the file do not exist initialize an empty test space
         TS0.resize(my_nrows_11,0);
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "reading TS time [s] = " << T_start << endl;
         fprintf(logfile,"reading TS time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // --- Create Reference Solution and Right Hand Side (RHS) -------------------------

      // NOTES: it creates a random reference solution as a DDMat


      if(rank==0){
         cout << endl << "Create Reference Solution and RHS" << endl;
         fprintf(logfile,"\nCreate Reference Solution and RHS\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      iReg nRHS;
      DDMat RefSol;
      DDMat RHS0;

      if ( RHS_Flag == 0 || RHS_Flag == 2){

         // Unitary or random solution

         nRHS = 1;
         RefSol.resize(my_nrows,nRHS);
         RHS0.resize(my_nrows,nRHS);

         // Set the reference solution
         rExt* ptr_RefSol = RefSol.get_ptr_coef_data();
         rExt add = 1.e+1*static_cast<rExt>(rank);
         iReg k = 0;
         if (RHS_Flag == 2){
            for (iReg i = 0; i < my_nrows; i++){
               for (iReg j = 0; j < nRHS; j++){
                  ptr_RefSol[k++] = static_cast<rExt>(i) + add;
               }
            }
         } else {
            for (iReg i = 0; i < my_nrows*nRHS; i++){
               ptr_RefSol[k++] = 1.0;
            }
         }

         // Create the corresponding RHS
         A0_SP.Prepare_MxV(nRHS);
         A0_SP.MxV(RefSol,RHS0,true);
         A0_SP.UndoPrepare_MxV();

      } else if ( RHS_Flag == 1 ) {

         // Unitary RHS
         nRHS = 1;
         RHS0.resize(my_nrows,nRHS);
         rExt *ptr_RHS = RHS0.get_ptr_coef_data();
         for (iReg i = 0; i < my_nrows*nRHS; i++){
            ptr_RHS[i] = 1.0;
         }

      } else if (RHS_Flag == 3 ){

         // Read the RHS from an external file
         if (!check_file_exist(RHS0_fname)){
            linsol_error("Driver","External RHS file is not present");
            MPI_Finalize();
            return 0;
         }
         // If the file exists read it
         //DDMat tmp = read_DDMat_ascii(my_nrows,RHS0_fname);
         DDMat tmp = read_DDMat(my_nrows,RHS0_fname);
         RHS0 = tmp;

      } else {
         linsol_error("Driver","Wrong value of RHS_Flag");
         MPI_Finalize();
         return 0;
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "rhs time [s] = " << T_start << endl;
         fprintf(logfile,"rhs time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // --- Partitioning Matrix ---------------------------------------------------------
      // NOTES: the function 'PartKway' creates through METIS a partition of the DSMat object.
      // Then, the function 'Redist' redistributes the entries currently stored locally
      // into this MPI rank to all the other ranks according to the partition generated by METIS.
      if(rank==0){
         cout << endl << "Partitioning Matrix" << endl;
         fprintf(logfile,"\nPartitioning Matrix\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      //DSMat *A;
      SadPointMat * A_SP;
      VEC<iGlo> vec_part_K, vec_inv_part_K;
      VEC<iGlo> vec_iperm_K, vec_inv_iperm_K;
      VEC<iGlo> vec_part_Ct, vec_inv_part_Ct;
      VEC<iGlo> vec_iperm_Ct, vec_inv_iperm_Ct;
      if (PART){
         if(rank==0){
            cout << endl << "   Making partitioning" << endl;
            fprintf(logfile,"\n   Making partitioning\n");
         }
         try{
            A0_SP.MatK.PartKway(rank,size,vec_part_K,vec_iperm_K);
            A0_SP.MatCt.PartRevAug(rank,size,vec_part_Ct,vec_iperm_Ct);
         } catch (linsol_error) {
            linsol_error("Driver","partitioning A0");
            MPI_Finalize();
            return 0;
         } catch (mpi_error) {
            linsol_error("Driver","MPI error partitioning A0");
            MPI_Finalize();
            return 0;
         }

         try{
            Create_InvPart(rank,size,my_nrows_11,vec_part_K,vec_iperm_K,vec_inv_part_K,
                           vec_inv_iperm_K);
            Create_InvPart(rank,size,my_nrows-my_nrows_11,vec_part_Ct,vec_iperm_Ct,vec_inv_part_Ct,
                           vec_inv_iperm_Ct);
         } catch (linsol_error) {
            linsol_error("Driver","creating inverse partitioning");
            MPI_Finalize();
            return 0;
         } catch (mpi_error) {
            linsol_error("Driver","MPI error creating inverse partitioning");
            MPI_Finalize();
            return 0;
         }

         try{
            A_SP = new SadPointMat();
            if(rank==0){
               cout << "   Partitioning MatK" << endl;
               fprintf(logfile,"   Partitioning MatK\n");
            }
            A0_SP.MatK.Redist(rank,size,vec_part_K,vec_iperm_K,A_SP->MatK);

            if(rank==0){
               cout << "   Partitioning MatCt" << endl;;
               fprintf(logfile,"   Partitioning MatCt\n");
            }
            A0_SP.MatCt.Redist_Rec(rank,size,vec_part_Ct,vec_iperm_Ct,vec_part_K,vec_iperm_K,A_SP->MatCt);

            if(rank==0){
               cout << "   Partitioning MatC" << endl << endl ;
               fprintf(logfile,"   Partitioning MatC\n\n");
            }
            A0_SP.MatC.Redist_Rec(rank,size,vec_part_K,vec_iperm_K,vec_part_Ct,vec_iperm_Ct,A_SP->MatC);

            A0_SP.~SadPointMat();
         } catch (linsol_error) {
            linsol_error("Driver","redistributing A0");
            MPI_Finalize();
            return 0;
         } catch (mpi_error) {
            linsol_error("Driver","MPI error redistributing A0");
            MPI_Finalize();
            return 0;
         }
      } else {
         A_SP = &A0_SP;
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "redist time [s] = " << T_start << endl;
         fprintf(logfile,"redist time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // --- Partitioning TS0 ------------------------------------------------------------
      // NOTES: the function 'Redist' redistributes the Test Space entries according to
      // the partition generated by METIS
      if(rank==0){
         cout << endl << "Partitioning TS0" << endl;
         fprintf(logfile,"\nPartitioning TS0\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //-----------------------------
      // Start timing              //
      MPI_Barrier(MPI_COMM_WORLD); //
      T_start = MPI_Wtime();       //
      //-----------------------------

      DDMat *TS;
      if (PART){
         try{
            TS = new DDMat();
            TS0.Redist(rank,size,vec_part_K,vec_iperm_K,*TS);
         } catch (linsol_error) {
            linsol_error("Driver","redistributing TS0");
            MPI_Finalize();
            return 0;
         } catch (mpi_error) {
            linsol_error("Driver","MPI error redistributing TS0");
            MPI_Finalize();
            return 0;
         }
      } else {
         TS = &TS0;
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "redist TS time [s] = " << T_start << endl;
         fprintf(logfile,"redist TS time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // --- Partitioning RHS ------------------------------------------------------------
      // NOTES: the function 'Redist' redistributes the RHS entries according to
      // the partition generated by METIS
      if(rank==0){
         cout << endl << "Partitioning RHS" << endl;
         fprintf(logfile,"\nPartitioning RHS\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      DDMat *RHS;
      RHS = new DDMat();
      if (PART){
         try{
            iReg n1,n2;
            DDMat RHS0_1,RHS0_2;
            DDMat RHS_1,RHS_2;
            n1 = my_nrows_11;
            n2 = my_nrows - my_nrows_11;
            rExt* ptr_RHS0 = RHS0.get_ptr_coef_data();
            RHS0_1.overlap(n1,nRHS,ptr_RHS0);
            RHS0_2.overlap(n2,nRHS,ptr_RHS0+n1);
            RHS0_1.Redist(rank,size,vec_part_K,vec_iperm_K,RHS_1);
            RHS0_2.Redist(rank,size,vec_part_Ct,vec_iperm_Ct,RHS_2);
            n1 = RHS_1.get_nrows();
            n2 = RHS_2.get_nrows();
            rExt* ptr_RHS_1 = RHS_1.get_ptr_coef_data();
            rExt* ptr_RHS_2 = RHS_2.get_ptr_coef_data();
            RHS->resize(n1+n2,nRHS);
            rExt* ptr_RHS = RHS->get_ptr_coef_data();
            for (iReg i=0; i<n1; i++)
               ptr_RHS[i] = ptr_RHS_1[i];
            for (iReg i=0; i<n2; i++)
               ptr_RHS[n1+i] = ptr_RHS_2[i];
            RHS0_1.assign_null();
            RHS0_2.assign_null();
         } catch (linsol_error) {
            linsol_error("Driver","redistributing RHS0");
            MPI_Finalize();
            return 0;
         } catch (mpi_error) {
            linsol_error("Driver","MPI error redistributing RHS0");
            MPI_Finalize();
            return 0;
         }
      } else {
         RHS = &RHS0;
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "redist rhs time [s] = " << T_start << endl;
         fprintf(logfile,"redist rhs time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      // Compute Reverse Augmented Preconditioner

      if(rank==0){
         cout << endl << "Reverse Augmented Computation" << endl;
         fprintf(logfile,"\nReverse Augmented Computation\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      try {
         //PREC_RA.Compute(A0_SP,&TS0);
         PREC_RA.Compute(*A_SP,TS);
      } catch (linsol_error){
         linsol_error("Driver","computing Reverse Augmented");
         MPI_Finalize();
         return 0;
      }

      //--------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "setup time [s] = " << T_start << endl;
         fprintf(logfile,"setup time [s] = %12.4e\n",T_start);
      }
      //--------------------------------------------------

      rExt const RA_time = T_start;

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      // Declare a dummy preconditioner
      IDMat dummy_PREC = IDMat(A_SP->get_nrows());

      // set solver
      try
      {
        solver->Set_Solver(*A_SP,*RHS,dummy_PREC,PREC_RA);
      }
      catch( linsol_error )
      {
        linsol_error( "Driver", "setting solver" );
        MPI_Finalize();
        return 0;
      }

      //---------------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "set solver time [s] = " << T_start << endl;
         fprintf(logfile,"set solver time [s] = %12.4e\n",T_start);
      }
      //---------------------------------------------------------

      // --- Solve -----------------------------------------------------------------------
      if(rank==0){
         cout << endl << "Solve" << endl;
         fprintf(logfile,"\nSolve\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //--------------------------
      // Start timing
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime();
      //--------------------------

      DDMat SOL(A_SP->get_nrows(),nRHS);

      // set solver
      try
      {
         string ConvProf = "./SOLVER_convProf";
         solver->Solve(*A_SP,*RHS,SOL,dummy_PREC,PREC_RA, &ConvProf);
      }
      catch( linsol_error )
      {
        linsol_error( "Driver", "solving solver" );
        MPI_Finalize();
        return 0;
      }

      //---------------------------------------------------------
      // Get time
      MPI_Barrier(MPI_COMM_WORLD);
      T_start = MPI_Wtime() - T_start;
      if(rank==0){
         cout << "solve time [s] = " << T_start << endl;
         fprintf(logfile,"solve time [s] = %12.4e\n",T_start);
      }
      //---------------------------------------------------------

      rExt SOLVER_time = T_start;

      iReg ITER = solver->get_ITER();
      rExt const * ptr_normRES = solver->get_ptr_normRES();
      if(rank==0){
         cout << "+ n Iter   = " << ITER << endl;
         fprintf(logfile,"+ n Iter   = %12d\n",ITER);
         for(iReg i = 0; i < nRHS; i++){
            cout << "  iRHS     = " << i << endl;
            cout << "  normRES  = " << ptr_normRES[i] << endl;
            fprintf(logfile,"  iRHS     = %12d\n",i);
            fprintf(logfile,"  normRES  = %12.4e\n",ptr_normRES[i]);
         }
      }

      //---------------------------------------------------------*/

      // Print preconditioner info
      PREC_RA.print_OperStats(nullptr);
      PREC_RA.print_OperStats(logfile);

      // Print preconditioner times
      PREC_RA.print_Times(nullptr);
      PREC_RA.print_Times(logfile);

      // --- Elaboration summary ---------------------------------------------------------
      if(rank==0){
         cout << endl << "Elaboration summary" << endl;
         fprintf(logfile,"\nElaboration summary\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      rExt const cost = static_cast< rExt >( PREC_RA.Get_FlopEst() ) / static_cast< rExt >( globalNterm );

      if(rank==0) {

         // total preconditioner time
         cout << "SETUP time [s]  = " << RA_time << endl;
         cout << "SOLVER time [s] = " << SOLVER_time << endl;
         cout << "SOLVER n Iter   = " << ITER        << endl;
         printf( "Total cost      = %8.6f\n", cost );
         cout << endl;

         fprintf(logfile,"SETUP time [s]  = %12.4e\n",RA_time);
         fprintf(logfile,"SOLVER time [s] = %12.4e\n",SOLVER_time);
         fprintf(logfile,"SOLVER n Iter   = %12i\n",  ITER);
         fprintf( logfile, "Total cost      = %8.6f\n", cost );
         fprintf(logfile,"\n");

         for(iReg i = 0; i < nRHS; i++){
            cout << "iRHS            = " << i                 << endl;
            cout << "normRES         = " << ptr_normRES[i] << endl;
            fprintf(logfile,"iRHS            = %12i\n",i);
            fprintf(logfile,"normRES         = %12.4e\n",ptr_normRES[i]);
         }

      }

      if( PART )
      {
        delete A_SP;
        delete TS;
        delete RHS;
      }

   // +++ END execution SCOPE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   }

   // Delete the objects
   delete solver;
   PREC_RA.~RevAugm();

   // --- Close Log file -----------------------------------------------------------------
   fclose(logfile);

   // --- MPI Finalize -------------------------------------------------------------------
   if(rank==0){ cout << endl << "MPI Finalize" << endl; }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;

}
