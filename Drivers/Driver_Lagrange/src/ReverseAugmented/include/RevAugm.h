#pragma once

#include "Preconditioner.h"          // to use: Preconditioner

#include "DSMat.h"                   // to use: DSMat
#include "aFSAI.h"                   // to use: aFSAI
#include "aAMG.h"                   // to use: aFSAI

#include "RevAugm_params.h"

/**
 * class RevAugm
 * @brief This class is used to manage the Reverse Augmented preconditioner.
 */
class RevAugm: public Preconditioner{

   //-------------------------------------------------------------------------------------

   // private members
   private:

   /**
    * @brief general aAMG parameters.
    */
   RevAugm_params params;

   /**
    * @brief Diagonal entries of the block[1,1] preconditioner.
    */
   VEC<rExt> D_coef;

   /**
    * @brief Coupling matrix (block [0,1])
    */
   DSMat *MatC = nullptr;

   /**
    * @brief Coupling matrix (block [1,0])
    */
   DSMat *MatCt = nullptr;

   /**
    * @brief Schur complement
    */
   DSMat Schur;

   /**
    * @brief Preconditioner for Stiffness matrix (block [0,0]).
    */
   Preconditioner *prec_K = nullptr;

   //-------------------------------------------------------------------------------------
   // private functions
   private:

   /**
    * @brief Compute the coefficient of the diagonal matrix.
    */
   void cpt_D_coef( DSMat const & MatK,
                    DSMat const * const MatC,
                    DSMat const * const MatCt );

   //-------------------------------------------------------------------------------------

   // public members
   public:

   /**
    * @brief Creates an empty object.
    */
   RevAugm ();

   /**
    * @brief Deletes the object.
    */
   ~RevAugm();

   /**
    * @brief Computes the Reverse Augmented preconditioner in a distributed memory environment.
    * @param [in] mat_A system matrix.
    */
   void Compute ( DSMat& __restrict__ mat_A_in , DDMat *const __restrict__ V0 = nullptr );


   /**
    * @brief Prepares the data structure necessary for MxV.
    */
   void Prepare_MxV(const iReg nRHS);

   /**
    * @brief Undo the data structure necessary for MxV.
    */
   void UndoPrepare_MxV();


   /**
    * @brief Computes the Saddle Point Matrix by vector Product.
    * @param [in] x vector that is multiplied by the matrix.
    * @param [out] b vector where result is stored.
    * @param [in] Barrier_flag if true global Barrier is called before return.
    */
   void MxV(DDMat& __restrict__ x, DDMat& __restrict__ b,
            bool Barrier_flag);

   /**
    * @brief Get nrows.
    */
   iReg get_nrows() const { return MatC->get_nrows() + MatCt->get_nrows();};

   private:

   iGlo Get_Globalnrows_inner() const{ return MatC->Get_Globalnrows() + MatCt->Get_Globalnrows();};

   bool Check_inner() const { return true; };

   public:

   void fPrintAll_ASCII(const string &filename) const {};

   void MxM ( DSMat& __restrict__ Mat, DSMat& __restrict__ PrecMat ) {};

   void Set_Parameters ( const void *params );

   /**
    * @brief Returns an estimate of the number of non-zeroes of the preconditioner
    *        on this process
    */
   iExt Get_nnz () const { return -1; };

   /**
    * @brief Returns an estimate of the number of FLOP needed to perform the
    *        preconditioner application on this process
    */
   iExt Get_FlopEst () const;

   /**
    * @brief Set the type and default parameters for the block [0,0]
    *
    */
   void PrecK_set_type(const iReg type);

   /**
    * @brief Set the type and default parameters for the block [0,0]
    *
    */
   void set_gamma(const rExt gamma_in);

   /**
    * @brief Prints the operator statistics for the whole hierarchy.
    * @param [in] *ofile output file.
    */
   void print_OperStats ( FILE *ofile );

   /**
    * @brief Prints the times needed to compute the whole hierarchy.
    * @param [in] *ofile output file.
    */
   void print_Times ( FILE *ofile );

   /**
    * @brief Set default parameters
    * @param[out] paramsOut Reverse Augmented parameters
    */
   void setDefaultParameters( RevAugm_params & paramsOut );

   //-------------------------------------------------------------------------------------

   private:

   /**
    * @brief Set default parameters for block 11 preconditioner
    * @param[out] paramsOut Reverse Augmented parameters
    */
   void setDefaultParametersBlock11( RevAugm_params & paramsOut );

   /**
    * @brief Set parameters for block 11 preconditioner
    * @param[in] paramsIn Reverse Augmented parameters
    */
   void setParametersBlock11( RevAugm_params const & paramsIn );

   //-------------------------------------------------------------------------------------

   // public setting functions NOT useful for external users
   public:

   /**
    * @brief Sets the multigrid cycle type.
    */
   void Set_CycleType(const iReg __restrict__ input_CycleType){
      params.precK_aAMG_params.gen_params.CycleType = input_CycleType;
   }

   /**
    * @brief Sets the maximum number of levels
    */
   void Set_maxNlevel(const iReg __restrict__ input_maxNlevel){
      params.precK_aAMG_params.gen_params.maxNlevel = input_maxNlevel;
   }

   /**
    * @brief Sets the maximum size of the coarsen grid
    */
   void Set_maxCoarseSize(const iReg __restrict__ input_maxCoarseSize){
      params.precK_aAMG_params.gen_params.maxCoarseSize = input_maxCoarseSize;
   }

   //-------------------------------------------------------------------------------------

   // public functions to set TSPACE's parameters (NOT useful for external users)
   public:

   /**
    * @brief Sets the Tspace's number of test vectors.
    */
   void TSPACE_Set_ntv(const iReg __restrict__ input_ntv){
      params.precK_aAMG_params.TSP_params.ntv = input_ntv;
   }

   /**
    * @brief Sets the Eigensolver used to compute test vectors.
    */
   void TSPACE_Set_method(const iReg __restrict__ input_method){
      params.precK_aAMG_params.TSP_params.method = input_method;
   }

   /**
    * @brief Sets the Tspace's Maximum number of iterations.
    */
   void TSPACE_Set_maxITER(const iReg __restrict__ input_maxITER){
      params.precK_aAMG_params.TSP_params.maxITER = input_maxITER;
   }

   /**
    * @brief Sets the Tspace's exit tolerance.
    */
   void TSPACE_Set_exitTOL(const rExt __restrict__ input_exitTOL){
      params.precK_aAMG_params.TSP_params.exitTOL = input_exitTOL;
   }

   /**
    * @brief Sets the Tspace's orth_freq.
    */
   void TSPACE_Set_orth_freq(const iReg __restrict__ input_orth_freq){
      params.precK_aAMG_params.TSP_params.orth_freq = input_orth_freq;
   }

   /**
    * @brief Sets the Tspace's ritz_freq.
    */
   void TSPACE_Set_ritz_freq(const iReg __restrict__ input_ritz_freq){
      params.precK_aAMG_params.TSP_params.ritz_freq = input_ritz_freq;
   }

   //-------------------------------------------------------------------------------------

   // public functions to set COARSEN's parameters (NOT useful for external users)
   public:

   /**
    * @brief Sets the Coarsen's type.
    */
   void COARSEN_Set_type(const iReg __restrict__ input_type){
      params.precK_aAMG_params.COARS_params.type = input_type;
   }

   /**
    * @brief Sets the Coarsen's tolerance for filtering.
    */
   void COARSEN_Set_tau(const rExt __restrict__ input_tau){
      params.precK_aAMG_params.COARS_params.tau = input_tau;
   }

   /**
    * @brief Sets the Coarsen's tolerance for filtering.
    */
   void COARSEN_Set_tau_jump(const rExt __restrict__ input_tau_jump){
      params.precK_aAMG_params.COARS_params.tau_jump = input_tau_jump;
   }

   //-------------------------------------------------------------------------------------

   // public functions to set PROLONGATION's parameters (NOT useful for external users)
   public:

   /**
    * @brief Sets the Prolongation's type.
    */
   void PROLONGATION_Set_type(const iReg input_type){
      params.precK_aAMG_params.PROLO_params.type = input_type;
   }

   /**
    * @brief Sets the Prolongation's power pattern.
    */
   void PROLONGATION_Set_nPow(const iReg input_nPow){
      params.precK_aAMG_params.PROLO_params.nPow = input_nPow;
   }

   /**
    * @brief Sets the Prolongation's smoothing.
    */
   void PROLONGATION_Set_smooth(iReg __restrict__ input_smooth){
      params.precK_aAMG_params.PROLO_params.smooth = input_smooth;
   }

   /**
    * @brief Sets the Prolongartion's maxrownrm.
    */
   void PROLONGATION_Set_maxrownrm(const rExt __restrict__ input_maxrownrm){
      params.precK_aAMG_params.PROLO_params.maxrownrm = input_maxrownrm;
   }

   /**
    * @brief Sets the Prolongartion's maxcond.
    */
   void PROLONGATION_Set_maxcond(const rExt __restrict__ input_maxcond){
      params.precK_aAMG_params.PROLO_params.maxcond = input_maxcond;
   }

   /**
    * @brief Sets the Prolongartion's itmax_vol.
    */
   void PROLONGATION_Set_itmax_vol(const iReg __restrict__ input_itmax_vol){
      params.precK_aAMG_params.PROLO_params.itmax_vol = input_itmax_vol;
   }

   /**
    * @brief Sets the Prolongartion's tol_vol.
    */
   void PROLONGATION_Set_tol_vol(const rExt __restrict__ input_tol_vol){
      params.precK_aAMG_params.PROLO_params.tol_vol = input_tol_vol;
   }

   /**
    * @brief Sets the Prolongartion's eps.
    */
   void PROLONGATION_Set_eps(const rExt __restrict__ input_eps){
      params.precK_aAMG_params.PROLO_params.eps = input_eps;
   }

   /**
    * @brief Sets the Prolongartion's dist_min.
    */
   void PROLONGATION_Set_dist_min(const iReg __restrict__ input_dist_min){
      params.precK_aAMG_params.PROLO_params.dist_min = input_dist_min;
   }

   /**
    * @brief Sets the Prolongartion's dist_max.
    */
   void PROLONGATION_Set_dist_max(const iReg __restrict__ input_dist_max){
      params.precK_aAMG_params.PROLO_params.dist_max = input_dist_max;
   }

   /**
    * @brief Sets the Prolongartion's mmax.
    */
   void PROLONGATION_Set_mmax(const iReg __restrict__ input_mmax){
      params.precK_aAMG_params.PROLO_params.mmax = input_mmax;
   }

   //-------------------------------------------------------------------------------------

   // public functions to set FILTERING's parameters (NOT useful for external users)
   public:

   /**
    * @brief Sets Prol_Weight
    */
   void FILTER_Set_Prol_Weight(const rExt input_Prol_Weight){
      params.precK_aAMG_params.FILT_params.Prol_Weight = input_Prol_Weight;
   }

   /**
    * @brief Sets Prol_Tol
    */
   void FILTER_Set_Prol_Tol(const rExt input_Prol_Tol){
      params.precK_aAMG_params.FILT_params.Prol_Tol = input_Prol_Tol;
   }

   /**
    * @brief Sets Oper_Weight
    */
   void FILTER_Set_Oper_Weight(const rExt input_Oper_Weight){
      params.precK_aAMG_params.FILT_params.Oper_Weight = input_Oper_Weight;
   }

   //-------------------------------------------------------------------------------------

   // public functions to set SMOOTHER's parameters (NOT useful for external users)
   public:

   // general SMOOTHER's parameters ------------------------------------------------------

   /**
    * @brief Sets the Smoother's type.
    */
   void SMOOTHER_Set_type(const iReg __restrict__ input_type){
      params.precK_aAMG_params.SMOOTH_params.type = input_type;
   }

   /**
    * @brief Sets the Smoother's weight.
    */
   void SMOOTHER_Set_omega(const rExt __restrict__ input_omega){
      params.precK_aAMG_params.SMOOTH_params.omega = input_omega;
   }

   /**
    * @brief Sets the Maximum number of iterations to compute max eigenvalue (estimate omega).
    */
   void SMOOTHER_Set_itmax_maxeig(const iReg __restrict__ input_itmax_maxeig){
      params.precK_aAMG_params.SMOOTH_params.itmax_maxeig = input_itmax_maxeig;
   }

   /**
    * @brief Sets the Exit tolerance to compute max eigenvalue (estimate omega).
    */
   void SMOOTHER_Set_tol_maxeig(const rExt __restrict__ input_tol_maxeig){
      params.precK_aAMG_params.SMOOTH_params.tol_maxeig = input_tol_maxeig;
   }

   /**
    * @brief Sets the Number of pre-smoothing steps.
    */
   void SMOOTHER_Set_nupre(const iReg __restrict__ input_nupre){
      params.precK_aAMG_params.SMOOTH_params.nupre = input_nupre;
   }

   /**
    * @brief Sets the Number of post-smoothing steps.
    */
   void SMOOTHER_Set_nupost(const iReg __restrict__ input_nupost){
      params.precK_aAMG_params.SMOOTH_params.nupost = input_nupost;
   }

   // aFSAI SMOOTHER's parameters --------------------------------------------------------

   /**
    * @brief Sets the Requested power for the communication Graph.
    */
   void SMOOTHER_Set_aFSAI_TargetPower(const iReg __restrict__ input_afsai_TargetPower){
      params.precK_aAMG_params.SMOOTH_params.aFSAI_SMOOTHER_params.TargetPower = input_afsai_TargetPower;
   }

   /**
    * @brief Sets the Number of steps for computing the adaptive FSAI Smoother.
    */
   void SMOOTHER_Set_aFSAI_nstep(const iReg __restrict__ input_afsai_nstep){
      params.precK_aAMG_params.SMOOTH_params.aFSAI_SMOOTHER_params.nstep = input_afsai_nstep;
   }

   /**
    * @brief Sets the Step size for computing the adaptive FSAI Smoother.
    */
   void SMOOTHER_Set_aFSAI_step_size(const iReg __restrict__ input_afsai_step_size){
      params.precK_aAMG_params.SMOOTH_params.aFSAI_SMOOTHER_params.step_size = input_afsai_step_size;
   }

   /**
    * @brief Sets the Exit tolerance for computing the adaptive FSAI Smoother.
    */
   void SMOOTHER_Set_aFSAI_eps(const rExt __restrict__ input_afsai_eps){
      params.precK_aAMG_params.SMOOTH_params.aFSAI_SMOOTHER_params.eps = input_afsai_eps;
   }

   // aFSAI's parameters --------------------------------------------------------

   /**
    * @brief Sets the Requested power for the communication Graph.
    */
   void Set_TargetPower(const iReg __restrict__ input_afsai_TargetPower){
      params.precK_aFSAI_params.TargetPower = input_afsai_TargetPower;
   }

   /**
    * @brief Sets the Number of steps for computing the adaptive FSAI.
    */
   void Set_nstep(const iReg __restrict__ input_afsai_nstep){
      params.precK_aFSAI_params.nstep = input_afsai_nstep;
   }

   /**
    * @brief Sets the Step size for computing the adaptive FSAI.
    */
   void Set_step_size(const iReg __restrict__ input_afsai_step_size){
      params.precK_aFSAI_params.step_size = input_afsai_step_size;
   }

   /**
    * @brief Sets the Exit tolerance for computing the adaptive FSAI.
    */
   void Set_eps(const rExt __restrict__ input_afsai_eps){
      params.precK_aFSAI_params.eps = input_afsai_eps;
   }

};
