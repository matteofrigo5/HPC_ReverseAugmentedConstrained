#include "cpt_D_coef_loc.h"

#if defined (LAPACKE)
   #include "cblas.h"    // to use: dnrm2,ddot
   #include "lapacke.h"  // to use: dpotrf,dpotrs
#endif
#if defined (MKL)
   #include <mkl_lapacke.h>  // to use: dpotrf,dpotrs
#endif

#include <cmath>          // to use: abs,sqrt

#include "linsol_error.h" // to throw linsol errors
#include "VEC.h"          // to use: VEC
#include "inl_blas1.h"    // to use: inl_ddot

void cpt_D_coef_loc(const iReg firstrow, const iReg nrows_L, const iExt * const iat_K, const iReg * const ja_K,
                    const rExt* const coef_K, const iExt * const iat_Ct, const iReg * const ja_Ct,
                    const rExt * coef_Ct, rExt* ptr_D_loc)
{
   VEC<rExt> vec_K_loc;
   VEC<rExt> vec_C_col;
   VEC<rExt> vec_SVD,vec_Work;

   for (iReg k=firstrow; k<firstrow+nrows_L; k++){

      // Gather row/col indeces in ind_ebe and Bt_row/B_col in Bt_row
      iReg nnz_row = iat_Ct[k+1]-iat_Ct[k];

      if (nnz_row != 0) {
         // Initialize matrix to be extracted and vector
         try {
             vec_K_loc.assign(nnz_row*nnz_row,0.0);
             vec_C_col.assign(nnz_row,0.0);
             vec_SVD.assign(nnz_row,0.0);
             vec_Work.assign(nnz_row,0.0);
          } catch (linsol_error) {
             throw linsol_error("cpt_D_coef_loc","Error allocating");
          }

         rExt* K_loc = vec_K_loc.data();
         rExt* C_col = vec_C_col.data();
         rExt* svd = vec_SVD.data();
         rExt* work = vec_Work.data();
         const rExt* Ct_row;
         const iReg*  Ct_ind;

         // Gather row/col indeces in ind_ebe and Bt_row/B_col in Bt_row
         Ct_ind = ja_Ct+iat_Ct[k];
         Ct_row = coef_Ct + iat_Ct[k];

         // Extract K_ebe
         for (iReg  i=0; i<nnz_row; i++){
            iReg ii = Ct_ind[i];
            C_col[i] = Ct_row[i];
            for (iReg j=i; j<nnz_row; j++){
               iReg jj = Ct_ind[j];
               iExt ind_loc = iat_K[ii];
               while (ind_loc<iat_K[ii+1]) {
                  if (ja_K[ind_loc] == jj){
                     K_loc[i+j*nnz_row] = coef_K[ind_loc];
                     K_loc[j+i*nnz_row] = coef_K[ind_loc];
                     break;
                  }
                  ind_loc++;
               }
            }
         }

         lapack_int info;

         // Compute max singular value eigenvalue
         info = LAPACKE_dgesvd(LAPACK_COL_MAJOR,'N','N',nnz_row,nnz_row,K_loc,nnz_row,svd,nullptr,nnz_row,nullptr,nnz_row,work);
         if(info != 0){throw linsol_error ("cpt_D_coef_loc","error in LAPACKE_dgeev");}

         rExt Knorm2 = 0.0;
         for(iReg i=0; i<nnz_row; i++){
            if(Knorm2<svd[i]) Knorm2 = svd[i];
         }

	 rExt den = inl_ddot(nnz_row,Ct_row,1,C_col,1);
	 if (den >= std::numeric_limits<rExt>::denorm_min())
            ptr_D_loc[k-firstrow] = Knorm2 / den;
         else
            ptr_D_loc[k-firstrow] = 1.0;

      } else {
         ptr_D_loc[k-firstrow] = 1.0;
      }
   }
}
