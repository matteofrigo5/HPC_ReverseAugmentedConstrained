#include "precision.h"

void cpt_D_coef_loc(const iReg firstrow, const iReg nrows_L, const iExt * const iat_K, const iReg * const ja_K,
                    const rExt* const coef_K, const iExt * const iat_Ct, const iReg * const ja_Ct,
                    const rExt * coef_Ct, rExt* ptr_D_loc);
