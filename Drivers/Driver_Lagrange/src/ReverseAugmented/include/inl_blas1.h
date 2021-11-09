/**
 * @file inl_blas1.h
 * @brief This header is used to manage the LEVEL 1 BLAS.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <math.h>    // to use: sqrt
#include <algorithm> // to use: min,max

#include "precision.h" // to use: iReg, rExt

/**
 * @brief DCOPY
 */
inline void inl_dcopy(const iReg n, const rExt *const v1, const iReg k1,
                      rExt *const v2, const iReg k2){
   iReg j = 0;
   for (iReg i = 0; i < n; i += k1){
      v2[i] = v1[j];
      j += k2;
   }
   return;
}

/**
 * @brief DDOT
 */
inline rExt inl_ddot(const iReg n, const rExt *const v1, const iReg k1,
                     const rExt *const v2, const iReg k2){
   rExt ddot = 0.0;
   iReg j = 0;
   for (iReg i = 0; i < n; i += k1){
      ddot += v1[i]*v2[j];
      j += k2;
   }
   return ddot;
}

/**
 * @brief DDNRM2
 */
inline rExt inl_dnrm2(const iReg n, const rExt *const v1, const iReg k1){
   rExt dnrm2 = 0.0;
   for (iReg i = 0; i < n; i += k1) dnrm2 += (v1[i])*(v1[i]);
   return sqrt(dnrm2);
}

/**
 * @brief DAXPY
 */
inline void inl_daxpy(const iReg n, const rExt alpha, const rExt *const x,
                      const iReg k1, rExt *const y, const iReg k2){
   iReg j = 0;
   for (iReg i = 0; i < n; i += k1){
      y[j] = y[j] + alpha*x[i];
      j += k2;
   }
   return;
}


