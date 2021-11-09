/**
 * @file RevAugm_params.h
 * @brief This header is used to manage the Reverse Augmented parameters.
 * @date April 2021
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @author Matteo Frigo
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "precision.h"    // to use: iReg, rExt
#include "aAMG_params.h"  // to use: aAMG_params
#include "aFSAI_params.h" // to use: aFSAI_params

typedef aAMG_params aAMG_block11_params;

typedef aFSAI_params aFSAI_block11_params;

/**
 * @brief Structure to store Reverse Augmented parameters.
 */
struct RevAugm_params {

   /**
    * @brief Requested Power for the communication Graph.
    */
   iReg type_precK;

   /**
    * @brief Relaxation parameter.
    */
   rExt gamma;

   /**
    * @brief aFSAI parameters for block 11
    */
   aFSAI_block11_params precK_aFSAI_params;

   /**
    * @brief aAMG parameters for block 11
    */
   aAMG_block11_params precK_aAMG_params;

};

