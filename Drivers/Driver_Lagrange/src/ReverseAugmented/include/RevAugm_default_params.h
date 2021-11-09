/**
 * @file RevAugm_default_params.h
 * @brief This header is used to manage the default Reverse Augmented parameters.
 * @date April 2021
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @author Matteo Frigo
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent
 */

#pragma once

/**
 * @brief default preconditioner type for block [0,0].
 * 0 aFSAI
 * 1 aAMG mech
 * 2 aAMG cfd
 */
#define RevAug_default_type_precK 1

/**
 * @brief default Relaxation parameter.
 */
#define RevAugm_default_gamma  1.0
