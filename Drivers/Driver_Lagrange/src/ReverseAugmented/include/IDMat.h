/**
 * @file IDMat.h
 * @brief This header is used to manage the Identity Matrix.
 * @date October 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent
 */

#pragma once

#include <mpi.h>
#include <omp.h>
#include <cstring>  // to use: memcpy
#include <iostream> // to use: sizeof

#include "Chronos.h"    // to use: Chronos
#include "precision.h"  // to use: iReg, iExt, rExt, type_MPI_iReg
#include "MatrixProd.h" // to use: MatrixProd

#include "linsol_error.h" // to throw erros
#include "mpi_error.h"    // to throw MPI_errors

using namespace std;

/**
 * class IDMat.
 * @brief This class is used to represent Identity Matrix.
 */
class IDMat: public MatrixProd{

   //-------------------------------------------------------------------------------------

   // Private members
   private:

  /**
   * @brief Number of rows.
   */
   iReg nrows = 0;

   //-------------------------------------------------------------------------------------

   // private functions
   private:

   /**
    * @brief Copy members.
    */
   void copy_IDMat ( IDMat& other ){
      copy_MatrixProd(other);
      this->nrows = other.nrows;
   }

   /**
    * @brief Evaluates the total number of rows.
    */
   iGlo Get_Globalnrows_inner() const;

   /**
    * Internal check of consistency of the IDMat instance.
    */
   bool Check_inner() const {return true;};

   //-------------------------------------------------------------------------------------

   // Public functions
   public:

   /**
    * @brief Copy operator.
    */
   IDMat& operator=( IDMat& other ){
      copy_IDMat(other);
      return *this;
   }

  /**
   * @brief Creates an empty object.
   */
   IDMat (): nrows{0} {}

   /**
    * @brief Creates an identity matrix of given size
    */
   IDMat (const iReg input_nrows){

      // Check errors
      if(input_nrows < 0) throw linsol_error("IDMat::IDMat","nrows < 0");

      // Set the number of rows
      nrows = input_nrows;

   }

   /**
    * @brief Deletes the object.
    */
   ~IDMat() {nrows = 0;}

   /**
    * @brief Retrieve the number of rows.
    */
   iReg get_nrows(){return nrows;}
   iReg get_nrows() const {return nrows;}

   /**
    * @brief Create all the data structure needed for MxV. This function has no
    *        effect, it is inserted only for compatibility with the parent class.
    * @param [in] nRHS number of Right Hand Sides.
    */
   void Prepare_MxV(const iReg nRHS) { Prepared_flag = true; nRHS_prep = nRHS; }

   /**
    * @brief Undo the data structure necessary for MxV. This function has no
    *        effect, it is inserted only for compatibility with the parent class.
    */
   void UndoPrepare_MxV() { Prepared_flag = false; nRHS_prep = 0; }

   /**
     * @brief Computes the Matrix by vector Product.
     * @param [in] x vector that is multiplied by the matrix.
     * @param [out] b vector where result is stored.
     * @param [in] Barrier_flag if true global Barrier is called before return.
     */
   void MxV(DDMat& __restrict__ x, DDMat& __restrict__ b, bool Barrier_flag);

};
