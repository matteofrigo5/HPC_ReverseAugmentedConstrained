/**
 * @file SadPointMat.h
 * @brief This header is used to manage the Saddle Point Matrix.
 * @date February 2021
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @author Matteo Frigo
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent
 */

#pragma once

#include "DSMat.h"

#include "Chronos.h"      // to use: Chronos
#include "linsol_error.h" // to throw erros
#include "mpi_error.h"    // to throw MPI_errors

/**
 * class SadPointMat.
 * @brief This class is used to represent Distributed Sparse Saddle Point problems.
 */
class SadPointMat: public DSMat{

   //-------------------------------------------------------------------------------------

   // private members
   private:

   // public members
   public:

   /**
    * @brief Stiffness matrix (block [0,0])
    */
   DSMat MatK;

   /**
    * @brief Coupling matrix (block [0,1])
    */
   DSMat MatC;

   /**
    * @brief Coupling matrix (block [1,0])
    */
   DSMat MatCt;

   //-------------------------------------------------------------------------------------

   // private functions
   private:
   /**
    * @brief Copy members.
    */
   void copy_SadPointMat ( SadPointMat& other ){
      try {
         // copy MatrixProd members
         copy_MatrixProd(other);

         // copy DSMat members
         this->MatK  = other.MatK;
         this->MatC  = other.MatC;
         this->MatCt = other.MatCt;
      } catch (linsol_error) {
         throw linsol_error ("SadPointMat::copy_SadPointMat","error in allocating memory");
      }
   }

   // public members
   public:

   /**
    * @brief Copy operator.
    */
   SadPointMat& operator=( SadPointMat& other ){
      copy_SadPointMat(other);
      return *this;
   }

   /**
    * @brief Creates an empty object.
    */
   SadPointMat ();

   /**
    * @brief Deletes the object.
    */
   ~SadPointMat();

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
   iReg get_nrows() const { return MatK.get_nrows() + MatCt.get_nrows(); }

   /**
    * Evaluates the total number of rows of the Saddle Point matrix.
    */
   iGlo Get_Globalnrows_inner() const override { return MatK.get_nrows() + MatCt.get_nrows();  }

   /**
    * Evaluates the total number of terms of the Saddle Point matrix.
    */
   iGlo Get_Globalnterm() const { return MatK.Get_Globalnterm()+MatC.Get_Globalnterm()+MatCt.Get_Globalnterm(); }

};

