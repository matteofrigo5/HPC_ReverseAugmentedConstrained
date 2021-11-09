#include "read_fnames.h"
#include "utilities.hpp"

iReg read_fnames(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                 string &K0_fname, string &C0_fname, string &Ct0_fname, string &RHS0_fname, string &TS0_fname){

   pugi::xml_node const node = input_xml.child( "Chronos" ).child( "Files" );

   std::string str;
   if( get_string( node, "matrixK", str ) == 0 )
   {
      K0_fname = str;
   }
   else
   {
      linsol_error( "Driver", "wrong XML tags - matrixK attribute is mandatory" );
      MPI_Finalize();
      return 1;
   }
   if( get_string( node, "matrixC", str ) == 0 )
   {
      C0_fname = str;
   }
   else
   {
      linsol_error( "Driver", "wrong XML tags - matrixC attribute is mandatory" );
      MPI_Finalize();
      return 1;
   }
   if( get_string( node, "matrixCt", str ) == 0 )
   {
      Ct0_fname = str;
   }
   else
   {
      linsol_error( "Driver", "wrong XML tags - matrixCt attribute is mandatory" );
      MPI_Finalize();
      return 1;
   }
   if( get_string( node, "testSpace", str ) == 0 )
   {
      TS0_fname = str;
   }
   else
   {
      errorMsg( "Warning: Initial test space file not provided" );
   }
   if( get_string( node, "rhs", str ) == 0 )
   {
      TS0_fname = str;
   }
   else
   {
      errorMsg( "Warning: rhs file not provided" );
   }

   return 0;
}
