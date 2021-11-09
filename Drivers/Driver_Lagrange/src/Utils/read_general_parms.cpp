#include "read_general_parms.h"
#include "utilities.hpp"

iReg read_general_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                        type_OMP_iReg &nthreads, iReg &verbosity, bool &PART){

   // Set default value
   nthreads = 1;
   verbosity = 0;
   PART = false;

   pugi::xml_node const node = input_xml.child( "Chronos" ).child( "GeneralParameters" );

   warning warn( "GeneralParameters" );

   if( get_int_value( node, "numThreads", nthreads ) != 0 )
   {
      warn.print( "numThreads" );
   }
   if( get_int_value( node, "verbosity", verbosity ) != 0 )
   {
      warn.print( "verbosity" );
   }
   if( get_lgl_value( node, "metisPart", PART ) != 0 )
   {
      warn.print( "metisPart" );
   }

   // Check correctness of parameters before exit
   bool Error = false;
   if( nthreads < 0 )
   {
      errorMsg( "Error: wrong nthreads value" );
      Error = true;
   }
   if( verbosity < 0 || verbosity > 3 )
   {
      errorMsg( "Error: wrong verbosity value" );
      Error = true;
   }
   if( Error )
   {
      linsol_error( "Driver", "wrong general parameters in input" );
      MPI_Finalize();
      return 1;
   }
   return 0;
}
