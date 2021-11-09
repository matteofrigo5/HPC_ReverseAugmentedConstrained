#include "get_solver_kind.h"
#include "utilities.hpp"

iReg get_solver_kind( pugi::xml_document const & input_xml )
{
  pugi::xml_node const node = input_xml.child( "Chronos" ).child( "Solver" );

  warning warnLvl0( "Solver" );

  std::string str;
  iReg solverID = 1;
  if( get_string( node, "kind", str ) == 0 )
  {
    if( stringCompare( str, "BiCGStab" ) )
    {
      solverID = 1;
    }
    else if( stringCompare( str, "GMRES" ) )
    {
      solverID = 2;
    }
  }
  else
  {
    warning( "kind" );
  }
  return solverID;
}
