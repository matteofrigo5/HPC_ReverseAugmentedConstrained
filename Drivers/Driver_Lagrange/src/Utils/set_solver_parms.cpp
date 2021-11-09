#include "set_solver_parms.h"
#include "utilities.hpp"

iReg set_solver_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                      LinSolver &SOLVER, iReg &RHS_Flag ){

   iReg maxiter;
   rExt tolerance;
   iReg init_sol;

   // Set solver default parameters
   RHS_Flag = 0;

   pugi::xml_node const node = input_xml.child( "Chronos" ).child( "Solver" );

   warning warnLvl0( "Solver" );
   if( node.child( "BiCGStab" ) != NULL )
   {
      pugi::xml_node const node_solv = node.child( "BiCGStab" );

      warning warnLvl1( warnLvl0 );
      warnLvl1.push( "BiCGStab" );

      int ival;
      double rval;
      get_int_value( node_solv, "maxIter", ival ) == 0       ? SOLVER.Set_maxITER( ival )          : warnLvl1.print( "maxIter" ); maxiter = 1000;//maxiter = ival;
      get_dbl_value( node_solv, "exitTol", rval ) == 0       ? SOLVER.Set_exitTOL( rval )          : warnLvl1.print( "exitTol" ); tolerance =1.e-6;//tolerance = rval;
      get_int_value( node_solv, "initSol", ival ) == 0       ? SOLVER.Set_initial_sol( ival )      : warnLvl1.print( "initSol" ); init_sol = 0; //init_sol = ival;
      get_int_value( node_solv, "printConvProf", ival ) == 0 ? SOLVER.Set_printConvProfile( ival ) : warnLvl1.print( "printConvProf" );
      if( get_int_value( node_solv, "rhsFlag", ival ) == 0 )
      {
        RHS_Flag = iReg( ival );
      }
      else
      {
        warnLvl1.print( "rhsFlag" );
      }
   }
   else
   {
     warnLvl0.print( "BiCGStab" );
   }

   if( node.child( "GMRES" ) != NULL )
   {
      pugi::xml_node const node_solv = node.child( "GMRES" );

      warning warnLvl1( warnLvl0 );
      warnLvl1.push( "GMRES" );

      int ival;
      double rval;
      get_int_value( node_solv, "maxIter", ival ) == 0       ? SOLVER.Set_maxITER( ival )          : warnLvl1.print( "maxIter" ); maxiter = ival;
      get_dbl_value( node_solv, "exitTol", rval ) == 0       ? SOLVER.Set_exitTOL( rval )          : warnLvl1.print( "exitTol" ); tolerance = rval;
      get_int_value( node_solv, "initSol", ival ) == 0       ? SOLVER.Set_initial_sol( ival )      : warnLvl1.print( "initSol" ); init_sol = ival;
      get_int_value( node_solv, "restart", ival ) == 0       ? SOLVER.Set_restart( ival )          : warnLvl1.print( "restart" );
      get_int_value( node_solv, "printConvProf", ival ) == 0 ? SOLVER.Set_printConvProfile( ival ) : warnLvl1.print( "printConvProf" );
      if( get_int_value( node_solv, "rhsFlag", ival ) == 0 )
      {
        RHS_Flag = iReg( ival );
      }
      else
      {
        warnLvl1.print( "rhsFlag" );
      }
   }
   else
   {
     warnLvl0.print( "GMRES" );
   }

   // Check correctness of parameters before set
   bool Error = false;
   if( maxiter <= 0 )
   {
      errorMsg( "Error: wrong maxIter value" );
      Error = true;
   }
   if( tolerance <= 0.0 )
   {
      errorMsg( "Error: wrong exitTol value" );
      Error = true;
   }
   if( init_sol < 0 || init_sol > 2 )
   {
      errorMsg( "Error: wrong initSol value" );
      Error = true;
   }
   if( Error )
   {
      linsol_error( "Driver", "wrong solver parameters in input" );
      MPI_Finalize();
      return 1;
   }

   return 0;
}
