#include "set_prec_parms.h"
#include "utilities.hpp"

void set_prec_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                    RevAugm &PREC_RA){

   bool USE_aFSAI = false;
   bool USE_aAMG_mech = false;
   bool USE_aAMG_cfd = false;

   pugi::xml_node const node = input_xml.child( "Chronos" ).child( "Preconditioner" );

   warning warnLvl0( "Preconditioner" );

   double rval;
   std::string str;

   get_dbl_value( node, "gamma", rval ) == 0 ? PREC_RA.set_gamma( rval ) : warnLvl0.print( "gamma" );
   if( get_string( node, "precSchur", str ) == 0 )
   {
      if( stringCompare( str, "aFSAI" ) )
      {
         USE_aFSAI = true;
      }
      else if( stringCompare( str, "AMG_mech" ) )
      {
         USE_aAMG_mech = true;
      }
      else if( stringCompare( str, "AMG_cfd" ) )
      {
         USE_aAMG_cfd = true;
      }
   }
   else
   {
      warning( "precSchur" );
   }

   if( node.child( "aFSAI" ) != NULL )
   {
      pugi::xml_node const node_prec = node.child( "aFSAI" );

      warning warnLvl1( warnLvl0 );
      warnLvl1.push( "aFSAI" );

      int ival;
      double rval;
      get_int_value( node_prec, "targetPower", ival ) == 0 ? PREC_RA.Set_TargetPower( ival ) : warnLvl1.print( "targetPower" );
      get_int_value( node_prec, "numStep", ival ) == 0     ? PREC_RA.Set_nstep( ival )       : warnLvl1.print( "numStep" );
      get_int_value( node_prec, "stepSize", ival ) == 0    ? PREC_RA.Set_step_size( ival )   : warnLvl1.print( "stepSize" );
      get_dbl_value( node_prec, "eps", rval ) == 0         ? PREC_RA.Set_eps( rval )         : warnLvl1.print( "eps" );
   }
   else
   {
      warnLvl0.print( "aFSAI" );
   }

   if( node.child( "AMG_mech" ) != NULL )
   {
      pugi::xml_node const node_prec = node.child( "AMG_mech" );

      warning warnLvl1( warnLvl0 );
      warnLvl1.push( "AMG_mech" );

      int ival;
      double rval;
      pugi::xml_node local_node;

      get_int_value( node_prec, "cycleType", ival ) == 0     ? PREC_RA.Set_CycleType( ival )     : warnLvl1.print( "cycleType" );
      get_int_value( node_prec, "maxNumLevels", ival ) == 0  ? PREC_RA.Set_maxNlevel( ival )     : warnLvl1.print( "maxNumLevels" );
      get_int_value( node_prec, "maxCoarseSize", ival ) == 0 ? PREC_RA.Set_maxCoarseSize( ival ) : warnLvl1.print( "maxCoarseSize" );

      local_node = node_prec.child( "TestSpace" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "TestSpace" );

        get_int_value( local_node, "ntv", ival ) == 0      ? PREC_RA.TSPACE_Set_ntv( ival )       : warnLvl2.print( "ntv" );
        get_int_value( local_node, "method", ival ) == 0   ? PREC_RA.TSPACE_Set_method( ival )    : warnLvl2.print( "method" );
        get_int_value( local_node, "maxIter", ival ) == 0  ? PREC_RA.TSPACE_Set_maxITER( ival )   : warnLvl2.print( "maxIter" );
        get_dbl_value( local_node, "exitTol", rval ) == 0  ? PREC_RA.TSPACE_Set_exitTOL( rval )   : warnLvl2.print( "exitTol" );
        get_int_value( local_node, "orthFreq", ival ) == 0 ? PREC_RA.TSPACE_Set_orth_freq( ival ) : warnLvl2.print( "orthFreq" );
        get_dbl_value( local_node, "ritzFreq", rval ) == 0 ? PREC_RA.TSPACE_Set_ritz_freq( rval ) : warnLvl2.print( "ritzFreq" );
      }
      else
      {
        warnLvl1.print( "TestSpace" );
      }

      local_node = node_prec.child( "Smoother" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Smoother" );

        get_int_value( local_node, "type", ival ) == 0             ? PREC_RA.SMOOTHER_Set_type( ival )              : warnLvl2.print( "type" );
        get_dbl_value( local_node, "omega", rval ) == 0            ? PREC_RA.SMOOTHER_Set_omega( rval )             : warnLvl2.print( "omega" );
        get_int_value( local_node, "maxIterMaxEig", ival ) == 0    ? PREC_RA.SMOOTHER_Set_itmax_maxeig( ival )      : warnLvl2.print( "maxIterMaxEig" );
        get_dbl_value( local_node, "tolMaxEig", rval ) == 0        ? PREC_RA.SMOOTHER_Set_tol_maxeig( rval )        : warnLvl2.print( "tolMaxEig" );
        get_int_value( local_node, "nuPre", ival ) == 0            ? PREC_RA.SMOOTHER_Set_nupre( ival )             : warnLvl2.print( "nuPre" );
        get_int_value( local_node, "nuPost", ival ) == 0           ? PREC_RA.SMOOTHER_Set_nupost( ival )            : warnLvl2.print( "nuPost" );
        get_int_value( local_node, "aFSAItargetPower", ival ) == 0 ? PREC_RA.SMOOTHER_Set_aFSAI_TargetPower( ival ) : warnLvl2.print( "aFSAItargetPower" );
        get_int_value( local_node, "aFSAInumStep", ival ) == 0     ? PREC_RA.SMOOTHER_Set_aFSAI_nstep( ival )       : warnLvl2.print( "aFSAInumStep" );
        get_int_value( local_node, "aFSAIstepSize", ival ) == 0    ? PREC_RA.SMOOTHER_Set_aFSAI_step_size( ival )   : warnLvl2.print( "aFSAIstepSize" );
        get_dbl_value( local_node, "aFSAIeps", rval ) == 0         ? PREC_RA.SMOOTHER_Set_aFSAI_eps( rval )         : warnLvl2.print( "aFSAIeps" );
      }
      else
      {
        warnLvl1.print( "Smoother" );
      }

      local_node = node_prec.child( "Coarsen" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Smoother" );

        get_int_value( local_node, "type", ival ) == 0    ? PREC_RA.COARSEN_Set_type( ival )     : warnLvl2.print( "type" );
        get_dbl_value( local_node, "tau", rval ) == 0     ? PREC_RA.COARSEN_Set_tau( rval )      : warnLvl2.print( "tau" );
        get_dbl_value( local_node, "tauJump", rval ) == 0 ? PREC_RA.COARSEN_Set_tau_jump( rval ) : warnLvl2.print( "tauJump" );
      }
      else
      {
        warnLvl1.print( "Coarsen" );
      }

      local_node = node_prec.child( "Prolongation" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Prolongation" );

        get_int_value( local_node, "type", ival ) == 0       ? PREC_RA.PROLONGATION_Set_type( ival )      : warnLvl2.print( "type" );
        get_int_value( local_node, "power", ival ) == 0      ? PREC_RA.PROLONGATION_Set_nPow( ival )      : warnLvl2.print( "power" );
        get_int_value( local_node, "smooth", ival ) == 0     ? PREC_RA.PROLONGATION_Set_smooth( ival )    : warnLvl2.print( "smooth" );
        get_dbl_value( local_node, "maxRowNorm", rval ) == 0 ? PREC_RA.PROLONGATION_Set_maxrownrm( rval ) : warnLvl2.print( "maxRowNorm" );
        get_dbl_value( local_node, "maxCond", rval ) == 0    ? PREC_RA.PROLONGATION_Set_maxcond( rval )   : warnLvl2.print( "maxCond" );
        get_int_value( local_node, "maxIterVol", ival ) == 0 ? PREC_RA.PROLONGATION_Set_itmax_vol( ival ) : warnLvl2.print( "maxIterVol" );
        get_dbl_value( local_node, "tolVol", rval ) == 0     ? PREC_RA.PROLONGATION_Set_tol_vol( rval )   : warnLvl2.print( "tolVol" );
        get_dbl_value( local_node, "eps", rval ) == 0        ? PREC_RA.PROLONGATION_Set_eps( rval )       : warnLvl2.print( "eps" );
        get_int_value( local_node, "minDist", ival ) == 0    ? PREC_RA.PROLONGATION_Set_dist_min( ival )  : warnLvl2.print( "minDist" );
        get_int_value( local_node, "maxDist", ival ) == 0    ? PREC_RA.PROLONGATION_Set_dist_max( ival )  : warnLvl2.print( "maxDist" );
        get_int_value( local_node, "mMax", ival ) == 0       ? PREC_RA.PROLONGATION_Set_mmax( ival )      : warnLvl2.print( "mMax" );
      }
      else
      {
        warnLvl1.print( "Prolongation" );
      }

      local_node = node_prec.child( "Filter" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Filter" );

        get_dbl_value( local_node, "prolWeight", rval ) == 0 ? PREC_RA.FILTER_Set_Prol_Weight( rval ) : warnLvl2.print( "prolWeight" );
        get_dbl_value( local_node, "prolTol", rval ) == 0    ? PREC_RA.FILTER_Set_Prol_Tol( rval )    : warnLvl2.print( "prolTol" );
        get_dbl_value( local_node, "operWeight", rval ) == 0 ? PREC_RA.FILTER_Set_Oper_Weight( rval ) : warnLvl2.print( "operWeight" );
      }
      else
      {
        warnLvl1.print( "Filter" );
      }
   }
   else
   {
      warnLvl0.print( "AMG_mech" );
   }

   if( node.child( "AMG_cfd" ) != NULL )
   {
      pugi::xml_node const node_prec = node.child( "AMG_cfd" );

      warning warnLvl1( warnLvl0 );
      warnLvl1.push( "AMG_cfd" );

      int ival;
      double rval;
      pugi::xml_node local_node;

      get_int_value( node_prec, "cycleType", ival ) == 0     ? PREC_RA.Set_CycleType( ival )     : warnLvl1.print( "cycleType" );
      get_int_value( node_prec, "maxNumLevels", ival ) == 0  ? PREC_RA.Set_maxNlevel( ival )     : warnLvl1.print( "maxNumLevels" );
      get_int_value( node_prec, "maxCoarseSize", ival ) == 0 ? PREC_RA.Set_maxCoarseSize( ival ) : warnLvl1.print( "maxCoarseSize" );

      local_node = node_prec.child( "TestSpace" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "TestSpace" );

        get_int_value( local_node, "ntv", ival ) == 0      ? PREC_RA.TSPACE_Set_ntv( ival )       : warnLvl2.print( "ntv" );
        get_int_value( local_node, "method", ival ) == 0   ? PREC_RA.TSPACE_Set_method( ival )    : warnLvl2.print( "method" );
        get_int_value( local_node, "maxIter", ival ) == 0  ? PREC_RA.TSPACE_Set_maxITER( ival )   : warnLvl2.print( "maxIter" );
        get_dbl_value( local_node, "exitTol", rval ) == 0  ? PREC_RA.TSPACE_Set_exitTOL( rval )   : warnLvl2.print( "exitTol" );
        get_int_value( local_node, "orthFreq", ival ) == 0 ? PREC_RA.TSPACE_Set_orth_freq( ival ) : warnLvl2.print( "orthFreq" );
        get_dbl_value( local_node, "ritzFreq", rval ) == 0 ? PREC_RA.TSPACE_Set_ritz_freq( rval ) : warnLvl2.print( "ritzFreq" );
      }
      else
      {
        warnLvl1.print( "TestSpace" );
      }

      local_node = node_prec.child( "Smoother" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Smoother" );

        get_int_value( local_node, "type", ival ) == 0             ? PREC_RA.SMOOTHER_Set_type( ival )              : warnLvl2.print( "type" );
        get_dbl_value( local_node, "omega", rval ) == 0            ? PREC_RA.SMOOTHER_Set_omega( rval )             : warnLvl2.print( "omega" );
        get_int_value( local_node, "maxIterMaxEig", ival ) == 0    ? PREC_RA.SMOOTHER_Set_itmax_maxeig( ival )      : warnLvl2.print( "maxIterMaxEig" );
        get_dbl_value( local_node, "tolMaxEig", rval ) == 0        ? PREC_RA.SMOOTHER_Set_tol_maxeig( rval )        : warnLvl2.print( "tolMaxEig" );
        get_int_value( local_node, "nuPre", ival ) == 0            ? PREC_RA.SMOOTHER_Set_nupre( ival )             : warnLvl2.print( "nuPre" );
        get_int_value( local_node, "nuPost", ival ) == 0           ? PREC_RA.SMOOTHER_Set_nupost( ival )            : warnLvl2.print( "nuPost" );
        get_int_value( local_node, "aFSAItargetPower", ival ) == 0 ? PREC_RA.SMOOTHER_Set_aFSAI_TargetPower( ival ) : warnLvl2.print( "aFSAItargetPower" );
        get_int_value( local_node, "aFSAInumStep", ival ) == 0     ? PREC_RA.SMOOTHER_Set_aFSAI_nstep( ival )       : warnLvl2.print( "aFSAInumStep" );
        get_int_value( local_node, "aFSAIstepSize", ival ) == 0    ? PREC_RA.SMOOTHER_Set_aFSAI_step_size( ival )   : warnLvl2.print( "aFSAIstepSize" );
        get_dbl_value( local_node, "aFSAIeps", rval ) == 0         ? PREC_RA.SMOOTHER_Set_aFSAI_eps( rval )         : warnLvl2.print( "aFSAIeps" );
      }
      else
      {
        warnLvl1.print( "Smoother" );
      }

      local_node = node_prec.child( "Coarsen" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Smoother" );

        get_int_value( local_node, "type", ival ) == 0    ? PREC_RA.COARSEN_Set_type( ival )     : warnLvl2.print( "type" );
        get_dbl_value( local_node, "tau", rval ) == 0     ? PREC_RA.COARSEN_Set_tau( rval )      : warnLvl2.print( "tau" );
        get_dbl_value( local_node, "tauJump", rval ) == 0 ? PREC_RA.COARSEN_Set_tau_jump( rval ) : warnLvl2.print( "tauJump" );
      }
      else
      {
        warnLvl1.print( "Coarsen" );
      }

      local_node = node_prec.child( "Prolongation" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Prolongation" );

        get_int_value( local_node, "type", ival ) == 0       ? PREC_RA.PROLONGATION_Set_type( ival )      : warnLvl2.print( "type" );
        get_int_value( local_node, "power", ival ) == 0      ? PREC_RA.PROLONGATION_Set_nPow( ival )      : warnLvl2.print( "power" );
        get_int_value( local_node, "smooth", ival ) == 0     ? PREC_RA.PROLONGATION_Set_smooth( ival )    : warnLvl2.print( "smooth" );
        get_dbl_value( local_node, "maxRowNorm", rval ) == 0 ? PREC_RA.PROLONGATION_Set_maxrownrm( rval ) : warnLvl2.print( "maxRowNorm" );
        get_dbl_value( local_node, "maxCond", rval ) == 0    ? PREC_RA.PROLONGATION_Set_maxcond( rval )   : warnLvl2.print( "maxCond" );
        get_int_value( local_node, "maxIterVol", ival ) == 0 ? PREC_RA.PROLONGATION_Set_itmax_vol( ival ) : warnLvl2.print( "maxIterVol" );
        get_dbl_value( local_node, "tolVol", rval ) == 0     ? PREC_RA.PROLONGATION_Set_tol_vol( rval )   : warnLvl2.print( "tolVol" );
        get_dbl_value( local_node, "eps", rval ) == 0        ? PREC_RA.PROLONGATION_Set_eps( rval )       : warnLvl2.print( "eps" );
        get_int_value( local_node, "minDist", ival ) == 0    ? PREC_RA.PROLONGATION_Set_dist_min( ival )  : warnLvl2.print( "minDist" );
        get_int_value( local_node, "maxDist", ival ) == 0    ? PREC_RA.PROLONGATION_Set_dist_max( ival )  : warnLvl2.print( "maxDist" );
        get_int_value( local_node, "mMax", ival ) == 0       ? PREC_RA.PROLONGATION_Set_mmax( ival )      : warnLvl2.print( "mMax" );
      }
      else
      {
        warnLvl1.print( "Prolongation" );
      }

      local_node = node_prec.child( "Filter" );
      if( local_node != NULL )
      {
        warning warnLvl2( warnLvl1 );
        warnLvl2.push( "Filter" );

        get_dbl_value( local_node, "prolWeight", rval ) == 0 ? PREC_RA.FILTER_Set_Prol_Weight( rval ) : warnLvl2.print( "prolWeight" );
        get_dbl_value( local_node, "prolTol", rval ) == 0    ? PREC_RA.FILTER_Set_Prol_Tol( rval )    : warnLvl2.print( "prolTol" );
        get_dbl_value( local_node, "operWeight", rval ) == 0 ? PREC_RA.FILTER_Set_Oper_Weight( rval ) : warnLvl2.print( "operWeight" );
      }
      else
      {
        warnLvl1.print( "Filter" );
      }
   }
   else
   {
      warnLvl0.print( "AMG_cfg" );
   }

   if( USE_aFSAI )
   {
      try
      {
         PREC_RA.PrecK_set_type(0);
      }
      catch( linsol_error )
      {
         throw linsol_error( "Set_prec_parms", "errorMsg in PrecK_set_type" );
      }
   }

   if( USE_aAMG_mech )
   {
      try
      {
         PREC_RA.PrecK_set_type(1);
      }
      catch( linsol_error )
      {
         throw linsol_error( "Set_prec_parms", "errorMsg in PrecK_set_type" );
      }
   }

   if( USE_aAMG_cfd )
   {
      try
      {
         PREC_RA.PrecK_set_type(2);
      }
      catch( linsol_error )
      {
         throw linsol_error( "Set_prec_parms", "errorMsg in PrecK_set_type" );
      }
   }
}
