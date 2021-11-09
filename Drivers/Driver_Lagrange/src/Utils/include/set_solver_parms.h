#include "precision.h"
#include "linsol_error.h"
#include "LinSolver.h"

// include pugixml for XML parsing
#include "pugixml.hpp"

iReg set_solver_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                      LinSolver &SOLVER, iReg &RHS_Flag );
