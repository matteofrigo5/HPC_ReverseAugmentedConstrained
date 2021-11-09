#include "precision.h"
#include "linsol_error.h"

// include pugixml for XML parsing
#include "pugixml.hpp"

iReg read_general_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                        type_OMP_iReg &nthreads, iReg &verbosity, bool &PART);
