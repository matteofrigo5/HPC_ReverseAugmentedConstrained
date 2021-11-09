#include "precision.h"
#include "linsol_error.h"
#include "RevAugm.h"

// include pugixml for XML parsing
#include "pugixml.hpp"

void set_prec_parms(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                    RevAugm &PREC_RA);
