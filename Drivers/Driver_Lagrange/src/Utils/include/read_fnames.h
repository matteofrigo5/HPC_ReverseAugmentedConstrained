#include "precision.h"
#include "linsol_error.h"

// include pugixml for XML parsing
#include "pugixml.hpp"

iReg read_fnames(const type_MPI_iReg rank, const pugi::xml_document &input_xml,
                 string &K0_fname, string &C0_fname, string &Ct0_fname, string &RHS0_fname, string &TS0_fname);
