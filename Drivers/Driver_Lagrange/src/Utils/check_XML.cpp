#include "check_XML.h"

iReg check_XML(const pugi::xml_document &input_xml){

   // Check if the XML is consisistent and the mandatory parameters are provided
   pugi::xml_node node;
   pugi::xml_attribute attr;
   node = input_xml.child("Chronos").child("Preconditioner");
   if (node == NULL){
      linsol_error("Driver","wrong XML tags - Preconditioner tag is mandatory");
      MPI_Finalize();
      return 1;
   }
   if (distance(node.begin(), node.end()) < 1) {
      linsol_error("Driver","wrong XML tags - Preconditioner tag must have at least one child");
      MPI_Finalize();
      return 1;
   }

   node = input_xml.child("Chronos").child("Solver");
   if (node == NULL){
      linsol_error("Driver","wrong XML tags - Solver tag is mandatory");
      MPI_Finalize();
      return 1;
   }

   // Read mandatory fnames
   attr = input_xml.child("Chronos").child("Solver").attribute("rhsFlag");
   if (attr != NULL){
      if (attr.as_int()==3){
         pugi::xml_attribute attr_names = input_xml.child("Chronos").child("fnames").attribute("rhs");
         if (attr_names == NULL){
            linsol_error("Driver","wrong XML tags - rhs attribute is mandatory if rhs_Flag is equal to 3");
            MPI_Finalize();
            return 1;
         }
      }
   }

   return 0;

}
