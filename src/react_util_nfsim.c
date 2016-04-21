#include "mcell_reactions.h"
#include "react.h"
#include <stdlib.h>
#include <string.h>

void calculate_reactant_orientation(struct abstract_molecule* reac, struct abstract_molecule* reac2, bool* orientation_flag1, bool* orientation_flag2, int* reactantOrientation1, int* reactantOrientation2){
  //calculate orientation information
  if(reac2 !=NULL){
    const char* compartment1 = extractSpeciesCompartmentFromNauty_c(reac->graph_data->graph_pattern);
    const char* compartment2 = extractSpeciesCompartmentFromNauty_c(reac2->graph_data->graph_pattern);
    compartmentStruct reactantCompartmentInfo1;
    compartmentStruct reactantCompartmentInfo2;

    // what compartment is the species now in
    reactantCompartmentInfo1 = getCompartmentInformation_c(compartment1);
    reactantCompartmentInfo2 = getCompartmentInformation_c(compartment2);

    //
    //reac is volume, reac2 is surface
    if((reac->flags & ON_GRID) == 0 &&  (reac2->flags & ON_GRID) != 0){
      *orientation_flag1 = 1;
      *orientation_flag2 = 1;
      *reactantOrientation2 = 1;
      if (strcmp(reactantCompartmentInfo1.outside,reactantCompartmentInfo2.name) == 0){
        
        *reactantOrientation1 = -1;
      }
      else{
        *reactantOrientation1 = 1;
      }
    }
    //reac2 is volume, reac is surface
    else if((reac2->flags & ON_GRID) == 0 &&  (reac->flags & ON_GRID) != 0){
      *orientation_flag1 = 1;
      *orientation_flag2 = 1;
      *reactantOrientation1 = 1;
      if (strcmp(reactantCompartmentInfo2.outside,reactantCompartmentInfo1.name) == 0){
        *reactantOrientation2 = -1;
      }
      else{
        *reactantOrientation2 = 1;
      }
    }
    //both surface molecules
    else if((reac2->flags & ON_GRID) != 0 &&  (reac->flags & ON_GRID) != 0){
      *orientation_flag1 = 1;
      *orientation_flag2 = 1;
      *reactantOrientation1 = 1;
      *reactantOrientation2 = 1;
    }
    //both volume molecules
    else{
      *orientation_flag1 =0;
      *orientation_flag2 =0;
      *reactantOrientation1 = 0;
      *reactantOrientation2 = 0;
    }

    free((char*)compartment1);
    free((char*)compartment2);
  }
  else if((reac->flags & ON_GRID) != 0){
    *orientation_flag1 = 1;
    *reactantOrientation1 = 1;
    *orientation_flag2 = 1;

  }
}