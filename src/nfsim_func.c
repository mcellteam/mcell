#include "nfsim_func.h"

double get_standard_diffusion(struct abstract_molecule* self){
    return self->properties->D;
}


double get_nfsim_diffusion(struct abstract_molecule* self){
    return 0.0;
}
