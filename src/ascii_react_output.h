#ifndef ASCII_REACT_OUTPUT_H
#define ASCII_REACT_OUTPUT_H

#include "mcell_structs.h"


int init_ascii_reaction_data(struct output_block *block_data,
                             struct volume* world);
int check_reaction_output_file(struct output_set *os);
int open_reaction_output_file(struct output_set *os, 
                              struct output_block *block_data,
                              struct volume *world);

int write_reaction_output(struct output_set *set,
                          struct volume* world);
#endif
