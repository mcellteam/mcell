#ifndef VOLUME_OUTPUT_H
#define VOLUME_OUTPUT_H

#include "mcell_structs.h"

int update_volume_output(struct volume *wrld, struct volume_output_item *vo);
int output_volume_output_item(struct volume *wrld,
                              char const *filename,
                              struct volume_output_item *vo);

#endif
