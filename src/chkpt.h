#ifndef CHKPT_H
#define CHKPT_H

/* header file for chkpt.c, MCell checkpointing functions */

int create_chkpt(struct volume *world, char const *filename);
int write_chkpt(struct volume *world, FILE *fs);
int read_chkpt(struct volume *world, FILE *fs);
void chkpt_signal_handler(struct volume *world, int signo);

#endif
