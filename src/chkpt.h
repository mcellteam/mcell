#ifndef CHKPT_H
#define CHKPT_H

/* header file for chkpt.c, MCell checkpointing functions */

int create_chkpt(char const *filename);
int write_chkpt(FILE *fs);
int read_chkpt(FILE *fs);
#ifndef _WIN32
void chkpt_signal_handler(int signo);
#endif

#endif
