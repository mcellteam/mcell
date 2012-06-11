#ifndef CHKPT_H
#define CHKPT_H

/* header file for chkpt.c, MCell checkpointing functions */

#ifdef MCELL_WITH_CHECKPOINTING /* only declare the functions if checkpointing is enabled */

int create_chkpt(char const *filename);
int write_chkpt(FILE *fs);
int read_chkpt(FILE *fs);
void chkpt_signal_handler(int signo);

#else

#warning included chkpt.h but checkpointing is not enabled

#endif

#endif
