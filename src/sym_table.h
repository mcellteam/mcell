#ifndef SYM_TABLE_H
#define SYM_TABLE_H

#include "mcell_structs.h"

#ifndef ISAAC64_H
/* These guys come in for free if we're using Jerkins' random numbers */
typedef uint32_t ub4;   /* unsigned 4-byte quantities */
typedef unsigned char ub1;        /* unsigned 1-byte quantities */
#endif

ub4 jerkins_hash(ub1 *sym);
unsigned long hash(char *sym);
struct sym_table *retrieve_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab);
struct sym_table *store_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab);
struct sym_table **init_symtab(int size);

#endif
