#ifndef SYM_TABLE_H
#define SYM_TABLE_H

#include "mcell_structs.h"

typedef unsigned long int ub4;   /* unsigned 4-byte quantities */
typedef unsigned char ub1;        /* unsigned 1-byte quantities */

ub4 jerkins_hash(ub1 *sym);
unsigned long hash(char *sym);
struct sym_table *retrieve_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab);
struct sym_table *store_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab);
struct sym_table **init_symtab(int size);

struct counter_hash_table **init_countertab(int size);
struct counter_hash_table *retrieve_counter(char *counter,
  struct reg_counter_ref *rcrp, struct counter_hash_table **countertab);
struct counter_hash_table *store_counter(char *counter,
  struct reg_counter_ref_list *rcrlp, struct counter_hash_table **countertab);

#endif
