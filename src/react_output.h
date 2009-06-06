#ifndef REACT_OUTPUT_H
#define REACT_OUTPUT_H

#include "mcell_structs.h"

/* Header file for reaction output routines */

extern int emergency_output_hook_enabled;
void install_emergency_output_hooks(void);

int truncate_output_file(char *name,double start_value);
void add_trigger_output(struct counter *c,struct output_request *ear,int n,short flags);
int flush_reaction_output(void);
int update_reaction_output(struct output_block *block);
int write_reaction_output(struct output_set *set,int final_chunk_flag);

struct output_expression* new_output_expr(struct mem_helper *oexpr_mem);
void set_oexpr_column(struct output_expression *oe,struct output_column *oc);
void learn_oexpr_flags(struct output_expression *oe);
struct output_expression* dupl_oexpr_tree(struct output_expression *root,struct mem_helper *oexpr_mem);
struct output_expression* first_oexpr_tree(struct output_expression *root);
struct output_expression* last_oexpr_tree(struct output_expression *root);
struct output_expression* next_oexpr_tree(struct output_expression *leaf);
struct output_expression* prev_oexpr_tree(struct output_expression *leaf);
void eval_oexpr_tree(struct output_expression *root,int skip_const);
void oexpr_flood_convert(struct output_expression *root,char old_oper,char new_oper);
char* oexpr_title(struct output_expression *root);

#endif
