#ifndef REACT_OUTPUT_H
#define REACT_OUTPUT_H

/* Header file for reaction output routines */

void update_reaction_output(struct output_list *olp);
void update_counter(struct counter_list *clp);
void eval_count_expr_tree(struct counter_list *clp);
int eval_count_expr(struct counter_list *operand1,
                    struct counter_list *operand2,
                    char oper,
                    struct counter_list *result);
double eval_double(double op1, double op2, char oper);


#endif
