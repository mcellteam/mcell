#ifndef REACT_OUTPUT_H
#define REACT_OUTPUT_H

/* Header file for reaction output routines */

int emergency_output();
int update_reaction_output(struct output_block *obp);
int write_reaction_output(struct output_block *obp,int final_chunk_flag);
int eval_count_expr_tree(struct output_evaluator *oep);
int eval_count_expr(struct output_evaluator *operand1,
                    struct output_evaluator *operand2,
                    char oper,
                    struct output_evaluator *result);
double eval_double(double op1, double op2, char oper);


#endif
