#ifndef MDLPARSE_UTIL_H
#define MDLPARSE_UTIL_H

#include "vector.h"
#include "mcell_structs.h"
#include "mdlparse.h"

void mdl_warning(struct mdlparse_vars *mpvp);
void swap_double(double *x, double *y);
void swap_int(int *x, int *y);
double *double_dup(double value);
struct name_list *concat_obj_name(struct name_list *name_list_end,char *name);
char *get_first_name(char *obj_name);
char *get_prefix_name(char *obj_name);
struct object *find_full_name(struct object *objp,char *full_name,
			      char *sub_name);
struct object *make_new_object(struct volume *volp,char *obj_name,
			       char *err_msg);
struct region *make_new_region(struct volume *volp,char *obj_name,
			       char *region_last_name,char *err_msg);
int copy_object(struct volume *volp,struct object *curr_objp,
		struct object *objp,struct object *objp2, char *err_msg);
char *concat_rx_name(char *name1, char *name2);

int equivalent_geometry(struct pathway *p1,struct pathway *p2,int n);
int prepare_reactions(struct mdlparse_vars *mpvp);

int make_cuboid(struct vector3 *p1,
                struct vector3 *p2,
                struct ordered_poly *opp);

int set_viz_state_value(struct object *objp, int viz_state);

void sort_num_expr_list(struct num_expr_list *head);

int my_sprintf(char *strp,char *format,struct arg *argp,u_int num_args);
int my_fprintf(FILE *outfile,char *format,struct arg *argp,u_int num_args);

struct counter *retrieve_reg_counter(struct volume *volp,
                                     void *vp,
                                     struct region *rp,
                                     byte counter_type);

struct counter *store_reg_counter(struct volume *volp,
                                  void *vp,
                                  struct region *rp,
                                  byte counter_type);

struct output_evaluator *init_counter(byte report_type,
                                      struct output_item *oip,
                                      struct counter *cp,
                                      u_int buffersize);
int insert_counter(byte report_type,
                       struct volume *volp,
                       struct output_item *oip,
                       struct output_evaluator *oep,
                       struct counter *cp,
                       u_int buffersize);
int build_count_tree(byte report_type,
                         struct volume *volp,
                         struct object *objp,
                         struct output_item *oip,
                         struct output_evaluator *oep,
                         void *vp,
                         u_int buffersize,
                         char *sub_name);

int handle_count_request(unsigned short sym_type,
                         void *value,
			 struct region *r,
			 struct object *obj,
			 byte report_type,
			 struct mdlparse_vars *mdlpvp);

#endif
