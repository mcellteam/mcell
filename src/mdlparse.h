#ifndef MDLPARSE_H
#define MDLPARSE_H

#include "mcell_structs.h"

struct mdlparse_vars {
  char *cval;
  char *cval_2;
  char *strval;
  int ival;
  double rval;

  struct sym_table *gp;
  struct sym_table *tp;
  char *sym_name;

  char *a_str;
  char str_buf[1024];
  char str_buf2[1024];
  char format_str[256];

  double val_1;
  double val_2;
  double tmp_dbl;
  double *dblp;

  struct species *specp;
  double mc_factor;
  double transport_mc_factor;
  double l_perp_bar;
  double l_perp_rms;
  double l_r_bar;
  double l_r_rms;

  int num_pos;
  struct num_expr_list *elp;
  struct num_expr_list *el_head;
  struct num_expr_list *el_tail;

  struct object *curr_obj;
  struct object *objp;
  struct object *objp2;
  struct object *top_objp;
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
  struct name_list *instance_name_list;
  struct name_list *instance_name_list_end;
  struct vector3 *pntp1,*pntp2;
  double tm[4][4];
  char *obj_name;
  char *prefix_name;
  char *full_name;

  struct release_site_obj *rsop;
  struct release_pattern *rpatp;
  struct frame_data_list *fdlp;
  struct viz_obj *vizp;
  int viz_state;
  int existing_state;

  FILE *file;

  char mdl_err_msg[1024];
  u_int line_num[MAX_INCLUDE_DEPTH];
  struct yy_buffer_state *include_stack[MAX_INCLUDE_DEPTH];
  char *include_filename[MAX_INCLUDE_DEPTH];
  char *cval_stack[MAX_INCLUDE_DEPTH];
  char *cval_2_stack[MAX_INCLUDE_DEPTH];
  u_int include_stack_ptr;
  byte include_flag;
  struct volume *vol;
};

int mdlparse(void *p);
void mdlerror(char *s,...);
int mdlparse_init(struct volume *vol);

#endif
