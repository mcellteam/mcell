#ifndef MDLPARSE_H
#define MDLPARSE_H

#include "mcell_structs.h"

#define RFLCT 0
#define TRANSP 1
#define SINK 2

struct mdlparse_vars {
  char *cval;
  char *cval_2;
  char *strval;
  int ival;
  double rval;

  struct sym_table *gp;
  struct sym_table *tp;
  struct sym_table *stp1;
  struct sym_table *stp2;
  struct sym_table *stp3;
  char *sym_name;

  char *a_str;
  char str_buf[1024];
  char str_buf2[1024];
  char temp_str[1024];
  char format_str[256];

  double val_1;
  double val_2;
  double tmp_dbl;
  double *dblp;
  int *intp;

  struct species *specp;
  double mc_factor;
  double transport_mc_factor;
  double l_perp_bar;
  double l_perp_rms;
  double l_r_bar;
  double l_r_rms;

  int num_pos;
  struct num_expr_list *elp;
  struct num_expr_list *elp_temp;
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
  struct vector3 *llf;
  struct vector3 *urb;
  struct vector3 *pntp1;
  struct vector3 *pntp2;
  double tm[4][4];
  char *obj_name;
  char *prefix_name;
  char full_name[1024];

  struct release_site_obj *rsop;
  struct release_pattern *rpatp;

  struct viz_obj *vizp;
  struct frame_data_list *fdlp;
  int viz_state;
  int existing_state;

  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct box_poly *bpp;
  struct element_data *edp;
  struct element_connection_list *connection_head;
  struct element_connection_list *connection_tail;
  struct element_connection_list *eclp;
  struct element_connection_list *eclp_temp;
  struct vertex_list *vertex_head;
  struct vertex_list *vertex_tail;
  struct vertex_list *vlp;
  struct vertex_list *vlp_temp;
  int n_walls;
  int n_verts;

  struct region *rp;
  struct region_list *region_list_head;
  struct region_list *rlp;
  struct element_list *element_list_head;
  struct element_list *elmlp;
  char *region_name;

  struct eff_dat *eff_dat_head;
  struct eff_dat *effdp;
  byte mol_quant_type;

  struct rxn *rxnp;
  struct mem_helper *path_mem;
  struct mem_helper *prod_mem;
  struct pathway *pathp;
  struct product *prodp;
  double fwd_km;
  double fwd_kcat;
  double back_km;
  double back_kcat;
  short orient_class;
  short orient_class1;
  short orient_class2;
  byte prod_all_3d;

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

struct element_connection_list {
  struct num_expr_list *connection_list;
  int n_verts;
  struct element_connection_list *next;
};

/**
 * A linked list used to store the coordinates of vertices and the   
 * corresponding normal vectors.
 */
struct vertex_list {
        struct vector3 *vertex;           /**< pointer to one polygon vertex */
        struct vector3 *normal;           /**< pointer to one polygon normal */
        struct vertex_list *next;         /**< pointer to next vertex_list */
};

int mdlparse(void *p);
void mdlerror(char *s,...);
int mdlparse_init(struct volume *vol);

#endif
