#ifndef MDLPARSE_H
#define MDLPARSE_H

#include "mcell_structs.h"


#define RFLCT 0
#define TRANSP 1
#define SINK 2
#define WINDW 3

#define WILDCARD_PRESENT   0x1
#define TRIGGER_PRESENT    0x2
#define COUNT_PRESENT      0x4
#define EXPRESSION_PRESENT 0x8

#define ORIENT_NOT_SET SHRT_MIN


struct arg { 
  byte arg_type; /* DBL, STR */
  void *arg_value;
};


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
  /* used when parsing wildcard names for viz_output */
  struct sym_table_list *sym_table_list_head;
  int count_flags; /* Keep track of wildcards and TRIGGER statements */
  char *sym_name;


  char str_buf[1024];
  char str_buf2[1024];
  char temp_str[1024];
  char format_str[256];
  char time_str[128];
  char *a_str;

  u_int num_args;
  struct arg arg_list[ARGSIZE];

  double val_1;
  double val_2;
  double tmp_dbl;
  double *dblp;
  int *intp;

  struct species *specp;
  /* Theoretical average diffusion distances for the molecule */
  double l_perp_bar;                 
  double l_perp_rms;
  double l_r_bar;
  double l_r_rms;

  long long num_pos;
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
  struct state_list *surf_state_head;
  struct state_list *eff_state_head;
  struct state_list *mol_state_head;
  struct state_list *slp;

  struct output_block *obp;
  long long n_output;
  long long output_freq;
  struct counter *cp;
  struct mem_helper *sym_list_mem;

  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct element_data *edp;
  struct tet_element_data *tedp;
  struct element_connection_list *connection_head;
  struct element_connection_list *connection_tail;
  struct element_connection_list *eclp;
  struct element_connection_list *eclp_temp;
  struct vertex_list *vertex_head;
  struct vertex_list *vertex_tail;
  struct vertex_list *vlp;
  struct vertex_list *vlp_temp;
  int n_walls;
  int n_walls_actual;
  int n_verts;
  double box_aspect_ratio;

  struct voxel_object *vop;
  struct ordered_voxel *ovp;
  int n_voxels;

  struct region *rp;
  struct region_list *region_list_head;
  struct region_list *rlp;
  struct element_list *element_list_head;
  struct element_list *elmlp;
  char *region_name;
  int exclude_me;
  
  struct eff_dat *eff_dat_head;
  struct eff_dat *effdp;
  byte mol_quant_type;

  struct rxn *rxnp;
  struct mem_helper *path_mem;
  struct mem_helper *prod_mem;
  struct pathway *pathp;
  struct rxn_pathname *rxpnp;
  struct product *prodp;
  double fwd_km;
  double fwd_kcat;
  double bkw_km;
  double bkw_kcat;
  char *fwd_rate_filename;
  char *bkw_rate_filename;
  int bidirectional_arrow;
  int catalytic_arrow;

  short orient_specified;  
  short orient_class;
  short orient_class1;
  short orient_class2;
  byte prod_all_3d;

  struct file_stream *filep;

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
  
  char *header_comment;
  byte exact_time_flag;
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
void mdlerror_nested(char *s);
int mdlparse_init(struct volume *vol);

#endif
