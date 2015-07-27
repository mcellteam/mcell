#ifndef MESHPARSE
#define MESHPARSE

#include <stdio.h>

struct num_expr_list {
  struct num_expr_list *next;
  double value; /* Value of one element of the expression */
};

struct num_expr_list_head {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

enum partition_axis_t {
  X_PARTS, /* X-axis partitions */
  Y_PARTS, /* Y-axis partitions */
  Z_PARTS  /* Z-axis partitions */
};

struct sym_table {
  struct sym_table *next; /* Chain to next symbol in this bin of the hash */
  int sym_type;           /* Symbol Type */
  char *name;             /* Name of symbol*/
  void *value;            /* Stored value, cast by sym_type */
};

struct object_list {
  struct object *obj_head;
  struct object *obj_tail;
};

struct vector3 {
  double x;
  double y;
  double z;
};

/* Object Type Flags */
enum object_type_t {
  META_OBJ,     /* Meta-object: aggregation of other objects */
  BOX_OBJ,      /* Box object: Polygonalized cuboid */
  POLY_OBJ,     /* Polygon list object: list of arbitrary triangles */
  REL_SITE_OBJ, /* Release site object */
  VOXEL_OBJ,    /* Voxel object (so-far unused) */
};

struct object {
  struct object *next;        /* Next sibling object */
  struct object *parent;      /* Parent meta object */
  struct object *first_child; /* First child object */
  struct object *last_child;  /* Last child object */
  struct sym_table *sym;      /* Symbol hash table entry for this object */
  char *last_name; /* Name of object without pre-pended parent object name */
  enum object_type_t object_type; /* Object Type Flags */
  void *contents;    /* Actual physical object, cast according to object_type */
  u_int num_regions; /* Number of regions defined on object */
  struct region_list *regions; /* List of regions for this object */
  int n_walls;                 /* Total number of walls in object */
  int n_walls_actual;          /* Number of non-null walls in object */
  struct wall *walls;          /* Array of walls in object */
  struct wall **wall_p; // Array of ptrs to walls in object (used at run-time)
  int n_verts;               /* Total number of vertices in object */
  struct vector3 **vertices; /* Array of pointers to vertices
                                (linked to "all_vertices" array) */
  double total_area;      /* Area of object in length units */
  u_int n_tiles;          /* Number of surface grid tiles on object */
  u_int n_occupied_tiles; /* Number of occupied tiles on object */
  double t_matrix[4][4];  /* Transformation matrix for object */
  short is_closed;              /* Flag that describes the geometry
                                   of the polygon object (e.g. for sphere
                                   is_closed = 1 and for plane is 0) */
};

void object_list_singleton(struct object_list *head, struct object *objp);
void add_object_to_list(struct object_list *head, struct object *objp);
struct vector3 *point_scalar(double val);
void mcell_free_numeric_list(struct num_expr_list *nel);
int advance_range(struct num_expr_list_head *list, double tmp_dbl);
int generate_range(struct num_expr_list_head *list, double start, double end, double step);

#endif
