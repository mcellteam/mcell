/* Container data structure for all physical objects */
struct object {
  struct object *next;        /* Next sibling object */
  struct object *parent;      /* Parent meta object */
  struct object *first_child; /* First child object */
  struct object *last_child;  /* Last child object */
  struct sym_entry *sym;      /* Symbol hash table entry for this object */
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

  bool periodic_x; // This flag only applies to box objects BOX_OBJ. If set
  bool periodic_y; // any volume molecules encountering the box surface in the x,
  bool periodic_z; // y or z direction are reflected back into the box as if they
                   // had entered the adjacent neighboring box */
};

/* Release Shape Flags */
enum release_shape_t {
  SHAPE_UNDEFINED = -1,  /* Not specified */
 SHAPE_SPHERICAL,       /* Volume enclosed by a sphere */
  SHAPE_CUBIC,           /* Volume enclosed by a cube */
  SHAPE_ELLIPTIC,        /* Volume enclosed by an ellipsoid */
  SHAPE_RECTANGULAR,     /* Volume enclosed by a rect. solid */
  SHAPE_SPHERICAL_SHELL, /* Surface of a sphere */
  SHAPE_REGION,          /* Inside/on the surface of an arbitrary region */
  SHAPE_LIST             /* Individiaul mol. placement by list */
};


