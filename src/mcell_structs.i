#define ORIENT_NOT_SET -100

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

/* Constants for notification levels */
enum notify_level_t {
  NOTIFY_NONE,  /* no output */
  NOTIFY_BRIEF, /* give a brief description (only used for a few types) */
  NOTIFY_FULL,  /* give a (possibly verbose) description */
};

/* Constants for warning levels */
enum warn_level_t {
  WARN_COPE, /* do something sensible and continue silently */
  WARN_WARN, /* do something sensible but emit a warning message */
  WARN_ERROR /* treat the warning and an error and stop */
};

enum overwrite_policy_t {
  FILE_UNDEFINED,  /* not specified */
  FILE_OVERWRITE,  /* always overwrite, even after a checkpoint */
  FILE_SUBSTITUTE, /* DEFAULT: append to entries earlier in time than "now", but
                      overwrite later entries */
  FILE_APPEND,        /* always append to file, even on a new run */
  FILE_APPEND_HEADER, /* always append to file, including the header, even on a
                         new run */
  FILE_CREATE, /* always create the file, or give an error if the file already
                  exists (to prevent overwriting) */
};

/* And finally we have some flags to say whether we're to count over */
/* the entire world or the volume enclosed by a region (set only one) */
#define REPORT_CONTENTS 1
#define REPORT_WORLD 0x20
#define REPORT_ENCLOSED 0x40
#define REPORT_TRIGGER 0x80

typedef unsigned char byte;

/* Timing pattern for molecule release from a release site. */
struct release_pattern {
  struct sym_entry *sym;   /* Symbol hash table entry for the pattern */
  double delay;            /* Delay between time 0 and first release event. */
  double release_interval; /* Time between release events within a train. */
  double train_interval; /* Time from the start of one train to the start of
                            the next one. */
  double train_duration; /* Length of the train. */
  int number_of_trains;  /* How many trains are produced. */
};

/* Data structure used to build boolean combinations of regions */
struct release_evaluator {
  byte op;    /* Region Expression Flags: the operation used */
  void *left; /* The left side of the expression--another evaluator or a region
                 object depending on bitmask of op */
  void *right; /* The right side--same thing */
};

/* Data that controls what output is written to a single file */
struct output_set {
  struct output_set *next;            /* Next data set in this block */
  struct output_block *block;         /* Which block do we belong to? */
  char *outfile_name;                 /* Filename */
  enum overwrite_policy_t file_flags; /* Overwrite Policy Flags: tells us how to
                                       * handle existing files */
  u_int chunk_count;    /* Number of buffered output chunks processed */
  char *header_comment; /* Comment character(s) for header */
  int exact_time_flag;  /* Boolean value; nonzero means print exact time in
                           TRIGGER statements */
  struct output_column *column_head; /* Data for one output column */
};

/* Data Output Timing Type */
/* Reaction and Viz data output timing */
enum output_timer_type_t {
  OUTPUT_BY_STEP,
  OUTPUT_BY_TIME_LIST,
  OUTPUT_BY_ITERATION_LIST,
};

/* Symbol hash table */
/* Used to parse and store user defined symbols from the MDL input file */
struct sym_table_head {
  struct sym_entry **entries;
  int n_entries;
  int n_bins;
};

/* Symbol hash table entry */
/* Used to parse and store user defined symbols from the MDL input file */
struct sym_entry {
  struct sym_entry *next; /* Chain to next symbol in this bin of the hash */
  int sym_type;           /* Symbol Type */
  char *name;             /* Name of symbol*/
  void *value;            /* Stored value, cast by sym_type */
};

/* Linked list of symbols */
/* Used to parse and retrieve user defined symbols having wildcards from the
 * MDL input file */
struct sym_table_list {
  struct sym_table_list *next;
  struct sym_entry *node; /* Symbol table entry matching a user input wildcard
                             string */
};

/* Linked list of numerical expressions */
/* Used for parsing MDL input file arithmetic expressions */
struct num_expr_list {
  struct num_expr_list *next;
  double value; /* Value of one element of the expression */
};
