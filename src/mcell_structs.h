#ifndef MCELL_STRUCTS
#define MCELL_STRUCTS

#include <limits.h>
#include <sys/types.h>
#include <stdio.h>

#include "rng.h"
#include "vector.h"
#include "mem_util.h"
#include "sched_util.h"

/*****************************************************/
/**  Brand new constants created for use in MCell3  **/
/*****************************************************/

/* Species flags */
   /* Walls have IS_SURFACE set, molecules do not. */
   /* Surface and grid molecules have ON_SURFACE set */
   /* Grid molecules have ON_GRID set */
   /* COUNT_ENCLOSED set if you count what happens inside closed region */
   /* (otherwise only count stuff happening at the surface) */
#define ON_SURFACE       0x01
#define ON_GRID          0x02
#define IS_SURFACE       0x04
#define NOT_FREE         0x07
#define TIME_VARY        0x08
#define CAN_MOLMOL       0x10
#define CAN_MOLGRID      0x20
#define CAN_MOLSURF      0x40
#define CAN_MOLWALL      0x80
#define CANT_INITIATE    0x100
#define COUNT_CONTENTS   0x1000
#define COUNT_HITS       0x2000
#define COUNT_SOME       0x3000
#define COUNT_RXNS       0x4000
#define COUNT_ENCLOSED   0x8000

/* rxn/mol/region counter report types */
/* Do not set both WORLD and ENCLOSED flags; ENCLOSED applies only to regions */
/* First set reports a single number */
#define REPORT_CONTENTS        0
#define REPORT_FRONT_HITS      1
#define REPORT_BACK_HITS       2
#define REPORT_FRONT_CROSSINGS 3
#define REPORT_BACK_CROSSINGS  4
#define REPORT_RXNS            5
/* Anything >= REPORT_MULTIPLE reports some combination of the above */
#define REPORT_MULTIPLE        6 
#define REPORT_ALL_HITS        7
#define REPORT_ALL_CROSSINGS   8
/* Concentration is kind of special. */
#define REPORT_CONCENTRATION   9
#define REPORT_ELAPSED_TIME    10
/* All basic report types can be masked with this value */
#define REPORT_TYPE_MASK       0x0F
/* And finally we have some flags to say whether we're to count over */
/* the entire world or the volume enclosed by a region (set only one) */
#define REPORT_WORLD           0x20
#define REPORT_ENCLOSED        0x40

/* rxn/mol/region counter flags */
/* Only set one of MOL_COUNTER or RXN_COUNTER */
/* Set ENCLOSING_COUNTER if the region is closed and counts inside itself */
#define MOL_COUNTER 1
#define RXN_COUNTER 2
#define ENCLOSING_COUNTER 4

#define MANIFOLD_UNCHECKED 0
#define NOT_MANIFOLD       1
#define IS_MANIFOLD        2

/* Reaction flags */
  /* RX_TRANSP signifies that a reaction is between a molecule and a TRANSPARENT wall */
  /* Any value equal to or less than RX_SPECIAL refers to a special wall type */
  /* RX_BLOCKED signals a reaction that cannot take place because the grid is full */
  /* Any value equal to or less than RX_NO_RX indicates that a reaction did not take place */
  /* RX_FLIP signals that a molecule flips its orientation (crosses a wall if it's free) */
  /* RX_DESTROY signals that the molecule no longer exists (so don't try to keep using it) */
  /* RX_A_OK signals that all is OK with a reaction, proceed as normal (reflect if you're free) */
  /* RX_NO_MEM signals a memory allocation error. */
#define RX_TRANSP  -3
#define RX_SPECIAL -3
#define RX_BLOCKED -2
#define RX_NO_RX   -2
#define RX_FLIP    -1
#define RX_LEAST_VALID_PATHWAY 0
#define RX_DESTROY  0
#define RX_A_OK     1
#define RX_NO_MEM   3


/* Flags for BSP trees to determine whether something is a node or a branch */
/* Will either have BRANCH_XN through _ZP, or _L, _R, _X, _Y, _Z. */
/* P is positive, N is negative. */
#define BRANCH_XN  0x01
#define BRANCH_XP  0x02
#define BRANCH_YN  0x04
#define BRANCH_YP  0x08
#define BRANCH_ZN  0x10
#define BRANCH_ZP  0x20

#define BRANCH_L  0x01
#define BRANCH_R  0x02
#define BRANCH_X  0x04
#define BRANCH_Y  0x08
#define BRANCH_Z  0x10


/* Constants equating integers with coordinates */
#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#define X_NEG 0
#define X_POS 1
#define Y_NEG 2
#define Y_POS 3
#define Z_NEG 4
#define Z_POS 5
#define NODIR 6

#define X_NEG_BIT 0x01
#define X_POS_BIT 0x02
#define Y_NEG_BIT 0x04
#define Y_POS_BIT 0x08
#define Z_NEG_BIT 0x10
#define Z_POS_BIT 0x20


/* Collision types for rays striking surfaces */
#define COLLIDE_REDO    -1
#define COLLIDE_MISS    0
#define COLLIDE_FRONT   1
#define COLLIDE_BACK    2

#define COLLIDE_MOL_M   3
#define COLLIDE_MOL_SP  4
#define COLLIDE_MOL_SN  5

#define COLLIDE_SV_NX   6
#define COLLIDE_SV_PX   7
#define COLLIDE_SV_NY   8
#define COLLIDE_SV_PY   9
#define COLLIDE_SV_NZ   10
#define COLLIDE_SV_PZ   11

#define COLLIDE_MASK    0x0F

#define COLLIDE_WALL    0x10
#define COLLIDE_MOL     0x20
#define COLLIDE_SUBVOL  0x40


/* Types for things we can hit */
#define VOL_COLLISION    1
#define WALL_COLLISION   2   
#define MOL_COLLISION    3


/* Flags for edges. */
/* BARE edges do not connect to anything. */
/* SHARED edges have walls with the same coordinate frame */
/* ROTONLY edges have walls whose coordinate frames are related by rotation */
/* TRANSROT edges require translation and rotation to move between walls */
#define EDGE_BARE     0
#define EDGE_SHARED   1
#define EDGE_ROTONLY  2
#define EDGE_TRANSROT 3


/* Size constants */
#define EPS_C 1e-12
#define GIGANTIC 1e140

/* Special rate constants used to define unusual reactions */
#define KCAT_RATE_TRANSPARENT -1.0


/* Abstract molecule flags */
/* RULES: only one of TYPE_GRID, TYPE_3D, TYPE_SURF set. */
/*   Only TYPE_3D and TYPE_SURF can ACT_DIFFUSE */
/*   ACT_NEWBIE beats ACT_INERT beats ACT_REACT */
/*   Can free up memory when nothing in IN_MASK */

/* Molecule type--grid molecule, 3D molecule, or surface molecule */
#define TYPE_GRID 0x001
#define TYPE_3D   0x002
#define TYPE_SURF 0x004
#define TYPE_MASK 0x007

/* NEWBIE molecules get scheduled before anything else happens to them. */
/* INERT molecules don't react, REACT molecules do */
/* CHANGE molecules have had their rate constant changed */
/* DIFFUSE molecules diffuse (duh!) */
#define ACT_DIFFUSE 0x008
#define ACT_INERT   0x010
#define ACT_REACT   0x020
#define ACT_NEWBIE  0x040
#define ACT_CHANGE  0x080

/* Flags telling us which linked lists the molecule appears in. */
#define IN_SCHEDULE 0x100
#define IN_SURFACE  0x200
#define IN_VOLUME   0x400
#define IN_MASK     0x700

/* Flags telling us what our counting status is */
#define COUNT_ME    0x800

/* End of abstract molecule flags. */


/* How big will we let the reaction table get? */
/* 0x100000 = 2 million */
#define MAX_RX_HASH 0x100000

/* How big will we let the count-by-region table get? */
/* 0x10000 = 128K */
#define MAX_COUNT_HASH 0x10000 

/* mask for count-by-region hash */
#define COUNT_HASHMASK 0xffff

/* What's the upper bound on the number of coarse partitions? */
#define MAX_COARSE_PER_AXIS 16
#define MIN_COARSE_PER_AXIS 6
#define MAX_TARGET_TIMESTEP 1.0e6
#define MIN_TARGET_TIMESTEP 10.0


/* Shapes for release sites */
#define SHAPE_SPHERICAL 0
#define SHAPE_CUBIC 1
#define SHAPE_ELLIPTIC 2
#define SHAPE_RECTANGULAR 3
#define SHAPE_SPHERICAL_SHELL 4


/* Flags for parser to indicate which axis we are partitioning */

#define X_PARTS 0
#define Y_PARTS 1
#define Z_PARTS 2


/* Flags for the exact disk computation */
#define EXD_HEAD  0
#define EXD_TAIL  1
#define EXD_CROSS 2
#define EXD_SPAN  3
#define EXD_OTHER 4


/*********************************************************/
/**  Constants used in MCell3 brought over from MCell2  **/
/*********************************************************/


/* Generic numerical constants */
#define EPSILON 1e-14
                                                                                
#define R_UINT_MAX 2.3283064365386962890625e-10
#define MY_PI 3.14159265358979323846
#define N_AV 6.0221415e23
#define ROUND_UP 0.5
                                                                                
/* Wall element shapes */
#define RECT_POLY 0
#define TRI_POLY 1
#define GEN_POLY 2
   /* RECT_POLY is rectangular */
   /* TRI_POLY is triangular */
   /* GEN_POLY is non-rectangular with >3 vertices */


/* Surface grid shapes */
#define RECTANGULAR 0
#define TRIANGULAR 1


/* Orientations */
/* Relative to a surface */
#define OUTWRD 1
#define INWRD -1

/* Relative to a molecule */
#define POS_POLE 1
#define NEG_POLE -1
#define POLE POS_POLE
                                                                                

/* Grid molecule site placement types */
/* Place either a certain density or an exact number of effectors */
#define EFFDENS 0
#define EFFNUM 1


/*******************************************************/
/**  Old constants copied from MCell2, may be broken  **/
/*******************************************************/

/* Parser parameters.  Probably need to be revisited. */
/* size of symbol hash table */
#define SYM_HASHSIZE 0x20000

/* mask for symbol table hash */
#define SYM_HASHMASK 0x1ffff

/* maximum number of args allowed in MDL "C"-style print statements */
#define ARGSIZE 255

/* maximum allowed nesting level of INCLUDE_FILE statements in MDL */
#define MAX_INCLUDE_DEPTH 16

/* default size of output count buffers */
#define COUNTBUFFERSIZE 10000
                                                                                

/* Data types to be stored in MDL parser symbol table: */

/* chemical reaction: */
#define RX 1

/* name of chemical reaction: */
#define RXPN 2

/* molecule type (i.e. species): */
#define MOL 3

/* polygon or box object: */
#define POLY 4

/* release site object: */
#define RSITE 5

/* meta-object: */
#define OBJ 6

/* release pattern: */
#define RPAT 7

/* object region: */
#define REG 8

/* integer type (used only for COUNT statements): */
#define INT 9

/* double: */
#define DBL 10

/* string: */
#define STR 11

/* array of doubles: */ 
#define ARRAY 12

/* file stream type for "C"-style file-io: */
#define FSTRM 13

/* expression type (used only in COUNT statements): */
#define EXPR 14

/* temporary place-holder type for assignment statements: */
#define TMP 15

                                                                                
/* Object types */
#define META_OBJ 0
#define BOX_OBJ 1
#define POLY_OBJ 2
#define REL_SITE_OBJ 3
                                                                                

/* Box sides */
#define TP 0
#define BOT 2
#define FRNT 4
#define BCK 6
#define LFT 8
#define RT 10
#define ALL_SIDES INT_MAX
                   
                                                             
/* Viz state values */
#define EXCLUDE_OBJ INT_MIN


/* output evaluator specifications. Broken until we finish rxn counting */
#define OVER_E 0
#define EACH_E 1
#define SPEC_E 2
                                                                                
#define OVER_L 0
#define EACH_L 1
#define SPEC_L 2
                                                                                
#define SUM 0
#define DT 1
#define CUM 2
                                                                                
#define A_EVENTS 0
#define INIT_EVENTS 1
#define INTER_EVENTS 2


/* Output evaluator index types. */
#define UNKNOWN 0
#define TIME_STAMP_VAL 1
#define INDEX_VAL 2


/* Reaction and Viz data output timing */
#define OUTPUT_BY_STEP 0 
#define OUTPUT_BY_TIME_LIST 1
#define OUTPUT_BY_ITERATION_LIST 2


/* Region counter type.  INIT probably broken. */
#define RX_STATE 0
#define INIT_TRANS 1
#define TRANSITIONS 2
#define MOL_TRANS_EACH 3
#define MOL_TRANS_ALL 4


/* Visualization stuff. */
/* Visualization modes. */
#define NO_VIZ_MODE 0
#define DX_MODE 1
#define DREAMM_V3_MODE 2
#define RK_MODE 3 
#define ASCII_MODE 4


/* Visualization frame data types. */
/* Used to select type of data to include in viz output files */
/* Will probably change significantly when we redesign DReAMM output format */
#define ALL_FRAME_DATA 0
#define EFF_POS 1
#define EFF_STATES 2
#define MOL_POS 3
#define MOL_STATES 4
#define MOL_POS_STATES 5
#define SURF_POS 6
#define SURF_STATES 7


/* release number methods */
#define CONSTNUM 0
#define GAUSSNUM 1
#define VOLNUM 2



/**********************************************/
/**  New/reworked structures used in MCell3  **/
/**********************************************/

typedef unsigned char byte;


/* If you don't include sys/types.h, #define SYS_TYPES_NOT_LOADED so */
/* you get the u_short/int/long set of types */
#ifdef SYS_TYPES_NOT_LOADED
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
#endif


/* Properties of one type of molecule or surface */
struct species
{
  u_int species_id;             /* Unique ID for this species */
  u_int hashval;                /* Hash value (may be nonunique) */
  struct sym_table *sym;        /* Symbol table entry (name) */
  struct eff_dat *eff_dat_head; /* if IS_SURFACE this points to head of
                                   effector data list associated with 
                                   surface class */
  
  u_int population;             /* How many of this species exist? */
  
  double D;                     /* Diffusion constant */
  double D_ref;                 /* Reference diffusion constant */
  double radius;                /* Molecular radius */
  double space_step;            /* Characteristic step length */
  double time_step;             /* Minimum (maximum?) sensible timestep */
/*short charge;*/               /* Electric charge. */
  u_short flags;                /* Free?  Membrane bound?  Membrane? */
  
  int viz_state;                /* Visualization state for output */
  byte checked;                 /* Bread crumb for graph traversal */
};


/* All pathways leading away from a given intermediate */
struct rxn
{
  struct rxn *next;          /* Next reaction with these reactants
                                but differing geometry */
  struct sym_table *sym;     /* ptr to symbol table entry for this rxn */
  
  u_int n_reactants;         /* How many reactants? (At least 1.) */
  int n_pathways;            /* How many pathways lead away? (Negative = special reaction)*/
  u_int *product_idx;        /* Index of 1st player for products of pathway */
  double *cum_probs;         /* Cumulative probabilities for (entering) all pathways */
  double *cat_probs;         /* Probabilities of leaving all pathways (<=0.0 is instant) */
  
  struct species **players;  /* Identities of reactants/products */
  short *geometries;         /* Geometries of reactants/products */

  struct t_func *prob_t;     /* List of probabilities changing over time */
  
  struct pathway *pathway_head; /* list/array of pathways built at parse-time */
};


/* Sets of reactions grouped by...? */
struct rxn_group
{
/* Someone else gets to fill this in. */
};


struct rxn_pathname {
  struct sym_table *sym;
  u_int hashval;
  struct pathway *path;
};


/* Parse-time structure for reaction pathways */
struct pathway {
  struct pathway *next;
  struct rxn_pathname *pathname; /* data for named reaction pathway or NULL */
  double count;                  /* How many times have we taken this path? */
  struct species *reactant1;     /* First reactant in reaction pathway */
  struct species *reactant2;     /* Second reactant (NULL if none) */
  struct species *reactant3;     /* Third reactant--surface type or NULL */
  double km;                     /* Rate constant */
  double kcat;                   /* Catalytic dead time */
  char* km_filename;             /* Filename for time-varying rates */
  short orientation1;            /* Orientation of first reactant */
  short orientation2;            /* Orientation of second reactant */
  short orientation3;            /* Orientation of third reactant */
  struct product *product_head;  /* Linked lists of species created */
  struct pathway_count_request *pcr;  /* Who is counting us? */
};

/* Parse-time structure for products of reaction pathways */
struct product {
  struct product *next;
  struct species *prod;          /* Molecule type to be created */
  short orientation;             /* Orientation to place molecule */
};


/* Piecewise constant function for time-varying reaction rates */
struct t_func
{
  struct t_func *next;
  double time;                    /* Time to switch to next rate */
  double value;                   /* Current rate */
  int path;                       /* Which rxn pathway is this for? */
};


/* Abstract structure that starts all molecule structures */
struct abstract_molecule
{
  struct abstract_molecule *next;  /* Next molecule in scheduling queue */
  double t;                        /* Scheduling time. */
  double t2;                       /* Dead time for catalysis & such. */
  short flags;                     /* Who am I, what am I doing, etc. */
  struct species *properties;      /* What type of molecule are we? */
  struct mem_helper *birthplace;   /* What was I allocated from? */
};


/* Freely diffusing or fixed molecules in solution */
struct molecule
{
  struct abstract_molecule *next;
  double t;
  double t2;
  short flags;
  struct species *properties;
  struct mem_helper *birthplace;
  
  struct vector3 pos;             /* Position in space */
  struct subvolume *subvol;       /* Partition we are in */
  
  struct cmprt_data *curr_cmprt;  /* Compartment we are in (for counting) */
  double birthday;                /* Time at which this particle was born */
  
  struct surface_grid *previous_grid;   /* Wall we were released from */
  int index;                            /* Index on that wall (don't rebind) */
  
  struct molecule *next_v;        /* Next molecule in this subvolume */
};


/* Freely diffusing or fixed molecules on a surface */
/* Same as struct molecule with extra stuff tacked on end */
struct surface_molecule
{
  struct abstract_molecule *next;
  double t;
  double t2;
  short flags;
  struct species *properties;
  struct mem_helper *birthplace;
  
  struct vector3 pos;               /* Position in the world */
  struct subvolume *subvol;         /* Partition we're in */

  struct region_data *curr_region;  /* Region we are in (for counting) */

  struct molecule *next_v;          /* Next molecule in this volume */

  short orient;                     /* Facing up or down? */
  struct vector2 s_pos;             /* Position in surface coordinates */
  struct wall *curr_wall;           /* The surface element we are on */
  
  struct surface_molecule *next_s;  /* Next molecule on this surface */
};


/* Fixed molecule on a grid on a surface */
struct grid_molecule
{
  struct abstract_molecule *next;
  double t;
  double t2;
  short flags;
  struct species *properties;
  struct mem_helper *birthplace;
  
  int grid_index;              /* Which gridpoint do we occupy? */
  short orient;                /* Which way do we point? */
  struct surface_grid *grid;   /* Our grid (which tells us our surface) */
};


struct edge
{
  struct wall *forward;     /* For which wall is this a forwards transform? */
  struct wall *backward;    /* For which wall is this a reverse transform? */
  
  struct vector2 translate;  /* Translation vector */
  double cos_theta;          /* Cosine of angle between bases */
  double sin_theta;          /* Sine of angle between bases */
  
  double length;             /* Length of the edge */
  double length_1;           /* Reciprocal of length of edge */
  
  int flags;
};


struct wall
{
  struct wall *next;              /* Next wall in the universe */
  
  struct species *surf_class;      /* Surface class for this wall */

  int side;                       /* index of this wall in its parent object */

/*
  int projection;
*/

  struct vector3 *vert[3];        /* Array of pointers to vertices */
  
  double uv_vert1_u;              /* Surface u-coord of 2nd corner (v=0) */
  struct vector2 uv_vert2;        /* Surface coords of third corner */

  struct edge *edges[3];          /* Array of pointers to each edge. */
  struct wall *nb_walls[3];       /* Array of pointers to neighboring walls */

  double area;                    /* Area of this element */
  
  struct vector3 normal;          /* Normal vector for this wall */
  struct vector3 unit_u;          /* U basis vector for this wall */
  struct vector3 unit_v;          /* V basis vector for this wall */
  double d;                       /* Distance to origin (point normal form) */
  
  struct surface_molecule *mol;   /* Head of list of surface molecules */
  int mol_count;                  /* How many surface molecules? */
  
  struct surface_grid *effectors; /* Grid of effectors for this wall */
  
  int viz_state;                  /* For display purposes--is short enough? */
  u_short flags;                  /* Flags for whether we need to count */

  struct object *parent_object;   /* The object we are a part of */
  struct storage *birthplace;     /* Where we live in memory */
  
  struct region_list *regions;    /* Regions that contain this wall */
};


/* Linked list of walls (for subvolumes) */
struct wall_list
{
  struct wall_list *next;   /* The next entry in the list */
  
  struct wall *this_wall;        /* The wall in this entry */
};


/* 3D vector of short integers */
struct short3D
{
  short x;
  short y;
  short z;
};


/* Point in space that will tell us which compartments we're in */

struct waypoint
{
  struct vector3 loc;              /* This is where the waypoint is */
  struct region_list *regions;     /* We are inside these regions */
  struct region_list *antiregions; /* We are out of these regions */
};


/* Vertices are stored in unbalanced ternary trees ordered by z-coord */
/* Why do we need to bother?  Also, this is a weird tree structure. */
struct vertex_tree
{
  struct vertex_tree *next;         /* Vertices with same z (this is weird) */
  struct vertex_tree *above;        /* Vertices with z larger than ours */
  struct vertex_tree *below;        /* Vertices with z smaller than ours */
  
  struct vector3 loc;               /* Vertex coordinate itself */
};


/* Contains space for molcules, walls, wall_lists, etc. */
struct storage
{
  struct mem_helper *list;  /* Wall lists */
  struct mem_helper *mol;   /* Molecules */
  struct mem_helper *smol;  /* Surface molecules */
  struct mem_helper *gmol;  /* Grid molecules */
  struct mem_helper *face;  /* Walls */
  struct mem_helper *join;  /* Edges */
  struct mem_helper *tree;  /* Vertices */
  struct mem_helper *effs;  /* Effector grids */
  struct mem_helper *coll;  /* Collision list */
  struct mem_helper *regl;  /* Region lists */
  struct mem_helper *exdv;  /* Vertex lists for interaction disk area */
  
  struct wall *wall_head;              /* Locally stored walls */
  int wall_count;                      /* How many local walls? */
  struct vertex_tree *vert_head;       /* Locally stored vertices */
  int vert_count;                      /* How many vertices? */
  
  struct schedule_helper *timer;       /* Local scheduler */
  double current_time;                 /* Local time */
  double max_timestep;                 /* Local maximum timestep */
};

/* Linked list of storage areas. */
struct storage_list
{
  struct storage_list *next;
  struct storage *store;
};


/* Walls and molecules in a spatial subvolume */
struct subvolume
{
  struct wall_list *wall_head; /* Head of linked list of intersecting walls */
  struct wall_list *wall_tail; /* Tail of list of walls */
  int wall_count;              /* How many walls intersect? */
  
  struct molecule *mol_head;   /* Head of linked list of molecules */
  int mol_count;               /* How many molecules are here? */
  
  int index;                   /* Index of subvolume (parallelization?) */
  
  struct short3D llf;          /* Indices of left lower front corner */
  struct short3D urb;          /* Indices of upper right back corner */
  
  short is_bsp;                /* Flags saying what the void pointers are */

  void *neighbor[6];           /* Subvolume or bsp_tree across each face */
  
  struct storage *local_storage;         /* Local storage */
};


/* Binary space partitioning tree for subvolume connections */
struct bsp_tree
{
  void *left;        /* The tree below the partition */
  void *right;       /* The tree above the partition */
  short partition;   /* The index of the partition */
  short flags;       /* Coordinate that is split, plus terminal node flags */
};


struct rxn_counter_data
{
  struct rxn_pathname *rxn_type;    /* named rxn we are counting */
  double n_rxn_at;                  /* # rxn occurrance on surface */
  double n_rxn_enclosed;            /* # rxn occurrance inside closed region */
};


struct move_counter_data
{
  struct species *mol_type;         /* species we are counting */
  double front_hits;               /* # hits on front of region (normal up) */
  double back_hits;                /* # hits on back of region */
  double front_to_back;            /* # crossings from front to back */
  double back_to_front;            /* # crossings from back to front */
  double scaled_hits;              /* To determine integrated concentration */
  int n_at;                        /* # molecules on region surface */
  int n_enclosed;                  /* # molecules inside closed region */
};


union counter_data
{
  struct rxn_counter_data rx;
  struct move_counter_data move;
};


/* Struct to count rxns or molecules within regions (where "within" includes */
/* on the inside of a fully closed surface, for 3D molecules) */
struct counter
{
  struct counter *next;
  byte counter_type;               /* MOL_COUNTER or RXN_COUNTER */
  struct region *reg_type;         /* Region we are counting on */
  union counter_data data;         /* data we are counting:
                                      reference data.move for move counter
                                      reference data.rx for rxn counter */
};


/* All data about the world */
struct volume
{
/*  struct vector3 corner[8];*/  /* Corners of the world */
#if 0
  struct vector3 llf;           /* left lower front corner of world */
  struct vector3 urb;           /* upper right back corner of world */
#endif
  
  int nx_parts;                 /* Number of coarse X partition boundaries */
  int ny_parts;                 /* Number of coarse Y partition boundaries */
  int nz_parts;                 /* Number of coarse Z partition boundaries */
  double *x_partitions;         /* Coarse X partition boundaries */
  double *y_partitions;         /* Coarse Y partition boundaries */
  double *z_partitions;         /* Coarse Z partition boundaries */
  
  int n_fineparts;        /* Number of fine partition boundaries (multiple of n_axis_partitions) */
  double *x_fineparts;           /* Fine X partition boundaries */
  double *y_fineparts;           /* Fine Y partition boundaries */
  double *z_fineparts;           /* Fine Z partition boundaries */
  
  int n_waypoints;              /* How many of these = (n_axis_p-3)^3 */
  struct waypoint *waypoints;   /* Contains compartment information */
  byte place_waypoints_flag;         /* We need to place waypoints
                                   if we count 3D diffusing molecules
                                   in regions */
  
  int n_subvols;                /* How many coarse subvolumes? */
  struct subvolume *subvol;     /* Array containing all subvolumes */
   
  int binning;                  /* How many real partitions per one in lookup? */
  struct subvolume **lookup;     /* 3D lookup array pointing at subvolumes */
  
  int n_walls;                  /* Total number of walls */
  int n_verts;
  
  int rx_hashsize;                 /* How many entries in our reaction 
                  			hash table? */
  int n_reactions;              /* How many reactions are there, total? */
  struct rxn **reaction_hash;   /* A hash table of all reactions. */
  struct mem_helper *rxn_mem;   /* Memory to store time-varying reactions */
  
  int count_hashmask;         /* Mask for looking up count hash table */
  struct counter **count_hash;/* count hash table */
  
  int n_species;                /* How many different species (molecules)? */
  struct species **species_list; /* Array of all species (molecules). */
  
  struct schedule_helper *releaser;
  
  struct mem_helper *storage_allocator;      /* Storage for storage list */
  struct storage_list *storage_head;   /* Linked list of all storage */
  
  double speed_limit;           /* How far can the fastest particle get in one timestep? */
  /*counts the number of times the molecules have had their positions updated */
  double diffusion_number;
  /* counts the number of timesteps that molecules have diffused for */
  double diffusion_cumsteps;

  struct schedule_helper *count_scheduler;
  
  /* Simulation initialization parameters  */
  struct sym_table **main_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct release_pattern *default_release_pattern;
  struct release_event_queue *release_event_queue_head;
  struct output_evaluator *count_zero;
  struct output_block *output_block_head;
  struct viz_obj *viz_obj_head;
  struct frame_data_list *frame_data_head;
  struct mem_helper *pathway_requester;
  struct species *g_mol;
  struct species *g_surf;
  double time_unit;
  double time_step_max;
  double space_step;
  double length_unit;
  double effector_grid_density;
  double rx_radius_3d;
  double *r_step;
  double *d_step;
  double r_num_directions;
  double sim_elapsed_time;
  double chkpt_elapsed_time;
  double chkpt_elapsed_time_start;
  double current_time;
  double current_start_time;
  double max_diffusion_step;
  double random_number_use;
  double ray_voxel_tests;
  double ray_polygon_tests;
  double ray_polygon_colls;
  double mol_mol_colls;
  double diffusion_steps;
  struct vector3 bb_min;	/**< bounding box minimum size */
  struct vector3 bb_max;	/**< bounding box maximum size */
  u_int tot_mols;
  struct rng_state *rng;
  u_int init_seed;
  long long it_time;
  double elapsed_time;
  long long start_time;
  u_int radial_directions;
  u_int radial_subdivisions;
  u_int num_directions;
  int directions_mask;
  int fully_random;
  int procnum;			/**< procedure number */
  int viz_mode;
  struct rk_mode_data *rk_mode_var;
  byte voxel_image_mode;
  byte voxel_volume_mode;
  char *molecule_prefix_name;
  double my_counter;

  /* MCell startup command line arguments */
  byte info_opt;
  u_int seed_seq;
  long long iterations;
  char *log_file_name;
  FILE *log_file;
  FILE *err_file;
  u_int log_freq;
  u_int chkpt_init;
  u_int chkpt_flag;
  long long chkpt_iterations;
  u_int chkpt_seq_num;
  char *chkpt_infile;
  char *chkpt_outfile;
  FILE *chkpt_infs;
  FILE *chkpt_outfs;
  FILE *chkpt_signal_file_tmp;
  char *mdl_infile_name;
  char *curr_file;
};


/* Grid over a surface containing grid_molecules */
struct surface_grid
{
  int n;                   /* Number of slots along each axis */

  double inv_strip_wid;    /* Reciprocal of the width of one strip */
  double vert2_slope;      /* Slope from vertex 0 to vertex 2 */
  double fullslope;        /* Slope of full width of triangle */
  struct vector2 vert0;    /* 2D coordinates of vertex zero */
  
  double binding_factor;   /* Binding probability correction factor for surface area */
  
  u_int n_tiles;           /* Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2) */
  u_int n_occupied;        /* Number of tiles occupied by grid_molecules */
  struct grid_molecule **mol;  /* Array of pointers to grid_molecule for each tile */
  
  int index;               /* Unique index into effector_table */
  
  struct subvolume *subvol;/* Best match for which subvolume we're in */
  struct wall *surface;    /* The wall that we are in */
};


/* Data structure to store information about collisions. */
struct collision
{
  struct collision *next;
  double t;                     /* Time of collision (may be slightly early) */
  
  void *target;                 /* Thing that we hit */
  int what;                     /* What kind of thing did we hit? */
  struct rxn *intermediate;     /* Reaction that told us we could hit it */
  struct vector3 loc;           /* Location of impact */
};


/* Data structures to store information about interaction disk geometry */
struct exd_vertex
{
  struct exd_vertex *next;
  double u,v;              /* x,y style coordinates */
  double r2,zeta;          /* r,theta style coordinates */
  struct exd_vertex *e;    /* Edge to next vertex */
  struct exd_vertex *span; /* List of edges spanning this point */
  int role;                /* Head, tail, whatever */
};



/* Data structures to describe release events */
struct release_event_queue {
  struct release_event_queue *next;
  double event_time;			/**< time of the release */

  struct release_site_obj *release_site;
  struct vector3 location;		/**< location of the release */
  int train_counter;			/**< counts executed trains */
  double train_high_time;		/**< time of the train's start */
};


struct release_site_obj {
	struct vector3 *location;	/**< location of release site */
	struct species *mol_type;	/**< species to be released */
	byte release_number_method;
	byte release_shape;
	int release_number;
	int mean_number;
	double mean_diameter;
	double concentration;
        double standard_deviation;
	struct vector3 *diameter;	

	double release_prob;
	struct release_pattern *pattern;
};

/* Timing pattern for molecule release from a release site. */
struct release_pattern {
        struct sym_table *sym;
	double delay;			/**< delay between time 0 
                                           and first release event. */
	double release_interval;	/**< time between release events 
                                           within a train. */
	double train_interval;		/**< time from the start of one train 
                                           to the start of the next one. */
	double train_duration;		/**< length of the train. */
	int number_of_trains;		/**< how many trains are produced. */
};

/* Holds information about a box with rectangular patches on it. */
struct subdivided_box
{
  int nx;
  int ny;
  int nz;
  double *x;
  double *y;
  double *z;
};

/******************************************************************/
/**  Everything below this line has been copied from MCell 2.69  **/
/******************************************************************/

/**
 * Linked list of all reaction data output blocks
 */
struct output_block {
	struct output_block *next;   /**< next reaction data output block*/
        double t;                   /**< scheduled time to update counters
                                         associated with this output block*/
	byte timer_type;    /**< output timer type: OUTPUT_BY_STEP,
                                     OUTPUT_BY_TIME_LIST,
                                     OUTPUT_BY_ITERATION_LIST */
	double step_time;                     /**< output frequency (secs)*/
        struct num_expr_list *time_list_head; /**< list of output times (secs)
                                                   or list of iterations*/
        struct num_expr_list *curr_time_ptr; /**< ptr to current position in
                                                   output time list*/
        u_int buffersize;                     /**< output chunk size*/
        u_int curr_buf_index;              /**< index to current position
                                                   in output buffers*/
        double *time_array;                 /**< array of output times
                                                 for current output chunk */
        u_int chunk_count;                  /**< number of chunks processed*/
	struct output_item *output_item_head; /**< list of count output
                                                     statements associated with
                                                     this output block*/
};


/**
 * Linked list of output statements associated with a given output block 
 */
struct output_item {
	struct output_item *next;  
	char *outfile_name;               /**< name of file to contain output*/
	struct output_evaluator *output_evaluator_head;  /**< list of counters 
                                                  associated with this
                                                  count output statement*/
	struct output_evaluator *count_expr;  /**< root of count expression tree
                                               to be evaluated for this
                                               count output statement*/
	struct output_item *next_column;
};


/**
 * Linked list of output evaluators to be evaluated at update time
 */
struct output_evaluator {
	struct output_evaluator *next;    /**< next item in count list*/
	byte update_flag;           /**< counter update necessary?*/
	byte reset_flag;            /**< reset temp_data to 0 on each iteration?*/
	byte index_type;            /**< flag indicating final_data is to be
	                              indexed by either
	                              TIME_STAMP_VAL or INDEX_VAL*/
	byte data_type;             /**< type of data to track:
                                      EXPR INT DBL*/
	u_int n_data;               /** buffer size */
	void *temp_data;            /**< ptr to intermediate data
                                         specified by type*/
	void *final_data;           /**< ptr to final outputable data
                                         specified by type*/
	struct output_evaluator *operand1;
	struct output_evaluator *operand2;
	char oper;
};


struct lig_output_evaluator {
	struct lig_output_evaluator *next;
        struct output_evaluator *output_evaluator;
};


struct lig_count_ref {
	struct lig_count_ref *next;
	unsigned short type;
        char *full_name;
        struct output_evaluator *output_evaluator;
};


struct viz_state_ref {
	struct viz_state_ref *next;
	int viz_state;
        char *full_name;
};

struct pathway_count_request
{
  struct pathway_count_request *next;
  struct output_evaluator *requester;
};




/**
 * Compartment. [\todo need more info]
 */
struct cmprt_data {
	struct cmprt_data *next;
	struct sym_table *sym;
        char *full_name;
        byte fully_closed;
	int instance;
	int *lig_count;
	double *conc;
	double volume;
	double vm;
        int n_corners;
        int n_walls;
	struct vector3 *corner;
	struct vector3 *vertex_normal;
	struct vector3 *normal;
	struct wall_list *wall_list;
	struct wall **wall;
	struct cmprt_data **neighbor;
};

struct cmprt_data_list {
	struct cmprt_data_list *next;
	struct cmprt_data *cmprt_data;
};

/**
 * A polygon list object, part of a surface.
 */
struct polygon_object {
	struct lig_count_ref *lig_count_ref;
					/**< ptr to list of lig_count_ref
	                                   structures: one for each time polygon
	                                   object is referenced in count stmt */
	struct viz_state_ref *viz_state_ref;
					/**< ptr to list of viz_state_ref
	                                   structures: one for each time polygon
	                                   object is referenced in
					   STATE_VALUES block */
        struct ordered_poly *polygon_data; /**< pointer to data structure
                                                holding polygon vertices etc... */
	struct subdivided_box *sb;      /**< Holds corners of box if necessary */
	int n_walls;			/**< Number of polygons in
                                             polygon object */
        int n_verts;                    /**< Number of vertices in
                                             polygon object */
	byte fully_closed;		/**< flag indicating closure of object */
        struct species **surf_class;    /** array of pointers to surface class, 
                                            one for each polygon */
	struct bit_array *side_removed; /**< Bit array; if bit is on, side is removed */
/*        struct eff_dat **eff_prop;*/	/**< array of ptrs to eff_dat data
					   structures, one for each polygon. */
};


/**
 * A general ordered polyhedron. 
 * That is, the vertices of each polygonal face are ordered according to
 * the right hand rule.
 */
struct ordered_poly {
	struct vector3 *vertex;         /**< Array of polygon vertices */
	struct vector3 *normal;         /**< Array of polygon normals */
	struct element_data *element; /**< Array element_data
                                              data structures */
	int n_verts;                    /**< Number of vertices in polyhedron */
	int n_walls;                    /**< Number of polygons in polyhedron */
};

/**
 * Data structure used to build one polygon.
 * This data structure is used to store the data from the MDL file
 * and to contruct each polygon of a polygon object.
 */
struct element_data {
        int vertex_index[3];              /**< Array of vertex indices forming a
                                           polygon. */
	int n_verts;                    /**< Number of vertices in polygon (always 3). */
};

/**
 * A compartment.
 */
struct cmprt {
	struct sym_table *sym;
        unsigned short type;
	int inst_count;
	struct lig_output_evaluator **lig_output_evaluator;  /**< array of ptrs to lig_output_evaluator
	                                        structures: one for each
	                                        ligand type */
        byte a_zone_lig;
        unsigned short side_stat[6];
        struct vector3 *vert1;
        struct vector3 *vert2;
        struct vector3 *a_zone_loc;
        byte *lig_prop[6];
        struct eff_dat *eff_prop[6];
	int *count_freq;
        int color[6];
};

/**
 * Linked list of surface effector placement data.
 */
struct eff_dat {
        struct eff_dat *next;
        struct species *eff; /* effector species to place on surface */
        byte quantity_type;  /* type is either EFFDENS or EFFNUM */
        double quantity; /* amount of effectors to place: density or number */
	short orientation;
};

/**
 * Linked list of elements.
 * [\todo what is this?]
 */
struct element_list {
        struct element_list *next;
        u_int begin;
        u_int end;
	struct element_special *special;
};

/* Elements can be patches on boxes or other regions */
struct element_special
{
  struct vector3 corner1;
  struct vector3 corner2;
  struct region *referent;
  byte exclude;
};

/**
 * Region of an object
 * If region is a manifold then it can be used as a volume and surface region.
 * Otherwise it can only be used as a surface region.
 */
struct region {
	struct sym_table *sym;
	u_int hashval;
        char *region_last_name;
	struct object *parent;
	struct element_list *element_list_head;
	struct bit_array *membership;
	struct reg_counter_ref_list *reg_counter_ref_list;
	struct eff_dat *eff_dat_head;
        struct species *surf_class;
	double area;
        u_short flags;
        byte manifold_flag;
};

/**
 * A list of regions
 */
struct region_list {
	struct region_list *next;
	struct region *reg;
};

/*
 *region counter reference, store all the info for counting reaction on regions
 */
struct reg_counter_ref {
	struct reg_counter_ref *next;  
	unsigned int counter;
	byte count_type; 	/*Three possible types:RX_STATE, TRANSITION, INIT_TRANS, MOL_TRANS*/
	byte count_method;	/* Three types:DT, SUM, CUM */
	struct region *parent;
	struct rx *state;
	struct rx *next_state;
	struct lig_transition_count **transition_count_each; /**< array of pointers to transition counter structures on region. One array element per rx mechanism in simulation. Indexed by parent_rx.rx_index */
};

struct reg_counter_ref_list {
	struct reg_counter_ref_list *next;
	struct reg_counter_ref *reg_counter_ref;
};

/*
 *counter hash table for counting on regions in each object
 */
struct counter_hash_table {
	struct counter_hash_table *next;
	char *name;
	void *value;
};


/**
 * Container data structure for all physical objects.
 */
struct object {
        struct object *next;		/**< ptr to next sibling object */
        struct object *parent;		/**< ptr to parent meta object */
        struct object *first_child;	/**< ptr to first child object */
        struct object *last_child;	/**< ptr to last child object */
        struct sym_table *sym;
        char *last_name;
        byte object_type;
        void *contents;			/**< ptr to actual physical object */
	struct lig_count_ref *lig_count_ref;
					/**< ptr to list of lig_count_ref
	                                   structures: one for each time meta-
	                                   object is referenced in count stmt */
        unsigned int num_regions;	/**< number of regions defined */
	struct region_list *regions;    /**< ptr to list of regions for 
					      this object */
	struct counter_hash_table **counter_hash_table;	/**<hash table for region counter in object*/
        int n_walls;                  /**< Total number of walls in object */
        int n_walls_actual;           /**< number of non-null walls in object */
        struct wall *walls;           /**< array of walls in object */
        struct wall **wall_p;         /**< array of ptrs to walls in object */
        int n_verts;                  /**< Total number of vertices in object */
        struct vector3 *verts;        /**< array of vertices in object */
        struct vector3 **vert_p;      /**< array of ptrs to verts in object */
        double total_area;            /**< area of object in length units */
        u_int n_tiles;                /**< number of tiles on object */
        u_int n_occupied_tiles;       /**< number of occupied tiles on object */
        struct mem_helper *edgemem;   /**< Storage for edges of object */
        struct viz_obj *viz_obj;
        int *viz_state;			/**< array of viz state values.
					   One for each element of object. */
        double t_matrix[4][4];		/**< transformation matrix for object */
};


/**
 * Doubly linked list of names.
 */
struct name_list {
        struct name_list *next;
        char *name;
        struct name_list *prev;
};

/**
 * Visualization objects.
 */
struct viz_obj {
	struct viz_obj *next;
	char *name;
        char *full_name;
	struct object *obj;
	struct viz_child *viz_child_head;
};

/**
 * Linked list of pointers to objects.
 * Used to point to child polygon or box objects to be visualized.
 */
struct viz_child {
  struct viz_child *next;
  struct object *obj;
};

/**
 * Geometric transformation data for a physical object.
 * Used to instantiate an object.
 */
struct transformation {
        struct vector3 translate;
        struct vector3 scale;
        struct vector3 rot_axis;
        double rot_angle;
};

/**
 * Molecule release pattern data.
 */


/**
 * Linked list of data to be output.
 * [\todo better description.]
 */
struct frame_data_list {
	struct frame_data_list *next;
        byte list_type;		/* data output timing type (OUTPUT_BY_TIME_LIST, etc) */
	int type;               /* visualization frame data type 
					(ALL_FRAME_DATA, etc.) */ 
	long long viz_iterationll;	/* value of the current iteration step. */
	long long n_viz_iterations;	/* number of iterations in the 
					iteration_list. */
	struct num_expr_list *iteration_list;   /* linked list of iteration 
							steps values */
	struct num_expr_list *curr_viz_iteration; /* points to the current
						 iteration in the linked list */
};


/**
 * Linked list of unique viz states.
 * required by certain viz output modes e.g. Renderman
 */
struct state_list {
  int state;
  char *name;
  struct state_list *next;
};


/**
 * A pointer to filehandle and it's real name.
 * Used for user defined file IO operations.
 */
struct file_stream {
	char *name;
	FILE *stream;
};

/**
 * Linked list of symbols.
 * Used to parse and store user defined symbols from the MDL input file.
 */
struct sym_table {
        struct sym_table *next;  /**< next symbol in symbol table*/
        unsigned short sym_type; /**< type of symbol stored -
                                   OBJ, RX, MOL, DBL, PNT ...*/
        char *name;              /**< name of symbol*/
        void *value;             /**< ptr to stored value*/
#ifdef KELP
		byte keep_alive;	/**< flag to indicate continued use of
							  this symbol table entry during computation */
		byte ref_count;		/**< number of times referenced in MDL file */
#endif
};

/**
 * Linked list of numerical expressions.
 * Used for parsing MDL input file arithmetic expressions.
 */
struct num_expr_list {
  struct num_expr_list *next;
  double value;
};
 

struct rk_mode_data
{
  int n_bins;
  int* bins;
  double *parts;
  struct vector3 *direction;
  int n_written;
};

#ifdef DEBUG
#define no_printf printf
#else
void no_printf(const char *,...);
#endif


#endif
