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


/* MCell version */
#define MCELL_VERSION "3.001" 

/* Species flags */
   /* Walls have IS_SURFACE set, molecules do not. */
   /* Grid molecules have ON_GRID set */
   /* COUNT_ENCLOSED set if you count what happens inside closed region */
   /* (otherwise only count stuff happening at the surface) */
#define ON_GRID          0x01
#define IS_SURFACE       0x02
#define NOT_FREE         0x03
#define TIME_VARY        0x04
#define CAN_MOLMOL       0x10
#define CAN_MOLGRID      0x20
#define CAN_MOLWALL      0x40
#define CAN_GRIDGRID     0x80
#define CAN_GRIDWALL     0x100
#define CANT_INITIATE    0x200
#define COUNT_CONTENTS   0x1000
#define COUNT_HITS       0x2000
#define COUNT_RXNS       0x4000
#define COUNT_ENCLOSED   0x8000
#define COUNT_SOME       0xF000

/* rxn/mol/region counter report types */
/* Do not set both WORLD and ENCLOSED flags; ENCLOSED applies only to regions */
/* First set reports a single number */
#define REPORT_NOTHING         0
#define REPORT_CONTENTS        1
#define REPORT_RXNS            2
#define REPORT_FRONT_HITS      3
#define REPORT_BACK_HITS       4
#define REPORT_FRONT_CROSSINGS 5
#define REPORT_BACK_CROSSINGS  6
/* Anything >= REPORT_MULTIPLE reports some combination of the above */
#define REPORT_MULTIPLE        7
#define REPORT_ALL_HITS        8
#define REPORT_ALL_CROSSINGS   9
/* Concentration is kind of special. */
#define REPORT_CONCENTRATION   10
#define REPORT_ELAPSED_TIME    11
/* All basic report types can be masked with this value */
#define REPORT_TYPE_MASK       0x0F
/* And finally we have some flags to say whether we're to count over */
/* the entire world or the volume enclosed by a region (set only one) */
#define REPORT_WORLD           0x20
#define REPORT_ENCLOSED        0x40
#define REPORT_TRIGGER         0x80

/* rxn/mol/region counter flags */
/* Only set one of MOL_COUNTER or RXN_COUNTER */
/* Set ENCLOSING_COUNTER if the region is closed and counts inside itself */
#define MOL_COUNTER 1
#define RXN_COUNTER 2
#define ENCLOSING_COUNTER 4
#define TRIG_COUNTER 8

#define MANIFOLD_UNCHECKED 0
#define NOT_MANIFOLD       1
#define IS_MANIFOLD        2

#define COUNT_RX_CONTENTS 1
#define COUNT_RX_ENCLOSED 2
#define COUNT_RX_SOME     3

/* Reaction flags */
  /* RX_TRANSP signifies that a reaction is between a molecule and a TRANSPARENT wall */
  /* Any value equal to or less than RX_SPECIAL refers to a special wall type */
  /* RX_BLOCKED signals a reaction that cannot take place because the grid is full */
  /* Any value equal to or less than RX_NO_RX indicates that a reaction did not take place */
  /* RX_FLIP signals that a molecule flips its orientation (crosses a wall if it's free) */
  /* RX_DESTROY signals that the molecule no longer exists (so don't try to keep using it) */
  /* RX_A_OK signals that all is OK with a reaction, proceed as normal (reflect if you're free) */
  /* RX_NO_MEM signals a memory allocation error. */
#define RX_REFLEC  -4
#define RX_TRANSP  -3
#define RX_SPECIAL -3
#define RX_BLOCKED -2
#define RX_NO_RX   -2
#define RX_FLIP    -1
#define RX_LEAST_VALID_PATHWAY 0
#define RX_DESTROY  0
#define RX_A_OK     1
#define RX_NO_MEM   3
#define RX_PATHWAY_BITS 32
#define RX_GREATEST_VALID_PATHWAY ((1LL << RX_PATHWAY_BITS)-1) 


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
/* EPS_C is the fractional difference between two values that is considered meaningful */
/* GIGANTIC is a distance that is larger than any possible simulation */
/* FOREVER is a time that cannot be reached within one simulation (too many timesteps) */
#define EPS_C 1e-12
#define GIGANTIC 1e140
#define FOREVER 1e20

/* Special rate constants used to define unusual reactions */
#define KCAT_RATE_TRANSPARENT -1.0
#define KCAT_RATE_REFLECTIVE -2.0


/* Abstract molecule flags */
/* RULES: only one of TYPE_GRID, TYPE_3D set. */
/*   ACT_NEWBIE beats ACT_INERT beats ACT_REACT */
/*   Can free up memory when nothing in IN_MASK */

/* Molecule type--grid molecule, 3D molecule, or surface molecule */
#define TYPE_GRID   0x001
#define TYPE_3D   0x002
#define TYPE_MASK 0x003

/* NEWBIE molecules get scheduled before anything else happens to them. */
/* INERT molecules don't react, REACT molecules do */
/* CHANGE molecules have had their rate constant changed */
/* DIFFUSE molecules diffuse (duh!) */
/* CLAMPED molecules diffuse for part of a timestep and don't react with surfaces */
#define ACT_DIFFUSE 0x008
#define ACT_INERT   0x010
#define ACT_REACT   0x020
#define ACT_NEWBIE  0x040
#define ACT_CHANGE  0x080
#define ACT_CLAMPED 0x1000

/* Flags telling us which linked lists the molecule appears in. */
#define IN_SCHEDULE 0x100
#define IN_SURFACE  0x200
#define IN_VOLUME   0x400
#define IN_MASK     0x700

/* Flags telling us what our counting status is */
#define COUNT_ME    0x800

/* End of abstract molecule flags. */


/* How big will we let the reaction table get? */
/* 0x400000 = 8 million */
#define MAX_RX_HASH 0x400000

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
#define SHAPE_UNDEFINED -1
#define SHAPE_SPHERICAL 0
#define SHAPE_CUBIC 1
#define SHAPE_ELLIPTIC 2
#define SHAPE_RECTANGULAR 3
#define SHAPE_SPHERICAL_SHELL 4
#define SHAPE_REGION 5
#define SHAPE_LIST 6


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

/* Negative numbers used for reaction disks */
/* Note: TARGET_OCCLUDED is assumed for any negative number not defined here */
#define TARGET_OCCLUDED    -1
#define EXD_OUT_OF_MEMORY  -2

/* Boolean set operations for releases on regions */
/* Set only one of NO_OP, UNION, INTERSECTION, SUBTRACTION */
#define REXP_NO_OP        1
#define REXP_UNION        2
#define REXP_INTERSECTION 4
#define REXP_SUBTRACTION  8
#define REXP_MASK         0x0F
#define REXP_LEFT_REGION  0x10
#define REXP_RIGHT_REGION 0x20

/* Distance in length units to search for a new site for a grid molecule */
/* after checkpointing.  Current site might be full, so a value >1 is */
/* advisable.  Being a little generous here. */
#define CHKPT_GRID_TOLERANCE 2.0


/* Constants for garbage collection of defunct molecules */
/* (those that were consumed when hit by another molecule) */
#define MIN_DEFUNCT_FOR_GC 1024
#define MAX_DEFUNCT_FRAC 0.2

/* Constants for notification levels */
#define NOTIFY_NONE 0
#define NOTIFY_BRIEF 1
#define NOTIFY_FULL 2
#define NOTIFY_CUSTOM 3

/* Constants for warning levels */
#define WARN_COPE 0
#define WARN_WARN 1
#define WARN_ERROR 2

/* Stuff to set surface diffusion behavior */
#define SURFACE_DIFFUSION_RETRIES 10

/* Flags for different types of file output */
#define FILE_UNDEFINED 0
#define FILE_OVERWRITE 1
#define FILE_SUBSTITUTE 2
#define FILE_APPEND 3
#define FILE_APPEND_HEADER 4
#define FILE_CREATE 5

/* Flags for output expressions */

#define OEXPR_TYPE_UNDEF 0x0
#define OEXPR_TYPE_INT 0x1
#define OEXPR_TYPE_DBL 0x2
#define OEXPR_TYPE_TRIG 0x3
#define OEXPR_TYPE_MASK 0x7
#define OEXPR_TYPE_CONST 0x8

#define OEXPR_LEFT_INT 0x10
#define OEXPR_LEFT_DBL 0x20
#define OEXPR_LEFT_TRIG 0x30
#define OEXPR_LEFT_REQUEST 0x40
#define OEXPR_LEFT_OEXPR 0x50
#define OEXPR_LEFT_MASK 0x70
#define OEXPR_LEFT_CONST 0x80

#define OEXPR_RIGHT_INT 0x100
#define OEXPR_RIGHT_DBL 0x200
#define OEXPR_RIGHT_TRIG 0x300
#define OEXPR_RIGHT_REQUEST 0x400
#define OEXPR_RIGHT_OEXPR 0x500
#define OEXPR_RIGHT_MASK 0x700
#define OEXPR_RIGHT_CONST 0x800



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

/* Viz output options */
#define VIZ_ALL_MOLECULES 0x01
#define VIZ_MOLECULES_STATES 0x02
#define VIZ_SURFACE_STATES 0x04

/*******************************************************/
/**  Old constants copied from MCell2, may be broken  **/
/*******************************************************/

/* Parser parameters.  Probably need to be revisited. */
/* size of symbol hash table 0x100000 = 1M */
#define SYM_HASHSIZE 0x100000

/* mask for symbol table hash */
#define SYM_HASHMASK 0x0FFFFF

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

/* Used only for TRIGGER statements */
#define TRIG_STRUCT 16

                                                                                
/* Object types */
#define META_OBJ 0
#define BOX_OBJ 1
#define POLY_OBJ 2
#define REL_SITE_OBJ 3
#define VOXEL_OBJ 4
                                                                                

/* Box sides */
#define TP 0
#define BOT 2
#define FRNT 4
#define BCK 6
#define LFT 8
#define RT 10
#define ALL_SIDES INT_MAX
                   
                                                             
/* Viz state values */
#define EXCLUDE_OBJ INT_MIN /*object is not visualized */
#define INCLUDE_OBJ INT_MAX /*object is visualized but state value is not set*/


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
#define TRIGGER_VAL 3


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
#define DREAMM_V3_GROUPED_MODE 3
#define RK_MODE 4
#define ASCII_MODE 5


/* Visualization frame data types. */
/* Used to select type of data to include in viz output files */
/* Will probably change significantly when we redesign DReAMM output format */
#define ALL_FRAME_DATA 0
#define EFF_POS 1
#define EFF_STATES 2
#define MOL_POS 3
#define MOL_ORIENT 4
#define MOL_STATES 5
#define SURF_POS 6 
#define SURF_STATES 7
#define MESH_GEOMETRY 8
#define REG_DATA 9
#define ALL_MOL_DATA  10
#define ALL_MESH_DATA 11

/* release number methods */
#define CONSTNUM 0
#define GAUSSNUM 1
#define VOLNUM 2
#define CCNNUM 3



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
  u_int chkpt_species_id;        /* Unique ID for this species from the 
                                   checkpoint file */
  u_int hashval;                /* Hash value (may be nonunique) */
  struct sym_table *sym;        /* Symbol table entry (name) */
  struct eff_dat *eff_dat_head; /* if IS_SURFACE this points to head of
                                   effector data list associated with 
                                   surface class */
  
  u_int population;             /* How many of this species exist? */
  
  double D;                     /* Diffusion constant */
  double D_ref;                 /* Reference diffusion constant */
  double radius;                /* Molecular radius */
  double area;                  /* Surface area consumed */
  double space_step;            /* Characteristic step length */
  double time_step;             /* Minimum (maximum?) sensible timestep */
/*short charge;*/               /* Electric charge. */
  u_short flags;                /* Free?  Membrane bound?  Membrane? */
  
  long long n_deceased;         /* Total number that have been destroyed. */
  double cum_lifetime;          /* Timesteps lived by destroyed molecules */
 
  int viz_state;                /* Visualization state for output */
  int region_viz_value;         /* Visualization state for surface class 
                                   for output */
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
  
  long long n_occurred;      /* How many times has this reaction occurred? */
  double n_skipped;          /* How many reactions were skipped due to probability overflow? */

  struct t_func *prob_t;     /* List of probabilities changing over time */
  
  struct pathway *pathway_head; /* list of pathways built at parse-time */
  struct pathway_info *info; /* Counts and names for each pathway */
};


/* Sets of reactions grouped by...? */
struct rxn_group
{
/* Someone else gets to fill this in. */
};


/* Named reaction.  Pathway may be NULL during running of simulation; it
is only needed during parsing. */
struct rxn_pathname {
  struct sym_table *sym;
  u_int hashval;
  u_int path_num;
  struct rxn *rx;
};


/* Parse-time structure for reaction pathways */
/* Everything except pathname can be deallocated after prepare_reactions */
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
};

/* Parse-time structure for products of reaction pathways */
struct product {
  struct product *next;
  struct species *prod;          /* Molecule type to be created */
  short orientation;             /* Orientation to place molecule */
};

/* Run-time info for each pathway */
/* Always do basic counts--do more sophisticated stuff if pathname!=NULL */
struct pathway_info
{
  double count;
  short count_flags;
  struct rxn_pathname *pathname;
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
  double birthday;                 /* Time at which this particle was born */
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
  double birthday;
  
  struct vector3 pos;             /* Position in space */
  struct subvolume *subvol;       /* Partition we are in */
  
  struct cmprt_data *curr_cmprt;  /* Compartment we are in (for counting) */
  
  struct wall *previous_wall;     /* Wall we were released from */
  int index;                      /* Index on that wall (don't rebind) */
  
  struct molecule *next_v;        /* Next molecule in this subvolume */
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
  double birthday;
  
  int grid_index;              /* Which gridpoint do we occupy? */
  short orient;                /* Which way do we point? */
  struct surface_grid *grid;   /* Our grid (which tells us our surface) */
  struct vector2 s_pos;        /* Where are we in surface coordinates? */
  
  struct grid_molecule *prev_g; /* Doubly linked list of molecules */
  struct grid_molecule *next_g; /* at this grid index */
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
  
  struct surface_grid *effectors; /* Grid of effectors for this wall */
  
  int viz_state;                  /* For display purposes--is short enough? */
  u_short flags;                  /* Flags for whether we need to count */

  struct object *parent_object;   /* The object we are a part of */
  struct storage *birthplace;     /* Where we live in memory */
  
  struct region_list *counting_regions; /* Counted-on regions containing this wall */
};


/* Linked list of walls (for subvolumes) */
struct wall_list
{
  struct wall_list *next;   /* The next entry in the list */
  
  struct wall *this_wall;        /* The wall in this entry */
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
  double n_rxn_at;                  /* # rxn occurrance on surface */
  double n_rxn_enclosed;            /* # rxn occurrance inside closed region */
};


struct move_counter_data
{
  double front_hits;               /* # hits on front of region (normal up) */
  double back_hits;                /* # hits on back of region */
  double front_to_back;            /* # crossings from front to back */
  double back_to_front;            /* # crossings from back to front */
  double scaled_hits;              /* To determine integrated concentration */
  int n_at;                        /* # molecules on region surface */
  int n_enclosed;                  /* # molecules inside closed region */
};

struct trig_counter_data
{
  double t_event;                    /* Event time (exact) */
  struct vector3 loc;                /* Real position of event */
  struct trigger_request *listeners; /* Places waiting to be notified */
};


union counter_data
{
  struct rxn_counter_data rx;
  struct move_counter_data move;
  struct trig_counter_data trig;
};


/* Struct to count rxns or molecules within regions (where "within" includes */
/* on the inside of a fully closed surface, for 3D molecules) */
struct counter
{
  struct counter *next;
  byte counter_type;               /* MOL_COUNTER or RXN_COUNTER, plus flags */
  struct region *reg_type;         /* Region we are counting on */
  void *target;                    /* Mol or rxn pathname we're counting */
  union counter_data data;         /* data for the count
                                      reference data.move for move counter
                                      reference data.rx for rxn counter 
                                      reference data.trig for trigger */
};


struct trigger_request
{
  struct trigger_request *next; /* Next request */
  struct output_request *ear;   /* Who wants to hear about the trigger */
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
  byte releases_on_regions_flag; /* Triggers special release site initialization */
  
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
  
  struct output_block *output_block_head;
  struct output_request *output_request_head;
  struct mem_helper *oexpr_mem;
  struct mem_helper *outp_request_mem;
  struct mem_helper *counter_mem;
  struct mem_helper *trig_request_mem;
  
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
  double *r_step_surface;
  double r_num_directions;
  double sim_elapsed_time;
  double chkpt_elapsed_real_time;    /** elapsed simulation time (in sec) for new 
                                    checkpoint */
  double chkpt_elapsed_real_time_start;  /**< start of the simulation time (in sec)
                                        for new checkpoint */
  double current_real_time;          /**< current simulation time in seconds */
  double current_start_real_time;    /**< simulation start time (in seconds) */
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
  u_int init_seed;              /**<  initial seed value for random function */
  u_int chkpt_byte_order_mismatch;   /**< flag that defines whether mismatch
                                      in byte order exists between machines
                                      that writes and reads checkpoint file.*/

  long long it_time;
  double elapsed_time;    /**< number of iterations after simulation starts */
  long long start_time;  /**< starting iteration number for the current run */
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
  char *file_prefix_name;
  char *mcell_version; 
  u_short viz_output_flag; /*takes  VIZ_ALL_MOLECULES  */
 
  /* Optional stuff */
  int use_expanded_list;
  int randomize_gmol_pos;
  double vacancy_search_dist2;
  byte surface_reversibility;
  byte volume_reversibility;

  /* MCell startup command line arguments */
  byte info_opt;
  u_int seed_seq;            /**< index in the seed_array */
  long long iterations;
  char *log_file_name;
  FILE *log_file;
  FILE *err_file;
  u_int log_freq;
  /* flag set to 0 if CHECKPOINT_INFILE can be opened for reading,
     otherwise - if there is no such file or it can't be read it is set to 1 */
  u_int chkpt_init;
  /* flag set to 1 if there are any CHECKPOINT statements in "mdl" file */
  u_int chkpt_flag;
  /* value of the CHECKPOINT_ITERATIONS keyword */
  long long chkpt_iterations;
  u_int chkpt_seq_num;
  char *chkpt_infile;
  char *chkpt_outfile;
  /* Handle for the file named "chkpt_infile" */
  FILE *chkpt_infs;
  FILE *chkpt_outfs;
  FILE *chkpt_signal_file_tmp;
  char *mdl_infile_name;
  char *curr_file;
  
  /* Notification/warning/output stuff */
  struct notifications *notify;
  
  /* Concentration clamp at surfaces */
  struct ccn_clamp_data *clamp_list;
  
  /* Nifty pointers for debugging */
  struct output_request *watch_orq;
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
  double t_matrix[4][4];                /**< transformation matrix */
  int train_counter;			/**< counts executed trains */
  double train_high_time;		/**< time of the train's start */
};


struct release_site_obj {
	struct vector3 *location;	/**< location of release site */
	struct species *mol_type;	/**< species to be released */
	byte release_number_method;
	byte release_shape;
	short orientation;
	int release_number;
	int mean_number;
	double mean_diameter;
	double concentration;
        double standard_deviation;
	struct vector3 *diameter;
	struct release_region_data *region_data;
	struct release_single_molecule *mol_list;

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

/* Extended data for complex releases on regions */
/* Not all fields are used for all release types */
struct release_region_data
{
  struct vector3 llf;   /* One corner of bounding box for release area */
  struct vector3 urb;   /* Opposite corner */

  int n_walls_included;  /* How many walls total */
  double *cum_area_list; /* Cumulative area of all walls */
  int *wall_index;       /* Indices of each wall (by object) */
  int *obj_index;        /* Indices for objects (in world list) */
  
  int n_objects;                  /* How many objects are there total */
  struct object **owners;         /* Array of pointers to each object */
  struct bit_array **in_release;  /* Bit array saying which walls are in release for each object */
  int *walls_per_obj;             /* Number of walls in release for each object */

  struct object *self;            /* A pointer to our own release object */
  struct release_evaluator *expression;  /* A set construction expression combining regions to form this release site */
};

/* Data structure used to store LIST releases */
struct release_single_molecule
{
  struct release_single_molecule *next;  /* Next in the list */
  struct species *mol_type;              /* Species to release */
  struct vector3 loc;                    /* Position to release it */
  short orient;                          /* Orientation (for 2D species) */
};

/* Data structure used to hold combinations of regions */
struct release_evaluator
{
  byte op;       /* The operation used (values are #def'ed as REXP_... */
  void *left;    /* The left side of the expression--another evaluator or a region object depending on bitmask of op */
  void *right;   /* The right side--same thing */
};

/* Holds information about a box with rectangular patches on it. */
struct subdivided_box
{
  int nx; /* number of subdivisions including box corners in X-direction */
  int ny; /* number of subdivisions including box corners in Y-direction */
  int nz; /* number of subdivisions including box corners in Z-direction */
  double *x;  /* array of X-coordinates of subdivisions */
  double *y;  /* array of Y-coordinates of subdivisions */
  double *z;  /* array of Z-coordinates of subdivisions */
};

/* Holds information about what we want dumped to the screen */
struct notifications
{
  /* Informational stuff */
  byte progress_report;
  byte diffusion_constants;
  byte reaction_probabilities;
  double reaction_prob_notify;
  byte partition_location;
  byte box_triangulation;
  byte custom_iterations;
  long long custom_iteration_value;
  byte release_events;
  byte file_writes;
  byte final_summary;
  
  /* Warning stuff */
  byte neg_diffusion;
  byte neg_reaction;
  byte high_reaction_prob;
  double reaction_prob_warn;
  byte close_partitions;
  byte degenerate_polys;
  byte overwritten_file;
  byte short_lifetime;
  long long short_lifetime_value;
  byte missed_reactions;
  double missed_reaction_value;
  byte missed_surf_orient;
  byte useless_vol_orient;
};

struct ccn_clamp_data
{
  struct ccn_clamp_data *next;     /* The next concentration clamp */
  struct species *surf_class;   /* Which surface class clamps? */
  struct species *mol;             /* Which molecule does it clamp? */
  double concentration;            /* At which concentration? */
  short orient;                    /* On which side? */
  struct object *objp;             /* Which object are we clamping? */
  struct bit_array *sides;         /* Which sides on that object? */
  int n_sides;                     /* How many are set */
  int *side_idx;                   /* Indices of the sides that are set */
  double *cum_area;                /* Cumulative area of all the sides */
  double scaling_factor;           /* Used to predict #mols/timestep */
  struct ccn_clamp_data *next_mol; /* Next molecule for this class */
  struct ccn_clamp_data *next_obj; /* Next object for this class */
};

struct output_block
{
  struct output_block *next;            /* Next in world or scheduler */
  double t;                             /* Scheduled time to update counters */
  
  byte timer_type;                      /* OUTPUT_BY_STEP, ...BY_TIME_LIST, ...BY_ITERATION_LIST */
  
  double step_time;                     /* Output interval (seconds) */
  struct num_expr_list *time_list_head; /* List of output times/iteration numbers */
  struct num_expr_list *time_now;       /* Current entry in list */
  
  u_int buffersize;                     /* Size of output buffer */
  u_int buf_index;                      /* Index into buffer (for non-triggers) */
  
  double *time_array;                   /* Array of output times (for non-triggers) */
  
  struct output_set *data_set_head;     /* Linked list of data sets (separate files) */
};

struct output_set
{
  struct output_set *next;             /* Next data set */
  struct output_block *block;          /* Which block do we belong to? */
  char *outfile_name;                  /* Filename */
  int file_flags;                      /* Tells us how to handle existing files */
  u_int chunk_count;                    /* Number of chunks processed */  
  char *header_comment;                /* Comment character(s) for header */
  struct output_column *column_head;   /* Data for one output column */
};

struct output_column
{
  struct output_column *next;       /* Next column */
  struct output_set *set;           /* Which set do we belong to? */
  byte data_type;                   /* INT, DBL, TRIG_STRUCT */
  double initial_value;             /* To continue existing cumulative counts--not implemented yet--and keep track of triggered data */
  void *buffer;                     /* Output buffer array (cast based on column_type) */
  struct output_expression *expr;   /* Evaluate this to calculate our value (NULL if trigger) */
};

struct output_expression
{
  struct output_column *column;     /* Which column are we going to? */
  int expr_flags;                   /* What kinds of things are to the left and right? */
  struct output_expression *up;     /* Parent output expression */
  void *left;                       /* Item on the left */
  void *right;                      /* Item on the right */
  char oper;                        /* Operation to apply to items */
  double value;                     /* Resulting value from operation */
  char *title;                      /* String describing what we've got */
};

struct output_request
{
  struct output_request *next;          /* Next request in global list */
  struct output_expression *requester;  /* Expression in which we appear */
  struct sym_table *count_target;       /* Mol/rxn we're supposed to count */
  struct sym_table *count_location;     /* Place we're supposed to count it */
  byte report_type;                     /* Flags telling us how to count */
};

struct output_trigger_data
{
  double t_iteration;          /* Time of the iteration of triggering event */
  double t_delta;              /* Offset of event time from iteration time */
  struct vector3 loc;          /* Position of event */
  int how_many;                /* Number of events */
  char *name;                  /* Name to give event */
};

/******************************************************************/
/**  Everything below this line has been copied from MCell 2.69  **/
/******************************************************************/

#if 0
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
	int file_flags;            /* Append, overwrite, etc. */
	int first_write;           /* Only write header on first write */
	char *header_comment;             /**< comment character(s) for header */
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
	                              indexed by either TIME_STAMP_VAL or INDEX_VAL
				      or is output when triggered (TRIGGER_VAL) */
	byte data_type;             /**< type of data to track:
                                      EXPR INT DBL TRIG_STRUCT*/
	u_int n_data;               /** buffer size */
	void *temp_data;            /**< ptr to intermediate data
                                         specified by type*/
	void *final_data;           /**< ptr to final outputable data
                                         specified by type*/
	struct output_evaluator *operand1;
	struct output_evaluator *operand2;
	char oper;
	char *column_title; /* Column title in output file */
};
#endif

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
 * A voxel list object, part of a volume.
 */
struct voxel_object {
	struct lig_count_ref *lig_count_ref;
					/**< ptr to list of lig_count_ref
	                                   structures: one for each time voxel 
	                                   object is referenced in count stmt */
	struct viz_state_ref *viz_state_ref;
					/**< ptr to list of viz_state_ref
	                                   structures: one for each time voxel 
	                                   object is referenced in
					   STATE_VALUES block */
        struct ordered_voxel *voxel_data; /**< pointer to data structure
                                                holding voxel vertices etc... */
	int n_voxels;			/**< Number of voxels in
                                             voxel object */
        int n_verts;                    /**< Number of vertices in
                                             voxel object */
	byte fully_closed;		/**< flag indicating closure of object */
};

/**
 * A general ordered polyhedron consisting from tetrahedrons (voxels). 
 * That is, the vertices of each polygonal face are ordered according to
 * the right hand rule.
 */
struct ordered_voxel {
	struct vector3 *vertex;         /**< Array of tetrahedron vertices */
	struct tet_element_data *element; /**< Array tet_element_data
                                              data structures */
	struct tet_neighbors_data *neighbor; /**< Array tet_neighbors_data
                                              data structures */
	int n_verts;                  /**< Number of vertices in polyhedron */
	int n_voxels;                 /**< Number of voxels in polyhedron */
};

/**
 * Data structure used to build one tetrahedron.
 * This data structure is used to store the data from the MDL file
 * and to contruct each tetrahedron of a voxel object.
 */
struct tet_element_data {
        int vertex_index[4];              /**< Array of vertex indices forming a
                                           tetrahedron. */
	int n_verts;                    /**< Number of vertices in tetrahedron (always 4). */
};

/**
 * This data structure is used to store the data about neighbors
 * of each tetrahedron of a voxel object.
 */
struct tet_neighbors_data {
        int neighbors_index[4];        /**< Array of indices pointing 
                                           to the neighbors of tetrahedron. */
	int n_neighbors;               /**< Number of neighbors of tetrahedron (always 4). */
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
        u_int begin; /*first number in the list of numbers */
        u_int end;   /* last number in the list of numbers */
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
        struct vector3 *bbox;  /* Vector of length 2, may be null 
                                 - NOT IMPLEMENTED yet*/
        int region_viz_value; /* used for visualization */
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
struct viz_obj 
{
	struct viz_obj *next;
	char *name;            /* name taken from OBJECT_FILE_PREFIXES
                  or FILENAME_PREFIXES or FILENAME assignment  */
        char *full_name;       /* full name of the object, like A.B.C */
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
 * Linked list of symbols.
 * Used to parse and store user defined symbols having wildcards
   from the MDL input file.
 */
struct sym_table_list {
  struct sym_table_list *next;
  struct sym_table *node;

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
