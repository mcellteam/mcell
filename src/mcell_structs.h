#ifndef MCELL_STRUCTS
#define MCELL_STRUCTS

#include <limits.h>
#include <sys/types.h>
#include <stdio.h>

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
   /* IS_ACTIVE is set if this molecule can do anything on its own */
#define ON_SURFACE  0x01
#define ON_GRID     0x02
#define IS_SURFACE  0x04
#define IS_ACTIVE   0x08
#define CAN_MOLMOL  0x10
#define CAN_MOLGRID 0x20
#define CAN_MOLSURF 0x40
#define CAN_MOLWALL 0x80


/* Reaction flags */
#if 0
#define RX_DESTROY   0x001
#define RX_FLIP      0x002
#define RX_PROD      0x004
#define RX_REFL      0x008
#define RX_2DESTROY  0x010
#define RX_2FLIP     0x020
#define RX_2PROD     0x040
#endif


/* Flags for BSP trees to determine whether something is a node or a branch */
/* Will either have BRANCH_XN through _ZP, or _L, _R, _X, _Y, _Z. */
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

#define X_NEG_BIT 0x01
#define X_POS_BIT 0x02
#define Y_NEG_BIT 0x04
#define Y_POS_BIT 0x08
#define Z_NEG_BIT 0x10
#define Z_POS_BIT 0x20


/* Collision types for rays striking surfaces */
#define COLLIDE_MISS    0
#define COLLIDE_FRONT   1
#define COLLIDE_BACK    2
#define COLLIDE_REDO    -1

#define COLLIDE_MOL_M   3
#define COLLIDE_MOL_SP  4
#define COLLIDE_MOL_SN  5

#define COLLIDE_SV_NX   6
#define COLLIDE_SV_PX   7
#define COLLIDE_SV_NY   8
#define COLLIDE_SV_PY   9
#define COLLIDE_SV_NZ   10
#define COLLIDE_SV_PZ   11

#define COLLIDE_WALL    0x10
#define COLLIDE_MOL     0x20
#define COLLIDE_SUBVOL  0x40
#define COLLIDE_MASK    0x0F


/* Types for things we can hit */
#define VOL_COLLISION    1
#define WALL_COLLISION   2   
#define MOL_COLLISION    3


/* Flags for edges. */
#define EDGE_BARE     0
#define EDGE_SHARED   1
#define EDGE_ROTONLY  2
#define EDGE_TRANSROT 3


/* Size constants */
#define EPS_C 1e-12
#define GIGANTIC 1e140


/* Abstract molecule flags */
/* RULES: only one of TYPE_GRID, TYPE_3D, TYPE_SURF set. */
/*   Only TYPE_3D and TYPE_SURF can ACT_DIFFUSE */
/*   ACT_NEWBIE beats ACT_INERT beats ACT_REACT beats ACT_BORED */
/*   Can free up memory when nothing in IN_MASK */

#define TYPE_GRID 0x001
#define TYPE_3D   0x002
#define TYPE_SURF 0x004
#define TYPE_MASK 0x007

#define ACT_DIFFUSE 0x008
#define ACT_INERT   0x010
#define ACT_REACT   0x020
#define ACT_NEWBIE  0x040

#define IN_SCHEDULE 0x100
#define IN_SURFACE  0x200
#define IN_VOLUME   0x400
#define IN_MASK     0x700


/* How big will we let the reaction table get? */
/* 0x100000 = 2 million */
#define MAX_RX_HASH 0x100000


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

/*********************************************************/
/**  Constants used in MCell3 brought over from MCell2  **/
/*********************************************************/


/* Generic numerical constants */
#define EPSILON 1e-14
                                                                                
#define R_UINT_MAX 2.3283064365386963e-10
#define MY_PI 3.14159265358979323846
#define N_AV 6.0221415e23
#define ROUND_UP 0.5
                                                                                
                                                                                                                                                                
/* Polygon list types */
#define BOX_POLY 0  
#define ORDERED_POLY 1  
#define UNORDERED_POLY 2  
   /* BOX_POLY is a right parallelepiped polyhedron */
   /* ORDERED_POLY is a SGI Inventor-style polyhedron with shared vertices */
   /* UNORDERED_POLY is a free-form collection of distinct polygons */


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
#define EFFDENS 0
#define EFFNUM 1
   /* Place either a certain density or an exact number of effectors */


/*******************************************************/
/**  Old constants copied from MCell2, may be broken  **/
/*******************************************************/

/* Parser parameters.  Probably need to be revisited. */
#define HASHSIZE 128
#define HASHMASK 0x7f
#define COUNTER_HASH 16 
#define COUNTER_HASHMASK 0xf
   /*LWW 6/13/03 Hashtable size for region counters in each OBJECT*/

#define PATHWAYSIZE 64
#define RXSIZE 2048
#define ARGSIZE 255
#define NUM_ADD_EFFECTORS 1024
#define MAX_INCLUDE_DEPTH 16
#define COUNTBUFFERSIZE 10000
                                                                                

/* Data types to be stored in sym_table, probably broken! */
#define RX 1
#define MOL 2
#define PNT 3
#define CMP 4
#define POLY 5
#define RSITE 6
#define OBJ 7
#define RPAT 8
#define REG 9
#define INT 10
#define DBL 11
#define STR 12
#define ARRAY 13
#define FSTRM 14
#define EXPR 15
#define TMP 16

                                                                                
/* Object types, probably okay. */
#define META_OBJ 0
#define BOX_OBJ 1
#define POLY_OBJ 2
#define REL_SITE_OBJ 3
                                                                                

/* Box sides.  This is a weird way to do it.  Why not bitmasks? */
#define TP 0
#define BOT 2
#define FRNT 4
#define BCK 6
#define LFT 8
#define RT 10
#define ALL_SIDES INT_MAX
                   
                                                             
/* Viz state values */
#define EXCLUDE_OBJ INT_MIN


/* Count list specifications.  INIT stuff is probably broken. */
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


/* Count list value types.  Probably okay. */
#define UNKNOWN 0
#define TIME_STAMP_VAL 1
#define INDEX_VAL 2


/* Reaction data output type */
#define FRAME_DATA 0
#define FREQ_DATA 1
   /* Use either frame data or frequency to control when to output */


/* Output list type */
#define FRAME_NUMBER 0
#define REAL_TIME 1


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
#define IRIT_MODE 2
#define RADIANCE_MODE 3
#define RAYSHADE_MODE 4
#define RENDERMAN_MODE 5
#define POVRAY_MODE 6
#define MCELL_MODE 7

/* Visualization frame data types. */
#define ALL_FRAME_DATA 0
#define EFF_POS 1
#define EFF_STATES 2
#define MOL_POS 3
#define MOL_STATES 4
#define SURF_POS 5
#define SURF_STATES 6


/* event types for release event queue */
#define TRAIN_HIGH_EVENT 0
#define TRAIN_LOW_EVENT 1
#define RELEASE_EVENT 2
                                                                                
/* release number methods */
#define CONSTNUM 0
#define GAUSSNUM 1
#define VOLNUM 2


/* Stimulus motion types.  What is this? */
#define FXD 0



/**********************************************/
/**  New/reworked structures used in MCell3  **/
/**********************************************/

typedef unsigned char byte;

#ifdef NOINCLUDE_SYS_TYPES
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
/*double time_step;*/           /* Minimum (maximum?) sensible timestep */
  short charge;                 /* Electric charge. */
  u_short flags;                /* Free?  Membrane bound?  Membrane? */
  
  int viz_state;                /* Visualization state for output */
  byte checked;                 /* Bread crumb for graph traversal */
};


/* All pathways leading away from a given intermediate */
struct rxn
{
  struct rxn *next;          /* Next reaction with these reactants */
  struct sym_table *sym;     /* ptr to symbol table entry for this rxn */
  
  u_int n_reactants;         /* How many reactants? (At least 1.) */
  u_int n_pathways;          /* How many pathways lead away? */
  u_int *product_idx;        /* Index of 1st player for products of pathway */
  double *cum_rates;         /* Cumulative rates for (entering) all pathways */
  double *cat_rates;         /* Rate of leaving all pathways (<=0.0 is instant) */
  
  struct species **players;  /* Identities of reactants/products */
  short *geometries;         /* Geometries of reactants/products */

  int n_rate_t_rxns;         /* How many pathways have varying rates? */
  int *rate_t_rxn_map;       /* Indices of pathways with varying rates */
  struct t_func *rate_t;     /* Rate over time for each varying pathway */
  struct t_func *jump_t;     /* Summary of transition times */
   
  u_int last_update;         /* When did we last update rates/counts? */
  
  u_int *rxn_count_dt;       /* How many times this timestep? */
  u_int *rxn_count_cum;      /* How many times ever? */
  struct pathway *pathway_head; /* list of pathways built at parse-time */
};


/* Sets of reactions grouped by...? */
struct rxn_group
{
/* Someone else gets to fill this in. */
};

/* Parse-time structure for reaction pathways */
struct pathway {
  struct pathway *next;
  struct species *reactant1;
  struct species *reactant2;
  struct species *reactant3;
  double km;
  double kcat;
  short orientation1;
  short orientation2;
  short orientation3;
  struct product *product_head;
};


/* Parse-time structure for products of reaction pathways */
struct product {
  struct product *next;
  struct species *prod;
  short orientation;
};


/* Piecewise constant function for time-varying reaction rates */
struct t_func
{
  int index;        /* Which constant part are we in? */
  int n;            /* How many pieces do we have? */
  u_int *time;      /* When are we done (units=timesteps)? */
  double *value;    /* What is the value now? */
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
  double path_length;
  int collisions;
  
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
  
  int viz_state;                  /* For display purposes */

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
  struct vector3 loc;
  struct cmprt_data **owners;
  int n_owners;
};


struct vertex_tree
{
  struct vertex_tree *next;
  struct vertex_tree *above;
  struct vertex_tree *below;
  
  struct vector3 loc;
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
  
  struct wall *wall_head;
  int wall_count;
  struct vertex_tree *vert_head;
  int vert_count;
  
  struct schedule_helper *timer;
  double current_time;
  double max_timestep;
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
  
  struct storage *mem;         /* Local storage */
};


/* Binary space partitioning tree for subvolume connections */
struct bsp_tree
{
  void *left;        /* The tree below the partition */
  void *right;       /* The tree above the partition */
  short partition;   /* The index of the partition */
  short flags;       /* Coordinate that is split, plus terminal node flags */
};


struct counter
{
  struct counter *next;
  int wall_id;
  int mol_id;
  int crossings;
  int impacts;
};


/* All data about the world */
struct volume
{
/*  struct vector3 corner[8];*/  /* Corners of the world */
#if 0
  struct vector3 llf;           /* left lower front corner of world */
  struct vector3 urb;           /* upper right back corner of world */
#endif
  
  int n_parts;                  /* Number of coarse partition boundaries */
  double *x_partitions;         /* Coarse X partition boundaries */
  double *y_partitions;         /* Coarse Y partition boundaries */
  double *z_partitions;         /* Coarse Z partition boundaries */
  
  int n_fineparts;        /* Number of fine partition boundaries (multiple of n_axis_partitions) */
  double *x_fineparts;           /* Fine X partition boundaries */
  double *y_fineparts;           /* Fine Y partition boundaries */
  double *z_fineparts;           /* Fine Z partition boundaries */
  
  int n_waypoints;              /* How many of these = (n_axis_p-3)^3 */
  struct waypoint *waypoints;   /* Contains compartment information */
  
  int n_subvols;                /* How many coarse subvolumes? */
  struct subvolume *subvol;     /* Array containing all subvolumes */
   
  int binning;                  /* How many real partitions per one in lookup? */
  struct subvolume **lookup;     /* 3D lookup array pointing at subvolumes */
  
  int n_walls;                  /* Total number of walls */
  int n_verts;
  
  int hashsize;                 /* How many entries in our hash table? */
  int n_reactions;              /* How many reactions are there, total? */
  struct rxn **reaction_hash;   /* A hash table of all reactions. */
  
  int collide_hashmask;         /* Mask for looking up collision hash table */
  struct counter **collide_hash;/* Collision hash table */
  
  int n_species;                /* How many different species? */
  struct species **species_list; /* Array of all species. */
  
  struct schedule_helper *releaser;
  
  struct mem_helper *storage_mem;      /* Storage for storage list */
  struct storage_list *storage_head;   /* Linked list of all storage */
  
  double speed_limit;           /* How far can the fastest particle get in one timestep? */
  double diffusion_number;
  double diffusion_cumsteps;
  
  /* Simulation initialization parameters  */
  struct sym_table **main_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct release_pattern *default_release_pattern;
  struct release_event_queue *release_event_queue_head;
  struct count_list *count_list;
  struct count_list *count_zero;
  struct output_list *output_list;
  struct viz_obj *viz_obj_head;
  struct frame_data_list *frame_data_head;
  struct species *g_mol;
  struct species *g_surf;
  double time_unit;
  double time_step_max;
  double length_unit;
  double effector_grid_density;
  double rx_radius_3d;
  double *r_step;
  double *d_step;
  double *factorial_r;
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
  double diffusion_steps;
  double sim_elapse_time;
  struct vector3 bb_min;
  struct vector3 bb_max;
  u_int tot_mols;
  u_int seed;
  u_int init_seed;
  u_int it_time;
  u_int start_time;
  u_int n_release_events;
  u_int radial_directions;
  u_int radial_subdivisions;
  u_int num_directions;
  int directions_mask;
  int fully_random;
  int procnum;
  int viz_mode;
  byte voxel_image_mode;
  byte voxel_volume_mode;
  char *molecule_prefix_name;

  /* MCell startup command line arguments */
  byte info_opt;
  u_int seed_seq;
  u_int iterations;
  char *log_file_name;
  FILE *log_file;
  u_int log_freq;
  u_int chkpt_init;
  u_int chkpt_flag;
  u_int chkpt_iterations;
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


/* Temporary data structure to store information about collisions. */
struct collision
{
  struct collision *next;
  double t;
  
  void *target;
  struct rxn *intermediate;
  struct vector3 loc;
  int what;
};


struct release_event_queue {
  struct release_event_queue *next;
  double event_time;

  struct release_site_obj *release_site;
  struct vector3 location;
  byte event_type;
  int event_counter;
  double train_high_time;
  int index;                   /**< unique index of this release_event */
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


struct release_pattern {
        struct sym_table *sym;
	double delay;
	double release_interval;
	double train_interval;
	double train_duration;
	int number_of_trains;
};

/******************************************************************/
/**  Everything below this line has been copied from MCell 2.69  **/
/******************************************************************/


/**
 * Node data structure. 
 * \todo why only left and right walls?
 */
struct node_dat {
	struct vector3 corner[8];
	struct wall *right_wall;
	struct wall *left_wall;
	byte left_subvol;
	byte right_subvol;
	int left_node;
	int right_node;
};

/**
 * Linked list for all counter data being output.
 * Modified by Lin-Wei 5/21/02
 */
struct output_list {
	struct output_list *next;   /**< next item in output list*/
	byte out_type;
	unsigned int id;		/**<unique id number for each REACTION_OUTPUT_DATA command*/
	int counter;
	int freq;
	int n_output;
	int index;
	struct counter_info *counter_info;
	struct reaction_list *reaction_list;
};

struct counter_info {
	struct counter_info *next;  
	char outfile_name[1024];              /**< name of variable being tracked*/
	struct count_list *count_list;
};

/**
 * Reaction output data list
 *
 */
struct reaction_list {
  struct reaction_list *next;
  byte list_type;
  int n_reac_iterations;
  int reac_iteration;
  int *array;
  struct num_expr_list *iteration_list;
  struct num_expr_list *curr_reac_iteration;
};

/**
 * Linked list of counters.
 */
struct count_list {
	struct count_list *next;    /**< next item in count list*/
	int n_output;			/**< total output time steps.*/
	int freq;
	unsigned int frame_index;
	byte reset_flag;            /**< reset temp_data to 0 on each iteration?*/
	byte update_flag;           /**< counter update necessary?*/
	byte data_type;             /**< type of data to track:
	                              MOL RX CMP EXPR INT DBL*/
	byte index_type;            /**< flag indicating final_data is to be
	                              indexed by either
	                              TIME_STAMP_VAL or INDEX_VAL*/
	int n_data;
	void *temp_data;            /**< ptr to data specified by type*/
	void *final_data;           /**< ptr to data specified by type*/
	struct count_list *operand1;
	struct count_list *operand2;
	char oper;
};

struct lig_count_list {
	struct lig_count_list *next;
        struct count_list *count_list;
};

struct lig_count_ref {
	struct lig_count_ref *next;
	unsigned short type;
        char *full_name;
        struct count_list *count_list;
};

struct viz_state_ref {
	struct viz_state_ref *next;
	int viz_state;
        char *full_name;
};

/**
 * Compartment. [\todo need more info]
 */
struct cmprt_data {
	struct cmprt_data *next;
	struct sym_table *sym;
        char *full_name;
        byte cmprt_type;            /**< type of compartment:
                                       BOX_POLY, ORDERED_POLY, UNORDERED_POLY */
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
	byte list_type;			/**< type of polygon list:
					   BOX_POLY, ORDERED_POLY,
                                           or UNORDERED_POLY */
        void *polygon_data;             /**< pointer to appropriate data structure
                                           holding polygon vertices etc...
                                           either:
                                           struct box_poly 
                                             for list_type BOX_POLY
                                           -or-
                                           struct ordered_poly 
                                             for list_type ORDERED_POLY
                                           -or-
                                           struct unordered_poly
                                             for list_type UNORDERED_POLY */
	int n_walls;			/**< Number of polygons in
                                             polygon object */
        int n_verts;                    /**< Number of vertices in
                                             polygon object */
	byte fully_closed;		/**< flag indicating closure of object */
        struct species **surf_class;    /** array of pointers to surface class, 
                                            one for each polygon */
        unsigned short *side_stat;	/**< array of side status values:
	                                   one for each polygon in object.
					   0 indicates removed
					   1 indicates exists. */
/*        struct eff_dat **eff_prop;*/	/**< array of ptrs to eff_dat data
					   structures, one for each polygon. */
};

/**
 * A polyhedron that is a box.
 */
struct box_poly {
        struct vector3 *llf;             /**< ptr to LLF of box polyhedron */
        struct vector3 *urb;             /**< prt to URB of box polyhedron */
};

/**
 * A general ordered polyhedron. 
 * That is, the vertices of each polygonal face are ordered according to
 * the right hand rule.
 */
struct ordered_poly {
	struct vector3 *vertex;         /**< Array of polygon vertices */
	struct vector3 *normal;         /**< Array of polygon normals */
	struct element_data *element_data; /**< Array element_data
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
        int *vertex_index;              /**< Array of vertex indices forming a
                                           polygon. */
	int n_verts;                    /**< Number of vertices in polygon. */
};

/**
 * A compartment.
 */
struct cmprt {
	struct sym_table *sym;
        unsigned short type;
	int inst_count;
	struct lig_count_list **lig_count_list;  /**< array of ptrs to lig_count_list
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
        unsigned int begin;
        unsigned int end;   
};

/**
 * A region.
 * [\todo what is this?]
 */
struct region {
	struct sym_table *sym;
	int hashval;
        char *region_last_name;
	struct object *parent;
	struct element_list *element_list_head;
	struct reg_counter_ref_list *reg_counter_ref_list;
	struct eff_dat *eff_dat_head;
        struct species *surf_class;
};

/**
 * Linked list of regions.
 * [\todo what is this?]
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
/*        struct eff_dat **eff_prop;*/	/**<  if this object is a
					   BOX_OBJ or POLY_OBJ this will be an
					   array of ptrs to eff_dat data
					   structures, one for each polygon. */
	struct cmprt_data *cmprt_data;	/**< if this object is a 
					   BOX_OBJ or POLY_OBJ this will
					   point to the cmprt_data struct
					   containing the instantiated object */
        int n_walls;                  /**< Total number of walls in object */
        struct wall *walls;           /**< array of walls in object */
        struct wall **wall_p;         /**< array of ptrs to walls in object */
        int n_verts;                  /**< Total number of vertices in object */
        struct vector3 *verts;        /**< array of vertices in object */
        struct vector3 **vert_p;      /**< array of ptrs to verts in object */
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
        byte list_type;
	int type;
	int viz_iteration;
	int n_viz_iterations;
	struct num_expr_list *iteration_list;
	struct num_expr_list *curr_viz_iteration;
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

#ifdef DEBUG
#define no_printf printf
#else
void no_printf(const char *,...);
#endif


#endif
