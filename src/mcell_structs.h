#ifndef MCELL_STRUCTS
#define MCELL_STRUCTS

#include <limits.h>
#include <sys/types.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include "config.h"

#include "rng.h"
#include "vector.h"
#include "mem_util.h"
#include "sched_util.h"
#include "util.h"

/* Macro for eliminating "unused variable" or "unused parameter" warnings. */
#define UNUSED(p) ((void) (p))

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

#if __GNUC__ < 3
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

/*****************************************************/
/**  Brand new constants created for use in MCell3  **/
/*****************************************************/

#define ORIENT_NOT_SET SHRT_MIN

/* Species flags */
   /* Surface classes have IS_SURFACE set, molecules do not. */
   /* Grid molecules have ON_GRID set */
   /* Volume molecules have NOT_FREE clear (i.e. flags&NOT_FREE==0) */
   /* Macromolecular complexes have IS_COMPLEX set */
   /* CAN_ flags specify what types of reactions this molecule can undergo */
   /* CANT_INITIATE means that this molecule may not trigger a reaction with another molecule */
   /* COUNT_TRIGGER means that someone wants to output a TRIGGER statement when something happens to this molecule */
   /* COUNT_CONTENTS is set if you're counting numbers of molecules in/on regions */
   /* COUNT_HITS is set if you're counting when the molecules hit regions */
   /* COUNT_RXNS is set if you're counting reactions involving this molecule */
   /* COUNT_ENCLOSED set if you count what happens inside closed region (vol molecules, or surface mols treated as if they were vol mols) */
   /* COUNT_SOME_MASK is a bitmask which is nonzero if any counting happens */
   /* CAN_REGION_BORDER is set when grid molecule can interact with region
      border that is declared REFLECTIVE/TRANSPARENT/ABSORPTIVE for that
      molecule */
   /* REGION_PRESENT set for the grid molecule when it is part of the
      SURFACE_CLASS definition and there are regions defined with
      this SURFACE_CLASS assigned */
#define ON_GRID          0x01
#define IS_SURFACE       0x02
#define NOT_FREE         0x03
#define TIME_VARY        0x04
#define CAN_MOLMOLMOL    0x08
#define CAN_MOLMOL       0x10
#define CAN_MOLGRID      0x20
#define CAN_MOLWALL      0x40
#define CAN_GRIDGRID     0x80
#define CAN_GRIDWALL     0x100
#define CAN_MOLMOLGRID   0x200
#define CANT_INITIATE    0x400
#define COUNT_TRIGGER    0x0800
#define COUNT_CONTENTS   0x1000
#define COUNT_HITS       0x2000
#define COUNT_RXNS       0x4000
#define COUNT_ENCLOSED   0x8000
#define COUNT_SOME_MASK  0xF800
#define CAN_MOLGRIDGRID  0x10000
#define CAN_GRIDGRIDGRID 0x20000
#define IS_COMPLEX       0x40000
#define SET_MAX_STEP_LENGTH  0x80000
#define CAN_REGION_BORDER 0x100000
#define REGION_PRESENT   0x200000


/* Abstract Molecule Flags */

/* RULES: only one of TYPE_GRID, TYPE_3D set. */
/*   ACT_NEWBIE beats ACT_INERT beats ACT_REACT */
/*   Can free up memory when nothing in IN_MASK */

/* Molecule type--grid molecule, 3D molecule, or surface molecule */
#define TYPE_GRID   0x001
#define TYPE_3D   0x002
#define TYPE_MASK 0x003

/* If this flag is set, this mol is a subunit in a complex */
#define COMPLEX_MEMBER 0x0004
#define COMPLEX_MASTER 0x2000

/* NEWBIE molecules get scheduled before anything else happens to them. */
/* INERT molecules don't react. */
/* ACT_REACT is set for molecules taking part in unimolecular reaction, or
   reaction with a surface */
/* ACT_INERT represent molecules in a catalytic dead-time */
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
/* And a mask to pick off all three IN_ flags */
#define IN_MASK     0x700

/* Flags telling us what our counting status is */
#define COUNT_ME    0x800

/* Flag indicating that a molecule is old enough to take the maximum timestep */
#define MATURE_MOLECULE 0x2000

/* End of Abstract Molecule Flags. */


/* Output Report Flags */
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


/* Manifold Flags */
enum manifold_flag_t
{
  MANIFOLD_UNCHECKED,   /* Manifold status is unknown */
  NOT_MANIFOLD,         /* Known to be non-manifold */
  IS_MANIFOLD           /* Known to be manifold */
};

/* Reaction flags */
  /* RX_ABSORB_REGION_BORDER signifies that a reaction is between a surface molecule and an ABSORPTIVE  region border */
  /* RX_REFLEC signifies that a reaction is between a molecule and a REFLECTIVE  wall */
  /* RX_TRANSP signifies that a reaction is between a molecule and a TRANSPARENT wall */
  /* Any value equal to or less than RX_SPECIAL refers to a special wall type */
  /* RX_BLOCKED signals a reaction that cannot take place because the grid is full */
  /* Any value equal to or less than RX_NO_RX indicates that a reaction did not take place */
  /* RX_FLIP signals that a molecule flips its orientation (crosses a wall if it's free) */
  /* RX_DESTROY signals that the molecule no longer exists (so don't try to keep using it) */
  /* RX_A_OK signals that all is OK with a reaction, proceed as normal (reflect if you're free) */
#define RX_ABSORB_REGION_BORDER  -5
#define RX_REFLEC  -4
#define RX_TRANSP  -3
#define RX_SPECIAL -3
#define RX_BLOCKED -2
#define RX_NO_RX   -2
#define RX_FLIP    -1
#define RX_LEAST_VALID_PATHWAY 0
#define RX_DESTROY  0
#define RX_A_OK     1
#define MAX_MATCHING_RXNS 64


/* Pathway flags */
/* TRANSPARENT means surface reaction between the molecule and TRANSPARENT wall */
/* REFLECTIVE means surface reaction between the molecule and REFLECTIVE wall */
/* CLAMP_CONC means surface reaction of CLAMP_CONCENTRATION type */
#define PATHW_TRANSP      0x0001
#define PATHW_REFLEC      0x0002
#define PATHW_ABSORP      0x0004
#define PATHW_CLAMP_CONC  0x0008   


/* BSP Flags */
/* (Currently unused--these are for self-subdividing subvolumes.) */
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


/* Direction Values */
#define X_NEG 0
#define X_POS 1
#define Y_NEG 2
#define Y_POS 3
#define Z_NEG 4
#define Z_POS 5
#define NODIR 6

/* Direction Bit Flags */
#define X_NEG_BIT 0x01
#define X_POS_BIT 0x02
#define Y_NEG_BIT 0x04
#define Y_POS_BIT 0x08
#define Z_NEG_BIT 0x10
#define Z_POS_BIT 0x20


/* Collision types for rays striking surfaces */
/* First a bunch of target types */
/* REDO happens if you can't tell whether you hit or not (hit near an edge, for example */
#define COLLIDE_REDO    -1
/* MISS means we hit nothing */
#define COLLIDE_MISS    0
/* FRONT and BACK are for surfaces */
#define COLLIDE_FRONT   1
#define COLLIDE_BACK    2
/* MOL_M is collision with a volume molecule */
#define COLLIDE_MOL_M   3
/* SV_?? is for collisions with subvolumes (negative and positive for each coordinate axis */
#define COLLIDE_SV_NX   4 
#define COLLIDE_SV_PX   5
#define COLLIDE_SV_NY   6
#define COLLIDE_SV_PY   7
#define COLLIDE_SV_NZ   8
#define COLLIDE_SV_PZ   9
/* A mask to pick off all of the collision target types */
#define COLLIDE_MASK    0x0F
/* Bitmasks for each of the major types of collision */
#define COLLIDE_WALL    0x10
#define COLLIDE_MOL     0x20 /* collision between 2 volume molecules */
#define COLLIDE_SUBVOL  0x40
#define COLLIDE_MOL_MOL 0x80  /* collision between 3 volume molecules */
#define COLLIDE_MOL_GRID 0x100 /* collision between 2 volume molecules
                                  and 1 grid molecule taken in the order
                                   mol-mol-grid */
#define COLLIDE_GRID_GRID 0x200 /* collision between 1 volume molecule
                                  and 2 grid molecules */
#define COLLIDE_GRID     0x400 /* bimolecular collision between moving
                                volume_molecule and grid_molecule */


/* Size constants */
/* EPS_C is the fractional difference between two values that is considered meaningful */
/* GIGANTIC is a distance that is larger than any possible simulation */
/* FOREVER is a time that cannot be reached within one simulation (too many timesteps) */
#define EPS_C 1e-12
#define SQRT_EPS_C 1e-6
#define GIGANTIC 1e140
#define FOREVER 1e20


/* How big will we let the reaction table get? */
/* 0x400000 = 8 million */
#define MAX_RX_HASH 0x400000

/* How big will we let the count-by-region table get? */
/* 0x10000 = 128K */
#define MAX_COUNT_HASH 0x10000 

/* mask for count-by-region hash */
#define COUNT_HASHMASK 0xffff

/* What's the upper bound on the number of coarse partitions? */
/* Not used for user-defined partitions */
#define MAX_COARSE_PER_AXIS 16
#define MIN_COARSE_PER_AXIS 6
#define MAX_TARGET_TIMESTEP 1.0e6
#define MIN_TARGET_TIMESTEP 10.0


/* Release Shape Flags */
enum release_shape_t
{
  SHAPE_UNDEFINED = -1,   /* Not specified */
  SHAPE_SPHERICAL,        /* Volume enclosed by a sphere */
  SHAPE_CUBIC,            /* Volume enclosed by a cube */
  SHAPE_ELLIPTIC,         /* Volume enclosed by an ellipsoid */
  SHAPE_RECTANGULAR,      /* Volume enclosed by a rect. solid */
  SHAPE_SPHERICAL_SHELL,  /* Surface of a sphere */
  SHAPE_REGION,           /* Inside/on the surface of an arbitrary region */
  SHAPE_LIST              /* Individiaul mol. placement by list */
};

/* Region Expression Flags */
/* Boolean set operations for releases on regions */
/* Set only one of NO_OP, UNION, INTERSECTION, SUBTRACTION, INCLUSION */
#define REXP_NO_OP        0x01
#define REXP_UNION        0x02
#define REXP_INTERSECTION 0x04
#define REXP_SUBTRACTION  0x08
#define REXP_INCLUSION    0x10
#define REXP_MASK         0x1F
#define REXP_LEFT_REGION  0x20
#define REXP_RIGHT_REGION 0x40

#ifdef MCELL_WITH_CHECKPOINTING
/* Distance in length units to search for a new site for a grid molecule */
/* after checkpointing.  Current site might be full, so a value >1 is */
/* advisable.  Being a little generous here. */
#define CHKPT_GRID_TOLERANCE 2.0
#endif


/* Constants for garbage collection of defunct molecules */
/* (those that were consumed when hit by another molecule) */
#define MIN_DEFUNCT_FOR_GC 1024
#define MAX_DEFUNCT_FRAC 0.2


/* Constants for notification levels */
enum notify_level_t
{
  NOTIFY_NONE,        /* no output */
  NOTIFY_BRIEF,       /* give a brief description (only used for a few types) */
  NOTIFY_FULL,        /* give a (possibly verbose) description */
};

/* Constants for warning levels */
enum warn_level_t
{
  WARN_COPE,          /* do something sensible and continue silently */
  WARN_WARN,          /* do something sensible but emit a warning message */
  WARN_ERROR          /* treat the warning and an error and stop */
};

/* Number of times to try diffusing on a surface before we give up (we might fail if the target grid is full) */
#define SURFACE_DIFFUSION_RETRIES 10

/* Overwrite Policy Flags */
/* Flags for different types of file output */
enum overwrite_policy_t
{
  FILE_UNDEFINED,     /* not specified */
  FILE_OVERWRITE,     /* always overwrite, even after a checkpoint */
  FILE_SUBSTITUTE,    /* DEFAULT: append to entries earlier in time than "now", but overwrite later entries */
  FILE_APPEND,        /* always append to file, even on a new run */
  FILE_APPEND_HEADER, /* always append to file, including the header, even on a new run */
  FILE_CREATE,        /* always create the file, or give an error if the file already exists (to prevent overwriting) */
};

/* Output Expression Flags */
/* INT means that this expression is an integer */
/* DBL means that this expression is a double */
/* TRIG means that this expression will actually be handled by a triggering event */
/* MASK lets us pick off the INT/DBL/TRIG flags */
/* CONST means that this expression will not change during runtime (compute at parse time and store) */
#define OEXPR_TYPE_UNDEF 0x0
#define OEXPR_TYPE_INT 0x1
#define OEXPR_TYPE_DBL 0x2
#define OEXPR_TYPE_TRIG 0x3
#define OEXPR_TYPE_MASK 0x7
#define OEXPR_TYPE_CONST 0x8

/* Same things for sub-expressions to the left, plus */
/* REQUEST means that the expression contains a request for a count statement, not real count data yet (needs to be initialized) */
/* OEXPR means that the expression is itself an expression that needs to be evaluated (not data) */
#define OEXPR_LEFT_INT 0x10
#define OEXPR_LEFT_DBL 0x20
#define OEXPR_LEFT_TRIG 0x30
#define OEXPR_LEFT_REQUEST 0x40
#define OEXPR_LEFT_OEXPR 0x50
#define OEXPR_LEFT_MACROREQUEST 0x60
#define OEXPR_LEFT_MASK 0x70
#define OEXPR_LEFT_CONST 0x80

/* Same things again for sub-expressions to the right */
#define OEXPR_RIGHT_INT 0x100
#define OEXPR_RIGHT_DBL 0x200
#define OEXPR_RIGHT_TRIG 0x300
#define OEXPR_RIGHT_REQUEST 0x400
#define OEXPR_RIGHT_OEXPR 0x500
#define OEXPR_RIGHT_MACROREQUEST 0x600
#define OEXPR_RIGHT_MASK 0x700
#define OEXPR_RIGHT_CONST 0x800


/* Magic value to indicate that a release pattern is actually a reaction */
/* Should be some number not between 0 and 1 that is also not -1 */
#define MAGIC_PATTERN_PROBABILITY 1.101001000100001

/* Output Trigger Flags */
/* Don't set both RXN and HIT flags */
#define TRIG_IS_RXN       0x1
#define TRIG_IS_HIT       0x2


/* Range of molecule indices used to avoid self-reactions for 3D unbinding */
#define DISSOCIATION_MAX -1000
#define DISSOCIATION_MIN -1000000000

#ifdef MCELL_WITH_CHECKPOINTING
/* Checkpoint related flags */
enum checkpoint_request_type_t
{
  CHKPT_NOT_REQUESTED,    /* No CP requested */
  CHKPT_SIGNAL_CONT,      /* CP requested via SIGUSR* signal, continue after CP */
  CHKPT_SIGNAL_EXIT,      /* CP requested via SIGUSR* signal, exit after CP  */
  CHKPT_ALARM_CONT,       /* CP requested via "alarm" signal, continue after CP */
  CHKPT_ALARM_EXIT,       /* CP requested via "alarm" signal, exit after CP */
  CHKPT_ITERATIONS_CONT,  /* CP requested due to iteration count, continue after CP */
  CHKPT_ITERATIONS_EXIT,  /* CP requested due to iteration count, exit after CP */
};
#endif

/*********************************************************/
/**  Constants used in MCell3 brought over from MCell2  **/
/*********************************************************/

/* 1/2^32 */
#define R_UINT_MAX 2.3283064365386962890625e-10

#define MY_PI 3.14159265358979323846
#define N_AV 6.0221417930e23
#define ROUND_UP 0.5
                                                                                
/* Placement Type Flags */
/* Place either a certain density or an exact number of surface molecules */
#define EFFDENS 0
#define EFFNUM 1

/* Viz output options */
#define VIZ_ALL_MOLECULES 0x01
#define VIZ_MOLECULES_STATES 0x02
#define VIZ_SURFACE_STATES 0x04
/* the formats below are valid only for the DREAMM_V3_MODE */
#define VIZ_MOLECULE_FORMAT_ASCII 0x08
#define VIZ_MOLECULE_FORMAT_BINARY 0x10
#define VIZ_MESH_FORMAT_ASCII 0x20
#define VIZ_MESH_FORMAT_BINARY 0x40
#define VIZ_ALL_MESHES 0x80

/************************************************************/
/**  Old constants copied from MCell2, some may be broken  **/
/************************************************************/

/* maximum allowed nesting level of INCLUDE_FILE statements in MDL */
#define MAX_INCLUDE_DEPTH 16

/* default size of output count buffers */
#define COUNTBUFFERSIZE 10000
                                                                                
/* Symbol types */
/* Data types for items in MDL parser symbol tables. */
enum symbol_type_t
{
  RX,               /* chemical reaction */
  RXPN,             /* name of chemical reaction */
  MOL,              /* molecule or surface class type (i.e. species) */
  OBJ,              /* meta-object */
  RPAT,             /* release pattern */
  REG,              /* object region */
  DBL,              /* double (numeric variable in MDL file) */
  STR,              /* string (text variable in MDL file) */
  ARRAY,            /* numeric array (array variable in MDL file) */
  FSTRM,            /* file stream type for "C"-style file-io in MDL file */
  TMP,              /* temporary place-holder type for assignment statements */
  VIZ_CHILD,        /* viz_child structures (in viz_output_block sym tables) */
};

/* Count column data types */
enum count_type_t
{
  COUNT_UNSET=-1,         /* no value specified */
  COUNT_DBL,              /* double */
  COUNT_INT,              /* integer type */
  COUNT_TRIG_STRUCT,      /* trigger_struct data type (for TRIGGER statements) */
};

/* Object Type Flags */
enum object_type_t
{
  META_OBJ,     /* Meta-object: aggregation of other objects */
  BOX_OBJ,      /* Box object: Polygonalized cuboid */
  POLY_OBJ,     /* Polygon list object: list of arbitrary triangles */
  REL_SITE_OBJ, /* Release site object */
  VOXEL_OBJ,    /* Voxel object (so-far unused) */
};

/* Box sides */
/* Note that there are two triangles per side, so we count up by two */
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

/* Data Output Timing Type */
/* Reaction and Viz data output timing */
enum output_timer_type_t
{
  OUTPUT_BY_STEP,
  OUTPUT_BY_TIME_LIST,
  OUTPUT_BY_ITERATION_LIST,
};

/* Visualization modes. */
enum viz_mode_t
{
  NO_VIZ_MODE,
  DX_MODE,
  DREAMM_V3_MODE,
  DREAMM_V3_GROUPED_MODE,
  RK_MODE,
  ASCII_MODE,
  CELLBLENDER_MODE,
};

/* Visualization Frame Data Type */
/* Used to select type of data to include in viz output files */
/* Will probably change significantly when we redesign DReAMM output format */
enum viz_frame_type_t
{
  ALL_FRAME_DATA,
  EFF_POS,
  EFF_STATES,
  MOL_POS,
  MOL_ORIENT,
  MOL_STATES,
  SURF_POS,
  SURF_STATES,
  MESH_GEOMETRY,
  REG_DATA,
  ALL_MOL_DATA,
  ALL_MESH_DATA,
  NUM_FRAME_TYPES,
};

/* Release Number Flags */
enum release_number_type_t
{
  CONSTNUM,
  GAUSSNUM,
  VOLNUM,
  CCNNUM,
  DENSITYNUM
};

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

/* Linked list used to separate molecules by species */
struct per_species_list
{
  struct per_species_list    *next;       /* pointer to next p-s-l */
  struct species             *properties; /* species for items in this bin */
  struct volume_molecule     *head;       /* linked list of mols */
};

/* Properties of one type of molecule or surface */
struct species
{
  u_int species_id;             /* Unique ID for this species */
#ifdef MCELL_WITH_CHECKPOINTING
  u_int chkpt_species_id;       /* Unique ID for this species from the 
                                   checkpoint file */
#endif
  u_int hashval;                /* Hash value (may be nonunique) */
  struct sym_table *sym;        /* Symbol table entry (name) */
  struct eff_dat *eff_dat_head; /* If IS_SURFACE this points to head of
                                   effector data list associated with 
                                   surface class */
  
  u_int population;             /* How many of this species exist? */
  
  double D;                     /* Diffusion constant */
  double D_ref;                 /* Reference diffusion constant. 
                                   For now - a placeholder for the future use */
  double space_step;            /* Characteristic step length */
  double time_step;             /* Minimum (maximum?) sensible timestep */
  double max_step_length;       /* maximum allowed random walk step */
  u_int flags;                /* Species Flags:  Vol Molecule? Surface Molecule?
                                   Surface Class? Counting stuff, etc... */
  
  long long n_deceased;         /* Total number that have been destroyed. */
  double cum_lifetime;          /* Timesteps lived by now-destroyed molecules */
 
  int region_viz_value;         /* Visualization state for surface class 
                                   for output */
  /* if species s a surface_class (IS_SURFACE) below there are linked lists of molecule names/orientations that may be present in special reactions for this surface class */
  struct name_orient *refl_mols;  /* names of the mols that REFLECT from surface */
  struct name_orient *transp_mols; /* names of the mols that are TRANSPARENT for surface */
  struct name_orient *absorb_mols; /* names of the mols that ABSORB at surface */
  struct name_orient *clamp_conc_mols; /* names of mols that CLAMP_CONC at surface */
};


/* All pathways leading away from a given intermediate */
struct rxn
{
  struct rxn *next;          /* next node in the reaction linked list where
                                each node contains only pathways with
                                equivalent geometry  */
  struct sym_table *sym;     /* Ptr to symbol table entry for this rxn */
  
  u_int n_reactants;         /* How many reactants? (At least 1.) */
  int n_pathways;            /* How many pathways lead away?
                                (Negative = special reaction, i.e. transparent etc...)*/
  double *cum_probs;         /* Cumulative probabilities for (entering) all pathways */
#if 0
  /* This is for Michaelis-Menten reaction kinetics, which are currently unimplemented. */
  double *cat_probs;         /* Probabilities of leaving all pathways (<=0.0 is instant) */
#endif
  struct complex_rate **rates; /* Rates for cooperative macromol subunit rxns */
  double max_fixed_p;        /* Maximum 'p' for region of p-space for all non-cooperative pathways */
  double min_noreaction_p;   /* Minimum 'p' for region of p-space which is always in the non-reacting "pathway". (note that cooperativity may mean that some values of p less than this still do not produce a reaction) */
  double pb_factor;          /* Conversion factor from rxn rate to rxn probability (used for cooperativity) */
  
  u_int *product_idx;        /* Index of 1st player for products of each pathway */
  struct species **players;  /* Identities of reactants/products */
  short *geometries;         /* Geometries of reactants/products */
  unsigned char *is_complex; /* Flags indicating which reactants/products are subunits of a complex.  array is NULL if not a macro-rxn. */
  
  long long n_occurred;      /* How many times has this reaction occurred? */
  double n_skipped;          /* How many reactions were skipped due to probability overflow? */

  struct t_func *prob_t;     /* List of probabilities changing over time, by pathway */
  
  struct pathway *pathway_head; /* List of pathways built at parse-time */
  struct pathway_info *info; /* Counts and names for each pathway */
};



/* User-defined name of a reaction pathway */
struct rxn_pathname {
  struct sym_table *sym;     /* Ptr to symbol table entry for this rxn name */
  u_int hashval;             /* Hash value for counting named rxns on regions */
  u_int path_num;            /* Pathway number in rxn */
  struct rxn *rx;            /* The rxn associated with this name */
  struct magic_list *magic;  /* A list of stuff that magically happens when the reaction happens */
};


/* Parse-time structure for reaction pathways */
/* Everything except pathname can be deallocated after prepare_reactions */
struct pathway {
  struct pathway *next;          /* Next pathway for this reaction */
  struct rxn_pathname *pathname; /* Data for named reaction pathway or NULL */
  struct species *reactant1;     /* First reactant in reaction pathway */
  struct species *reactant2;     /* Second reactant (NULL if none) */
  struct species *reactant3;     /* Third reactant (NULL if none) */
  unsigned char is_complex[3];   /* flag indicating whether each reactant must be a subunit in a complex */
  double km;                     /* Rate constant */
  char* km_filename;             /* Filename for time-varying rates */
  struct complex_rate *km_complex; /* Rate "constant" for cooperative subunit rxns */
  short orientation1;            /* Orientation of first reactant */
  short orientation2;            /* Orientation of second reactant */
  short orientation3;            /* Orientation of third reactant */
  struct product *product_head;  /* Linked lists of species created */
  char *prod_signature;         /* string created from the names of
                                   products put in alphabetical order */
  short flags;                  /* flags describing special reactions -
                                REFLECTIVE, TRANSPARENT, CLAMP_CONCENTRATION */

};

/* Parse-time structure for products of reaction pathways */
struct product {
  struct product *next;
  struct species *prod;          /* Molecule type to be created */
  short orientation;             /* Orientation to place molecule */
  unsigned char is_complex;      /* flag indicating whether product is to be a subunit in a complex */
};

/* Run-time info for each pathway */
/* Always do basic counts--do more sophisticated stuff if pathname!=NULL */
struct pathway_info
{
  double count;                 /* How many times the pathway has been taken */
  struct rxn_pathname *pathname;  /* The name of the pathway or NULL */
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
/* Used to make C structs act like C++ objects */
struct abstract_molecule
{
  struct abstract_molecule *next;  /* Next molecule in scheduling queue */
  double t;                        /* Scheduling time. */
  double t2;                       /* Time of next unimolecular reaction */
  short flags;                     /* Abstract Molecule Flags: Who am I, what am I doing, etc. */
  struct species *properties;      /* What type of molecule are we? */
  struct mem_helper *birthplace;   /* What was I allocated from? */
  double birthday;                 /* Time at which this particle was born */
  struct abstract_molecule **cmplx;  /* Other molecules forming this complex, if we're part of a complex (0: master, 1...n subunits) */
};


/* Volume molecules: freely diffusing or fixed in solution */
struct volume_molecule
{
  struct abstract_molecule *next;
  double t;
  double t2;
  short flags;
  struct species *properties;
  struct mem_helper *birthplace;
  double birthday;
  struct volume_molecule **cmplx; /* Other molecules forming this complex, if we're part of a complex (0: master, 1...n subunits) */
  
  struct vector3 pos;             /* Position in space */
  struct subvolume *subvol;       /* Partition we are in */
  
  struct wall *previous_wall;     /* Wall we were released from */
  int index;                      /* Index on that wall (don't rebind) */
  
  struct volume_molecule **prev_v;       /* Previous molecule in this subvolume */
  struct volume_molecule *next_v;        /* Next molecule in this subvolume */
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
  struct grid_molecule **cmplx; /* Other molecules forming this complex, if we're part of a complex (0: master, 1...n subunits) */
  
  unsigned int grid_index;     /* Which gridpoint do we occupy? */
  short orient;                /* Which way do we point? */
  struct surface_grid *grid;   /* Our grid (which tells us our surface) */
  struct vector2 s_pos;        /* Where are we in surface coordinates? */
};


/* Used to transform coordinates of grid molecules diffusing between adjacent walls */
struct edge
{
  struct wall *forward;     /* For which wall is this a forwards transform? */
  struct wall *backward;    /* For which wall is this a reverse transform? */
  
  struct vector2 translate;  /* Translation vector between coordinate systems */
  double cos_theta;          /* Cosine of angle between coordinate systems */
  double sin_theta;          /* Sine of angle between coordinate systems */
  
  double length;             /* Length of the shared edge */
  double length_1;           /* Reciprocal of length of shared edge */
};


struct wall
{
  struct wall *next;              /* Next wall in the universe */
  
  struct surf_class_list *surf_class_head;      /* linked list of surface classes for this wall (multiple surface classes may come from the overlapping regions */
  int num_surf_classes;           /* number of attached surface classes */

  int side;                       /* index of this wall in its parent object */

  struct vector3 *vert[3];        /* Array of pointers to vertices */

  double uv_vert1_u;              /* Surface u-coord of 2nd corner (v=0) */
  struct vector2 uv_vert2;        /* Surface coords of third corner */

  struct edge *edges[3];          /* Array of pointers to each edge. */
  struct wall *nb_walls[3];       /* Array of pointers to walls that share an edge*/

  double area;                    /* Area of this element */
  
  struct vector3 normal;          /* Normal vector for this wall */
  struct vector3 unit_u;          /* U basis vector for this wall */
  struct vector3 unit_v;          /* V basis vector for this wall */
  double d;                       /* Distance to origin (point normal form) */
  
  struct surface_grid *grid; /* Grid of effectors for this wall */
  
  u_short flags;                  /* Count Flags: flags for whether and what we need to count */

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


/* A linked list used to store the coordinates of vertices and the 
   corresponding normal vectors */
struct vertex_list
{
   struct vector3 *vertex;   		/* pointer to one polygon vertex */
   struct vertex_list *next;		/* pointer to next vertex list */
};


/* Grid over a surface containing grid_molecules */
struct surface_grid
{
  int n;                   /* Number of slots along each axis */

  double inv_strip_wid;    /* Reciprocal of the width of one strip */
  double vert2_slope;      /* Slope from vertex 0 to vertex 2 */
  double fullslope;        /* Slope of full width of triangle */
  struct vector2 vert0;    /* Projection of vertex zero onto unit_u and unit_v of wall */
  
  double binding_factor;   /* Binding probability correction factor for surface area */
  
  u_int n_tiles;           /* Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2) */
  u_int n_occupied;        /* Number of tiles occupied by grid_molecules */
  struct grid_molecule **mol;  /* Array of pointers to grid_molecule for each tile */
  
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


/* Point in space that will tell us which compartments we're in
   as determined by tracing from infinity */
struct waypoint
{
  struct vector3 loc;              /* This is where the waypoint is */
  struct region_list *regions;     /* We are inside these regions */
  struct region_list *antiregions; /* We are outside of (but hit) these regions */
};


/* Contains local memory and scheduler for molecules, walls, wall_lists, etc. */
struct storage
{
  struct mem_helper *list;  /* Wall lists */
  struct mem_helper *mol;   /* Molecules */
  struct mem_helper *gmol;  /* Grid molecules */
  struct mem_helper *face;  /* Walls */
  struct mem_helper *join;  /* Edges */
  struct mem_helper *grids;  /* Effector grids */
  struct mem_helper *coll;  /* Collision list */
  struct mem_helper *sp_coll;  /* Collision list - helps in trimolecular reactions*/
  struct mem_helper *tri_coll;  /* Collision list for trimolecular collisions */
  struct mem_helper *regl;  /* Region lists */
  struct mem_helper *exdv;  /* Vertex lists for exact interaction disk area */
  struct mem_helper *pslv;  /* Per-species-lists for vol mols */
  
  struct wall *wall_head;              /* Locally stored walls */
  int wall_count;                      /* How many local walls? */
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
  
  struct pointer_hash mol_by_species;   /* table of species->molecule list */
  struct per_species_list *species_head;
  int mol_count;               /* How many molecules are here? */
  
  struct short3D llf;          /* Indices of left lower front corner */
  struct short3D urb;          /* Indices of upper right back corner */
  
  short world_edge;            /* Direction Bit Flags that are set for SSVs at edge of world */ 
/*  short is_bsp;   */             /* Flags saying what the void pointers are */

  /* void *neighbor[6];   */        /* Subvolume or bsp_tree across each face */
  
  struct storage *local_storage; /* Local memory and scheduler */
};


/* Binary space partitioning tree for subvolume connections */
/* Not used in the current release but may be useful in future. */
/*
struct bsp_tree
{
  void *left;        // The tree below the partition 
  void *right;       // The tree above the partition 
  short partition;   // The index of the partition 
  short flags;       // BSP Flags: coordinate that is split, plus terminal node flags 
};
*/

/* Count data specific to named reaction pathways */
struct rxn_counter_data
{
  double n_rxn_at;                  /* # rxn occurrance on surface */
  double n_rxn_enclosed;            /* # rxn occurrance inside closed region */
};


/* Count data specific to molecules */
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


/* Counter data specific to trigger events */
struct trig_counter_data
{
  double t_event;                    /* Event time (exact) */
  struct vector3 loc;                /* Real position of event */
  short orient;                      /* For MOL_COUNTER: molecule orientation */
  struct trigger_request *listeners; /* Places waiting to be notified */
};


/* List of output items that need to know about this specific trigger event */
struct trigger_request
{
  struct trigger_request *next; /* Next request */
  struct output_request *ear;   /* Who wants to hear about the trigger */
};


/* Shared memory for appropriate counts */
union counter_data
{
  struct rxn_counter_data rx;
  struct move_counter_data move;
  struct trig_counter_data trig;
};


/* Struct to count rxns or molecules within regions (where "within" includes */
/* on the inside of a fully closed surface) */
struct counter
{
  struct counter *next;
  byte counter_type;               /* Counter Type Flags (MOL_COUNTER etc.) */
  struct region *reg_type;         /* Region we are counting on */
  void *target;                    /* Mol or rxn pathname we're counting (as indicated by counter_type) */
  short orientation;               /* requested grid molecule orientation */
  union counter_data data;         /* data for the count:
                                      reference data.move for move counter
                                      reference data.rx for rxn counter 
                                      reference data.trig for trigger */
};


enum magic_types { magic_undefined,magic_release };

struct magic_list
{
  struct magic_list *next;
  void *data;
  enum magic_types type;
};

/* All data about the world */
struct volume
{
  /* Coarse partitions are input by the user */
  /* They may also be generated automagically */
  /* They mark the positions of initial partition boundaries */
  int nx_parts;                 /* Number of coarse X partition boundaries */
  int ny_parts;                 /* Number of coarse Y partition boundaries */
  int nz_parts;                 /* Number of coarse Z partition boundaries */
  double *x_partitions;         /* Coarse X partition boundaries */
  double *y_partitions;         /* Coarse Y partition boundaries */
  double *z_partitions;         /* Coarse Z partition boundaries */
  int mem_part_x;               /* Granularity of memory-partition binning for the X-axis */
  int mem_part_y;               /* Granularity of memory-partition binning for the Y-axis */
  int mem_part_z;               /* Granularity of memory-partition binning for the Z-axis */
  int mem_part_pool;            /* Scaling factor for sizes of memory pools in each storage. */
  
  /* Fine partitions are intended to allow subdivision of coarse partitions */
  /* Subdivision is not yet implemented */
  int n_fineparts;              /* Number of fine partition boundaries */
  double *x_fineparts;          /* Fine X partition boundaries */
  double *y_fineparts;          /* Fine Y partition boundaries */
  double *z_fineparts;          /* Fine Z partition boundaries */
  
  int n_waypoints;              /* How many waypoints (one per subvol) */
  struct waypoint *waypoints;   /* Waypoints contain fully-closed region information */
  byte place_waypoints_flag;    /* Used to save memory if waypoints not needed */

  int n_subvols;                /* How many coarse subvolumes? */
  struct subvolume *subvol;     /* Array containing all subvolumes */
   
  int n_walls;                  /* Total number of walls */
  int n_verts;                  /* Total number of vertices */
  struct vector3 *all_vertices;  /* Central repository of vertices with a
                                    partial order imposed by natural ordering
                                    of "storages" */
  struct wall_list **walls_using_vertex; /* Array of linked lists of walls using                                         a vertex (has the size of
                                         "all_vertices" array */ 
  int rx_hashsize;                 /* How many slots in our reaction hash table? */
  int n_reactions;                 /* How many reactions are there, total? */
  struct rxn **reaction_hash;      /* A hash table of all reactions. */
  struct mem_helper *tv_rxn_mem;   /* Memory to store time-varying reactions */
  
  int count_hashmask;                      /* Mask for looking up count hash table */
  struct counter **count_hash;             /* Count hash table */
  struct schedule_helper *count_scheduler; /* When to generate reaction output */

  struct schedule_helper *volume_output_scheduler; /* When to generate volume output */
  
  int n_species;                 /* How many different species (molecules)? */
  struct species **species_list; /* Array of all species (molecules). */
  
  struct schedule_helper *releaser; /* Scheduler for release events */
  
  struct mem_helper *storage_allocator; /* Memory for storage list */
  struct storage_list *storage_head;    /* Linked list of all local memory/schedulers */
  
  double speed_limit;           /* How far can the fastest particle get in one timestep? */

  struct sym_table_head *fstream_sym_table;   /* Global MDL file stream symbol hash table */
  struct sym_table_head *var_sym_table;   /* Global MDL variables symbol hash table */
  struct sym_table_head *rxn_sym_table;   /* RXN symbol hash table */
  struct sym_table_head *obj_sym_table;   /* Objects symbol hash table */
  struct sym_table_head *reg_sym_table;   /* Regions symbol hash table */
  struct sym_table_head *mol_sym_table;   /* Molecule type symbol hash table */
  struct sym_table_head *rpat_sym_table;  /* Release pattern hash table */
  struct sym_table_head *rxpn_sym_table;  /* Named reaction pathway hash table */

  struct object *root_object;         /* Root of the object template tree */
  struct object *root_instance;       /* Root of the instantiated object tree */

  struct release_pattern *default_release_pattern;  /* release once at t=0 */

  struct volume_output_item *volume_output_head; /* List of all volume data output items */
  
  struct macro_count_request *macro_count_request_head;
  struct output_block *output_block_head;     /* Global list of reaction data output blocks */
  struct output_request *output_request_head; /* Global list linking COUNT statements to internal variables */
  struct mem_helper *oexpr_mem;               /* Memory to store output_expressions */
  struct mem_helper *outp_request_mem;        /* Memory to store output_requests */
  struct mem_helper *counter_mem;             /* Memory to store counters (for counting molecules/reactions on regions) */
  struct mem_helper *trig_request_mem;        /* Memory to store listeners for trigger events */
  struct mem_helper *magic_mem;               /* Memory used to store magic lists for reaction-triggered releases and such */
  double elapsed_time;                        /* Used for concentration measurement */

  /* Visualization state */
  struct viz_output_block *viz_blocks;        /* VIZ_OUTPUT blocks from file */

  struct species *all_mols;   /* Refers to ALL_MOLECULES keyword */
  struct species *all_volume_mols; /* Refers to ALL_VOLUME_MOLECULES keyword */
  struct species *all_surface_mols; /* Refers to ALL_SURFACE_MOLECULES keyword */

  double time_unit;        /* Duration of one global time step in real time */
                           /* Used to convert between real time and internal time */
  double time_step_max;    /* Maximum internal time that a molecule may diffuse */

  double grid_density;     /* Density of grid for surface molecules, number per um^2 */
  double length_unit;      /* Internal unit of distance, 1/sqrt(grid_density), in microns */
  double r_length_unit;    /* Reciprocal of length_unit to avoid division */
  double rx_radius_3d;     /* Interaction radius for reactions between volume molecules */

  double space_step;       /* User-supplied desired average diffusion distance for volume molecules */

  double *r_step;          /* Lookup table of 3D diffusion step lengths */
  double *d_step;          /* Lookup table of 3D diffusion direction vectors */
  double *r_step_surface;  /* Lookup table of 2D diffusion step lengths */
  double *r_step_release;  /* Lookup table of diffusion lengths for 3D release */
  u_int radial_subdivisions; /* Size of 2D and 3D step length lookup tables */
  u_int radial_directions; /* Requested size of 3D direction lookup table */
  u_int num_directions;    /* Actual size of 3D direction lookup table */
  int directions_mask;     /* Mask to obtain RNG bits for direction lookup */
  int fully_random;        /* If set, generate directions with trig functions instead of lookup table */
  int dissociation_index;  /* Used to keep 3D products from reacting with each other too soon */

  int complex_placement_attempts;   /* How many times will we try to place each complex before giving up? */

#ifdef MCELL_WITH_CHECKPOINTING
  long long chkpt_iterations; /* Number of iterations to advance before checkpointing */
  u_int chkpt_init; /* Set if this is the initial run of a simulation with no previous checkpoints */
  u_int chkpt_flag; /* Set if there are any CHECKPOINT statements in "mdl" file */
  u_int chkpt_seq_num;        /* Number of current run in checkpoint sequence */

  char *chkpt_infile;         /* Name of checkpoint file to read from */
  char *chkpt_outfile;        /* Name of checkpoint file to write to */
  u_int chkpt_byte_order_mismatch;   /* Flag that defines whether mismatch
                                      in byte order exists between the saved
                                      checkpoint file and the machine reading it */

  double chkpt_elapsed_real_time;    /* elapsed simulation time (in sec) for new checkpoint */
  double chkpt_elapsed_real_time_start;  /* start of the simulation time (in sec) for new checkpoint */
#endif
  double current_real_time;          /* current simulation time in seconds */
  double current_start_real_time;    /* simulation start time (in seconds) */

  long long diffusion_number;      /* Total number of times molecules have had their positions updated */
  double diffusion_cumtime;     /* Total time spent diffusing by all molecules */
  long long ray_voxel_tests;       /* How many ray-subvolume intersection tests have we performed */
  long long ray_polygon_tests;     /* How many ray-polygon intersection tests have we performed */
  long long ray_polygon_colls;     /* How many ray-polygon intersections have occured */
  /* below "mol" means volume molecule, "grid" means surface molecule */
  long long mol_mol_colls;         /* How many mol-mol collisions have occured */
  long long mol_grid_colls;     /* How many mol-grid collisions have occured */
  long long grid_grid_colls;     /* How many grid-grid collisions have occured */
  long long mol_wall_colls;     /* How many mol-wall collisions have occured */
  long long mol_mol_mol_colls;     /* How many mol-mol-mol collisions have occured */
  long long mol_mol_grid_colls;     /* How many mol-mol-grid collisions have occured */
  long long mol_grid_grid_colls;     /* How many mol-grid-grid collisions have occured */
  long long grid_grid_grid_colls;     /* How many grid-grid-grid collisions have occured  */

  struct vector3 bb_llf;	/* llf corner of world bounding box */
  struct vector3 bb_urb;	/* urb corner of world bounding box */

  struct rng_state *rng;        /* State of the random number generator (currently isaac64) */
  u_int init_seed;              /* Initial seed value for random number generator */

  long long it_time;      /* How many iterations have been run so far */
  long long start_time;   /* Starting iteration number for the current run */

  int procnum;            /* Processor number for a parallel run */
  int quiet_flag;               /* Quiet mode */
  int with_checks_flag;  /* Check geometry for overlapped walls? */

  struct mem_helper *coll_mem;  /* Collision list */
  struct mem_helper *sp_coll_mem;  /* Collision list (trimol) */
  struct mem_helper *tri_coll_mem;  /* Collision list (trimol) */
  struct mem_helper *exdv_mem;  /* Vertex lists for exact interaction disk area */

  char const *mcell_version;     /* Current version number.
                              Format is "3.XX.YY" where XX is major release number (for new features)
                              and YY is minor release number (for patches) */
 
  int use_expanded_list;       /* If set, check neighboring subvolumes for mol-mol interactions */
  int randomize_gmol_pos;      /* If set, always place surface molecule at random location instead of center of grid */
  double vacancy_search_dist2; /* Square of distance to search for free grid location to place surface product */
  byte surface_reversibility;  /* If set, match unbinding diffusion distribution to binding distribution at surface */
  byte volume_reversibility;   /* If set, match unbinding diffusion distribution to binding distribution in volume */

  /* MCell startup command line arguments */
  u_int seed_seq;            /* Seed for random number generator */
  long long iterations;      /* How many iterations to run */
  unsigned long log_freq;    /* Interval between simulation progress reports, default scales as sqrt(iterations) */
  char *mdl_infile_name;     /* Name of MDL file specified on command line */
  char const *curr_file;     /* Name of MDL file currently being parsed */
  
  /* XXX: Why do we allocate this on the heap rather than including it inline? */
  struct notifications *notify; /* Notification/warning/output flags */
  
  struct ccn_clamp_data *clamp_list;  /* List of objects at which volume molecule concentrations should be clamped */
  
#ifdef MCELL_WITH_CHECKPOINTING
  /* Flags for asynchronously-triggered checkpoints */
  enum checkpoint_request_type_t checkpoint_requested; /* Flag indicating whether a checkpoint has been requested. */
  unsigned int checkpoint_alarm_time;   /* number of seconds between checkpoints */
  int continue_after_checkpoint;        /* 0: exit after chkpt, 1: continue after chkpt */
  long long last_checkpoint_iteration;  /* Last iteration when chkpt was created */
#endif
  time_t begin_timestamp;               /* Time since epoch at beginning of 'main' */
  char *initialization_state;           /* NULL after initialization completes */
  /* flags that tells whether reactions of certain types are present
     in the simulation (used for the molecule collision report, also see above
     the corresponding counters) */
  int mol_mol_reaction_flag;
  int mol_grid_reaction_flag;
  int grid_grid_reaction_flag;
  int mol_wall_reaction_flag;
  int mol_mol_mol_reaction_flag;
  int mol_mol_grid_reaction_flag;
  int mol_grid_grid_reaction_flag;
  int grid_grid_grid_reaction_flag;
  /* shared walls information per mesh vertex is created when there are
     reactions present with more than one surface reactant
     or more than one surface product */
  int create_shared_walls_info_flag;
  /* resource usage during initialization */
  struct timeval u_init_time; /* user time */ 
  struct timeval s_init_time; /* system time */
  time_t t_start; /* global start time */
  byte reaction_prob_limit_flag; /* checks whether there is at least one 
                                    reaction with probability greater 
                                    than 1 including variable rate reactions */ 
};


/* Data structure to store information about collisions. */
struct collision
{
  struct collision *next;
  double t;                     /* Time of collision (may be slightly early) */
  
  void *target;                 /* Thing that we hit: wall, molecule, subvol etc */
  int what;                     /* Target-type Flags: what kind of thing did we hit? */
  struct rxn *intermediate;     /* Reaction that told us we could hit it */
  struct vector3 loc;           /* Location of impact */
};

/* Special type of collision - used when moving molecule
   can engage in tri-molecular  collisions */
struct sp_collision
{
  struct sp_collision *next;
  double t;                     /* Time of collision (may be slightly early) */
  double t_start;                /* Start time of random walk */
  struct vector3 pos_start;      /* Start position of random walk */
  struct subvolume *sv_start;    /* Start subvolume */
  
  struct species *moving;       /* Species of the moving molecule */
  void *target;                 /* Thing that we hit: wall, molecule, subvol etc */
  int what;                     /* Target-type Flags: what kind of thing did we hit? */
  struct vector3 disp;          /* Random walk displacement for the moving molecule */
  struct vector3 loc;           /* Location of impact */
};



/* Data structure to store information about trimolecular and bimolecular 
   collisions. */
struct tri_collision
{
  struct tri_collision *next;
  double t;                     /* Time of collision (may be slightly early) */
  
  void *target1;                 /* First thing that we hit: wall, molecule, subvol etc */
  void *target2;                 /* Second thing that we hit: wall, molecule, subvol etc - always the furthest from the moving molecule */
  short orient;                  /* orientation of the moving volume_molecule
                                    when it hits the surface_molecule */
  int what;                     /* Target-type Flags: what kind of thing did we hit? */
  struct rxn *intermediate;     /* Reaction that told us we could hit                                    target1 and/or target2  */
  struct vector3 loc;           /* Assumed location of impact */
  struct vector3 loc1;           /* Location of impact with first target */
  struct vector3 loc2;           /* Location of impact with second target */
  struct vector3 last_walk_from; /* Location of mol. before last step before final collision */
  double factor;                /* Result of "exact_disk()" with both targets
                                   or scaling coef. for MOL_WALL interaction */
  double local_prob_factor;    /* coefficient depending on the number of 
                                  nearest neighbors for MOL_GRID_GRID 
                                  interaction */
  struct wall *wall;          /* pointer to the wall in the collision if
                                 such exists  */
};

/* Data structures to store information about exact interaction disk geometry */
struct exd_vertex
{
  struct exd_vertex *next;
  double u,v;              /* x,y style coordinates */
  double r2,zeta;          /* r,theta style coordinates */
  struct exd_vertex *e;    /* Edge to next vertex */
  struct exd_vertex *span; /* List of edges spanning this point */
  int role;                /* Exact Disk Flags: Head, tail, whatever */
};


/* Data structures to describe release events */
struct release_event_queue {
  struct release_event_queue *next;
  double event_time;			 /* Time of the release */
  struct release_site_obj *release_site; /* What to release, where to release it, etc */
  double t_matrix[4][4];                 /* transformation matrix for location of release site */
  int train_counter;			 /* counts executed trains */
  double train_high_time;		 /* time of the train's start */
};


/* Release site information  */
struct release_site_obj {
  struct vector3 *location;	  /* location of release site */
  struct species *mol_type;	  /* species to be released */
  byte release_number_method;     /* Release Number Flags: controls how release_number is used (enum release_number_type_t) */
  int8_t release_shape;           /* Release Shape Flags: controls shape over which to release (enum release_shape_t) */
  short orientation;              /* Orientation of released surface molecules */
  double release_number;             /* Number to release */
  double mean_diameter;           /* Diameter for symmetric releases */
  double concentration;           /* Concentration of molecules to release.
                     Units are Molar for volume molecules, and number per um^2 for surface molecules. */
  double standard_deviation;      /* Standard deviation of release_number for GAUSSNUM,
                     or of mean_diameter for VOLNUM */
  struct vector3 *diameter;       /* x,y,z diameter for geometrical release shapes */
  struct release_region_data *region_data; /* Information related to release on regions */
  struct release_single_molecule *mol_list; /* Information related to release by list */

  double release_prob;            /* Probability of releasing at scheduled time */
  struct release_pattern *pattern;  /* Timing of releases by virtual function generator */
  char *name; /* Fully referenced name of the instantiated release_site */
};


/* Timing pattern for molecule release from a release site. */
struct release_pattern {
  struct sym_table *sym;        /* Symbol hash table entry for the pattern */
  double delay;			/* Delay between time 0 and first release event. */
  double release_interval;	/* Time between release events within a train. */
  double train_interval;	/* Time from the start of one train to the start of the next one. */
  double train_duration;	/* Length of the train. */
  int number_of_trains;		/* How many trains are produced. */
};


/* Extended data for complex releases on regions,
   including boolean combinations of regions.
   Not all fields are used for all release types. */
struct release_region_data
{
  struct vector3 llf;   /* One corner of bounding box for release volume */
  struct vector3 urb;   /* Opposite corner */

  int n_walls_included;  /* How many walls total */
  double *cum_area_list; /* Cumulative area of all walls */
  int *wall_index;       /* Indices of each wall (by object) */
  int *obj_index;        /* Indices for objects (in owners array) */
  
  int n_objects;                 /* How many objects are there total */
  struct object **owners;        /* Array of pointers to each object */
  struct bit_array **in_release; /* Array of bit arrays; each bit array says which walls are in release for an object */
  int refinement;                /* Number of volume-inclusion checks to do (2d releases only) */
  int *walls_per_obj;            /* Number of walls in release for each object */

  struct object *self;           /* A pointer to our own release site object */
  struct release_evaluator *expression;  /* A set-construction expression combining regions to form this release site */
};


/* Data structure used to build boolean combinations of regions */
struct release_evaluator
{
  byte op;       /* Region Expression Flags: the operation used */
  void *left;    /* The left side of the expression--another evaluator or a region object depending on bitmask of op */
  void *right;   /* The right side--same thing */
};


/* Data structure used to store LIST releases */
struct release_single_molecule
{
  struct release_single_molecule *next;
  struct species *mol_type;              /* Species to release */
  struct vector3 loc;                    /* Position to release it */
  short orient;                          /* Orientation (for 2D species) */
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
  /* Informational stuff, most possible values NOTIFY_FULL or NOTIFY_NONE */
  /* see corresponding keywords */
  enum notify_level_t progress_report;              /* PROGRESS_REPORT */
  enum notify_level_t diffusion_constants;          /* DIFFUSION_CONSTANT_REPORT */
  enum notify_level_t reaction_probabilities;       /* PROBABILITY_REPORT */
  enum notify_level_t time_varying_reactions;       /* VARYING_PROBABILITY_REPORT */
  double reaction_prob_notify;                      /* PROBABILITY_REPORT_THRESHOLD */
  enum notify_level_t partition_location;           /* PARTITION_LOCATION_REPORT */
  enum notify_level_t box_triangulation;            /* BOX_TRIANGULATION_REPORT */
  enum notify_level_t iteration_report;             /* ITERATION_REPORT */
  long long custom_iteration_value;                 /* ITERATION_REPORT */
  enum notify_level_t throughput_report;            /* THROUGHPUT_REPORT */
  enum notify_level_t checkpoint_report;            /* CHECKPOINT_REPORT */
  enum notify_level_t release_events;               /* RELEASE_EVENT_REPORT */
  enum notify_level_t file_writes;                  /* FILE_OUTPUT_REPORT */
  enum notify_level_t final_summary;                /* FINAL_SUMMARY */
  enum notify_level_t reaction_output_report;       /* REACTION_OUTPUT_REPORT */
  enum notify_level_t volume_output_report;         /* VOLUME_OUTPUT_REPORT */
  enum notify_level_t viz_output_report;            /* VIZ_OUTPUT_REPORT */
  enum notify_level_t molecule_collision_report;    /* MOLECULE_COLLISION_REPORT */
  
  /* Warning stuff, possible values IGNORED, WARNING, ERROR */
  /* see corresponding keywords */
  enum warn_level_t neg_diffusion;                  /* NEGATIVE_DIFFUSION_CONSTANT */
  enum warn_level_t neg_reaction;                   /* NEGATIVE_REACTION_RATE */
  enum warn_level_t high_reaction_prob;             /* HIGH_REACTION_PROBABILITY */
  double reaction_prob_warn;                        /* HIGH_PROBABILITY_THRESHOLD */
  enum warn_level_t close_partitions;               /* CLOSE_PARTITION_SPACING */
  enum warn_level_t degenerate_polys;               /* DEGENERATE_POLYGONS */
  enum warn_level_t overwritten_file;               /* OVERWRITTEN_OUTPUT_FILE */
  enum warn_level_t short_lifetime;                 /* LIFETIME_TOO_SHORT */
  long long short_lifetime_value;                   /* LIFETIME_THRESHOLD */
  enum warn_level_t missed_reactions;               /* MISSED_REACTIONS */
  double missed_reaction_value;                     /* MISSED_REACTION_THRESHOLD */
  enum warn_level_t missed_surf_orient;             /* MISSING_SURFACE_ORIENTATION */
  enum warn_level_t useless_vol_orient;             /* USELESS_VOLUME_ORIENTATION */
  enum warn_level_t complex_placement_failure;      /* COMPLEX_PLACEMENT_FAILURE */
  long long complex_placement_failure_threshold;    /* COMPLEX_PLACEMENT_FAILURE_THRESHOLD */
  enum warn_level_t mol_placement_failure;          /* MOLECULE_PLACEMENT_FAILURE */
  enum warn_level_t invalid_output_step_time;       /* INVALID_OUTPUT_STEP_TIME */
};

/* Information related to concentration clamp surfaces, by object */
struct ccn_clamp_data
{
  struct ccn_clamp_data *next;     /* The next concentration clamp, by surface class */
  struct species *surf_class;      /* Which surface class clamps? */
  struct species *mol;             /* Which molecule does it clamp? */
  double concentration;            /* At which concentration? */
  short orient;                    /* On which side? */
  struct object *objp;             /* Which object are we clamping? */
  struct bit_array *sides;         /* Which walls in that object? */
  int n_sides;                     /* How many walls? */
  int *side_idx;                   /* Indices of the walls that are clamped */
  double *cum_area;                /* Cumulative area of all the clamped walls */
  double scaling_factor;           /* Used to predict #mols/timestep */
  struct ccn_clamp_data *next_mol; /* Next clamp, by molecule, for this class */
  struct ccn_clamp_data *next_obj; /* Next clamp, by object, for this class */
};

/* Structure for a VOLUME_DATA_OUTPUT item */
struct volume_output_item
{
  /* Do not move or reorder these 2 items.  scheduler depends upon them */
  struct volume_output_item    *next;
  double                        t;

  char                         *filename_prefix;

  /* what? */
  int                           num_molecules;
  struct species              **molecules;          /* sorted by address */

  /* where? */
  struct vector3                location;
  struct vector3                voxel_size;
  int                           nvoxels_x;
  int                           nvoxels_y;
  int                           nvoxels_z;

  /* when? */
  enum output_timer_type_t      timer_type;
  double                        step_time;
  int                           num_times;
  double                       *times;              /* in numeric order  */
  double                       *next_time;          /* points into times */
};

/* Data for a single REACTION_DATA_OUTPUT block */
struct output_block
{
  struct output_block *next;            /* Next in world or scheduler */
  double t;                             /* Scheduled time to update counters */
  
  enum output_timer_type_t timer_type;  /* Data Output Timing Type (OUTPUT_BY_STEP, etc) */
  
  double step_time;                     /* Output interval (seconds) */
  struct num_expr_list *time_list_head; /* List of output times/iteration numbers */
  struct num_expr_list *time_now;       /* Current entry in list */
  
  u_int buffersize;                     /* Size of output buffer */
  u_int trig_bufsize;                   /* Size of output buffer for triggers */
  u_int buf_index;                      /* Index into buffer (for non-triggers) */
  
  double *time_array;                   /* Array of output times (for non-triggers) */
  
  struct output_set *data_set_head;     /* Linked list of data sets (separate files) */
};


/* Data that controls what output is written to a single file */
struct output_set
{
  struct output_set *next;             /* Next data set in this block */
  struct output_block *block;          /* Which block do we belong to? */
  char *outfile_name;                  /* Filename */
  enum overwrite_policy_t file_flags;  /* Overwrite Policy Flags: tells us how to handle existing files */
  u_int chunk_count;                   /* Number of buffered output chunks processed */  
  char *header_comment;                /* Comment character(s) for header */
  int exact_time_flag;                 /* Boolean value; nonzero means print exact time in TRIGGER statements */
  struct output_column *column_head;   /* Data for one output column */
};


/* Data that controls what data is written to one column of output file */
struct output_column
{
  struct output_column *next;       /* Next column in this set */
  struct output_set *set;           /* Which set do we belong to? */
  enum count_type_t data_type;      /* Type of data in this column. */
  double initial_value;             /* To continue existing cumulative counts--not implemented yet--and keep track of triggered data */
  void *buffer;                     /* Output buffer array (cast based on data_type) */
  struct output_expression *expr;   /* Evaluate this to calculate our value (NULL if trigger) */
};


/* Expression evaluation tree to compute output value for one column */
struct output_expression
{
  struct output_column *column;     /* Which column are we going to? */
  int expr_flags;                   /* Output Expression Flags: what are we and what is left and right? */
  struct output_expression *up;     /* Parent output expression */
  void *left;                       /* Item on the left */
  void *right;                      /* Item on the right */
  char oper;                        /* Operation to apply to items */
  double value;                     /* Resulting value from operation */
  char *title;                      /* String describing this expression for user */
};


/* Information about a requested COUNT or TRIGGER statement */
/* Used during initialization to link output expressions with appropriate target,
   and instruct MCell3 to collect appropriate statistics at the target. */
struct output_request
{
  struct output_request *next;          /* Next request in global list */
  struct output_expression *requester;  /* Expression in which we appear */
  struct sym_table *count_target;       /* Mol/rxn we're supposed to count */
  short count_orientation;              /* orientation of the molecule
                                           we are supposed to count */
  struct sym_table *count_location;     /* Object or region on which we're supposed to count it */
  byte report_type;                     /* Output Report Flags telling us how to count */
};


/* Data stored when a trigger event happens */
struct output_trigger_data
{
  double t_iteration;      /* Iteration time of the triggering event (in sec) */
  double event_time;          /* Exact time of the  event */
  struct vector3 loc;      /* Position of event */
  int how_many;            /* Number of events */
  short orient;            /* Orientation information */
  short flags;             /* Output Trigger Flags */
  char *name;              /* Name to give event */
};



/******************************************************************/
/**  Everything below this line has been copied from MCell 2.69  **/
/******************************************************************/


/* A polygon list object, part of a surface. */
struct polygon_object {
  int n_verts;                       /* Number of vertices in polyhedron */
  struct vertex_list *parsed_vertices;   /* Temporary linked list */
  int n_walls;                       /* Number of triangles in polyhedron */
  struct element_data *element;      /* Array specifying the vertex connectivity of each triangle */
  struct subdivided_box *sb;         /* Holds corners of box if necessary */
  struct bit_array *side_removed;    /* Bit array; if bit is set, side is removed */
};


/* Data structure used to build one triangular polygon
  according to the connectivity in the MDL file. */
struct element_data {
  int vertex_index[3];  /* Array of vertex indices forming a triangle */
};


/* A voxel list object, part of a volume */
struct voxel_object
{
  struct vector3 *vertex;               /* Array of tetrahedron vertices */
  struct tet_element_data *element;     /* Array tet_element_data data structures */
  struct tet_neighbors_data *neighbor;  /* Array tet_neighbors_data data structures */
  int n_verts;                          /* Number of vertices in polyhedron */
  int n_voxels;                         /* Number of voxels in polyhedron */
};

/**
 * Data structure used to build one tetrahedron.
 * This data structure is used to store the data from the MDL file
 * and to contruct each tetrahedron of a voxel object.
 */
struct tet_element_data
{
  int vertex_index[4];                  /* Array of vertex indices forming a tetrahedron. */
};


/**
 * This data structure is used to store the data about neighbors
 * of each tetrahedron of a voxel object.
 */
struct tet_neighbors_data
{
  int neighbors_index[4];               /* Array of indices pointing to the neighbors of tetrahedron. */
};


/* Surface molecule placement data */
struct eff_dat {
  struct eff_dat *next;
  struct species *eff; /* Species to place on surface */
  byte quantity_type;  /* Placement Type Flags: either EFFDENS or EFFNUM */
  double quantity;     /* Amount of surface molecules to place by density or number */
  short orientation;   /* Orientation of molecules to place */
};


/* Linked list of wall index ranges for regions */
struct element_list {
  struct element_list *next;
  u_int begin;                      /* First number in the range */
  u_int end;                        /* Last number in the range */
  struct element_special *special;  /* Pre-existing region or patch on box */
};


/* Elements can be patches on boxes or pre-existing regions */
struct element_special
{
  struct vector3 corner1;  /* Corner of patch on box */
  struct vector3 corner2;  /* Opposite corner of patch on box */
  struct region *referent; /* Points to pre-existing region on object */
  byte exclude;            /* If set, remove elements rather than include them */
};


/* Region of an object */
/* If region is a manifold then it can be used as both a volume and surface region.
   Otherwise it can only be used as a surface region. */
struct region {
  struct sym_table *sym;                  /* Symbol hash table entry for this region */
  u_int hashval;                          /* Hash value for counter hash table */
  char *region_last_name;                 /* Name of region without prepended object name */
  struct object *parent;                  /* Parent of this region */
  struct element_list *element_list_head; /* List of element ranges comprising this region (used at parse time) */
  struct bit_array *membership;           /* Each bit indicates whether the corresponding wall is in the region */
  struct eff_dat *eff_dat_head;           /* List of surface molecules to add to region */ 
  struct species *surf_class;             /* Surface class of this region */
  struct vector3 *bbox;                   /* Array of length 2 to hold corners of region bounding box
                                             (used for release in region) */
  int region_viz_value;                   /* Used for visualization */
  double area;                            /* Area of region */
  u_short flags;                          /* Counting subset of Species Flags */
  byte manifold_flag;                     /* Manifold Flags: If IS_MANIFOLD, region is a closed manifold and thus defines a volume */
  struct pointer_hash *boundaries;       /* hash table of edges that constitute
                                            external boundary of the region */
  int region_has_all_elements;           /* flag that tells whether the region
                                            contains ALL_ELEMENTS (effectively
                                            comprises the whole object) */
};


/* A list of regions */
struct region_list {
  struct region_list *next;
  struct region *reg;        /* A region */
};


/* Container data structure for all physical objects */
struct object {
  struct object *next;          /* Next sibling object */
  struct object *parent;        /* Parent meta object */
  struct object *first_child;	/* First child object */
  struct object *last_child;	/* Last child object */
  struct sym_table *sym;        /* Symbol hash table entry for this object */
  char *last_name;              /* Name of object without pre-pended parent object name */
  enum object_type_t  object_type; /* Object Type Flags */
  void *contents;		/* Actual physical object, cast according to object_type */
  u_int num_regions;	        /* Number of regions defined on object */
  struct region_list *regions;  /* List of regions for this object */
  int n_walls;                  /* Total number of walls in object */
  int n_walls_actual;           /* Number of non-null walls in object */
  struct wall *walls;           /* Array of walls in object */
  struct wall **wall_p;         /* Array of ptrs to walls in object (used at run-time) */
  int n_verts;                  /* Total number of vertices in object */
  struct vector3 **vertices;    /* Array of pointers to vertices 
                                         (linked to "all_vertices" array) */
  double total_area;            /* Area of object in length units */
  u_int n_tiles;                /* Number of surface grid tiles on object */
  u_int n_occupied_tiles;       /* Number of occupied tiles on object */
  double t_matrix[4][4];	/* Transformation matrix for object */
};


/* Doubly linked list of object names */
struct name_list {
  struct name_list *next;
  struct name_list *prev;
  char *name;              /* An object name */
};

/* Linked list of names-orientations. Used in printing special reactions report
   for surface classes */
struct name_orient{
   struct name_orient *next;
   char *name; /* molecule name */
   int orient; /* molecule orientation */
};



/* Visualization objects */
struct viz_dx_obj 
{
  struct viz_dx_obj *next;
  char *name;                        /* Name taken from OBJECT_FILE_PREFIXES
                                    or FILENAME_PREFIXES or FILENAME assignment */
  char *full_name;                   /* Full name of the object, like A.B.C */
  struct object *obj;                /* The object being visualized */
  struct viz_child *viz_child_head;  /* List of child objects to visualize */
  struct viz_output_block *parent;   /* Parent block to whom we belong */
  struct viz_child **actual_objects; /* Pointers to actual objects to visualize */
  int n_actual_objects;              /* Number of actual objects to visualize */
};


/* Tree of pointers to objects (mirrors standard geometry hierarchy). */
/* Contains child polygon or box objects to be visualized. */
struct viz_child
{
  struct viz_child *next;
  struct viz_child *parent;
  struct viz_child *children;
  struct object *obj;      /* An object to visualize*/
  int *viz_state;          /* Array of viz state values, one for each element of object. */
};

struct visualization_state
{
  /* Iteration numbers */
  long long last_meshes_iteration;
  long long last_mols_iteration;

  /* Tokenized filename prefix */
  char *filename_prefix_basename;
  char *filename_prefix_dirname;

  /* All visualized volume molecule species */
  int n_vol_species;
  struct species **vol_species;

  /* All visualized grid molecule species */
  int n_grid_species;
  struct species **grid_species;

  /* Iteration numbers and times of outputs */
  struct iteration_counter output_times;
  struct iteration_counter mesh_output_iterations;
  struct iteration_counter vol_mol_output_iterations;
  struct iteration_counter grid_mol_output_iterations;

  /* For DREAMM V3 Grouped output, combined group member strings */
  struct string_buffer combined_group_members;

  /* For DREAMM V3 Grouped output, the current object number */
  int dx_main_object_index;

  /* For DREAMM V3 Grouped output, the last iteration for certain outputs */
  long long dreamm_last_iteration_meshes;
  long long dreamm_last_iteration_vol_mols;
  long long dreamm_last_iteration_surf_mols;

  /* For DREAMM V3 Ungrouped output, path of 'frame data dir' and iter dir */
  char *frame_data_dir;
  char *iteration_number_dir;
};

struct viz_output_block
{
  struct viz_output_block *next;                /* Link to next block */
  struct frame_data_list *frame_data_head;      /* head of the linked list of viz frames to output */
  struct visualization_state viz_state_info;    /* miscellaneous state for viz_output code */
  enum viz_mode_t viz_mode;
  char *molecule_prefix_name;
  char *file_prefix_name;
  u_short viz_output_flag; /* Takes  VIZ_ALL_MOLECULES, VIZ_MOLECULES_STATES, etc.  */
  int *species_viz_states;

  int default_mesh_state;                       /* Only set if (viz_output_flag & VIZ_ALL_MESHES) */
  int default_mol_state;                        /* Only set if (viz_output_flag & VIZ_ALL_MOLECULES) */

  /* DREAMM-mode only. */
  struct viz_child **dreamm_object_info; /* Pointers to actual objects to visualize */
  struct object **dreamm_objects;
  int n_dreamm_objects;              /* Number of actual objects to visualize */

  /* DX-mode only: head of linked list of OBJECT_FILE_PREFIXES. */
  struct viz_dx_obj *dx_obj_head;

  /* RK-mode only. */
  struct rk_mode_data *rk_mode_var;

  /* Parse-time only: Tables to hold temporary information. */
  struct sym_table_head *viz_children;
  struct pointer_hash parser_species_viz_states;
};

/* Geometric transformation data for a physical object */
struct transformation {
  struct vector3 translate;  /* X,Y,Z translation vector */
  struct vector3 scale;      /* X,Y,Z scaling factors */
  struct vector3 rot_axis;   /* Vector defining an axis of rotation */
  double rot_angle;          /* Rotation angle in degrees */
};


/* Linked list of viz data to be output */
struct frame_data_list
{
  struct frame_data_list *next;
  enum output_timer_type_t list_type;       /* Data Output Timing Type (OUTPUT_BY_TIME_LIST, etc) */
  enum viz_frame_type_t type;               /* Visualization Frame Data Type (ALL_FRAME_DATA, etc) */ 
  long long viz_iteration;	            /* Value of the current iteration step. */
  long long n_viz_iterations;	            /* Number of iterations in the iteration_list. */
  struct num_expr_list *iteration_list;     /* Linked list of iteration steps values */
  struct num_expr_list *curr_viz_iteration; /* Points to the current iteration in the linked list */
};


/* A pointer to filehandle and its real name */
/* Used for user defined file IO operations */
struct file_stream {
  char *name;   /* File name */
  FILE *stream; /* File handle structure */
};

/* Symbol hash table */
/* Used to parse and store user defined symbols from the MDL input file */
struct sym_table_head {
  struct sym_table **entries; 
  int n_entries;
  int n_bins;
};

/* Symbol hash table entry */
/* Used to parse and store user defined symbols from the MDL input file */
/* XXX: This is a poorly named structure.  Maybe "sym_entry" or even just "symbol"? */
struct sym_table {
  struct sym_table *next;       /* Chain to next symbol in this bin of the hash */
  int sym_type;                 /* Symbol Type */
  char *name;                   /* Name of symbol*/
  void *value;                  /* Stored value, cast by sym_type */
};


/* Linked list of symbols */
/* Used to parse and retrieve user defined symbols having wildcards from the MDL input file */
struct sym_table_list {
  struct sym_table_list *next;
  struct sym_table *node;      /* Symbol table entry matching a user input wildcard string */
};


/* Linked list of numerical expressions */
/* Used for parsing MDL input file arithmetic expressions */
struct num_expr_list {
  struct num_expr_list *next;
  double value;                /* Value of one element of the expression */
};
 

/* Histogram visualization output mode used by Rex */
struct rk_mode_data
{
  int n_bins;
  int* bins;
  double *parts;
  struct vector3 *direction;
  int n_written;
};

/* Linked list of surface classes */
struct surf_class_list {
  struct surf_class_list *next;
  struct species *surf_class;
};

/* Linked list of edges - used for REGION borders */
struct edge_list
{
  struct edge_list *next;
  struct edge *ed;
};

/* Data about hits/crossing on the region border */
struct hit_data
{
  struct hit_data *next;
  struct region_list *count_regions; /* list of regions we counting on */
  int direction; /* 1 - INSIDE_OUT, 0 - OUTSIDE_IN */
  int crossed;  /* 1 - if crossed, 0 - if not */
  short orientation;  /* orientation of the grid molecule */
  struct vector3 loc;  /* location of the hit */
  double t;            /* time of the hit */
};

#endif
