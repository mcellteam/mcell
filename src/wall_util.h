#ifndef MCELL_WALL_UTIL
#define MCELL_WALL_UTIL

#include "mcell_structs.h"

/* Temporary data stored about an edge of a polygon */  
struct poly_edge
{
  struct poly_edge *next; /* Next edge in a hash table. */

  double v1x;             /* X coord of starting point */
  double v1y;             /* Y coord of starting point */
  double v1z;             /* Z coord of starting point */
  double v2x;             /* X coord of ending point */
  double v2y;             /* Y coord of ending point */
  double v2z;             /* Z coord of ending point */
  
  int face1;              /* Index of wall on one side of edge */
  int face2;              /* Index of wall on other side of edge */
  int edge1;              /* Which edge of wall1 are we? */
  int edge2;              /* Which edge of wall2 are we? */
  int n;                  /* How many walls share this edge? */
};


/* Hash table for rapid order-invariant lookup of edges. */
struct edge_hashtable
{
  struct poly_edge *data; /* Array of polygon edges */   
  
  int nkeys;              /* Length of array */
  int stored;            /* How many things do we have in the table? */
  int distinct;           /* How many of those are distinct? */
};


int edge_equals(struct poly_edge *e1,struct poly_edge *e2);
int edge_hash(struct poly_edge *pe,int nkeys);

int ehtable_init(struct edge_hashtable *eht,int nkeys);
int ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe);
void ehtable_kill(struct edge_hashtable *eht);

int surface_net( struct wall **facelist, int nfaces );
void init_edge_transform(struct edge *e,int edgenum);
int sharpen_object(struct object *parent);
int sharpen_world();

void jump_away_line(struct vector3 *p,struct vector3 *v,double k,
                    struct vector3 *A,struct vector3 *B,struct vector3 *n);
double touch_wall(struct vector3 *point,struct vector3 *move,struct wall *face);
int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt);
int collide_mol(struct vector3 *point,struct vector3 *move,
                struct abstract_molecule *a,double *t,struct vector3 *hitpt);

int intersect_box(struct vector3 *llf,struct vector3 *urb,struct wall *w);

void init_tri_wall(struct object *objp,int side,
                   struct vector3 *v0,struct vector3 *v1,struct vector3 *v2);

struct wall_list* wall_to_vol(struct wall *w, struct subvolume *sv);
struct vector3* localize_vertex(struct vector3 *p, struct storage *stor);
struct wall* localize_wall(struct wall *w, struct storage *stor);
struct wall* distrubute_wall(struct wall *w);
int distribute_object(struct object *parent);
int distribute_world();


#endif

