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

void ehtable_init(struct edge_hashtable *eht,int nkeys);
void ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe);
void ehtable_kill(struct edge_hashtable *eht);

int surface_net( struct wall **facelist, int nfaces );
void init_edge_transform(struct edge *e,int edgenum);

void jump_away_line(struct vector3 *p,struct vector3 *v,double k,
                    struct vector3 *A,struct vector3 *B,struct vector3 *n);
int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt);
int collide_mol(struct vector3 *point,struct vector3 *move,
                struct abstract_molecule *a,double *t,struct vector3 *hitpt);

#endif

