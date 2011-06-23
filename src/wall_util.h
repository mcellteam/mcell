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

/* This linked list node is used in walls overlap test */
struct wall_aux_list
{
  struct wall *this_wall;  /* wall */
  double d_prod;    /* dot product of wall's normal and random vector */
  struct wall_aux_list *next; 
};

struct plane {
	struct vector3 n;	/* Plane normal.  Points x on the plane satisfy
                                   dot_prod(n,x) = d */
        double d;                /* d = dot_prod(n,p) for a given point p on 
                                    the plane */
};


int edge_equals(struct poly_edge *e1,struct poly_edge *e2);
int edge_hash(struct poly_edge *pe,int nkeys);

int ehtable_init(struct edge_hashtable *eht,int nkeys);
int ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe);
void ehtable_kill(struct edge_hashtable *eht);

int surface_net( struct wall **facelist, int nfaces );
void init_edge_transform(struct edge *e,int edgenum);
int sharpen_object(struct object *parent);
int sharpen_world(void);

double closest_interior_point(struct vector3 *pt,struct wall *w,struct vector2 *ip,double r2);

int find_edge_point(struct wall *here,struct vector2 *loc,struct vector2 *disp,struct vector2 *edgept);
struct wall* traverse_surface(struct wall *here,struct vector2 *loc,int which,struct vector2 *newloc);
int is_manifold(struct region *r);

void jump_away_line(struct vector3 *p,struct vector3 *v,double k,
                    struct vector3 *A,struct vector3 *B,struct vector3 *n);
double touch_wall(struct vector3 *point,struct vector3 *move,struct wall *face);
int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt,int update_move);
int collide_mol(struct vector3 *point,struct vector3 *move,
                struct abstract_molecule *a,double *t,struct vector3 *hitpt);

int intersect_box(struct vector3 *llf,struct vector3 *urb,struct wall *w);

void init_tri_wall(struct object *objp,int side,
                   struct vector3 *v0,struct vector3 *v1,struct vector3 *v2);

struct wall_list* wall_to_vol(struct wall *w, struct subvolume *sv);
struct wall* localize_wall(struct wall *w, struct storage *stor);
struct wall* distrubute_wall(struct wall *w);
int distribute_object(struct object *parent);
int distribute_world(void);

void closest_pt_point_triangle(struct vector3 *p, struct vector3 *a, struct vector3 *b, struct vector3 *c, struct vector3 *final_result);
int test_sphere_triangle(struct vector3 *s, double radius, struct vector3 *a, struct vector3 *b, struct vector3 *c, struct vector3 *p);
void compute_plane(struct vector3 *a, struct vector3 *b, struct vector3 *c, struct plane *p);
int test_sphere_plane(struct vector3 *s, double radius, struct plane *p); 
int test_sphere_ray(struct vector3 *p, struct vector3 *d, struct vector3 *s,
           double radius, double *t, struct vector3 *q);
int test_segment_plane(struct vector3 *a, struct vector3 *b, struct plane *p, double *t, struct vector3 *q);
int test_bounding_boxes(struct vector3 *llf1, struct vector3 *urb1, struct vector3 *llf2, struct vector3 *urb2);

int surface_point_in_region(struct object *ob,int wall_n,struct vector3 *v,struct release_evaluator *expr);
int release_onto_regions(struct release_site_obj *rso,struct grid_molecule *g,int n);

void push_wall_to_list(struct wall_list **wall_nbr_head, struct wall *w);
void delete_wall_list(struct wall_list *wl_head);

struct wall_list* find_nbr_walls_shared_one_vertex(struct wall *origin, int  *shared_vert);
int wall_share_vertex(struct wall *w, struct vector3 *vert);
int walls_share_full_edge(struct wall *w1, struct wall *w2);
struct region_list * find_region_by_wall(struct object *parent, struct wall *this_wall);
int is_wall_edge_region_border(struct wall *this_wall, struct edge *this_edge);
int find_shared_edge_index_of_neighbor_wall(struct wall *orig_wall, struct wall  *nbr_wall);
void find_neighbor_wall_and_edge(struct wall *orig_wall, int orig_edge_ind,     struct wall **nbr_wall, int *nbr_edge_ind);
int wall_contains_both_vertices(struct wall *w, struct vector3 *vert_A, struct vector3 *vert_B);
int are_walls_coincident(struct wall *w1, struct wall *w2, double eps);
int are_walls_coplanar(struct wall *w1, struct wall *w2, double eps);
int overlap_coplanar_walls(struct wall *w1, struct wall *w2);

int overlap_tri_tri_3d(double  p1[3], double  q1[3], double  r1[3],
		       double  p2[3], double  q2[3], double  r2[3],
		       double  N1[3], double  N2[3]);
int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]);

void sorted_insert_wall_aux_list(struct wall_aux_list **headRef, struct wall_aux_list *newNode);
void delete_wall_aux_list(struct wall_aux_list *head);
int walls_belong_to_same_region(struct wall *w1, struct wall * w2);
  
#endif
