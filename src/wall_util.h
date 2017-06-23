/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#pragma once

#include "mcell_structs.h"

/* Temporary data stored about an edge of a polygon */
struct poly_edge {
  struct poly_edge *next; /* Next edge in a hash table. */

  double v1x; /* X coord of starting point */
  double v1y; /* Y coord of starting point */
  double v1z; /* Z coord of starting point */
  double v2x; /* X coord of ending point */
  double v2y; /* Y coord of ending point */
  double v2z; /* Z coord of ending point */

  int face[2]; /* wall indices on side of edge */
  int edge[2]; /* which edge of wall1/2 are we? */
  int n;     /* How many walls share this edge? */
};

/* Hash table for rapid order-invariant lookup of edges. */
struct edge_hashtable {
  struct poly_edge *data; /* Array of polygon edges */

  int nkeys;    /* Length of array */
  int stored;   /* How many things do we have in the table? */
  int distinct; /* How many of those are distinct? */
};

/* This linked list node is used in walls overlap test */
struct wall_aux_list {
  struct wall *this_wall; /* wall */
  double d_prod;          /* dot product of wall's normal and random vector */
  struct wall_aux_list *next;
};

struct plane {
  struct vector3 n; /* Plane normal.  Points x on the plane satisfy
                       dot_prod(n,x) = d */
  double d; /* d = dot_prod(n,p) for a given point p on
               the plane */
};

int edge_equals(struct poly_edge *e1, struct poly_edge *e2);
int edge_hash(struct poly_edge *pe, int nkeys);

int ehtable_init(struct edge_hashtable *eht, int nkeys);
int ehtable_add(struct edge_hashtable *eht, struct poly_edge *pe);
void ehtable_kill(struct edge_hashtable *eht);

int surface_net(struct wall **facelist, int nfaces);
void init_edge_transform(struct edge *e, int edgenum);
int sharpen_object(struct object *parent);

int sharpen_world(struct volume *world);

double closest_interior_point(struct vector3 *pt, struct wall *w,
                              struct vector2 *ip, double r2);

int find_edge_point(struct wall *here, struct vector2 *loc,
                    struct vector2 *disp, struct vector2 *edgept);
struct wall *traverse_surface(struct wall *here, struct vector2 *loc, int which,
                              struct vector2 *newloc);
int is_manifold(struct region *r, int count_regions_flag);

void jump_away_line(struct vector3 *p, struct vector3 *v, double k,
                    struct vector3 *A, struct vector3 *B, struct vector3 *n,
                    struct rng_state *rng);

int collide_wall(struct vector3 *point, struct vector3 *move, struct wall *face,
                 double *t, struct vector3 *hitpt, int update_move,
                 struct rng_state *rng, struct notifications *notify,
                 long long *polygon_tests);

int collide_mol(struct vector3 *point, struct vector3 *move,
                struct abstract_molecule *a, double *t, struct vector3 *hitpt,
                double rx_radius_3d);

int intersect_box(struct vector3 *llf, struct vector3 *urb, struct wall *w);

void init_tri_wall(struct object *objp, int side, struct vector3 *v0,
                   struct vector3 *v1, struct vector3 *v2);

struct wall_list *wall_to_vol(struct wall *w, struct subvolume *sv);

struct wall *localize_wall(struct wall *w, struct storage *stor);

int distribute_object(struct volume *world, struct object *parent);

int distribute_world(struct volume *world);

void closest_pt_point_triangle(struct vector3 *p, struct vector3 *a,
                               struct vector3 *b, struct vector3 *c,
                               struct vector3 *final_result);

int test_bounding_boxes(struct vector3 *llf1, struct vector3 *urb1,
                        struct vector3 *llf2, struct vector3 *urb2);

int release_onto_regions(struct volume *world, struct release_site_obj *rso,
                         struct surface_molecule *sm, int n);

struct surface_molecule *place_single_molecule(struct volume *state,
                                               struct wall *w,
                                               unsigned int grid_index,
                                               struct species *spec,
                                               short flags, short orientation,
                                               double t, double t2,
                                               double birthday,
                                               struct periodic_image *periodic_box,
                                               struct vector3 *pos3d);


void push_wall_to_list(struct wall_list **wall_nbr_head, struct wall *w);
void delete_wall_list(struct wall_list *wl_head);

struct wall_list *find_nbr_walls_shared_one_vertex(struct volume *world,
                                                   struct wall *origin,
                                                   long long int *shared_vert);

int walls_share_full_edge(struct wall *w1, struct wall *w2);

struct region_list *find_region_by_wall(struct wall *this_wall);

struct name_list *find_regions_names_by_wall(
    struct wall *w, struct string_buffer *ignore_regs);

struct region_list *
find_restricted_regions_by_wall(struct volume *world, struct wall *this_wall,
                                struct surface_molecule *sm);

struct region_list *
find_restricted_regions_by_object(struct volume *world, struct object *obj,
                                  struct surface_molecule *sm);

int are_restricted_regions_for_species_on_object(struct volume *world,
                                                 struct object *obj,
                                                 struct surface_molecule *sm);

int is_wall_edge_region_border(struct wall *this_wall, struct edge *this_edge);

int is_wall_edge_restricted_region_border(struct volume *world,
                                          struct wall *this_wall,
                                          struct edge *this_edge,
                                          struct surface_molecule *sm);

int find_shared_edge_index_of_neighbor_wall(struct wall *orig_wall,
                                            struct wall *nbr_wall);

void find_neighbor_wall_and_edge(struct wall *orig_wall, int orig_edge_ind,
                                 struct wall **nbr_wall, int *nbr_edge_ind);

int wall_contains_both_vertices(struct wall *w, struct vector3 *vert_A,
                                struct vector3 *vert_B);

int are_walls_coincident(struct wall *w1, struct wall *w2, double eps);

int are_walls_coplanar(struct wall *w1, struct wall *w2, double eps);

void sorted_insert_wall_aux_list(struct wall_aux_list **headRef,
                                 struct wall_aux_list *newNode);

void delete_wall_aux_list(struct wall_aux_list *head);

int walls_belong_to_at_least_one_different_restricted_region(
    struct volume *world, struct wall *w1, struct surface_molecule *sm1,
    struct wall *w2, struct surface_molecule *sm2);

int wall_belongs_to_all_regions_in_region_list(struct wall *w,
                                               struct region_list *rlp_head);

int wall_belongs_to_any_region_in_region_list(struct wall *w,
                                              struct region_list *rlp_head);

int region_belongs_to_region_list(struct region *rp, struct region_list *head);

void find_wall_center(struct wall *w, struct vector3 *center);
