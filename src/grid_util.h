#ifndef MCELL_GRID
#define MCELL_GRID

#include "mcell_structs.h"

#define TILE_CHECKED 0x01

/* contains information about the neigbors of the tile */
struct tile_neighbor{
  struct surface_grid *grid; /* surface grid the tile is on */
  unsigned int idx; /* index on that tile */
  short int flag;  /* flag */
  struct tile_neighbor *next;
};

void xyz2uv(struct vector3 *a,struct wall *w,struct vector2 *b);
void uv2xyz(struct vector2 *a,struct wall *w,struct vector3 *b);

int xyz2grid(struct vector3 *v,struct surface_grid *g);
void grid2xyz(struct surface_grid *g,int idx,struct vector3 *v);

int uv2grid(struct vector2 *v,struct surface_grid *g);
void grid2uv(struct surface_grid *g,int idx,struct vector2 *v);
void grid2uv_random(struct surface_grid *g,int idx,struct vector2 *v);

void init_grid_geometry(struct surface_grid *g);
int create_grid(struct wall *w,struct subvolume *guess);

void grid_neighbors(struct surface_grid *grid,int idx, int create_grid_flag, struct surface_grid **nb_grid,int *nb_idx);
int nearest_free(struct surface_grid *g,struct vector2 *v,double max_d2,double *found_dist2);
struct wall *search_nbhd_for_free(struct wall *origin,struct vector2 *point,double max_d2,int *found_idx,
                             int (*ok)(void*,struct wall*),void *context);

			     
int grid_release_check(struct release_region_data *rrd,int obj_n,int wall_n,int grid_n, struct release_evaluator *expr);

void delete_tile_neighbor_list(struct tile_neighbor *head);
void push_tile_neighbor_to_list(struct tile_neighbor **head, struct surface_grid *grid, int idx);
int push_tile_neighbor_to_list_with_checking(struct tile_neighbor **head, struct surface_grid *grid, int idx);
int add_more_tile_neighbors_to_list_fast(struct tile_neighbor **head, struct surface_grid *orig_grid, int orig_strip, int orig_stripe, int orig_flip, struct vector3 *start, struct vector3 *end, int edge_index, struct surface_grid *new_grid);
int is_neighbor_tile(struct surface_grid *orig_grid, int orig_idx, struct surface_grid *new_grid, int new_idx);
int get_tile_neighbor_from_list_of_vacant_neighbors(struct tile_neighbor *head, int index, struct surface_grid **grid, int *idx);
void find_closest_position(struct surface_grid *grid1, int idx1, struct surface_grid *grid2, int idx2, struct vector2 *p);                     
int is_inner_tile(struct surface_grid *g, int idx);
void find_neighbor_tiles(struct grid_molecule *g, int create_grid_flag, int search_for_reactant, struct tile_neighbor **tile_nbr_head, int *list_length);
void grid_all_neighbors_for_inner_tile(struct surface_grid *grid,int idx,struct vector2 *pos, struct tile_neighbor **tile_nbr_head, int *list_length);
void grid_all_neighbors_across_walls_through_edges(struct grid_molecule *g, int create_grid_flag, int search_for_reactant, struct tile_neighbor **tile_nbr_head, int *list_length);
void grid_all_neighbors_across_walls_through_vertices(struct grid_molecule *g, struct wall_list *wall_nbr_head, int create_grid_flag, int search_for_reactant,  struct tile_neighbor **tile_nbr_head, int *list_length);
void append_tile_neighbor_list(struct tile_neighbor **head1, struct tile_neighbor **head2);

void get_tile_vertices(struct surface_grid *sg, int idx, int *flp, struct vector2 *R, struct vector2 *S, struct vector2 *T);
int tile_orientation(struct vector2 *v, struct surface_grid *g);
double get_tile_area(struct vector2 *A, struct vector2 *B, struct vector2 *C);

int move_strip_up(struct surface_grid *grid, int idx);
int move_strip_down(struct surface_grid *grid, int idx);
void place_product_shared_segment(struct vector2 *R_shared, struct vector2 *S_shared, struct vector2 *T, struct vector2* prod, double k1, double k2);
void place_product_shared_vertex(struct vector2 *R_shared, struct vector2 *S, struct vector2 *T, struct vector2 *prod, double k1, double k2);
void place_product_close_to_segment_endpoint(struct vector2 *S, struct vector2 *E, struct vector2 *prod, double k1, double k2);
int is_corner_tile(struct surface_grid *g, int idx);
void find_shared_vertices_corner_tile_parent_wall(struct surface_grid *g, int idx, int *shared_vert);
void find_shared_vertices_for_neighbor_walls(struct wall *orig_wall, struct wall *nb_wall, int *shared_vert_1, int *shared_vert_2); 
int find_wall_vertex_for_corner_tile(struct surface_grid *grid, int idx);
int is_grid_molecule_behind_restrictive_boundary(struct grid_molecule *gm, struct wall *wall);
#endif
