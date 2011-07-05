#ifndef INIT_H
#define INIT_H

#include "mcell_structs.h"

int init_notifications(void);
int init_sim(void);

int init_species(void);
int init_geom(void);
int init_partitions(void);

int instance_obj(struct object *objp,
		 double (*im)[4]);

int instance_release_site(struct object *objp,
			  double (*im)[4]);

int instance_polygon_object(struct object *objp,
		double (*im)[4]);

int init_regions(void);
void init_clamp_lists(void);

int instance_obj_regions(struct object *objp);

int init_wall_regions(struct object *objp);

int init_effectors(void);
int instance_obj_effectors(struct object *objp);
int init_wall_effectors(struct object *objp);

int init_effectors_by_density(struct wall *w, struct eff_dat *eff_dat_head);

int init_effectors_by_number(struct object *objp, struct region_list *rlp);

void cube_corners(struct vector3 *p1,
                  struct vector3 *p2,
                  struct vector3 *corner);

void cube_face(struct vector3 *corner, struct vector3 **face, int i);

void cube_faces(struct vector3 *corner, struct vector3 *(*face)[4]);

void swap_double(double *x, double *y);

int init_releases(void);

void publish_special_reactions_report(struct species *sp, struct name_list *vol_species_name_list, struct name_list *surf_species_name_list);

int accumulate_vertex_counts_per_storage(struct object *objp, int *num_vertices_this_storage, double (*im)[4]);

int accumulate_vertex_counts_per_storage_polygon_object(struct object *objp, int *num_vertices_this_storage, double (*im)[4]);

int which_storage_contains_vertex(struct vector3 *v); 

int fill_world_vertices_array(struct object *objp, int *num_vertices_this_storage, double (*im)[4]);
int fill_world_vertices_array_polygon_object(struct object *objp, int *num_vertices_this_storage, double (*im)[4]);
void check_for_conflicting_surface_classes(struct wall *w);	
void check_for_conflicts_in_surface_class(struct species *sp);
struct species * get_species_by_name(char *name); 
void create_name_lists_of_volume_and_surface_mols(struct name_list **vol_species_name_list, struct name_list **surf_species_name_list);        			
void remove_molecules_name_list(struct name_list **nlist);    
int check_for_overlapped_walls(void);
    			
#endif
