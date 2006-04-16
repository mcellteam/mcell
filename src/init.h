#ifndef INIT_H
#define INIT_H

#include "mcell_structs.h"
#include "mdlparse.h"

void init_credits(void);

int init_notifications();
int init_sim(void);

int init_species(void);
int init_geom(void);
int init_partitions(void);

int instance_obj(struct object *objp,
		 double (*im)[4],
		 struct viz_obj *vizp,
		 struct lig_count_ref *lcrp,
		 char *sub_name);

int instance_release_site(struct object *objp,
			  double (*im)[4]);

int instance_polygon_object(struct object *objp,
		double (*im)[4],
		struct viz_obj *vizp,
		struct lig_count_ref *obj_lcrp,
		char *full_name);

int init_regions();
void init_clamp_lists();

int instance_obj_regions(struct object *objp, char *sub_name);

int init_wall_regions(struct object *objp, char *full_name);

int init_effectors();
int instance_obj_effectors(struct object *objp);
int init_wall_effectors(struct object *objp);

int init_effectors_by_density(struct wall *w, struct eff_dat *eff_dat_head);

int init_effectors_by_number(struct object *objp, struct region_list *rlp);

int compute_bb(struct object *objp,
	       double (*im)[4],
	       char *sub_name);

int compute_bb_release_site(struct object *objp,
			    double (*im)[4]);

int compute_bb_polygon_object(struct object *objp,
		double (*im)[4],
		char *full_name);

void cube_corners(struct vector3 *p1,
                  struct vector3 *p2,
                  struct vector3 *corner);

void cube_face(struct vector3 *corner, struct vector3 **face, int i);

void cube_faces(struct vector3 *corner, struct vector3 *(*face)[4]);

void swap_double(double *x, double *y);


/*
int init_nodes(void);

struct wall *init_wall(struct vector3 **face, struct vector3 **face_vertex_normal, int n_verts);

int init_region_counter(struct polygon_object *pop, struct cmprt_data *cdp, struct region_list *rlp);

int init_counter_hash_table (struct object *objp, struct region_list *region_list);

struct region_list *init_region_list_for_wall(struct region_list *rlp,int index);

int init_effector_grid(struct wall *wp);

int init_effectors_by_density(struct wall *wp, struct eff_dat *effdp_head);

int init_effector_table(struct wall *wp);

int init_rx_table(struct parent_rx *prxp);

int init_release_event_table(struct release_event_queue *reqp);

void destroy_sym_value(struct sym_table *sym);
*/

int init_releases();

#endif
