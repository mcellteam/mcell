
#include "dump_state.h"

#include <iostream>
#include <string>
#include <cassert>

#define MAX_ARRAY_ITEMS 16
#define MAX_SUBVOLUMES 4
#define DUMP_ARRAY_NEWLINE_COUNT 8

#define IND_ADD2(ind) (std::string(ind) + "  ").c_str()

using namespace std;

std::ostream & operator<<(std::ostream &out, const timeval &a) {
    out << a.tv_sec << "s, " << a.tv_usec << "us";
    return out;
}

std::ostream & operator<<(std::ostream &out, const vector3 &a) {
    out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
    return out;
}

std::ostream & operator<<(std::ostream &out, const int3D &a) {
    out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
    return out;
}


std::ostream & operator<<(std::ostream &out, const pointer_hash &a) {
    //out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
	out << "(pointer_hash - TODO)";
  return out;
}


std::ostream & operator<<(std::ostream &out, const reaction_flags &a) {

#define DUMP_ITEM(name) #name << "=" << a.name << ", "
	out
		<< DUMP_ITEM(vol_vol_reaction_flag)
		<< DUMP_ITEM(vol_surf_reaction_flag)
		<< DUMP_ITEM(surf_surf_reaction_flag)
		<< DUMP_ITEM(vol_wall_reaction_flag)
		<< DUMP_ITEM(vol_vol_vol_reaction_flag)
		<< DUMP_ITEM(vol_vol_surf_reaction_flag)
		<< DUMP_ITEM(vol_surf_surf_reaction_flag)
		<< DUMP_ITEM(vol_surf_surf_reaction_flag)
	;
#undef DUMP_ITEM
  return out;
}

void dump_double_array(int num, const char* num_name, double* values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [double[]] \t\t" << comment;
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
  	if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
  		cout << "\n" << ind << "  ";
  	}
  	cout << i << ":" << values[i] << ", ";
  }
  if (num >= MAX_ARRAY_ITEMS) {
  	cout << "...";
  }
	cout << "\n";
}

void dump_vector3_array(int num, const char* num_name, vector3* values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [vector3[]] \t\t" << comment;
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
  	if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
  		cout << "\n" << ind << "  ";
  	}
  	cout << i << ":" << values[i] << ", ";
  }
  if (num >= MAX_ARRAY_ITEMS) {
  	cout << "...";
  }
  cout << "\n";
}

void dump_molecule_flags(short flags, const char* ind) {
	cout << ind << "flags: \t\t" << flags << " [short] \t\t /* Abstract Molecule Flags: Who am I, what am I doing, etc. */\n";

	cout << ind << "  ";
#define DUMP_FLAG(f, mask) if (f & mask != 0) cout << #mask << "(0x" << hex << mask << "), ";
	DUMP_FLAG(flags, TYPE_SURF)
	DUMP_FLAG(flags, TYPE_VOL)
	DUMP_FLAG(flags, ACT_DIFFUSE)
	DUMP_FLAG(flags, ACT_REACT)
	DUMP_FLAG(flags, ACT_NEWBIE)
	DUMP_FLAG(flags, ACT_CHANGE)
	DUMP_FLAG(flags, ACT_CLAMPED)
	DUMP_FLAG(flags, IN_SCHEDULE)
	DUMP_FLAG(flags, IN_SURFACE)
	DUMP_FLAG(flags, IN_VOLUME)
#undef DUMP_FLAG
	cout << "\n";
}

void dump_abstract_molecule(abstract_molecule* amp, const char* ind) {
  cout << ind << "next: *\t\t" << (void*)amp->next << " [abstract_molecule] \t\t /* Next molecule in scheduling queue */ //???\n";
  cout << ind << "t: \t\t" << amp->t << " [double] \t\t                      /* Scheduling time. */\n";
  cout << ind << "t2: \t\t" << amp->t2 << " [double] \t\t                     /* Time of next unimolecular reaction */\n";
  dump_molecule_flags(amp->flags, "  ");

  cout << ind << "properties: *\t\t" << (void*)amp->properties << " [species] \t\t    /* What type of molecule are we? */\n";
  cout << ind << "birthplace: *\t\t" << (void*)amp->birthplace << " [mem_helper] \t\t /* What was I allocated from? */\n";
  cout << ind << "birthday: \t\t" << amp->birthday << " [double] \t\t               /* Real time at which this particle was born */\n";
  cout << ind << "id: \t\t" << amp->id << " [u_long] \t\t                     /* unique identifier of this molecule */\n";
  cout << ind << "periodic_box: *\t\t" << (void*)amp->periodic_box << " [periodic_image] \t\t  /* track the periodic box a molecule is in */\n";

  cout << ind << "graph_data: *\t\t" << (void*)amp->graph_data << " [graph_data] \t\t /* nfsim graph structure data */\n";
#if 0 // nfsim function pointers, not important for now
  u_int (*get_flags)(void *); /* (nfsim) returns the reactivity flags associated with this particle */
  double (*get_diffusion)(void *);        /* (nfsim) returns the diffusion value */
  double (*get_time_step)(void *);        /* (nfsim) function pointer to a method that returns the time step */
  double (*get_space_step)(void *);       /* (nfsim) function pointer to a method that returns the space step */
  /* end structs used by the nfsim integration */
#endif
  cout << ind << "mesh_name: *\t\t" << amp->mesh_name << " [char] \t\t                /* Name of mesh that molecule is either in (volume molecule) or on (surface molecule) */\n";
}

void dump_volume_molecule(volume_molecule* amp, const char* ind) {
	cout << "***TODO!\n";
}

void dump_surface_molecule(surface_molecule* amp, const char* ind) {
	cout << "***TODO!\n";
}

void dump_molecules(int num_all_molecules, molecule_info **all_molecules) {
  cout << "all_molecules: **\t\t" << (void**)all_molecules << " [molecule_info] \n";
  cout << "num_all_molecules: \t\t" << num_all_molecules << " [int] \n";
	for (int i = 0; i < num_all_molecules; i++) {
		molecule_info *curr = all_molecules[i];
		assert(curr != NULL);

		abstract_molecule *amp = curr->molecule;
		assert(amp != NULL);

    if ((amp->properties->flags & NOT_FREE) == 0) {
      volume_molecule *vmp = (volume_molecule *)amp;
      dump_volume_molecule(vmp, "  ");
    } else if ((amp->properties->flags & ON_GRID) != 0) {
      surface_molecule *smp = (surface_molecule *)amp;
      dump_surface_molecule(smp, "  ");
    } else {
    	dump_abstract_molecule(amp, "  ");
    }

	  cout << "reg_names: *\t\t" << curr->reg_names << " [string_buffer] \t\t   /* Region names */\n";
	  cout << "mesh_names: *\t\t" << curr->mesh_names << " [string_buffer] \t\t  /* Mesh names that molec is nested in */\n";
	  cout << "pos: \t\t" << curr->pos << " [vector3] \t\t /* Position in space */\n";
	  cout << "orient: \t\t" << curr->orient << " [short] \t\t /* Which way do we point? */\n";
	}

}


void dump_region(region* reg, const char* ind) {
	// TODO
	cout << ind << "reg: *\t\t" << (void*)reg << " [region] \t\t\n";
  cout << ind << "  " << "sym: *\t\t" << (void*)reg->sym << " [sym_entry] \t\t  /* Symbol hash table entry for this region */\n";
  cout << ind << "  " << "hashval: \t\t" << reg->hashval << " [u_int] \t\t          /* Hash value for counter hash table */\n";
  cout << ind << "  " << "region_last_name: *\t\t" << reg->region_last_name << " [char] \t\t /* Name of region without prepended object name */\n";
  cout << ind << "  " << "parent: *\t\t" << (void*)reg->parent << " [object] \t\t  /* Parent of this region */\n";
  cout << ind << "  " << "element_list_head: *\t\t" << (void*)reg->element_list_head << " [element_list] \t\t /* List of element ranges comprising this region (used at parse time) */\n";
  cout << ind << "  " << "membership: *\t\t" << (void*)reg->membership << " [bit_array] \t\t /* Each bit indicates whether the corresponding wall is in the region */\n";
  cout << ind << "  " << "sm_dat_head: *\t\t" << (void*)reg->sm_dat_head << " [sm_dat] \t\t /* List of surface molecules to add to region */\n";
  cout << ind << "  " << "surf_class: *\t\t" << (void*)reg->surf_class << " [species] \t\t /* Surface class of this region */\n";
  cout << ind << "  " << "bbox: *\t\t" << (void*)reg->bbox << " [vector3] \t\t /* Array of length 2 to hold corners of region bounding box (used for release in region) */\n";
  cout << ind << "  " << "area: \t\t" << reg->area << " [double] \t\t          /* Area of region */\n";
  cout << ind << "  " << "flags: \t\t" << reg->flags << " [u_short] \t\t        /* Counting subset of Species Flags */\n";
  cout << ind << "  " << "manifold_flag: \t\t" << reg->manifold_flag << " [byte] \t\t   /* Manifold Flags: If IS_MANIFOLD, region is a closed manifold and thus defines a volume */\n";
  cout << ind << "  " << "volume: \t\t" << reg->volume << " [double] \t\t                   /* volume of region for closed manifolds */\n";
  cout << ind << "  " << "boundaries: *\t\t" << (void*)reg->boundaries << " [pointer_hash] \t\t /* hash table of edges that constitute external boundary of the region */\n";
  cout << ind << "  " << "region_has_all_elements: \t\t" << reg->region_has_all_elements << " [int] \t\t /* flag that tells whether the region contains ALL_ELEMENTS (effectively comprises the whole object) */\n";
}

void dump_region_list(region_list *rl, const char* regions_name, const char* comment, const char* ind) {
  cout << ind << regions_name << ": *\t\t" << rl << " [region_list] \t\t " << comment << "\n";
  for (struct region_list *r = rl; rl != NULL; rl = rl->next) {
  	dump_region(r->reg, ind);
  }
}


void dump_one_waypoint(waypoint* wp, const char* ind) {
	assert(wp != NULL);
	cout << ind << "waypoint: *\t\t" << (void*)wp << " [waypoint] \t\t\n";
  cout << ind << "  " << "loc: \t\t" << wp->loc << " [vector3] \t\t          /* This is where the waypoint is */\n";
  dump_region_list(wp->regions, "regions", "/* We are inside these regions */", IND_ADD2(ind));
  dump_region_list(wp->antiregions, "antiregions", "/* We are outside of (but hit) these regions */", IND_ADD2(ind));
}


void dump_waypoints(int n_waypoints, waypoint* waypoints) {
  cout << "n_waypoints: \t\t" << n_waypoints << " [int] \t\t/* How many waypoints (one per subvol) */\n";
  cout << "waypoints: *\t\t" << (void*)waypoints << " [waypoint] \t\t/* Waypoints contain fully-closed region information */\n";

	for (int i = 0; i < n_waypoints; i++) {
		dump_one_waypoint(&waypoints[i], "  ");
	}
}


void dump_one_subvolume(subvolume* sv, const char* ind) {
	assert(sv != NULL);
	cout << ind << "subvolume: *\t\t" << (void*)sv << " [subvolume] \t\t\n";
  cout << ind << "  " << "wall_head: *\t\t" << (void*)sv->wall_head << " [wall_list] \t\t /* Head of linked list of intersecting walls */\n";
  cout << ind << "  " << "mol_by_species: \t\t" << sv->mol_by_species << " [pointer_hash] \t\t /* table of speciesv->molecule list */\n";
  cout << ind << "  " << "species_head: *\t\t" << (void*)sv->species_head << " [per_species_list] \t\t\n";
  cout << ind << "  " << "mol_count: \t\t" << sv->mol_count << " [int] \t\t /* How many molecules are here? */\n";
  cout << ind << "  " << "llf: \t\t" << sv->llf << " [int3D] \t\t /* Indices of left lower front corner */\n";
  cout << ind << "  " << "urb: \t\t" << sv->urb << " [int3D] \t\t /* Indices of upper right back corner */\n";
  cout << ind << "  " << "world_edge: \t\t" << sv->world_edge << " [short] \t\t /* Direction Bit Flags that are set for SSVs at edge of world */\n";
  cout << ind << "  " << "local_storage: *\t\t" << (void*)sv->local_storage << " [storage] \t\t /* Local memory and scheduler */\n";
}


void dump_subvolumes(int n_subvols, const char* num_name, const char* num_comment, subvolume* subvols, const char* subvols_name, const char* name_comment, const char* ind) {
	cout << ind << num_name << ": \t\t" << n_subvols << " [int] \t\t" << num_comment << "\n";
	cout << ind << subvols_name << ": *\t\t" << (void*)subvols << " [subvolume] \t\t" << name_comment << "\n";

	for (int i = 0; i < n_subvols && i < MAX_SUBVOLUMES; i++) {
		dump_one_subvolume(&subvols[i], "  ");
	}
	if (n_subvols >= MAX_SUBVOLUMES) {
		cout << ind << "  " << "... (total " << n_subvols << ")\n";
	}
}

extern "C" void dump_volume(struct volume* s) {

  cout << "********* volume dump **************\n";

  cout << "bond_angle: \t\t" << s->bond_angle << " [double] \t\t/* Defines the default bond rotation angle between molecules. Default is 0. */\n";

  // These are only used with dynamic geometry
  cout << "dg_parse: *\t\t" << (void*)s->dg_parse << " [dyngeom_parse_vars] \n";
  cout << "dynamic_geometry_filename: *\t\t" << (void*)s->dynamic_geometry_filename << " [char] \n";
  dump_molecules(s->num_all_molecules, s->all_molecules);

  cout << "names_to_ignore: *\t\t" << (void*)s->names_to_ignore << " [string_buffer] \n";

  /* Coarse partitions are input by the user */
  /* They may also be generated automagically */
  /* They mark the positions of initial partition boundaries */
  dump_double_array(s->nx_parts, "nx_parts", s->x_partitions, "x_partitions", "/* Coarse X partition boundaries */", "");
  dump_double_array(s->ny_parts, "ny_parts", s->y_partitions, "y_partitions", "/* Coarse Y partition boundaries */", "");
  dump_double_array(s->nz_parts, "nz_parts", s->z_partitions, "z_partitions", "/* Coarse Z partition boundaries */", "");

  cout << "mem_part_x: \t\t" << s->mem_part_x << " [int] \t\t/* Granularity of memory-partition binning for the X-axis */\n";
  cout << "mem_part_y: \t\t" << s->mem_part_y << " [int] \t\t/* Granularity of memory-partition binning for the Y-axis */\n";
  cout << "mem_part_z: \t\t" << s->mem_part_z << " [int] \t\t/* Granularity of memory-partition binning for the Z-axis */\n";
  cout << "mem_part_pool: \t\t" << s->mem_part_pool << " [int] \t\t/* Scaling factor for sizes of memory pools in each storage */\n";

  /* Fine partitions are intended to allow subdivision of coarse partitions */
  /* Subdivision is not yet implemented */
  dump_double_array(s->n_fineparts, "n_fineparts", s->x_fineparts, "x_fineparts", "/* Fine X partition boundaries */", "");
  dump_double_array(s->n_fineparts, "n_fineparts", s->y_fineparts, "y_fineparts", "/* Fine Y partition boundaries */", "");
  dump_double_array(s->n_fineparts, "n_fineparts", s->z_fineparts, "z_fineparts", "/* Fine Z partition boundaries */", "");

  cout << "periodic_traditional: \t\t" << s->periodic_traditional << " [bool] \n";

  dump_waypoints(s->n_waypoints, s->waypoints);

  cout << "place_waypoints_flag: \t\t" << (int)s->place_waypoints_flag << " [byte] \t\t/* Used to save memory if waypoints not needed */\n";

  dump_subvolumes(s->n_subvols, "n_subvols", "/* How many coarse subvolumes? */", s->subvol, "subvol", "/* Array containing all subvolumes */", "");

  cout << "n_walls: \t\t" << s->n_walls << " [int] \t\t/* Total number of walls */\n";
  dump_vector3_array(s->n_walls, "n_walls", s->all_vertices, "all_vertices", "/* Central repository of vertices with a partial order imposed by natural ordering of storages*/", "");
return;
  /* Array of linked lists of walls using a vertex (has the size of
   * "all_vertices" array */
  cout << "walls_using_vertex: **\t\t" << (void**)s->walls_using_vertex << " [wall_list] \n";
  cout << "rx_hashsize: \t\t" << s->rx_hashsize << " [int] \t\t/* How many slots in our reaction hash table? */\n";
  cout << "n_reactions: \t\t" << s->n_reactions << " [int] \t\t/* How many reactions are there, total? */\n";
  cout << "reaction_hash: **\t\t" << (void**)s->reaction_hash << " [rxn] \t\t/* A hash table of all reactions. */\n";
  cout << "tv_rxn_mem: *\t\t" << (void*)s->tv_rxn_mem << " [mem_helper] \t\t/* Memory to store time-varying reactions */\n";

  cout << "count_hashmask: \t\t" << s->count_hashmask << " [int] \t\t/* Mask for looking up count hash table */\n";
  cout << "count_hash: **\t\t" << (void**)s->count_hash << " [counter] \t\t/* Count hash table */\n";
  cout << "count_scheduler: *\t\t" << (void*)s->count_scheduler << " [schedule_helper] // When to generate reaction output\n";
  cout << "counter_by_name: *\t\t" << (void*)s->counter_by_name << " [sym_table_head] \n";

  cout << "volume_output_scheduler: *\t\t" << (void*)s->volume_output_scheduler << " [schedule_helper] \t\t/* When to generate volume output\n";

  cout << "n_species: \t\t" << s->n_species << " [int] \t\t/* How many different species (molecules)? */\n";

  //NFSim stats
  cout << "n_NFSimSpecies: \t\t" << s->n_NFSimSpecies << " [int] \t\t/* number of graph patterns encountered during the NFSim simulation */\n";
  cout << "n_NFSimReactions: \t\t" << s->n_NFSimReactions << " [int] \t\t/* number of reaction rules discovered through an NFSim simulation */\n";
  cout << "n_NFSimPReactions: \t\t" << s->n_NFSimPReactions << " [int] \t\t/* number of potential reactions found in the reaction network (not necessarely triggered) */\n";

  cout << "species_list: **\t\t" << (void**)s->species_list << " [species] \t\t/* Array of all species (molecules). */\n";

  // This is used to skip over certain sections in the parser when using
  // dynamic geometries.
  cout << "dynamic_geometry_flag: \t\t" << s->dynamic_geometry_flag << " [int] \n";
  cout << "disable_polygon_objects: \t\t" << s->disable_polygon_objects << " [int] \n";

  // List of all the dynamic geometry events that need to be scheduled
  cout << "dynamic_geometry_head: *\t\t" << (void*)s->dynamic_geometry_head << " [dg_time_filename] \n";

  // Memory to store time and MDL names for dynamic geometry
  cout << "dynamic_geometry_events_mem: *\t\t" << (void*)s->dynamic_geometry_events_mem << " [mem_helper] \n";

  // Scheduler for dynamic geometry
  cout << "dynamic_geometry_scheduler: *\t\t" << (void*)s->dynamic_geometry_scheduler << " [schedule_helper] \n";
  cout << "releaser: *\t\t" << (void*)s->releaser << " [schedule_helper] \t\t/* Scheduler for release events */\n";

  cout << "storage_allocator: *\t\t" << (void*)s->storage_allocator << " [mem_helper] \t\t/* Memory for storage list */\n";
  cout << "storage_head: *\t\t" << (void*)s->storage_head << " [storage_list] \t\t/* Linked list of all local  memory/schedulers\n";

  cout << "current_mol_id: \t\t" << s->current_mol_id << " [u_long] \t\t/* next unique molecule id to use*/\n";

  cout << "speed_limit: \t\t" << s->speed_limit << " [double] // How far can the fastest particle get in one timestep?\n";

  cout << "fstream_sym_table: *\t\t" << (void*)s->fstream_sym_table << " [sym_table_head] \t\t/* Global MDL file stream symbol hash table */\n";

  cout << "var_sym_table: *\t\t" << (void*)s->var_sym_table << " [sym_table_head] \t\t/* Global MDL variables symbol hash table */\n";

  cout << "rxn_sym_table: *\t\t" << (void*)s->rxn_sym_table << " [sym_table_head] \t\t/* RXN symbol hash table */\n";
  cout << "obj_sym_table: *\t\t" << (void*)s->obj_sym_table << " [sym_table_head] \t\t/* Objects symbol hash table */\n";
  cout << "reg_sym_table: *\t\t" << (void*)s->reg_sym_table << " [sym_table_head] \t\t/* Regions symbol hash table */\n";
  cout << "mol_sym_table: *\t\t" << (void*)s->mol_sym_table << " [sym_table_head] \t\t/* Molecule type symbol hash table */\n";
  cout << "rpat_sym_table: *\t\t" << (void*)s->rpat_sym_table << " [sym_table_head] \t\t/* Release pattern hash table */\n";
  cout << "rxpn_sym_table: *\t\t" << (void*)s->rxpn_sym_table << " [sym_table_head] \t\t/* Named reaction pathway hash table */\n";
  cout << "mol_ss_sym_table: *\t\t" << (void*)s->mol_ss_sym_table << " [sym_table_head] \t\t/* Spatially structured molecule symbol hash table */\n";

  cout << "root_object: *\t\t" << (void*)s->root_object << " [object] \t\t/* Root of the object template tree */\n";
  cout << "root_instance: *\t\t" << (void*)s->root_instance << " [object] \t\t/* Root of the instantiated object tree */\n";
  cout << "periodic_box_obj: *\t\t" << (void*)s->periodic_box_obj << " [object] \n";

  cout << "default_release_pattern: *\t\t" << (void*)s->default_release_pattern << " [release_pattern] \t\t/* release once at t=0 */\n";

  cout << "volume_output_head: *\t\t" << (void*)s->volume_output_head << " [volume_output_item] \t\t/* List of all volume data output items */\n";


  cout << "output_block_head: *\t\t" << (void*)s->output_block_head << " [output_block] \t\t/* Global list of reaction data output blocks */\n";
  cout << "output_request_head: *\t\t" << (void*)s->output_request_head << " [output_request] \t\t/* Global list linking COUNT statements to internal variables \n";


  cout << "oexpr_mem: *\t\t" << (void*)s->oexpr_mem << " [mem_helper] \t\t/* Memory to store output_expressions */\n";
  cout << "outp_request_mem: *\t\t" << (void*)s->outp_request_mem << " [mem_helper] \t\t/* Memory to store output_requests */\n";
  cout << "counter_mem: *\t\t" << (void*)s->counter_mem << " [mem_helper] \t\t/* Memory to store counters (for counting molecules/reactions on regions)\n";

  cout << "trig_request_mem: *\t\t" << (void*)s->trig_request_mem << " [mem_helper] \t\t/* Memory to store listeners for trigger events */\n";
  cout << "magic_mem: *\t\t" << (void*)s->magic_mem << " [mem_helper] \t\t/* Memory used to store magic lists for reaction-triggered releases and such \n";

  cout << "elapsed_time: \t\t" << s->elapsed_time << " [double] \t\t/* Used for concentration measurement */\n";

  /* Visualization state */
  cout << "viz_blocks: *\t\t" << (void*)s->viz_blocks << " [viz_output_block] \t\t/* VIZ_OUTPUT blocks from file */\n";

  cout << "all_mols: *\t\t" << (void*)s->all_mols << " [species] \t\t/* Refers to ALL_MOLECULES keyword */\n";
  cout << "all_volume_mols: *\t\t" << (void*)s->all_volume_mols << " [species] // Refers to ALL_VOLUME_MOLECULES keyword\n";
  cout << "all_surface_mols: *\t\t" << (void*)s->all_surface_mols << " [species] // Refers to ALL_SURFACE_MOLECULES keyword\n";

  cout << "time_unit: \t\t" << s->time_unit << " [double] \t\t/* Duration of one global time step in real time */\n";
					/* Used to convert between real time and internal time */
  cout << "time_step_max: \t\t" << s->time_step_max << " [double] \t\t/* Maximum internal time that a molecule may diffuse */\n";

  cout << "grid_density: \t\t" << s->grid_density << " [[double] \t\t/* Density of grid for surface molecules, number per um^2 */";
  cout << "length_unit: \t\t" << s->length_unit << " [double] \t\t/* Internal unit of distance, 1/sqrt(grid_density), in microns\n";

  cout << "r_length_unit: \t\t" << s->r_length_unit << " [double] \t\t/* Reciprocal of length_unit to avoid division */\n";
  cout << "rx_radius_3d: \t\t" << s->rx_radius_3d << " [double] \t\t/* Interaction radius for reactions between volume molecules \n";


  cout << "space_step: \t\t" << s->space_step << " [double] \t\t/* User-supplied desired average diffusion distance for volume molecules\n";

  cout << "r_step: *\t\t" << (void*)s->r_step << " [double] \t\t/* Lookup table of 3D diffusion step lengths */\n";
  cout << "d_step: *\t\t" << (void*)s->d_step << " [double] \t\t/* Lookup table of 3D diffusion direction vectors */\n";
  cout << "r_step_surface: *\t\t" << (void*)s->r_step_surface << " [double] \t\t/* Lookup table of 2D diffusion step lengths */\n";
  cout << "r_step_release: *\t\t" << (void*)s->r_step_release << " [double] \t\t/* Lookup table of diffusion lengths for 3D release */\n";
  cout << "radial_subdivisions: \t\t" << s->radial_subdivisions << " [u_int] \t\t/* Size of 2D and 3D step length lookup tables */\n";
  cout << "radial_directions: \t\t" << s->radial_directions << " [u_int] \t\t/* Requested size of 3D direction lookup table */\n";
  cout << "num_directions: \t\t" << s->num_directions << " [u_int] \t\t/* Actual size of 3D direction lookup table */\n";
  cout << "directions_mask: \t\t" << s->directions_mask << " [int] \t\t/* Mask to obtain RNG bits for direction lookup */\n";
  cout << "fully_random: \t\t" << s->fully_random << " [int] \t\t/* If set, generate directions with trig functions instead of lookup table \n";

  cout << "dissociation_index: \t\t" << s->dissociation_index << " [int] \t\t/* Used to keep 3D products from reacting with each  other too soon \n";

  cout << "chkpt_iterations: \t\t" << s->chkpt_iterations << " [long long] \t\t/* Number of iterations to advance before checkpointing\n";

  cout << "chkpt_init: \t\t" << s->chkpt_init << " [u_int] \t\t/* Set if this is the initial run of a simulation with no previous checkpoints\n";

  cout << "chkpt_flag: \t\t" << s->chkpt_flag << " [u_int] \t\t/* Set if there are any CHECKPOINT statements in mdl file */\n";
  cout << "chkpt_seq_num: \t\t" << s->chkpt_seq_num << " [u_int] \t\t/* Number of current run in checkpoint sequence */\n";
  cout << "keep_chkpts: \t\t" << s->keep_chkpts << " [int] \t\t/* flag to indicate if checkpoints should be kept */\n";

  cout << "chkpt_infile: *\t\t" << (void*)s->chkpt_infile << " [char] \t\t/* Name of checkpoint file to read from */\n";
  cout << "chkpt_outfile: *\t\t" << (void*)s->chkpt_outfile << " [char] \t\t/* Name of checkpoint file to write to */\n";
  cout << "chkpt_byte_order_mismatch: \t\t" << s->chkpt_byte_order_mismatch << " [u_int] \t\t/* Flag that defines whether mismatch in byte order exists between the saved checkpoint file and the machine reading it\n";
  cout << "chkpt_start_time_seconds: \t\t" << s->chkpt_start_time_seconds << " [double] \t\t/* start of the simulation time (in sec) for new checkpoint\n";

  cout << "current_time_seconds: \t\t" << s->current_time_seconds << " [double] \t\t/* current simulation time in seconds */\n";
  /* simulation start time (in seconds) or time of most recent checkpoint */
  cout << "simulation_start_seconds: \t\t" << s->simulation_start_seconds << " [double] \n";

  cout << "diffusion_number: \t\t" << s->diffusion_number << " [long long] \t\t/* Total number of times molecules have had their positions updated */\n";

  cout << "diffusion_cumtime: \t\t" << s->diffusion_cumtime << " [double] \t\t/* Total time spent diffusing by all molecules */\n";
  cout << "ray_voxel_tests: \t\t" << s->ray_voxel_tests << " [long long] \t\t/* How many ray-subvolume intersection tests have we performed */\n";

  cout << "ray_polygon_tests: \t\t" << s->ray_polygon_tests << " [long long] \t\t/* How many ray-polygon intersection tests have  we performed */\n";

  cout << "ray_polygon_colls: \t\t" << s->ray_polygon_colls << " [long long] \t\t/* How many ray-polygon intersections have occured */\n";

  cout << "dyngeom_molec_displacements: \t\t" << s->dyngeom_molec_displacements << " [long long] \t\t/* Total number of dynamic geometry molecule displacements */\n";

  /* below "vol" means volume molecule, "surf" means surface molecule */
  cout << "vol_vol_colls: \t\t" << s->vol_vol_colls << " [long long] \t\t/* How many vol-vol collisions have occured */\n";
  cout << "vol_surf_colls: \t\t" << s->vol_surf_colls << " [long long] \t\t/* How many vol-surf collisions have occured */\n";
  cout << "surf_surf_colls: \t\t" << s->surf_surf_colls << " [long long] \t\t/* How many surf-surf collisions have occured */\n";
  cout << "vol_wall_colls: \t\t" << s->vol_wall_colls << " [long long] \t\t/* How many vol-wall collisions have occured */\n";
  cout << "vol_vol_vol_colls: \t\t" << s->vol_vol_vol_colls << " [long long] // How many vol-vol-vol collisions have occured\n";
  cout << "vol_vol_surf_colls: \t\t" << s->vol_vol_surf_colls << " [long long] \t\t/* How many vol-vol-surf collisions have occured */\n";
  cout << "vol_surf_surf_colls: \t\t" << s->vol_surf_surf_colls << " [long long] \t\t/* How many vol-surf-surf collisions have occured */\n";

  cout << "surf_surf_surf_colls: \t\t" << s->surf_surf_surf_colls << " [long long] \t\t/* How many surf-surf-surf collisions have occured */\n";


  cout << "bb_llf: \t\t" << s->bb_llf << " [vector3] \t\t/* llf corner of world bounding box */\n";
  cout << "bb_urb: \t\t" << s->bb_urb << " [vector3] \t\t/* urb corner of world bounding box */\n";

  cout << "rng: *\t\t" << (void*)s->rng << " [rng_state] \t\t/* State of the random number generator (currently isaac64) */\n";

  cout << "init_seed: \t\t" << s->init_seed << " [u_int] \t\t/* Initial seed value for random number generator */\n";

  cout << "current_iterations: \t\t" << s->current_iterations << " [long long] \t\t/* How many iterations have been run so far */\n";
  // Starting iteration number for current run or iteration of most recent
  // checkpoint
  cout << "start_iterations: \t\t" << s->start_iterations << " [long long] \n";
  cout << "last_timing_time: \t\t" << s->last_timing_time << " [timeval] \t\t/* time and iteration of last timing event */\n";
  cout << "last_timing_iteration: \t\t" << s->last_timing_iteration << " [long long] \t\t/* during the main run_iteration loop */\n";

  cout << "procnum: \t\t" << s->procnum << " [int] \t\t/* Processor number for a parallel run */\n";
  cout << "quiet_flag: \t\t" << s->quiet_flag << " [int] \t\t/* Quiet mode */\n";
  cout << "with_checks_flag: \t\t" << s->with_checks_flag << " [int] \t\t/* Check geometry for overlapped walls? */\n";

  cout << "coll_mem: *\t\t" << (void*)s->coll_mem << " [mem_helper] \t\t/* Collision list */\n";
  cout << "sp_coll_mem: *\t\t" << (void*)s->sp_coll_mem << " [mem_helper] \t\t/* Collision list (trimol) */\n";
  cout << "tri_coll_mem: *\t\t" << (void*)s->tri_coll_mem << " [mem_helper] \t\t/* Collision list (trimol) */\n";
  cout << "exdv_mem: *\t\t" << (void*)s->exdv_mem << " [mem_helper] // Vertex lists f    ct interaction disk area\n";

  /* Current version number. Format is "3.XX.YY" where XX is major release
   * number (for new features) and YY is minor release number (for patches) */
  cout << "mcell_version: *\t\t" << (void*)s->mcell_version << " [const char] \n";

  cout << "use_expanded_list: \t\t" << s->use_expanded_list << " [int] \t\t/* If set, check neighboring subvolumes for mol-mol interactions */\n";

  cout << "randomize_smol_pos: \t\t" << s->randomize_smol_pos << " [int] \t\t/* If set, always place surface molecule at random location instead of center of grid */\n";

  cout << "vacancy_search_dist2: \t\t" << s->vacancy_search_dist2 << " [double] \t\t/* Square of distance to search for free grid  location to place surface product */\n";

  cout << "surface_reversibility: \t\t" << s->surface_reversibility << " [byte] \t\t/* If set, match unbinding diffusion distribution  to binding distribution at surface */\n";

  cout << "volume_reversibility: \t\t" << s->volume_reversibility << " [byte] \t\t/* If set, match unbinding diffusion distribution to binding distribution in volume */\n";


  /* If set to NEAREST_TRIANGLE, molecules are moved to a random location
   * slightly offset from the enclosing wall. If set to NEAREST_POINT, then
   * they are moved to the closest point on that wall (still slightly offset).
   * */
  cout << "dynamic_geometry_molecule_placement: \t\t" << s->dynamic_geometry_molecule_placement << " [int] \n";

  /* MCell startup command line arguments */
  cout << "seed_seq: \t\t" << s->seed_seq << " [u_int] \t\t/* Seed for random number generator */\n";
  cout << "iterations: \t\t" << s->iterations << " [long long] \t\t/* How many iterations to run */\n";
  cout << "log_freq: \t\t" << s->log_freq << " [unsigned long] \t\t/* Interval between simulation progress reports, default scales as sqrt(iterations) */\n";

  cout << "mdl_infile_name: *\t\t" << (void*)s->mdl_infile_name << " [char] \t\t/* Name of MDL file specified on command line */\n";
  cout << "curr_file: *\t\t" << (void*)s->curr_file << " [const char ] \t\t/* Name of MDL file currently being parsed */\n";

  // XXX: Why do we allocate this on the heap rather than including it inline?
  cout << "notify: *\t\t" << (void*)s->notify << " [notifications] \t\t/* Notification/warning/output flags */\n";

  cout << "clamp_list: *\t\t" << (void*)s->clamp_list << " [ccn_clamp_data] \t\t/* List of objects at which volume molecule concentrations should be clamped */\n";



  /* Flags for asynchronously-triggered checkpoints */

  /* Flag indicating whether a checkpoint has been requested. */
  cout << "checkpoint_requested: \t\t" << s->checkpoint_requested << " [checkpoint_request_type_t] \n";
  cout << "checkpoint_alarm_time: \t\t" << s->checkpoint_alarm_time << " [unsigned int] \t\t// number of seconds between checkpoints\n";
  cout << "continue_after_checkpoint: \t\t" << s->continue_after_checkpoint << " [int] \t\t/* 0: exit after chkpt, 1: continue after chkpt */\n";
  cout << "last_checkpoint_iteration: \t\t" << s->last_checkpoint_iteration << " [long long] \t\t/* Last iteration when chkpt was created */\n";
  cout << "begin_timestamp: \t\t" << s->begin_timestamp << " [time_t] \t\t/* Time since epoch at beginning of 'main' */\n";
  cout << "initialization_state: *\t\t" << (void*)s->initialization_state << " [char] \t\t/* NULL after initialization completes */\n";
  cout << "rxn_flags: \t\t" << s->rxn_flags << " [reaction_flags] \n";
  /* shared walls information per mesh vertex is created when there are
	 reactions present with more than one surface reactant or more than one
	 surface product */
  cout << "create_shared_walls_info_flag: \t\t" << s->create_shared_walls_info_flag << " [int] \n";
  /* resource usage during initialization */
	cout << "u_init_time: \t\t" << s->u_init_time << " [timeval] \t\t/* user time */\n";
  cout << "s_init_time: \t\t" << s->s_init_time << " [timeval] \t\t/* system time */\n";
  cout << "t_start: \t\t" << s->t_start << " [time_t] \t\t/* global start time */\n";
  cout << "reaction_prob_limit_flag: \t\t" << s->reaction_prob_limit_flag << " [byte] \t\t/* checks whether there is at least one reaction with probability greater than 1 including variable rate reactions */\n";



  //JJT: Checks if we will be communicating with nfsim
  cout << "nfsim_flag: \t\t" << s->nfsim_flag << " [int] \n";
  cout << "global_nfsim_volume: *\t\t" << (void*)s->global_nfsim_volume << " [species] \n";
  cout << "global_nfsim_surface: *\t\t" << (void*)s->global_nfsim_surface << " [species] \n";

  cout << "species_mesh_transp: *\t\t" << (void*)s->species_mesh_transp << " [pointer_hash] \n";
}
