struct mesh_region_string_buffs {
  struct string_buffer *old_inst_mesh_names;
  struct string_buffer *old_region_names;
};

int mcell_do_dg(struct volume *state, struct poly_object *poly_obj);

/*int mcell_update_geometry(struct volume *state);*/
/*int mcell_destroy_everything(*/
/*    struct volume *state,*/
/*    struct mesh_region_string_buffs *string_buffs);*/

/*int mcell_reinitialize(*/
/*    struct volume *state,*/
/*    struct mesh_region_string_buffs *string_buffs);*/
