struct mesh_region_string_buffs {
  struct string_buffer *old_inst_mesh_names;
  struct string_buffer *old_region_names;
};

int mcell_change_geometry(struct volume *state, struct poly_object_list *pobj_list);
