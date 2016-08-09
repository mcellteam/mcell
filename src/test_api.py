import mcellSwig as m

world = m.mcell_create()

m.mcell_init_state(world)

m.mcell_set_time_step(world,1e-6)
m.mcell_set_iterations(world,100)

species_def = m.mcell_species_spec()
species_def.name = "x"
species_def.D = 1e-6
species_def.is_2d = 0
species_def.custom_time_step = 0
species_def.target_only = 0
species_def.max_step_length = 0

sym = m.mcell_symbol();
mcell_sym = m.mcell_create_species(world, species_def, sym)
m.mcell_print_name(mcell_sym)

scene = m.object()
m.mcell_create_instance_object(world, "Scene", scene)

m.mcell_init_simulation(world)
m.mcell_run_simulation(world)
