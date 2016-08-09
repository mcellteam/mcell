import mcellSwig as m

world = m.mcell_create()

m.mcell_init_state(world)

m.mcell_set_time_step(world,1e-6)
m.mcell_set_iterations(world,100)

sp = m.mcell_species_spec()
sp.name = "x"
sp.D = 1e-6
sp.is_2d = 0
sp.custom_time_step = 0
sp.target_only = 0
sp.max_step_length = 0

sym = m.mcell_symbol;

m.mcell_create_species(world, sp, sym)

m.mcell_init_simulation(world)

m.mcell_run_simulation(world)
