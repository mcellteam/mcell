import mcellSwig as m

world = m.mcell_create()
print(world)
print(dir(world))

m.mcell_init_state(world)

m.mcell_set_time_step(world,1e-6)
m.mcell_set_iterations(world,100)

m.mcell_init_simulation(world)

m.mcell_run_simulation(world)



