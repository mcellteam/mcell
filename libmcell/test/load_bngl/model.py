#!/usr/bin/env python3

import sys
import os

MCELL_DIR = os.environ.get('MCELL_DIR', '')
if MCELL_DIR:
    sys.path.append(os.path.join(MCELL_DIR, 'lib'))
else:
    print("Error: variable MCELL_DIR that is used to find the mcell library was not set.")
    sys.exit(1)

import mcell as m

from parameters import *

if len(sys.argv) == 3 and sys.argv[1] == '-seed':
    # overwrite value SEED defined in module parameters
    SEED = int(sys.argv[2])


subsystem = m.Subsystem()

subsystem.load_bngl_molecule_types_and_reaction_rules('test.bngl')

#print(subsystem)

print(subsystem.find_elementary_molecule_type('A'))

cplx_inst_A = m.ComplexInstance([subsystem.find_elementary_molecule_type('A').inst()])

import geometry


rel_a = m.ReleaseSite(
    name = 'rel_a',
    complex_instance = cplx_inst_A,
    region = geometry.Cube,
    location = m.Vec3(0, 0, 0),
    site_diameter = 0,
    number_to_release = 10
)

#print(rel_a)

#print(subsystem)

#sys.exit(0)

# ---- surface classes assignment ----


# ---- instantiation data ----

instantiation = m.InstantiationData()
instantiation.add_geometry_object(geometry.Cube)
instantiation.add_release_site(rel_a)


viz_output = m.VizOutput(
    mode = m.VizMode.ASCII,
    filename_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene',
    every_n_timesteps = 1
)

observables = m.Observables()
observables.add_viz_output(viz_output)


model = m.Model()

# ---- configuration ----

model.config.time_step = TIME_STEP
model.config.seed = SEED
model.config.total_iterations_hint = ITERATIONS

model.config.partition_dimension = 20.02
model.config.subpartition_dimension = 2.002

# ---- default configuration overrides ----


# ---- add components ----

model.add_subsystem(subsystem)
model.add_instantiation_data(instantiation)
model.add_observables(observables)

# ---- initialization and execution ----

print(model)

model.initialize()

if DUMP:
    model.dump_internal_state()

if EXPORT_DATA_MODEL and model.viz_outputs:
    model.export_data_model()

model.run_iterations(ITERATIONS)
model.end_simulation()
