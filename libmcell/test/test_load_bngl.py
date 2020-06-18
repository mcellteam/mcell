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

subs = m.Subsystem()

subs.load_bngl_molecule_types_and_reaction_rules('test.bngl')

print(subs)

