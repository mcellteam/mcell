#!/usr/bin/env python3

import mcell as m

s = m.Species('a', diffusion_constant_2d=1, diffusion_constant_3d=2)

# should fail with ValueError: Only one of fields diffusion_constant_2d or diffusion_constant_3d can be set for simple species