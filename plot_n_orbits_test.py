import numpy as np
from math import sqrt

import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

# tspan = 3600 * 24 * 1.0
tspan = 60 * 60 * 24 * 1.0
dt = 20.0

if __name__ == '__main__':
  r_mag = pd.earth['radius'] + 400
  v_mag = np.sqrt(pd.earth['mu'] / r_mag)
  r0 = [r_mag, 0, 0]
  v0 = [0, v_mag, 0]

  r_mag = pd.earth['radius'] + 1000
  v_mag = sqrt(pd.earth['mu'] / r_mag) * 1.3
  r00 = [r_mag, 0, 0]
  v00 = [0, v_mag, 0]

  op0 = OP(r0, v0, tspan, dt)
  op00 = OP(r00, v00, tspan, dt)

  op0.propagate_orbit()
  op00.propagate_orbit()

  t.plot_n_orbits([op0.rs, op00.rs], labels=['0', '1'], cb=pd.earth, show_plot=True)