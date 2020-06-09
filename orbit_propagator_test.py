import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d.axes3d import Axes3D
plt.style.use('dark_background')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

cb = pd.earth

if __name__ == '__main__':
  # initial conditions of orbit parameters
  r_mag = cb['radius'] + 1500 # km -- Choose orbit height
  v_mag = np.sqrt(cb['mu'] / r_mag)

  # initial position and velocity vectors
  r0 = [0,r_mag,0] # initial position
  v0 = [0,0,v_mag] # intial velocity vector

  print(r0)
  print(v0)

  # timespan
  tspan = 2 * 24 * 60 * 60.0 # minutes -- Choose simulation timespan

  # timestep
  dt = 25.0

  op = OP(r0, v0, tspan, dt, cb=cb)
  op.propagate_orbit()
  op.plot_3d(show_plot=True)
