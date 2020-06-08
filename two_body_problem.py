import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d.axes3d import Axes3D
plt.style.use('dark_background')

def plot(r):
  # 3D plot
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_subplot(111, projection='3d')

  # plot trajectory and starting point
  ax.plot(r[:,0], r[:,1], r[:,2], 'b', label='Trajectory')
  ax.plot([r[0,0]], [r[0,1]], [r[0,2]], 'bo', label='Starting Position')

  r_plot = earth_radius

  # plot earth
  _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
  _x = r_plot * np.cos(_u) * np.sin(_v)
  _y = r_plot * np.sin(_u) * np.sin(_v)
  _z = r_plot * np.cos(_v)
  ax.plot_surface(_x, _y, _z, cmap='Greens')

  # plot X, Y, Z vectors (arrows)
  l = r_plot * 2.0
  x, y, z = [[0,0,0], [0,0,0], [0,0,0]] # where arrows start
  u, v, w = [[l,0,0], [0,l,0], [0,0,l]] # where arrows end
  ax.quiver(x, y, z, u, v, w, color='y')

  # check for custom axes limits
  max_val = np.max(np.abs(r))

  # set labels and title
  ax.set_xlim([-max_val, max_val])
  ax.set_ylim([-max_val, max_val])
  ax.set_zlim([-max_val, max_val])
  ax.set_xlabel('X (km)')
  ax.set_ylabel('Y (km)')
  ax.set_zlabel('Z (km)')
  ax.set_aspect('equal')
  #ax.plot_title()
  plt.legend()

  plt.show()

earth_radius = 6378.0 # km
earth_mu = 398600.0 # km^3 / seconds^2

def diffy_q(t, y, mu):
  # unpack state
  rx, ry, rz, vx, vy, vz = y
  r = np.array([rx, ry, rz])

  # norm of the radius vector
  norm_r = np.linalg.norm(r)

  # two body acceleration
  ax, ay, az = -r * mu / norm_r**3

  return [vx,vy,vz,ax,ay,az]

if __name__ == '__main__':
  # initial conditions of orbit parameters
  r_mag = earth_radius + 1000 # km -- Choose orbit height
  v_mag = np.sqrt(earth_mu / r_mag)

  # initial position and velocity vectors
  r0 = [r_mag,0,0]
  v0 = [0,v_mag,0]

  # timespan
  tspan = 106 * 60.0 # minutes -- Choose simulation timespan

  # timestep
  dt = 25.0

  # total number of steps
  n_steps = int(np.ceil(tspan / dt))

  # initialize arrays
  ys = np.zeros((n_steps, 6))
  ts = np.zeros((n_steps, 1))

  # initial conditions
  y0 = r0 + v0
  ys[0] = np.array(y0)
  step = 1

  # initiate solver
  solver = ode(diffy_q)
  solver.set_integrator('lsoda')
  solver.set_initial_value(y0, 0)
  solver.set_f_params(earth_mu)

  # propagate orbit
  while solver.successful() and step < n_steps:
    solver.integrate(solver.t + dt)
    ts[step] = solver.t
    ys[step] = solver.y
    step += 1
  
  rs = ys[:,:3]

  plot(rs)
