import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

import planetary_data as pd
import tools as t

def null_perts(): # no perturbations
  return {
    'J2': False,
    'aero': False,
    'moon_grav': False,
    'solar_gravity': False
  }

class OrbitPropagator:
  def __init__(self, state0, tspan, dt, coes=False, deg=True, cb=pd.earth, perts=null_perts()): # dt = timestep, cb = central body
    if coes:
      self.r0, self.v0, date = t.coes2rv(state0, deg=deg, mu=cb['mu'])
    else:
      self.r0 = state0[:3]
      self.v0 = state0[3:]

    self.y0 = list(self.r0) + list(self.v0)
    self.tspan = tspan
    self.dt = dt
    self.cb = cb

    # total number of steps
    self.n_steps = int(np.ceil(self.tspan / self.dt))

    # initialize variables
    self.ts = np.zeros((self.n_steps, 1))
    self.ys = np.zeros((self.n_steps, 6))
    self.ts[0] = 0
    self.ys[0,:] = self.y0
    self.step = 1

    # initiate solver
    self.solver = ode(self.diffy_q)
    self.solver.set_integrator('lsoda')
    self.solver.set_initial_value(self.y0, 0)

    # define perturbations dictionary
    self.perts = perts

    self.propagate_orbit()
  
  def propagate_orbit(self):
    # propagate orbit
    while self.solver.successful() and self.step < self.n_steps:
      self.solver.integrate(self.solver.t + self.dt)
      self.ts[self.step] = self.solver.t
      self.ys[self.step] = self.solver.y
      self.step += 1
    
    self.rs = self.ys[:,:3]
    self.vs = self.ys[:,3:]
  
  def diffy_q(self, t, y):
    # unpack state
    rx, ry, rz, vx, vy, vz = y
    r = np.array([rx, ry, rz])

    # norm of the radius vector
    norm_r = np.linalg.norm(r)

    # two body acceleration
    a = -r * self.cb['mu'] / norm_r**3

    # J2 perturbation
    if self.perts['J2']:
      z2 = r[2]**2
      r2 = norm_r**2
      tx = r[0] / norm_r * (5 * z2 / r2 - 1)
      ty = r[1] / norm_r * (5 * z2 / r2 - 1)
      tz = r[2] / norm_r * (5 * z2 / r2 - 3)

      a_j2 = 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2 / norm_r**4 * np.array([tx,ty,tz])

      a += a_j2

    return [vx,vy,vz,a[0],a[1],a[2]]
  
  def calculate_coes(self, degrees=True):
    print('Calculating COEs...')

    self.coes = np.zeros((self.n_steps, 6))

    for n in range(self.n_steps):
      self.coes[n,:] = t.rv2coes(self.rs[n,:], self.vs[n,:], mu=self.cb['mu'], degrees=degrees)

  def plot_3d(self, show_plot=False, save_plot=False, title='Some Plot'):
    # 3D plot
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111, projection='3d')

    # plot trajectory and starting point
    ax.plot(self.rs[:,0], self.rs[:,1], self.rs[:,2], 'b', label='Trajectory')
    ax.plot([self.rs[0,0]], [self.rs[0,1]], [self.rs[0,2]], 'bo', label='Starting Position')

    # plot earth
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = self.cb['radius'] * np.cos(_u) * np.sin(_v)
    _y = self.cb['radius'] * np.sin(_u) * np.sin(_v)
    _z = self.cb['radius'] * np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap='Greens')

    # plot X, Y, Z vectors (arrows)
    l = self.cb['radius'] * 2.0
    x, y, z = [[0,0,0], [0,0,0], [0,0,0]] # where arrows start
    u, v, w = [[l,0,0], [0,l,0], [0,0,l]] # where arrows end
    ax.quiver(x, y, z, u, v, w, color='y')

    # check for custom axes limits
    max_val = np.max(np.abs(self.rs))

    # set labels and title
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('equal')
    
    ax.set_title(title)
    plt.legend()

    if show_plot:
      plt.show()
    if save_plot:
      plt.savefig(title.replace(" ", "_") + '.png', dpi=300)