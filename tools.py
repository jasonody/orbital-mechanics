import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

import planetary_data as pd

d2r = np.pi / 180.0

def plot_n_orbits(rs, labels, cb=pd.earth, show_plot=False, save_plot=False, title='Many Orbits'):
    # 3D plot
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111, projection='3d')

    # plot trajectory and starting point
    n = 0
    for r in rs:
      ax.plot(r[:,0], r[:,1], r[:,2], label=labels[n])
      ax.plot([r[0,0]], [r[0,1]], [r[0,2]])
      n+=1

    # plot earth
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = cb['radius'] * np.cos(_u) * np.sin(_v)
    _y = cb['radius'] * np.sin(_u) * np.sin(_v)
    _z = cb['radius'] * np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap='Greens')

    # plot X, Y, Z vectors (arrows)
    l = cb['radius'] * 2.0
    x, y, z = [[0,0,0], [0,0,0], [0,0,0]] # where arrows start
    u, v, w = [[l,0,0], [0,l,0], [0,0,l]] # where arrows end
    ax.quiver(x, y, z, u, v, w, color='y')

    # check for custom axes limits
    max_val = np.max(np.abs(rs))

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