import math as m
import datetime
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
  
def coes2rv(coes, deg=False, mu=pd.earth['mu']):
  if deg:
    # conver to radians
    a, e, i, ta, aop, raan, date = coes
    i *= d2r
    ta *= d2r
    aop *= d2r
    raan *= d2r
  else:
    a, e, i, ta, aop, raan, date = coes
  
  E = ecc_anomonly([ta, e], 'tae')

  r_norm = a * (1 - e**2) / (1 + e * np.cos(ta))

  # calculate r and v vectors in perifocal fram
  r_perif = r_norm * np.array([m.cos(ta), m.sin(ta), 0])
  v_perif = m.sqrt(mu * a) / r_norm * np.array([-m.sin(E), m.cos(E) * m.sqrt(1 - e**2), 0])

  # rotation matrix from perifocal to ECI
  perif2eci = np.transpose(eci2perif(raan, aop, i))

  # calculate r and v vectors in inertial framea
  r = np.dot(perif2eci, r_perif)
  v = np.dot(perif2eci, v_perif)

  return r, v, date

def eci2perif(raan, aop, i):
  row0 = [-m.sin(raan) * m.cos(i) * m.sin(aop) + m.cos(raan) * m.cos(aop), m.cos(raan) * m.cos(i) * m.sin(aop) + m.sin(raan) * m.cos(aop), m.sin(i) * m.sin(aop)]
  row1 = [-m.sin(raan) * m.cos(i) * m.cos(aop) - m.cos(raan) * m.sin(aop), m.cos(raan) * m.cos(i) * m.cos(aop) - m.sin(raan) * m.sin(aop), m.sin(i) * m.cos(aop)]
  row2 = [m.sin(raan) * m.sin(i), -m.cos(raan) * m.sin(i), m.cos(i)]

  return np.array([row0,row1,row2])

def ecc_anomonly(arr, method, tol=1e-8):
  if method == 'newton':
    # Newton's method for iteratively finding E
    Me, e = arr
    E1 = 0.0
    if Me < np.pi / 2.0: 
      E0 = Me + e / 2.0
    else:
      E0 = Me - e
    
    for n in range (200): # arbitrary max number of steps
      ratio = (E0 - e * np.sin(E0) - Me) / (1 - e * np.cos(E0))
      if abs(ratio) < tol:
        if n == 0:
          return E0
        else:
          return E1
      else: 
        E1 = E0 - ratio
        E0 = E1
    # did not converge
    return False 
  elif method == 'tae':
    ta, e = arr
    
    return 2 * m.atan(m.sqrt((1 - e) / (1 + e)) * m.tan(ta / 2.0))
  else:
    print('Invalid method for accentric anomoly')

def tle2coes(tle_filename, mu=pd.earth['mu']):
  # read tle file
  with open(tle_filename, 'r') as f:
    lines = f.readlines()
  
  # separate into three lines
  line0 = lines[0].strip() # name of satellite
  line1 = lines[1].strip().split()
  line2 = lines[2].strip().split()

  # epoch (year and day)
  epoch = line1[3]
  year, month, day, hour = calc_epoch(epoch)

  # collect coes
  # inclination
  i = float(line2[2])
  # right ascension of ascending node
  raan = float(line2[3])
  # eccentricity
  e = float('0.' + line2[4])
  # arguement of perigee
  aop = float(line2[5])
  # mean anomaly
  Me = float(line2[6])
  # mean motion
  mean_motion = float(line2[7]) # revolutions per day
  # period
  T = 1 / mean_motion * 24 * 3600 # seconds
  # semi major axis
  a = (T**2 * mu / 4.0 / np.pi**2)**(1 / 3.0)

  # calculate eccentric anomaly
  E = ecc_anomonly([Me,e], 'newton')

  # calculate true anomaly
  ta = true_anomaly([E,e])

  return a, e, i, ta, aop, raan, [year, month, day, hour]

def calc_epoch(epoch):
  # epoch year
  year = int('20' + epoch[:2])

  epoch = epoch[2:].split('.')

  # day of year
  day_of_year = int(epoch[0]) - 1

  # decimal hour of day
  hour = float('.' + epoch[1]) * 24.0

  # get year-month-day
  date = datetime.date(year, 1, 1) + datetime.timedelta(day_of_year)

  # extract month and day
  month = float(date.month)
  day = float(date.day)

  return year, month, day, hour

def true_anomaly(arr):
  E, e = arr

  return 2 * np.arctan(np.sqrt((1+e) / (1-e)) * np.tan(E / 2.0))

def tle2rv(tle_filename):
  return coes2rv(tle2coes(tle_filename))