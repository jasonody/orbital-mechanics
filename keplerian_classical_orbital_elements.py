import planetary_data as pd
import tools as t
from OrbitPropagatorKep import OrbitPropagator as OP

tspan = 3600 * 24 * 1.0
dt = 10.0

cb = pd.earth

if __name__ == '__main__':
  # coes: a (pd radius + average of perigee and apogee), e (eccentricity), i (inclination), ta (true anomoly), aop (arguement of perigee), raan (right ascension of ascending node)

  # ISS (low Earth orbit)
  # Find orbital elements here: https://www.heavens-above.com/orbit.aspx?satid=25544
  c0 = [cb['radius'] + 414.0, 0.0002136, 51.6417, 0.0, 45.6438, 356.6431]

  # geostationary orbit
  c1 = [cb['radius'] + 35800.0, 0.0, 0.0, 0.0, 0.0, 0.0]

  # random orbit
  c2 = [cb['radius'] + 3000.0, 0.3, 20.0, 0.0, 15.0, 40.0]

  # Hubble space telescope
  # https://www.heavens-above.com/orbit.aspx?satid=20580
  c3 = [cb['radius'] + 538.0, 0.0002727, 28.4742, 0.0, 68.4651, 337.8887]

  # create orbit propagator
  op0 = OP(c0, tspan, dt, coes=True)
  op1 = OP(c1, tspan, dt, coes=True)
  op2 = OP(c2, tspan, dt, coes=True)
  op3 = OP(c3, tspan, dt, coes=True)

  op0.propagate_orbit()
  op1.propagate_orbit()
  op2.propagate_orbit()
  op3.propagate_orbit()

  t.plot_n_orbits([op0.rs, op1.rs, op2.rs, op3.rs], labels=['ISS', 'GEO', 'Random', 'Hubble'], show_plot=True)