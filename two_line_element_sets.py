import planetary_data as pd
import tools as t
from OrbitPropagatorKep import OrbitPropagator as OP

# time parameters
tspan = 3600 * 24 * 1.0
dt = 10.0

# central body
cb = pd.earth

if __name__ == '__main__':
  op0 = OP(t.tle2coes('TLEs/iss.tle.txt'), tspan, dt, coes=True) # https://celestrak.com/NORAD/elements/stations.txt
  op1 = OP(t.tle2coes('TLEs/delphini.tle.txt'), tspan, dt, coes=True)
  op2 = OP(t.tle2coes('TLEs/starlink-31.tle.txt'), tspan, dt, coes=True)
  
  t.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels=['ISS', 'Delphini', 'Starlink-31'], show_plot=True)