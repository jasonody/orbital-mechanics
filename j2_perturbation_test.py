import planetary_data as pd
import tools as t
from OrbitPropagatorKep import OrbitPropagator as OP
from OrbitPropagatorKep import null_perts

tspan = 3600 * 24 * 5.0 # 5 days
dt = 100.0

cb = pd.earth

if __name__ == '__main__':
  perts = null_perts()
  perts['J2'] = True
  
  op = OP(t.tle2coes('TLEs/iss1.tle.txt'), tspan, dt, coes=True, cb=cb, perts=perts)

  op.plot_3d(show_plot=True)