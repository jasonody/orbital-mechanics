import planetary_data as pd
import tools as t
from OrbitPropagatorKep import OrbitPropagator as OP
from OrbitPropagatorKep import null_perts

tspan = 3600 * 24 * 2
dt = 100.0

cb = pd.earth

if __name__ == '__main__':
  perts = null_perts()
  perts['J2'] = True

  op = OP(t.tle2coes('TLEs/iss.tle.txt'), tspan, dt, coes=True, perts=perts)

  op.calculate_coes()

  # if you see "Calculating COEs", all is good :)