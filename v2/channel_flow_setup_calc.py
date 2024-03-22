import numpy as np

Cf_coeff = 0.073
Cf_exp = -0.25

Re_bulk = 6900.0

Re_tau = np.sqrt(0.5*Cf_coeff*(2.0**(Cf_exp))*(Re_bulk**(2.0+Cf_exp)))

print("Bulk velocity Reynolds number: %.3f" % Re_bulk)
print("Friction velocity Reynolds number: %.3f" % Re_tau)