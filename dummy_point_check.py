import numpy as np
q = -0.76505532392946469

dh = np.float64(2.0)*np.pi/np.float64(8.0)

val = (q+np.float64(1.0))*dh
pstr = "%18.16e\n" % val
print(pstr)
# output: 1.8452511708580682e-01