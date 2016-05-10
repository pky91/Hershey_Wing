from airfoil_parameterization import AirfoilAnalysis
import numpy as np
import scipy.optimize as optimize


CST = np.asarray([0.25, 0.25, 0.25, 0.25, -0.25, -0.25, -0.25, -0.25])
airfoil_analysis_options = dict(AnalysisMethod='CFD', AirfoilParameterization='CST', Re=5e5, ComputeGradient=False,
                                CFDiterations=1000, CFDprocessors=2, cfdConfigFile='inv_NACA0012.cfg', ParallelAirfoils=False,
                                fd_step=1e-6, cs_step=1e-20)

ub = [1, 1, 1, 1, 0.2, 0.2, 0.2, 0.2]
lb = [0.05, 0.05, 0.05, -1, -1, -1, -1, -1]

af = AirfoilAnalysis(CST, airfoil_analysis_options)
x_original, y_original = af.getCoordinates()


bounds = []
for i in range(len(lb)):
    bounds.append((lb[i], ub[i]))

def obj(CST):
    af = AirfoilAnalysis(CST, airfoil_analysis_options)
    cl, cd = af.computeDirect(np.radians(5.0), 5e5)
    return -cl/cd

res = optimize.minimize(obj, CST, method='SLSQP', jac=None, bounds=bounds)
print res.message
print res.x
print res.fun

w = res.x
af = AirfoilAnalysis(res.x, airfoil_analysis_options)
x_new, y_new = af.getCoordinates()

import matplotlib.pylab as plt
plt.figure()
plt.plot(x_original, y_original, '^', label='Original')
plt.plot(x_new, y_new, 's', label='Optimized')
plt.legend(loc='best')
plt.show()




