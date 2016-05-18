from airfoil_parameterization import AirfoilAnalysis
import numpy as np
import scipy.optimize as optimize


CST = np.asarray([0.25, 0.25, 0.25, 0.25, -0.25, -0.25, -0.25, -0.25])
airfoil_analysis_options = dict(AnalysisMethod='CFD', AirfoilParameterization='CST', Re=5e5, ComputeGradient=True,
                                CFDiterations=500, CFDprocessors=16, cfdConfigFile='inv_NACA0012.cfg', ParallelAirfoils=False,
                                fd_step=1e-6, cs_step=1e-20, FreeFormDesign=True, BEMSpline='XFOIL', alphas=np.linspace(-15, 15, 30))

ub = [1, 1, 1, 1, 0.2, 0.2, 0.2, 0.2]
lb = [0.205, 0.205, 0.205, 0.205, -1, -1, -1, -1]

af = AirfoilAnalysis(CST, airfoil_analysis_options)
x_original, y_original = af.getCoordinates()

print x_original
print y_original


bounds = []
for i in range(len(lb)):
    bounds.append((lb[i], ub[i]))

def obj(CST):
    print CST
    af = AirfoilAnalysis(CST, airfoil_analysis_options)
    cl, cd, dcl_dalpha, dcd_dalpha, dcl_dRe, dcd_dRe, dcl_dafp, dcd_dafp, lexitflag = af.computeDirect(np.radians(5.0), 5e5)
    print cl, cd, dcl_dafp, dcd_dafp
    return -cl/cd

# res = optimize.minimize(obj, CST, method='SLSQP', jac=None, bounds=bounds, options={'eps':1e-3})
# print res.message
# print res.x
# print res.fun

# af = AirfoilAnalysis(res.x, airfoil_analysis_options)
# x_new, y_new = af.getCoordinates()

# print x_new
# print y_new

import matplotlib.pylab as plt
plt.figure()
plt.plot(x_original, y_original, '^', label='Original')
plt.plot(x_new, y_new, 's', label='Optimized')
plt.legend(loc='best')
plt.show()




