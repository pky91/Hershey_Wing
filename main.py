from airfoil_parameterization import AirfoilAnalysis
import numpy as np
import scipy.optimize as optimize



CST = np.asarray([-0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25])
airfoil_analysis_options = dict(AnalysisMethod='CFD', AirfoilParameterization='CST', Re=5e5, ComputeGradient=False,
                                CFDiterations=5000, CFDprocessors=16, cfdConfigFile='inv_NACA0012.cfg', ParallelAirfoils=False,
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
    if airfoil_analysis_options['ComputeGradient']:
        cl, cd, dcl_dalpha, dcd_dalpha, dcl_dRe, dcd_dRe, dcl_dafp, dcd_dafp, lexitflag = af.computeDirect(np.radians(5.0), 5e5)
        dld_dcst = np.array(dcl_dafp/cd-cl/cd**2*dcd_dafp)
    else:
        cl, cd = af.computeDirect(np.radians(5.0), 5e5)
        dld_dcst = 0.0
        dcd_dafp = 0.0
        dcl_dafp = 0.0
    return cl, cd, dcl_dafp, dcd_dafp, -cl/cd, -dld_dcst

h = 5e-3
cl, cd, dcl_dafp, dcd_dafp, ld, dld_dcst_zero = obj(CST)
fld = np.zeros(len(CST))
for i in range(len(CST)):
    CST[i] = CST[i]+h
    cl1, cd1, dcl_dafp1, dcd_dafp1, ld_1, dld_dcst1 = obj(CST)
    fld[i] = (ld_1 - ld)/h
    CST[i] = CST[i]-h

airfoil_analysis_options['ComputeGradient'] = True
cl, cd, dcl_dafp, dcd_dafp, ld, dld_dcst = obj(CST)

print "cl", cl
print "cd", cd
print "cl_afp", dcl_dafp
print "cd_afp", dcd_dafp
print"***"

print "Adjoint gradient", dld_dcst 
print "Finite difference", fld
print "Error", (dld_dcst - fld) / fld

# X_Axis = np.array([1, 2, 3, 4, 5, 6, 7, 8])
# import matplotlib.pylab as plt
# plt.figure()
# plt.plot(X_Axis, dld_dcst, linestyle='-', marker='^', color='r', label='Adjoint')
# plt.plot(X_Axis, fld, linestyle='-', marker='o', color='g', label='Finite-Diff.')
# plt.xlabel('CST parameters')
# plt.ylabel('-(Lift/Drag Ratio)')
# plt.legend(loc='best')
# plt.show()

res = optimize.minimize(obj, CST, method='SLSQP', jac=True, bounds=bounds, options={'eps':1e-3})
print res.message
print res.x
print res.fun

af = AirfoilAnalysis(res.x, airfoil_analysis_options)
x_new, y_new = af.getCoordinates()

print x_new
print y_new

import matplotlib.pylab as plt
plt.figure()
plt.plot(x_original, y_original, '^', label='Original')
plt.plot(x_new, y_new, 's', label='Optimized')
plt.legend(loc='best')
plt.show()

print "Done"



