import matplotlib.pylab as plt
from airfoil_parameterization import AirfoilAnalysis
import numpy as np
from openmdao.api import Problem, SqliteRecorder
from openmdao.drivers.pyoptsparse_driver import pyOptSparseDriver
from openmdao.api import IndepVarComp, Component, ExecComp, Group, ScipyGMRES

CST = np.asarray([0.25, 0.25, 0.25, 0.25, -0.25, -0.25, -0.25, -0.25])


airfoil_analysis_options = dict(AnalysisMethod='CFD', AirfoilParameterization='CST',
                                CFDiterations=10000, CFDprocessors=16, FreeFormDesign=True, BEMSpline='XFOIL', maxDirectAoA=1000, fd_step=1e-6, cs_step=1e-20,
                                alphas=np.linspace(-15, 15, 30), Re=5e5, ComputeGradient=False, cfdConfigFile='inv_NACA0012.cfg', ParallelAirfoils=True)
af = AirfoilAnalysis(CST, airfoil_analysis_options)
x_old, y_old = af.getCoordinates()


class Coefficients(Component):
    def __init__(self):
        super(Coefficients, self).__init__()

        self.add_param('CST', val=np.zeros(8))
        self.add_param('alpha', val=0.0)
        self.add_output('cl', val=0.0)
        self.add_output('cd', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        print params['CST']
        af = AirfoilAnalysis(params['CST'], airfoil_analysis_options)
        if airfoil_analysis_options['ComputeGradient']:
            cl, cd, self.dcl_dalpha, self.dcd_dalpha, dcl_dRe, dcd_dRe, self.dcl_dafp, self.dcd_dafp, lexitflag = af.computeDirect(np.radians(params['alpha']), 5e5)
            print self.dcl_dafp, self.dcd_dafp
        else:
            cl, cd = af.computeDirect(np.radians(params['alpha']), 5e5)
        unknowns['cl'] = cl
        unknowns['cd'] = cd

    def linearize(self, params, unknowns, resids):
        J = {}
        J['cl', 'alpha'] = self.dcl_dalpha*np.pi/180.
        J['cd', 'alpha'] = self.dcd_dalpha*np.pi/180.
        J['cl', 'CST'] = self.dcl_dafp.reshape(1,8)
        J['cd', 'CST'] = self.dcd_dafp.reshape(1,8)
        return J

class liftdrag(Component):
    def __init__(self):
        super(liftdrag, self).__init__()

        self.add_param('cl', val=0.0)
        self.add_param('cd', val=0.0)
        self.add_output('LD', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        print params['cl'], params['cd']
        unknowns['LD'] = -params['cl']/params['cd']

    def linearize(self, params, unknowns, resids):
        J = {}
        J['LD', 'cl'] = -1./params['cd']
        J['LD', 'cd'] = params['cl']/(params['cd']**2)
        return J

class G(Group):
    def __init__(self):
        super(G, self).__init__()
        self.add('CST', IndepVarComp('CST', np.zeros(8)), promotes=['*'])
        self.add('alpha', IndepVarComp('alpha', 0.0), promotes=['*'])
        self.add('coeff', Coefficients(), promotes=['*'])
        self.add('liftd', liftdrag(), promotes=['*'])

        self.fd_options['step_size'] = 1e-3
        self.fd_options['form'] =  'forward'
        self.fd_options['step_type'] = 'relative'

p = Problem(root=G())
#p.driver = pyOptSparseDriver()
#p.driver.options['optimizer'] = 'SNOPT'
#p.driver.add_desvar('alpha', lower=-10., upper=20.)
lower = np.ones((1,8))*[[-0.6, -0.76, -0.4, -0.25, 0.13, 0.16, 0.13, 0.1]]
upper = np.ones((1,8))*[[-0.13, -0.16, -0.13, 0.15, 0.55, 0.55, 0.4, 0.4]]
p.driver.add_desvar('CST', lower=lower, upper=upper)#, scaler=np.asarray([0.1, 0.005, 0.1, 0.1]))
# p.driver.add_objective('LD')
p.setup()
p['CST'] = CST
p['alpha'] = 5.0
p.run()

# p.check_partial_derivatives()
# p.check_total_derivatives()
#gradad = p.calc_gradient(['alpha', 'CST'], ['LD'], mode='auto')
airfoil_analysis_options['ComputeGradient'] = False
gradfd = p.calc_gradient(['alpha', 'CST'], ['LD'], mode='fd')
print "GRADIENT"
#print gradad
print gradfd

#plt.figure()
#plt.plot(range(len(gradfd[0])), (np.asarray(gradfd[0])), '*-', label='FD')
#plt.plot(range(len(gradad[0])), (np.asarray(gradad[0])), 's-', label='AD')
#plt.legend(loc='best')
#plt.show()

print p['LD']
print p['CST']
print p['alpha']

af_new = AirfoilAnalysis(p['CST'], airfoil_analysis_options)
x_new, y_new = af_new.getCoordinates()

xfoil_CST = [-0.13, -0.16, -0.13, 0.15, 0.210786, 0.289470, 0.4, 0.190282]
xfoil_alpha = 6.450669
xfoil_LD = -117.66

af_xfoil = AirfoilAnalysis(xfoil_CST, airfoil_analysis_options)
x_xfoil, y_xfoil = af_xfoil.getCoordinates()

plt.figure()
plt.plot(x_old, y_old, '*', label='Original')
plt.plot(x_xfoil, y_xfoil, '^', label='XFOIL')
plt.plot(x_new, y_new, 's', label='New')
plt.legend(loc='best')
plt.show()
