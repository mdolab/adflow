
import numpy as np
import pprint
from pygeo import DVGeometry, DVConstraints

DVGEO_CLASSES = (DVGeometry,)
try:
    from pygeo import DVGeometryVSP
    DVGEO_CLASSES = (DVGeometry, DVGeometryVSP)
except ImportError:
    pass

from baseclasses import AeroProblem

from adflow import ADFLOW

from openmdao.api import ImplicitComponent
from openmdao.core.analysis_error import AnalysisError

from .om_utils import get_dvs_and_cons

class OM_STATES_COMP(ImplicitComponent):
    """OpenMDAO component that wraps the flow solve"""

    def initialize(self):
        self.options.declare('ap', types=AeroProblem)
        self.options.declare('dvgeo', types=DVGEO_CLASSES, allow_none=True, default=None)
        self.options.declare('solver')
        self.options.declare('use_OM_KSP', default=False, types=bool,
            desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")

        # self.options.declare('max_procs', default=64, types=int)

        self.options['distributed'] = True

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

    def setup(self):
        solver = self.options['solver']
        ap = self.options['ap']
        geo = self.options['dvgeo']

        # flag to keep track of clean restarts. We need this because sometimes
        # the current solution fails even if it already was a clean restart because
        # the previous one also failed. In these cases, we dont want to retry
        # with a clean restart because there is no point.
        self.cleanRestart = True

        self.ap_vars,_ = get_dvs_and_cons(ap=ap)
        self.geo_vars,_ = get_dvs_and_cons(geo=geo)

        if self.comm.rank == 0:
            print('adding ap var inputs')
        for (args, kwargs) in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size, units=kwargs['units'])
            if self.comm.rank == 0:
                print(name)

        if self.comm.rank == 0:
            print('adding geo var inputs')

        for (args, kwargs) in self.geo_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)
            if self.comm.rank == 0:
                print(name)

        local_state_size = solver.getStateSize()
        self.add_output('states', shape=local_state_size)

        # apparently needed to initialize the state arrays
        # but also somehow borks up the restarts
        # solver.getResidual(ap)

        #self.declare_partials(of='states', wrt='*')

    def _set_ap(self, inputs):
        tmp = {}
        for (args, kwargs) in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]

        self.options['ap'].setDesignVars(tmp)

    def _set_geo(self, inputs, update_jacobian=True):
        dvgeo = self.options['dvgeo']
        if dvgeo is None:
            return

        tmp = {}
        for (args, kwargs) in self.geo_vars:
            name = args[0]
            tmp[name] = inputs[name]
        try:
            self.options['dvgeo'].setDesignVars(tmp, update_jacobian)
        except TypeError: # this is needed because dvGeo and dvGeoVSP have different APIs
            self.options['dvgeo'].setDesignVars(tmp)


    def _set_states(self, outputs):
        self.options['solver'].setStates(outputs['states'])


    def apply_nonlinear(self, inputs, outputs, residuals):

        self._set_states(outputs)
        self._set_ap(inputs)
        self._set_geo(inputs, update_jacobian=False)

        ap = self.options['ap']
        residuals['states'] = self.options['solver'].getResidual(ap)


    def solve_nonlinear(self, inputs, outputs):

        solver = self.options['solver']
        ap = self.options['ap']

        if self._do_solve:

            self._set_ap(inputs)
            self._set_geo(inputs, update_jacobian=False)
            ap.solveFailed = False # might need to clear this out?
            ap.fatalFail = False

            solver(ap)

            if ap.fatalFail:
                if self.comm.rank == 0:
                    print('###############################################################')
                    print('# Solve Fatal Fail. Analysis Error')
                    print('###############################################################')

                raise AnalysisError('ADFLOW Solver Fatal Fail')


            if ap.solveFailed: # the mesh was fine, but it didn't converge
                # if the previous iteration was already a clean restart, dont try again
                if self.cleanRestart:
                    print('###############################################################')
                    print('# This was a clean restart. Will not try another one.')
                    print('###############################################################')
                    solver.resetFlow(ap)
                    self.cleanRestart = True
                    raise AnalysisError('ADFLOW Solver Fatal Fail')

                # the previous iteration restarted from another solution, so we can try again
                # with a re-set flowfield for the initial guess.
                else:
                    if self.comm.rank == 0:
                        print('###############################################################')
                        print('# Solve Failed, attempting a clean restart!')
                        print('###############################################################')

                    ap.solveFailed = False
                    ap.fatalFail = False
                    solver.resetFlow(ap)
                    solver(ap)

                    if ap.solveFailed or ap.fatalFail: # we tried, but there was no saving it
                        print('###############################################################')
                        print('# Clean Restart failed. There is no saving this one!')
                        print('###############################################################')

                        # re-set the flow for the next iteration:
                        solver.resetFlow(ap)
                        # set the reset flow flag
                        self.cleanRestart = True
                        raise AnalysisError('ADFLOW Solver Fatal Fail')

                    # see comment for the same flag below
                    else:
                        self.cleanRestart = False

            # solve did not fail, therefore we will re-use this converged flowfield for the next iteration.
            # change the flag so that if the next iteration fails with current initial guess, it can retry
            # with a clean restart
            else:
                self.cleanRestart = False


        outputs['states'] = solver.getStates()


    def linearize(self, inputs, outputs, residuals):

        self.options['solver']._setupAdjoint()

        self._set_ap(inputs)
        self._set_geo(inputs, update_jacobian=False)
        self._set_states(outputs)

        #print('om_states linearize')


    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):

        solver = self.options['solver']
        ap = self.options['ap']
        geo = self.options['dvgeo']

        #self._set_ap(inputs)
        #self._set_geo(inputs)
        #self._set_states(outputs)

        if mode == 'fwd':
            if 'states' in d_residuals:
                xDvDot = {}
                for var_name in d_inputs:
                    xDvDot[var_name] = d_inputs[var_name]
                if 'states' in d_outputs:
                    wDot = d_outputs['states']
                else:
                    wDot = None

                dwdot = solver.computeJacobianVectorProductFwd(xDvDot=xDvDot,
                                                               wDot=wDot,
                                                               residualDeriv=True)
                d_residuals['states'] += dwdot

        elif mode == 'rev':
            if 'states' in d_residuals:
                resBar = d_residuals['states']

                wBar, xDVBar = solver.computeJacobianVectorProductBwd(
                    resBar=resBar,
                    wDeriv=True, xDvDeriv=True)

                if 'states' in d_outputs:
                    d_outputs['states'] += wBar

                for dv_name, dv_bar in xDVBar.items():
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += dv_bar.flatten()


    def solve_linear(self, d_outputs, d_residuals, mode):
        solver = self.options['solver']
        ap = self.options['ap']
        if self.options['use_OM_KSP']:
                if mode == 'fwd':
                    d_outputs['states'] = solver.globalNKPreCon(d_residuals['states'], d_outputs['states'])
                elif mode == 'rev':
                    d_residuals['states'] = solver.globalAdjointPreCon(d_outputs['states'], d_residuals['states'])
        else:
            if mode == 'fwd':
                d_outputs['states'] = solver.solveDirectForRHS(d_residuals['states'])
            elif mode == 'rev':
                #d_residuals['states'] = solver.solveAdjointForRHS(d_outputs['states'])
                solver.adflow.adjointapi.solveadjoint(d_outputs['states'], d_residuals['states'], True)

        return True, 0, 0
