from pygeo import DVGeometry, DVConstraints
from baseclasses import AeroProblem

from adflow import ADFLOW

from openmdao.api import ImplicitComponent
from openmdao.core.analysis_error import AnalysisError

from om_utils import get_dvs_and_cons




class OM_STATES_COMP(ImplicitComponent):
    """OpenMDAO component that wraps the flow solve"""

    def initialize(self):
        self.metadata.declare('ap', types=AeroProblem)
        self.metadata.declare('dvgeo', types=DVGeometry, allow_none=True, default=None)
        self.metadata.declare('solver', types=ADFLOW)
        self.metadata.declare('use_OM_KSP', default=False, types=bool, 
            desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")

        # self.metadata.declare('max_procs', default=64, types=int)

        self.distributed = True

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

    def setup(self):
        solver = self.metadata['solver']
        ap = self.metadata['ap']
        geo = self.metadata['dvgeo']

        self.ap_vars,_ = get_dvs_and_cons(ap=ap)
        self.geo_vars,_ = get_dvs_and_cons(geo=geo)

        for (args, kwargs) in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size, units=kwargs['units'])

        for (args, kwargs) in self.geo_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)

        local_state_size = solver.getStateSize()

        # apparently needed to initialize the state arrays
        solver.getResidual(ap)
        
        self.add_output('states', shape=local_state_size)

        self.declare_partials(of='states', wrt='*')

    def _set_ap(self, inputs):
        tmp = {}
        for (args, kwargs) in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]
        # self.metadata['ap'].setDesignVars(tmp)

    def _set_geo(self, inputs):
        dvgeo = self.metadata['dvgeo']
        if dvgeo is None: 
            return 

        tmp = {}
        for (args, kwargs) in self.geo_vars:
            name = args[0]
            tmp[name] = inputs[name]
        self.metadata['dvgeo'].setDesignVars(tmp)

    def _set_states(self, outputs):
        self.metadata['solver'].setStates(outputs['states'])

    def apply_nonlinear(self, inputs, outputs, residuals):
        
        self._set_ap(inputs)
        self._set_geo(inputs)
        self._set_states(outputs)
        
        ap = self.metadata['ap']
        residuals['states'] = self.metadata['solver'].getResidual(ap)

    def solve_nonlinear(self, inputs, outputs):
        solver = self.metadata['solver']
        ap = self.metadata['ap']

        if self._do_solve: 
            self._set_ap(inputs)
            self._set_geo(inputs)
            ap.solveFailed = False # might need to clear this out?
            ap.fatalFail = False
    
            solver(ap)

            if ap.solveFailed:
                if self.comm.rank == 0:
                    print('###############################################################')
                    print('#Solve Tolerance Not Reached, attempting a clean restart!')
                    print('###############################################################')

                solver.resetFlow(ap)
                solver(ap)

            if ap.fatalFail:
                print('###############################################################')
                print('#Solve Fatal Fail. Analysis Error')
                print('###############################################################')

                raise AnalysisError('ADFLOW Solver Fatal Fail')

        outputs['states'] = solver.getStates()

        #print('foobar fatalFail', ap.fatalFail)
        #print('foobar solveFailed', ap.solveFailed)

    def linearize(self, inputs, outputs, residuals):
        self.metadata['solver']._setupAdjoint()

    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):
        solver = self.metadata['solver']
        ap = self.metadata['ap']
        geo = self.metadata['dvgeo']


        self._set_ap(inputs)
        self._set_geo(inputs)
        self._set_states(outputs)


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
        solver = self.metadata['solver']
        ap = self.metadata['ap']
        if self.metadata['use_OM_KSP']:
                if mode == 'fwd':
                    d_outputs['states'] = solver.globalNKPreCon(d_residuals['states'], d_outputs['states'])
                elif mode == 'rev':
                    d_residuals['states'] = solver.globalAdjointPreCon(d_outputs['states'], d_residuals['states'])
        else:
            if mode == 'fwd':
                d_outputs['states'] = solver.solveDirectForRHS(d_residuals['states'])
            elif mode == 'rev':
                # import numpy as np
                # print('states seed', np.linalg.norm(d_outputs['states']))
                d_residuals['states'] = solver.solveAdjointForRHS(d_outputs['states'])


        return True, 0, 0
