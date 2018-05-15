
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
        self.metadata.declare('ap', types=AeroProblem)
        self.metadata.declare('dvgeo', types=DVGEO_CLASSES, allow_none=True, default=None)
        self.metadata.declare('solver')
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

    def _set_dvs(self, inputs, update_jacobian=True): 
        dvgeo = self.metadata['dvgeo']
        ap = self.metadata['ap']

        tmp = {}
        for name in inputs.keys():
            tmp[name] = inputs[name]

        try: 
            dvgeo.setDesignVars(tmp, update_jacobian)
        except TypeError: # this is needed because dvGeo and dvGeoVSP have different APIs
            dvgeo.setDesignVars(tmp)

        ap.setDesignVars(tmp) 

    def _set_states(self, outputs):
        self.metadata['solver'].setStates(outputs['states'])
        
    
    def apply_nonlinear(self, inputs, outputs, residuals):
        
        self._set_dvs(inputs)
        self._set_states(outputs)
        
        ap = self.metadata['ap']
        residuals['states'] = self.metadata['solver'].getResidual(ap)

    
    def solve_nonlinear(self, inputs, outputs):

        solver = self.metadata['solver']
        ap = self.metadata['ap']

        if self._do_solve: 
            self._set_dvs(inputs)
            ap.solveFailed = False # might need to clear this out?
            ap.fatalFail = False
        
            solver(ap)

            if ap.fatalFail:
                if self.comm.rank == 0:
                    print('###############################################################')
                    print('#Solve Fatal Fail. Analysis Error')
                    print('###############################################################')

                raise AnalysisError('ADFLOW Solver Fatal Fail')


            if ap.solveFailed: # the mesh was fine, but it didn't converge
                if self.comm.rank == 0:
                    print('###############################################################')
                    print('#Solve Failed, attempting a clean restart!')
                    print('###############################################################')

                ap.solveFailed = False
                ap.fatalFail = False
                solver.resetFlow(ap)
                solver(ap)

                if ap.solveFailed or ap.fatalFail: # we tried, but there was no saving it
                    print('###############################################################')
                    print('#Clean Restart failed. There is no saving this one!')
                    print('###############################################################')

                    raise AnalysisError('ADFLOW Solver Fatal Fail')


        outputs['states'] = solver.getStates()

    
    def linearize(self, inputs, outputs, residuals):

        self.metadata['solver']._setupAdjoint()

        self._set_dvs(inputs)
        self._set_states(outputs)

        #print('om_states linearize')
        
    
    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):

        solver = self.metadata['solver']
        ap = self.metadata['ap']
        geo = self.metadata['dvgeo']

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
                #d_residuals['states'] = solver.solveAdjointForRHS(d_outputs['states'])
                solver.adflow.adjointapi.solveadjoint(d_outputs['states'], d_residuals['states'], True)

        return True, 0, 0
