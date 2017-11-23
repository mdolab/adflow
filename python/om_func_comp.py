import numpy as np

from openmdao.api import ExplicitComponent

from baseclasses import AeroProblem
from pygeo import DVGeometry
from adflow import ADFLOW

from .om_utils import get_dvs_and_cons

FUNCS_UNITS={
    'mdot': 'kg/s', 
    'mavgptot': 'Pa', 
    'mavgps': 'Pa', 
    'mavgttot': 'degK',
    'drag': 'N',
    'dragpressure': 'N', 
    'dragviscous': 'N',
    'dragmomentum': 'N',

}

class OM_FUNC_COMP(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('ap', types=AeroProblem,)
        self.metadata.declare('dvgeo', types=DVGeometry, allow_none=True, default=None)
        self.metadata.declare('solver', types=ADFLOW)

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

    def setup(self):
        solver = self.metadata['solver']
        ap = self.metadata['ap']
        geo = self.metadata['dvgeo']

        self.ap_vars,_ = get_dvs_and_cons(ap=ap)
        self.geo_vars,_ = get_dvs_and_cons(geo=geo)

        for (args, kwargs) in self.geo_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)

        for (args, kwargs) in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size, units=kwargs['units'])

        local_state_size = self.metadata['solver'].getStateSize()
        local_state_sizes = self.comm.allgather(local_state_size)
        iproc = self.comm.rank

        ind1 = np.sum(local_state_sizes[:iproc])
        ind2 = np.sum(local_state_sizes[:iproc+1])

        self.add_input('states', src_indices=np.arange(ind1, ind2, dtype=int))

        for f_name, f_meta in solver.adflowCostFunctions.items():
            f_type = f_meta[1] 
            units = None
            if f_type in FUNCS_UNITS: 
                units = FUNCS_UNITS[f_type]
            self.add_output(f_name, shape=1, units=units)


            self.declare_partials(of=f_name, wrt='*')
                

    def _set_ap(self, inputs):
        tmp = {}
        for (args, kwargs) in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]
        self.metadata['ap'].setDesignVars(tmp)
        self.metadata['solver'].setAeroProblem(self.metadata['ap'])

    def _set_geo(self, inputs):
        dvgeo = self.metadata['dvgeo']
        if dvgeo is None: 
            return 

        tmp = {}
        for (args, kwargs) in self.geo_vars:
            name = args[0]
            tmp[name] = inputs[name]
        self.metadata['dvgeo'].setDesignVars(tmp)

    def _set_states(self, inputs):
        self.metadata['solver'].setStates(inputs['states'])

    def _get_func_name(self, name):
        return '%s_%s' % (self.metadata['ap'].name, name.lower())

    def compute(self, inputs, outputs):
        solver = self.metadata['solver']
        ap = self.metadata['ap']



        #actually setting things here triggers some kind of reset, so we only do it if you're actually solving
        if self._do_solve: 
            self._set_ap(inputs)
            self._set_geo(inputs)
            self._set_states(inputs)

        funcs = {}
        solver.evalFunctions(ap, funcs, ap.evalFuncs)

        for name in ap.evalFuncs:
            outputs[name] = funcs[self._get_func_name(name)]


    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        solver = self.metadata['solver']
        ap = self.metadata['ap']

        if self._do_solve: 
            self._set_ap(inputs)
            self._set_geo(inputs)
            self._set_states(inputs)
        
        if mode == 'fwd':
            xDvDot = {}
            for key in ap.DVs:
                if key in d_inputs:
                    mach_name = key.split('_')[0]
                    xDvDot[mach_name] = d_inputs[key]
            for (args, kwargs) in self.metadata['geo_vars'].variables:
                name = args[0]
                xDvDot[name] = d_inputs[name]

            if 'states' in d_inputs:
                wDot = d_inputs['states']
            else:
                wDot = None

            funcsdot = solver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot,
                wDot=wDot,
                funcDeriv=True)

            for name in ap.evalFuncs:
                func_name = name.lower()
                if name in d_outputs:
                    d_outputs[name] += funcsdot[func_name]

        elif mode == 'rev':
            funcsBar = {}
            for name in ap.evalFuncs:
                func_name = name.lower()
                if name in d_outputs:
                    funcsBar[func_name] = d_outputs[name] / self.comm.size

            # print('funcsBar', funcsBar)
            wBar, xDVBar = solver.computeJacobianVectorProductBwd(
                funcsBar=funcsBar,
                wDeriv=True, xDvDeriv=True)

            if 'states' in d_inputs:
                d_inputs['states'] += wBar
                # print('wBar', np.linalg.norm(wBar))

            for dv_name, dv_bar in xDVBar.items():
                if dv_name in d_inputs:
                    d_inputs[dv_name] += dv_bar.flatten()
