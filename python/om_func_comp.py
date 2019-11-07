
import numpy as np

from openmdao.api import ExplicitComponent


from baseclasses import AeroProblem
from adflow import ADFLOW

from pygeo import DVGeometry

from collections import OrderedDict

DVGEO_CLASSES = (DVGeometry,)
try:
    from pygeo import DVGeometryVSP
    DVGEO_CLASSES = (DVGeometry, DVGeometryVSP)
except ImportError:
    pass

from .om_utils import get_dvs_and_cons

FUNCS_UNITS={
    'mdot': 'kg/s',
    'mavgptot': 'Pa',
    'mavgps': 'Pa',
    'aavgptot': 'Pa',
    'aavgps': 'Pa',
    'mavgttot': 'degK',
    'mavgvx':'m/s',
    'mavgvy':'m/s',
    'mavgvz':'m/s',
    'drag': 'N',
    'lift': 'N',
    'dragpressure': 'N',
    'dragviscous': 'N',
    'dragmomentum': 'N',
    'fx': 'N',
    'fy': 'N',
    'fz': 'N',
    'forcexpressure': 'N',
    'forceypressure': 'N',
    'forcezpressure': 'N',
    'forcexviscous': 'N',
    'forceyviscous': 'N',
    'forcezviscous': 'N',
    'forcexmomentum': 'N',
    'forceymomentum': 'N',
    'forcezmomentum': 'N',
    'flowpower': 'W',
    'area':'m**2',
}

class OM_FUNC_COMP(ExplicitComponent):

    def initialize(self):
        self.options.declare('ap', types=AeroProblem,)
        self.options.declare('dvgeo', types=DVGEO_CLASSES, allow_none=True, default=None)
        self.options.declare('solver')

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

        #self.distributed=True

    def setup(self):
        solver = self.options['solver']
        ap = self.options['ap']
        geo = self.options['dvgeo']

        self.ap_vars,_ = get_dvs_and_cons(ap=ap)
        self.geo_vars,_ = get_dvs_and_cons(geo=geo)

        for (args, kwargs) in self.geo_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)

        for (args, kwargs) in self.ap_vars:
            name = args[0]
            size = args[1]
            value = 1.
            if 'value' in kwargs:
                value = kwargs['value']

            self.add_input(name, shape=size, units=kwargs['units'])

        local_state_size = self.options['solver'].getStateSize()
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

            if self.comm.rank == 0:
                print("adding adflow func as output: {}".format(f_name))
            self.add_output(f_name, shape=1, units=units)

            #self.declare_partials(of=f_name, wrt='*')

    def _set_ap(self, inputs):
        tmp = {}
        for (args, kwargs) in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name][0]

        self.options['ap'].setDesignVars(tmp)
        #self.options['solver'].setAeroProblem(self.options['ap'])

    def _set_geo(self, inputs, update_jacobian=True):
        dvgeo = self.options['dvgeo']
        if dvgeo is None:
            return

        tmp = {}
        for (args, kwargs) in self.geo_vars:
            name = args[0]
            tmp[name] = inputs[name]

        # if self.comm.rank == 0:
        #     import pprint
        #     pprint.pprint(tmp)
        try:
            self.options['dvgeo'].setDesignVars(tmp, update_jacobian)
        except TypeError: # this is needed because dvGeo and dvGeoVSP have different APIs
            self.options['dvgeo'].setDesignVars(tmp)

    def _set_states(self, inputs):
        self.options['solver'].setStates(inputs['states'])

    def _get_func_name(self, name):
        return '%s_%s' % (self.options['ap'].name, name.lower())

    def compute(self, inputs, outputs):
        solver = self.options['solver']
        ap = self.options['ap']
        #print('funcs compute')
        #actually setting things here triggers some kind of reset, so we only do it if you're actually solving
        if self._do_solve:
            self._set_ap(inputs)
            self._set_geo(inputs, update_jacobian=False)
            self._set_states(inputs)


        funcs = {}

        eval_funcs = [f_name for f_name, f_meta in solver.adflowCostFunctions.items()]
        solver.evalFunctions(ap, funcs, eval_funcs)
        #solver.evalFunctions(ap, funcs)

        #for name in ap.evalFuncs:
        for name in solver.adflowCostFunctions.keys():
            f_name = self._get_func_name(name)
            if f_name in funcs:
                outputs[name.lower()] = funcs[f_name]

    def _compute_partials(self, inputs, J):

        #self._set_ap(inputs)
        #self._set_geo(inputs)
        #self._set_states(inputs)

        #print('om_funcs linearize')

        pass


    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        solver = self.options['solver']
        ap = self.options['ap']


        #self.options['solver'].setAeroProblem(ap)
        #print('func matvec')

        #if self._do_solve:
        #    self._set_ap(inputs)
        #    self._set_geo(inputs)
        #    self._set_states(inputs)

        if mode == 'fwd':
            xDvDot = {}
            for key in ap.DVs:
                if key in d_inputs:
                    mach_name = key.split('_')[0]
                    xDvDot[mach_name] = d_inputs[key]
            for (args, kwargs) in self.geo_vars:
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

            for name, meta in ap.evalFuncs:
                func_name = name.lower()
                if name in d_outputs:
                    d_outputs[name] += funcsdot[func_name]

        elif mode == 'rev':
            funcsBar = {}

            for name, meta  in solver.adflowCostFunctions.items():
                func_name = name.lower()

                # we have to check for 0 here, so we don't include any unnecessary variables in funcsBar
                # becasue it causes ADflow to do extra work internally even if you give it extra variables, even if the seed is 0
                if func_name in d_outputs and d_outputs[func_name] != 0.:
                    funcsBar[func_name] = d_outputs[func_name][0] / self.comm.size

            # because of the 0 checking, the funcsBar is now only correct on the root proc,
            # so we need to broadcast it to everyone. its not actually imporant that the seeds are the same everywhere,
            # but the keys in the dictionary need to be the same.
            funcsBar = self.comm.bcast(funcsBar, root=0)


            #print(funcsBar, flush=True)

            d_input_vars = list(d_inputs.keys())
            n_input_vars = len(d_input_vars)

            wBar = None
            xDVBar = None

            if 'states' in d_inputs and n_input_vars==1 :
                wBar = solver.computeJacobianVectorProductBwd(
                    funcsBar=funcsBar,
                    wDeriv=True)
            elif ('states' not in d_input_vars) and n_input_vars:
                xDVBar = solver.computeJacobianVectorProductBwd(
                    funcsBar=funcsBar,
                    xDvDeriv=True)

            elif ('states' in d_input_vars) and n_input_vars:
                wBar, xDVBar = solver.computeJacobianVectorProductBwd(
                    funcsBar=funcsBar,
                    wDeriv=True, xDvDeriv=True)
            else: # nothing to do, so why is this being called?
                return

            if wBar is not None:
                d_inputs['states'] += wBar

            if xDVBar is not None:
                for dv_name, dv_bar in xDVBar.items():
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += dv_bar.flatten()


