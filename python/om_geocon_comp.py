
import numpy as np

from six import iteritems, itervalues

from openmdao.api import ExplicitComponent

from pygeo import DVGeometry, DVConstraints

DVGEO_CLASSES = (DVGeometry,)
try: 
    from pygeo import DVGeometryVSP
    DVGEO_CLASSES = (DVGeometry, DVGeometryVSP)
except ImportError: 
    pass


from .om_utils import get_dvs_and_cons

class OM_GEOCON_COMP(ExplicitComponent):
    """OpenMDAO Component that wraps the geometry constraints calculation"""

    def initialize(self):
        self.options.declare('dvgeo', types=DVGEO_CLASSES)
        self.options.declare('dvcon', types=DVConstraints)

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

    def setup(self):
        self.dvs, _ = get_dvs_and_cons(self.options['dvgeo'])
        for (args, kwargs) in self.dvs:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)


        funcsSens = {}
        self.options['dvcon'].evalFunctionsSens(funcsSens, includeLinear=True)

        for cons in itervalues(self.options['dvcon'].constraints):
            for name, con in iteritems(cons):
                if self.comm.rank == 0: 
                    print('nonlinear con', name)
                self.add_output(name, shape=con.nCon)
                jac = funcsSens[name]
                for wrt_var, subjac in iteritems(jac):
                    self.declare_partials(of=name, wrt=wrt_var)

        for name, con in iteritems(self.options['dvcon'].linearCon):
            if self.comm.rank == 0: 
                print('linear_con', name)
            self.add_output(name, shape=con.ncon)
            jac = funcsSens[name]
            for wrt_var, subjac in iteritems(jac):
                self.declare_partials(of=name, wrt=wrt_var, val=subjac)


    def _set_geo(self, inputs, updateJacobian=True):
        tmp = {}
        for (args, kwargs) in self.dvs:
            name = args[0]
            tmp[name] = inputs[name]
        
        try: 
            self.options['dvgeo'].setDesignVars(tmp, update_jacobian)
        except TypeError: # this is needed because dvGeo and dvGeoVSP have different APIs
            self.options['dvgeo'].setDesignVars(tmp)

    def compute(self, inputs, outputs):
    
        if self._do_solve: 
            # self._set_geo(inputs, updateJacobian=False)
            pass 

        funcs = {}
        self.options['dvcon'].evalFunctions(funcs, includeLinear=True)

        for cons in itervalues(self.options['dvcon'].constraints):
            for name, con in iteritems(cons):
                #print('nonlinear foobar', name, funcs[name])
                outputs[name] = funcs[name]

        for name, con in iteritems(self.options['dvcon'].linearCon):
            #print('linear foobar', name, funcs[name])
            outputs[name] = funcs[name]

    def compute_partials(self, inputs, J):
        #print('om_geocon linearize')
        if self._do_solve:
            #self._set_geo(inputs)
            pass 

        funcsSens = {}

        self.options['dvcon'].evalFunctionsSens(funcsSens)

        for of_var, jac in funcsSens.items():
            #print('of: ', of_var)
            for wrt_var, subjac in jac.items():
                #print('  wrt: ', wrt_var)
                J[of_var, wrt_var] = subjac

