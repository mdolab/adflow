import numpy as np

from six import iteritems, itervalues

from openmdao.api import ExplicitComponent

from pygeo import DVConstraints, DVGeometry

from .om_utils import get_dvs_and_cons

class OM_GEOCON_COMP(ExplicitComponent):
    """OpenMDAO Component that wraps the geometry constraints calculation"""

    def initialize(self):
        self.metadata.declare('dvgeo', type_=DVGeometry)
        self.metadata.declare('dvcon', type_=DVConstraints)

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

    def setup(self):
        self.dvs, _ = get_dvs_and_cons(self.metadata['dvgeo'])
        for (args, kwargs) in self.dvs:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)



        funcsSens = {}
        self.metadata['dvcon'].evalFunctionsSens(funcsSens, includeLinear=True)

        for cons in itervalues(self.metadata['dvcon'].constraints):
            for name, con in iteritems(cons):
                #print('nonlinear con', name)
                self.add_output(name, shape=con.nCon)
                jac = funcsSens[name]
                for wrt_var, subjac in iteritems(jac):
                    self.declare_partials(of=name, wrt=wrt_var, val=subjac)

        for name, con in iteritems(self.metadata['dvcon'].linearCon):
            #print('linear_con', name)
            self.add_output(name, shape=con.ncon)
            jac = funcsSens[name]
            for wrt_var, subjac in iteritems(jac):
                self.declare_partials(of=name, wrt=wrt_var, val=subjac)



    def _set_geo(self, inputs):
        tmp = {}
        for (args, kwargs) in self.dvs:
            name = args[0]
            tmp[name] = inputs[name]

        self.metadata['dvgeo'].setDesignVars(tmp)
        #self.metadata['dvcon'].setDesignVars(tmp)

    def compute(self, inputs, outputs):
        
        if self._do_solve: 
            self._set_geo(inputs)

        funcs = {}
        self.metadata['dvcon'].evalFunctions(funcs, includeLinear=True)

        for cons in itervalues(self.metadata['dvcon'].constraints):
            for name, con in iteritems(cons):
                #print('nonlinear foobar', con)
                outputs[name] = funcs[name]

        for name, con in iteritems(self.metadata['dvcon'].linearCon):
            #print('linear foobar', con)
            outputs[name] = funcs[name]

    def compute_partials(self, inputs, J):
        
        if self._do_solve:
            self._set_geo(inputs)

        funcsSens = {}

        self.metadata['dvcon'].evalFunctionsSens(funcsSens)

        for of_var, jac in funcsSens.items():
            #print('of: ', of_var)
            for wrt_var, subjac in jac.items():
                #print('  wrt: ', wrt_var)
                J[of_var, wrt_var] = subjac
