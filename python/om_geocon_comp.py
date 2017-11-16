import numpy as np

from six import iteritems, itervalues

from openmdao.api import ExplicitComponent

from pygeo import DVConstraints, DVGeometry

from mach_opt_prob_utils import get_dvs_and_cons

class GeometryConstraintsComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('dvgeo', type_=DVGeometry, required=True)
        self.metadata.declare('dvcon', type_=DVConstraints, required=True)

    def setup(self):
        self.dvs, _ = get_dvs_and_cons(self.metadata['dvgeo'])
        for (args, kwargs) in self.dvs:
            name = args[0]
            size = args[1]
            self.add_input(name, shape=size)
            
        for cons in itervalues(self.metadata['dvcon'].constraints):
            for name, con in iteritems(cons):
                #print('nonlinear con', name)
                self.add_output(name, shape=con.nCon)

        funcsSens = {}
        self.metadata['dvcon'].evalFunctionsSens(funcsSens, includeLinear=True)

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
        self._set_geo(inputs)
        
        funcsSens = {}

        self.metadata['dvcon'].evalFunctionsSens(funcsSens)
        
        for of_var, jac in funcsSens.items():
            #print('of: ', of_var)
            for wrt_var, subjac in jac.items(): 
                #print('  wrt: ', wrt_var)
                J[of_var, wrt_var] = subjac
