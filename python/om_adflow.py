import types

from pygeo import DVGeometry, DVConstraints
from baseclasses import AeroProblem

from pywarpustruct import USMesh
from adflow import ADFLOW

from openmdao.api import Group, PetscKSP, LinearRunOnce, LinearUserDefined, IndepVarComp

from .om_states_comp import OM_STATES_COMP
from .om_func_comp import OM_FUNC_COMP
from .om_geocon_comp import OM_GEOCON_COMP

from .om_utils import get_dvs_and_cons


class OM_ADFLOW(Group):
    """Group integrating the states, funcs, and geometry constraints into a single system"""

    def initialize(self):
        self.metadata.declare('ap', types=AeroProblem)
        self.metadata.declare('aero_options', types=dict)
        self.metadata.declare('mesh_options', types=dict, default={})   
        # self.metadata.declare('family_groups', types=dict, default={})
        self.metadata.declare('adflow_setup_cb', types=types.FunctionType, allow_none=True, default=None)
        self.metadata.declare('dvgeo', types=DVGeometry, allow_none=True, default=None)
        self.metadata.declare('dvcon', types=DVConstraints, allow_none=True, default=None)
        self.metadata.declare('use_OM_KSP', default=False, types=bool, 
            desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")
        self.metadata.declare('owns_indeps', default=False, 
            desc="when True will create IndepVarComp with outputs for all des vars")
        self.metadata.declare('debug', default=False)

        # This option is mainly used in testing when the solver object needs to be run by itself 
        #     and also used in the wrapper. You probably don't want to use this 
        #     unless you really know what you're doing!!!
        self.metadata.declare('solver', types=ADFLOW, allow_none=True, default=None,
            desc='optional argument to allow an existing solver instance to be passed in.')

        # self.metadata.declare('max_procs', default=64, types=int)

    def setup(self):
        ap = self.metadata['ap']
        dvgeo = self.metadata['dvgeo']
        dvcon = self.metadata['dvcon']

        solver = self.metadata['solver']
        if solver is None: 
            solver = ADFLOW(options=self.metadata['aero_options'], 
                            comm=self.comm, 
                            debug=self.metadata['debug'])
        print('foobar', solver)
        self.solver = solver

        # for fg in self.metadata['family_groups']:
        #     solver.addFamilyGroup(fg[0], fg[1])

        if self.metadata['adflow_setup_cb'] is not None: 
           self.metadata['adflow_setup_cb'](solver)  
        # for f_name, f_meta in self.metadata['functions'].items():
        #     fg = f_meta[0]
        #     f_type = f_meta[1]
        #     solver.addFunction(f_type, fg, f_name)

        # necessary to get some memory properly allocated
        solver.getResidual(ap)


        if self.metadata['mesh_options']: 
            mesh = USMesh(options=self.metadata['mesh_options'], comm=self.comm)
            solver.setMesh(mesh)

        if dvgeo is not None: 
            solver.setDVGeo(dvgeo)

        des_vars, constraints = get_dvs_and_cons(ap, dvgeo, dvcon)
        geo_vars, _ = get_dvs_and_cons(geo=dvgeo)

        des_var_names = [dv[0][0] for dv in des_vars]
        geo_var_names = [dv[0][0] for dv in geo_vars]

        if self.metadata['owns_indeps']: 
            indeps = self.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
            for (args, kwargs) in des_vars: 
                name = args[0]
                size = args[1]
                value= kwargs['value']
                if 'units' in kwargs: 
                    indeps.add_output(name, value, units=kwargs['units'])
                else: 
                    indeps.add_output(name, value)

        if dvcon is not None:

            geocon = OM_GEOCON_COMP(dvgeo=dvgeo, dvcon=dvcon)
            self.add_subsystem('geocon', geocon, promotes_inputs=geo_var_names)

            #for cons in dvcon.constraints.values():
            #    for name, con in cons.items():
            #        self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, scaler=con.scale)
            #for name, con in dvcon.linearCon.items():
            #    self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, linear=True)

        states = OM_STATES_COMP(ap=ap, dvgeo=dvgeo, solver=solver,
                            use_OM_KSP=self.metadata['use_OM_KSP'])

        self.add_subsystem('states', states, promotes_inputs=des_var_names)

        # this lets the OpenMDAO KSPsolver converge the adjoint using the ADflow PC
        if self.metadata['use_OM_KSP']:
            states.linear_solver = PetscKSP(iprint=2, atol=1e-8, rtol=1e-8, maxiter=300, ksp_type='gmres')
            states.linear_solver.precon = LinearUserDefined()

        functionals = OM_FUNC_COMP(ap=ap, dvgeo=dvgeo, solver=solver)
        self.add_subsystem('functionals', functionals,
                           promotes_inputs=des_var_names)
        self.connect('states.states', 'functionals.states')
