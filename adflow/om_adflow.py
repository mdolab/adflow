import types
from pygeo import DVGeometry, DVConstraints
from baseclasses import AeroProblem

from idwarp import USMesh
from adflow import ADFLOW

from openmdao.api import Group, PETScKrylov, LinearRunOnce, LinearUserDefined, IndepVarComp

from .om_states_comp import OM_STATES_COMP
from .om_func_comp import OM_FUNC_COMP
from .om_geocon_comp import OM_GEOCON_COMP

from .om_utils import get_dvs_and_cons


class OM_ADFLOW(Group):
    """Group integrating the states, funcs, and geometry constraints into a single system"""

    def initialize(self):
        self.options.declare('ap', types=AeroProblem)
        self.options.declare('setup_cb', types=types.FunctionType, allow_none=True, default=None)

        self.options.declare('use_OM_KSP', default=False, types=bool,
            desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")

        self.options.declare('owns_indeps', default=False,
            desc="when True will create IndepVarComp with outputs for all des vars")
        self.options.declare('debug', default=False)

    def setup(self):
        ap = self.options['ap']

        self.solver, self.mesh, self.dvgeo, self.dvcon = self.options['setup_cb'](self.comm)

        # necessary to get some memory properly allocated
        #self.solver.getResidual(ap)

        des_vars, constraints = get_dvs_and_cons(ap, self.dvgeo, self.dvcon)
        geo_vars, _ = get_dvs_and_cons(geo=self.dvgeo)

        self.des_var_names = des_vars
        self.geo_var_names = geo_vars
        self.constraint_names = constraints

        des_var_names = [dv[0][0] for dv in des_vars]
        geo_var_names = [dv[0][0] for dv in geo_vars]

        if self.options['owns_indeps']:
            indeps = self.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
            for (args, kwargs) in des_vars:
                name = args[0]
                size = args[1]
                value= kwargs['value']
                if 'units' in kwargs:
                    indeps.add_output(name, value, units=kwargs['units'])
                else:
                    indeps.add_output(name, value)

        states = OM_STATES_COMP(ap=ap, dvgeo=self.dvgeo, solver=self.solver,
                                use_OM_KSP=self.options['use_OM_KSP'])

        self.add_subsystem('states', states, promotes_inputs=des_var_names)

        # this lets the OpenMDAO KSPsolver converge the adjoint using the ADflow PC
        if self.options['use_OM_KSP']:
            states.linear_solver = PETScKrylov(iprint=2, atol=1e-8, rtol=1e-8, maxiter=300, ksp_type='gmres')
            states.linear_solver.precon = LinearUserDefined()

        functionals = OM_FUNC_COMP(ap=ap, dvgeo=self.dvgeo, solver=self.solver)
        self.add_subsystem('functionals', functionals,
                           promotes_inputs=des_var_names, max_procs=self.comm.size)
        self.connect('states.states', 'functionals.states')


        if self.dvcon is not None:

            geocon = OM_GEOCON_COMP(dvgeo=self.dvgeo, dvcon=self.dvcon)
            self.add_subsystem('geocon', geocon, promotes_inputs=geo_var_names)

            for cons in self.dvcon.constraints.values():
                for name, con in cons.items():
                    if self.comm.rank == 0:
                        print('adding nonlinear dvcon constarint: ', name)
                    self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, ref=1/con.scale, vectorize_derivs=True)
            for name, con in self.dvcon.linearCon.items():
                if self.comm.rank == 0:
                    print('adding linear dvcon constraint: ', name)
                self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, linear=True, vectorize_derivs=True)

