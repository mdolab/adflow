from pygeo import DVGeometry, DVConstraints
from baseclasses import AeroProblem

from pywarpustruct import USMesh
from adflow import ADFLOW

from openmdao.api import Group, PetscKSP, LinearRunOnce, LinearUserDefined

from comp_states import StatesComp
from comp_functionals import FunctionalsComp
from comp_geocon import GeometryConstraintsComp

from mach_opt_prob_utils import get_dvs_and_cons

class OM_ADFLOW(Group):

    def initialize(self):
        self.metadata.declare('aero_options', required=True)
        self.metadata.declare('mesh_options', required=True)
        self.metadata.declare('family_groups', default={})
        self.metadata.declare('functions', default={})
        self.metadata.declare('ap', type_=AeroProblem, required=True)
        self.metadata.declare('dvgeo', type_=DVGeometry, required=True)
        self.metadata.declare('dvcon', default=None) # TODO: type check this, but still allow None
        self.metadata.declare('use_OM_solver', default=False, type_=bool)
        self.metadata.declare('max_procs', default=64, type_=int)

    def setup(self):
        ap = self.metadata['ap']
        dvgeo = self.metadata['dvgeo']
        dvcon = self.metadata['dvcon']

        solver = ADFLOW(options=self.metadata['aero_options'], comm=self.comm)
        for fg in self.metadata['family_groups']:
            solver.addFamilyGroup(fg[0], fg[1])
        for f_name, f_meta in self.metadata['functions'].items():
            fg = f_meta[0]
            f_type = f_meta[1]
            solver.addFunction(f_type, fg, f_name)

        mesh = USMesh(options=self.metadata['mesh_options'], comm=self.comm)
        solver.setDVGeo(dvgeo)
        solver.setMesh(mesh)

        des_vars, constraints = get_dvs_and_cons(ap, dvgeo, dvcon)
        geo_vars, _ = get_dvs_and_cons(geo=dvgeo)

        des_var_names = [dv[0][0] for dv in des_vars]
        geo_var_names = [dv[0][0] for dv in geo_vars]

        if dvcon:

            geocon = GeometryConstraintsComp(dvgeo=dvgeo, dvcon=dvcon)
            self.add_subsystem('geocon', geocon, promotes_inputs=geo_var_names)

            #for cons in dvcon.constraints.values():
            #    for name, con in cons.items():
            #        self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, scaler=con.scale)
            #for name, con in dvcon.linearCon.items():
            #    self.add_constraint('geocon.%s' % name, lower=con.lower, upper=con.upper, linear=True)

        states = StatesComp(ap=ap, dvgeo=dvgeo, solver=solver,
                            max_procs=self.metadata['max_procs'], use_OM_solver=self.metadata['use_OM_solver'])

        self.add_subsystem('states', states, promotes_inputs=des_var_names)

        if self.metadata['use_OM_solver']:
            states.linear_solver = PetscKSP(iprint=2, atol=1e-8, rtol=1e-8, maxiter=300, ksp_type='gmres')
            states.linear_solver.precon = LinearUserDefined()

        functionals = FunctionalsComp(ap=ap, dvgeo=dvgeo, solver=solver, max_procs=self.metadata['max_procs'])
        self.add_subsystem('functionals', functionals,
                           promotes_inputs=des_var_names)
        self.connect('states.states', 'functionals.states')
