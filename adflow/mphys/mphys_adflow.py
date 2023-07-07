from pprint import pprint as pp
import numpy as np
from adflow import ADFLOW
from idwarp import USMesh, MultiUSMesh
from mphys.builder import Builder
from openmdao.api import AnalysisError, ExplicitComponent, Group, ImplicitComponent


from .om_utils import get_dvs_and_cons


class ADflowMesh(ExplicitComponent):
    """
    Component to get the partitioned initial surface mesh coordinates

    """

    def initialize(self):
        self.options.declare("aero_solver", recordable=False)

    def setup(self):
        self.aero_solver = self.options["aero_solver"]

        self.x_a0 = self.aero_solver.getSurfaceCoordinates(
            groupName=self.aero_solver.meshFamilyGroup, includeZipper=False
        ).flatten(order="C")

        coord_size = self.x_a0.size
        self.add_output(
            "x_aero0",
            distributed=True,
            shape=coord_size,
            desc="initial aerodynamic surface node coordinates",
            tags=["mphys_coordinates"],
        )

    def mphys_add_coordinate_input(self):
        self.add_input(
            "x_aero0_points", distributed=True, shape_by_conn=True, desc="aerodynamic surface with geom changes"
        )

        # return the promoted name and coordinates
        return "x_aero0_points", self.x_a0

    def mphys_get_surface_mesh(self):
        return self.x_a0

    def mphys_get_triangulated_surface(self, groupName=None):
        # this is a list of lists of 3 points
        # p0, v1, v2

        return self._getTriangulatedMeshSurface(groupName=groupName)

    def _getTriangulatedMeshSurface(self, groupName=None, **kwargs):
        """
        This function returns a trianguled verision of the surface
        mesh on all processors. The intent is to use this for doing
        constraints in DVConstraints.

        Returns
        -------
        surf : list
           List of points and vectors describing the surface. This may
           be passed directly to DVConstraint setSurface() function.
        """

        if groupName is None:
            groupName = self.aero_solver.allWallsGroup

        # Obtain the points and connectivity for the specified
        # groupName
        pts = self.aero_solver.comm.allgather(self.aero_solver.getSurfaceCoordinates(groupName, **kwargs))
        conn, faceSizes = self.aero_solver.getSurfaceConnectivity(groupName)
        conn = np.array(conn).flatten()
        conn = self.aero_solver.comm.allgather(conn)
        faceSizes = self.aero_solver.comm.allgather(faceSizes)
        # Triangle info...point and two vectors
        p0 = []
        v1 = []
        v2 = []

        # loop over the faces
        for iProc in range(len(faceSizes)):
            connCounter = 0
            for iFace in range(len(faceSizes[iProc])):
                # Get the number of nodes on this face
                faceSize = faceSizes[iProc][iFace]
                faceNodes = conn[iProc][connCounter : connCounter + faceSize]
                # Start by getting the centerpoint
                ptSum = [0, 0, 0]
                for i in range(faceSize):
                    # idx = ptCounter+i
                    idx = faceNodes[i]
                    ptSum += pts[iProc][idx]

                avgPt = ptSum / faceSize

                # Now go around the face and add a triangle for each adjacent pair
                # of points. This assumes an ordered connectivity from the
                # meshwarping
                for i in range(faceSize):
                    idx = faceNodes[i]
                    p0.append(avgPt)
                    v1.append(pts[iProc][idx] - avgPt)
                    if i < (faceSize - 1):
                        idxp1 = faceNodes[i + 1]
                        v2.append(pts[iProc][idxp1] - avgPt)
                    else:
                        # wrap back to the first point for the last element
                        idx0 = faceNodes[0]
                        v2.append(pts[iProc][idx0] - avgPt)

                # Now increment the connectivity
                connCounter += faceSize
        return [p0, v1, v2]

    def compute(self, inputs, outputs):
        if "x_aero0_points" in inputs:
            outputs["x_aero0"] = inputs["x_aero0_points"]
        else:
            outputs["x_aero0"] = self.x_a0

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            if "x_aero0_points" in d_inputs:
                d_outputs["x_aero0"] += d_inputs["x_aero0_points"]
        elif mode == "rev":
            if "x_aero0_points" in d_inputs:
                d_inputs["x_aero0_points"] += d_outputs["x_aero0"]


class ADflowWarper(ExplicitComponent):
    """
    OpenMDAO component that wraps the warping.

    """

    def initialize(self):
        self.options.declare("aero_solver", recordable=False)
        # self.options.declare('use_OM_KSP', default=False, types=bool,
        #    desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")

    def setup(self):
        # self.set_check_partial_options(wrt='*',directional=True)

        self.solver = self.options["aero_solver"]
        solver = self.solver

        # self.ap_vars,_ = get_dvs_and_cons(ap=ap)

        # state inputs and outputs
        local_volume_coord_size = solver.mesh.getSolverGrid().size

        self.add_input("x_aero", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])
        self.add_output("adflow_vol_coords", distributed=True, shape=local_volume_coord_size, tags=["mphys_coupling"])

        # self.declare_partials(of='adflow_vol_coords', wrt='x_aero')

    def compute(self, inputs, outputs):
        solver = self.solver

        x_a = inputs["x_aero"].reshape((-1, 3))
        solver.setSurfaceCoordinates(x_a, groupName=solver.meshFamilyGroup)
        solver.updateGeometryInfo()
        outputs["adflow_vol_coords"] = solver.mesh.getSolverGrid()

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            if "adflow_vol_coords" in d_outputs:
                if "x_aero" in d_inputs:
                    dxS = d_inputs["x_aero"]
                    dxV = self.solver.mesh.warpDerivFwd(dxS)
                    d_outputs["adflow_vol_coords"] += dxV

        elif mode == "rev":
            if "adflow_vol_coords" in d_outputs:
                if "x_aero" in d_inputs:
                    dxV = d_outputs["adflow_vol_coords"]
                    self.solver.mesh.warpDeriv(dxV)
                    dxS = self.solver.mesh.getdXs()
                    dxS = self.solver.mapVector(
                        dxS, self.solver.meshFamilyGroup, self.solver.designFamilyGroup, includeZipper=False
                    )
                    d_inputs["x_aero"] += dxS.flatten()


class ADflowSolver(ImplicitComponent):
    """
    OpenMDAO component that wraps the ADflow flow solver

    """

    def initialize(self):
        self.options.declare("aero_solver", recordable=False)
        # self.options.declare('use_OM_KSP', default=False, types=bool,
        #    desc="uses OpenMDAO's PestcKSP linear solver with ADflow's preconditioner to solve the adjoint.")
        self.options.declare("restart_failed_analysis", default=False)
        self.options.declare("err_on_convergence_fail", default=False)

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True
        self.analysis_error_on_failure = True

    def setup(self):
        # self.set_check_partial_options(wrt='*',directional=True)

        self.restart_failed_analysis = self.options["restart_failed_analysis"]
        self.err_on_convergence_fail = self.options["err_on_convergence_fail"]
        self.solver = self.options["aero_solver"]
        solver = self.solver

        # this is the solution counter for failed solution outputs.
        # the converged solutions are written by the adflow functionals group
        self.solution_counter = 0

        # flag to keep track if the current solution started from a clean restart,
        # or it was restarted from the previous converged state.
        self.cleanRestart = True

        # state inputs and outputs
        local_state_size = solver.getStateSize()

        self.add_input("adflow_vol_coords", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])
        self.add_output("adflow_states", distributed=True, shape=local_state_size, tags=["mphys_coupling"])

        # self.declare_partials(of='adflow_states', wrt='*')

    def _set_ap(self, inputs, print_dict=True):
        tmp = {}
        for args, _kwargs in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]

        # enable if you want to print all aero dv inputs
        # if self.comm.rank == 0:
        #     print('aero dv inputs:')
        #     pp(tmp)

        self.ap.setDesignVars(tmp)
        if self.comm.rank == 0 and print_dict:
            pp(tmp)

    def set_ap(self, ap):
        # this is the external function to set the ap to this component
        self.ap = ap

        self.ap_vars, _ = get_dvs_and_cons(ap=ap)

        # parameter inputs
        if self.comm.rank == 0:
            print("adding ap var inputs")
        for args, kwargs in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(
                name, distributed=False, shape=size, val=kwargs["value"], units=kwargs["units"], tags=["mphys_input"]
            )
            if self.comm.rank == 0:
                print("%s (%s)" % (name, kwargs["units"]))

    def _set_states(self, outputs):
        self.solver.setStates(outputs["adflow_states"])

    def apply_nonlinear(self, inputs, outputs, residuals):
        solver = self.solver

        self._set_states(outputs)
        self._set_ap(inputs, print_dict=False)

        ap = self.ap

        # Set the warped mesh
        # solver.mesh.setSolverGrid(inputs['adflow_vol_coords'])
        # ^ This call does not exist. Assume the mesh hasn't changed since the last call to the warping comp for now
        # TODO we must fix this. we need to put in the modified volume coordinates

        # flow residuals
        residuals["adflow_states"] = solver.getResidual(ap)

    def solve_nonlinear(self, inputs, outputs):
        solver = self.solver
        ap = self.ap
        if self._do_solve:
            # Set the warped mesh
            # solver.mesh.setSolverGrid(inputs['adflow_vol_coords'])
            # ^ This call does not exist. Assume the mesh hasn't changed since the last call to the warping comp for now

            self._set_ap(inputs)
            ap.solveFailed = False  # might need to clear this out?
            ap.fatalFail = False

            # do not write solution files inside the solver loop
            solver(ap, writeSolution=False)

            if ap.fatalFail:
                if self.comm.rank == 0:
                    print("###############################################################")
                    print("# Solve Fatal Fail. Analysis Error")
                    print("###############################################################")

                raise AnalysisError("ADFLOW Solver Fatal Fail")

            if ap.solveFailed:
                if self.restart_failed_analysis:  # the mesh was fine, but it didn't converge
                    # if the previous iteration was already a clean restart, dont try again
                    if self.cleanRestart:
                        if self.comm.rank == 0:
                            print("###############################################################")
                            print("# This was a clean restart. Will not try another one.")
                            print("###############################################################")

                        # write the solution so that we can diagnose
                        solver.writeSolution(baseName="analysis_fail", number=self.solution_counter)
                        self.solution_counter += 1

                        if self.analysis_error_on_failure:
                            solver.resetFlow(ap)
                            self.cleanRestart = True
                            raise AnalysisError("ADFLOW Solver Fatal Fail")

                    # the previous iteration restarted from another solution, so we can try again
                    # with a re-set flowfield for the initial guess.
                    else:
                        if self.comm.rank == 0:
                            print("###############################################################")
                            print("# Solve Failed, attempting a clean restart!")
                            print("###############################################################")

                        # write the solution so that we can diagnose
                        solver.writeSolution(baseName="analysis_fail", number=self.solution_counter)
                        self.solution_counter += 1

                        ap.solveFailed = False
                        ap.fatalFail = False
                        solver.resetFlow(ap)
                        solver(ap, writeSolution=False)

                        if ap.solveFailed or ap.fatalFail:  # we tried, but there was no saving it
                            if self.comm.rank == 0:
                                print("###############################################################")
                                print("# Clean Restart failed. There is no saving this one!")
                                print("###############################################################")

                            # write the solution so that we can diagnose
                            solver.writeSolution(baseName="analysis_fail", number=self.solution_counter)
                            self.solution_counter += 1

                            if self.analysis_error_on_failure:
                                # re-set the flow for the next iteration:
                                solver.resetFlow(ap)
                                # set the reset flow flag
                                self.cleanRestart = True
                                raise AnalysisError("ADFLOW Solver Fatal Fail")

                        # see comment for the same flag below
                        else:
                            self.cleanRestart = False

                elif self.err_on_convergence_fail:
                    # the solve failed but we dont want to re-try. We also want to raise an analysis error
                    if self.comm.rank == 0:
                        print("###############################################################")
                        print("# Solve Failed, not attempting a clean restart")
                        print("###############################################################")
                    raise AnalysisError("ADFLOW Solver Fatal Fail")

                else:
                    # the solve failed, but we dont want to restart or raise an error, we will just let this one pass
                    pass

            # solve did not fail, therefore we will re-use this converged flowfield for the next iteration.
            # change the flag so that if the next iteration fails with current initial guess, it can retry
            # with a clean restart
            else:
                self.cleanRestart = False

        outputs["adflow_states"] = solver.getStates()

    def linearize(self, inputs, outputs, residuals):
        solver = self.solver
        ap = self.ap

        self._set_ap(inputs, print_dict=False)
        self._set_states(outputs)

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):
        solver = self.solver
        ap = self.ap

        self._set_ap(inputs, print_dict=False)
        self._set_states(outputs)

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

        if mode == "fwd":
            if "adflow_states" in d_residuals:
                xDvDot = {}
                for var_name in d_inputs:
                    xDvDot[var_name] = d_inputs[var_name]
                if "adflow_vol_coords" in d_inputs:
                    xVDot = d_inputs["adflow_vol_coords"]
                else:
                    xVDot = None
                if "adflow_states" in d_outputs:
                    wDot = d_outputs["adflow_states"]
                else:
                    wDot = None

                dwdot = solver.computeJacobianVectorProductFwd(
                    xDvDot=xDvDot, xVDot=xVDot, wDot=wDot, residualDeriv=True
                )
                d_residuals["adflow_states"] += dwdot

        elif mode == "rev":
            if "adflow_states" in d_residuals:
                resBar = d_residuals["adflow_states"]

                wBar, xVBar, xDVBar = solver.computeJacobianVectorProductBwd(
                    resBar=resBar, wDeriv=True, xVDeriv=True, xDvDeriv=False, xDvDerivAero=True
                )

                if "adflow_states" in d_outputs:
                    d_outputs["adflow_states"] += wBar

                if "adflow_vol_coords" in d_inputs:
                    d_inputs["adflow_vol_coords"] += xVBar

                for dv_name, dv_bar in xDVBar.items():
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += dv_bar.flatten()

    def solve_linear(self, d_outputs, d_residuals, mode):
        solver = self.solver
        ap = self.ap

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

        # the adjoint might not be set up regardless if we changed APs
        # this is because the first call with any AP will not have this set up, so we have to check
        # if we changed APs, then we also freed adjoint memory,
        # and then again we would need to setup adjoint again
        # finally, we generally want to avoid extra calls here
        # because this routine can be call multiple times back to back in a LBGS solver.
        if not solver.adjointSetup:
            solver._setupAdjoint()

        if self.comm.rank == 0:
            print("Solving linear in mphys_adflow", flush=True)
        if mode == "fwd":
            d_outputs["adflow_states"] = solver.solveDirectForRHS(d_residuals["adflow_states"])
        elif mode == "rev":
            # d_residuals['adflow_states'] = solver.solveAdjointForRHS(d_outputs['adflow_states'])
            solver.adflow.adjointapi.solveadjoint(d_outputs["adflow_states"], d_residuals["adflow_states"], True)

        return True, 0, 0


class ADflowForces(ExplicitComponent):
    """
    OpenMDAO component that wraps force integration

    """

    def initialize(self):
        self.options.declare("aero_solver", recordable=False)

    def setup(self):
        # self.set_check_partial_options(wrt='*',directional=True)

        self.solver = self.options["aero_solver"]
        solver = self.solver

        self.add_input("adflow_vol_coords", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])
        self.add_input("adflow_states", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])

        local_surface_coord_size = solver.mesh.getSurfaceCoordinates().size
        self.add_output("f_aero", distributed=True, shape=local_surface_coord_size, tags=["mphys_coupling"])

        # self.declare_partials(of='f_aero', wrt='*')

    def _set_ap(self, inputs, print_dict=True):
        tmp = {}
        for args, _kwargs in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]

        self.ap.setDesignVars(tmp)
        if self.comm.rank == 0 and print_dict:
            pp(tmp)

    def set_ap(self, ap):
        # this is the external function to set the ap to this component
        self.ap = ap

        self.ap_vars, _ = get_dvs_and_cons(ap=ap)

        # parameter inputs
        # if self.comm.rank == 0:
        #     print('adding ap var inputs:')
        for args, kwargs in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, distributed=False, shape=size, units=kwargs["units"], tags=["mphys_input"])
            # if self.comm.rank == 0:
            #     print('%s (%s)'%(name, kwargs['units']))

    def _set_states(self, inputs):
        self.solver.setStates(inputs["adflow_states"])

    def compute(self, inputs, outputs):
        solver = self.solver

        self._set_ap(inputs)

        # Set the warped mesh
        # solver.mesh.setSolverGrid(inputs['adflow_vol_coords'])
        # ^ This call does not exist. Assume the mesh hasn't changed since the last call to the warping comp for now
        # TODO we must fix this. we need to put in the modified volume coordinates
        self._set_states(inputs)

        outputs["f_aero"] = solver.getForces().flatten(order="C")

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        solver = self.solver
        ap = self.ap

        self._set_ap(inputs, print_dict=False)

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

        if mode == "fwd":
            if "f_aero" in d_outputs:
                xDvDot = {}
                for var_name in d_inputs:
                    xDvDot[var_name] = d_inputs[var_name]
                if "adflow_states" in d_inputs:
                    wDot = d_inputs["adflow_states"]
                else:
                    wDot = None
                if "adflow_vol_coords" in d_inputs:
                    xVDot = d_inputs["adflow_vol_coords"]
                else:
                    xVDot = None
                if not (xVDot is None and wDot is None):
                    dfdot = solver.computeJacobianVectorProductFwd(xDvDot=xDvDot, xVDot=xVDot, wDot=wDot, fDeriv=True)
                    d_outputs["f_aero"] += dfdot.flatten()

        elif mode == "rev":
            if "f_aero" in d_outputs:
                fBar = d_outputs["f_aero"]

                wBar, xVBar, xDVBar = solver.computeJacobianVectorProductBwd(
                    fBar=fBar, wDeriv=True, xVDeriv=True, xDvDeriv=False, xDvDerivAero=True
                )

                if "adflow_vol_coords" in d_inputs:
                    d_inputs["adflow_vol_coords"] += xVBar
                if "adflow_states" in d_inputs:
                    d_inputs["adflow_states"] += wBar

                for dv_name, dv_bar in xDVBar.items():
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += dv_bar.flatten()


class AdflowHeatTransfer(ExplicitComponent):
    """
    OpenMDAO component that wraps heat transfer integration

    """

    def initialize(self):
        self.options.declare("aero_solver")

    def setup(self):
        # self.set_check_partial_options(wrt='*',directional=True)

        self.solver = self.options["aero_solver"]
        solver = self.solver

        local_nodes, _ = solver._getSurfaceSize(solver.allIsothermalWallsGroup)

        self.add_input("adflow_vol_coords", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])
        self.add_input("adflow_states", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])

        self.add_output(
            "q_convect",
            distributed=True,
            val=np.ones(local_nodes) * -499,
            shape=local_nodes,
            units="W/m**2",
            tags=["mphys_coupling"],
        )

        # self.declare_partials(of='f_aero', wrt='*')

    def _set_ap(self, inputs):
        tmp = {}
        for args, _kwargs in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]

        self.ap.setDesignVars(tmp)

    def set_ap(self, ap):
        # this is the external function to set the ap to this component
        self.ap = ap

        self.ap_vars, _ = get_dvs_and_cons(ap=ap)

        # parameter inputs
        if self.comm.rank == 0:
            print("adding ap var inputs")
        for args, kwargs in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, distributed=False, shape=size, units=kwargs["units"], tags=["mphys_input"])
            if self.comm.rank == 0:
                print(name)

    def _set_states(self, inputs):
        self.solver.setStates(inputs["adflow_states"])

    def compute(self, inputs, outputs):
        solver = self.solver

        ## already done by solver
        self._set_ap(inputs)

        # Set the warped mesh
        # solver.mesh.setSolverGrid(inputs['adflow_vol_coords'])
        # ^ This call does not exist. Assume the mesh hasn't changed since the last call to the warping comp for now
        # TODO we must fix this. we need to put in the modified volume coordinates

        #
        # self._set_states(inputs)

        outputs["q_convect"] = solver.getHeatFluxes().flatten(order="C")
        # print()

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        solver = self.solver
        ap = self.options["ap"]

        self._set_ap(inputs)

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

        if mode == "fwd":
            if "q_convect" in d_outputs:
                xDvDot = {}
                for var_name in d_inputs:
                    xDvDot[var_name] = d_inputs[var_name]
                if "adflow_states" in d_inputs:
                    wDot = d_inputs["adflow_states"]
                else:
                    wDot = None
                if "adflow_vol_coords" in d_inputs:
                    xVDot = d_inputs["adflow_vol_coords"]
                else:
                    xVDot = None
                if not (xVDot is None and wDot is None):
                    dhfdot = solver.computeJacobianVectorProductFwd(xDvDot=xDvDot, xVDot=xVDot, wDot=wDot, hfDeriv=True)
                    dhfdot_map = np.zeros((dhfdot.size, 3))
                    dhfdot_map[:, 0] = dhfdot.flatten()
                    dhfdot_map = self.solver.mapVector(
                        dhfdot_map, self.solver.allWallsGroup, self.solver.allIsothermalWallsGroup
                    )
                    dhfdot = dhfdot_map[:, 0]
                    d_outputs["q_convect"] += dhfdot

        elif mode == "rev":
            if "q_convect" in d_outputs:
                hfBar = d_outputs["q_convect"]

                hfBar_map = np.zeros((hfBar.size, 3))
                hfBar_map[:, 0] = hfBar.flatten()
                hfBar_map = self.solver.mapVector(
                    hfBar_map, self.solver.allIsothermalWallsGroup, self.solver.allWallsGroup
                )
                hfBar = hfBar_map[:, 0]

                wBar, xVBar, xDVBar = solver.computeJacobianVectorProductBwd(
                    hfBar=hfBar, wDeriv=True, xVDeriv=True, xDvDeriv=True
                )

                if "adflow_vol_coords" in d_inputs:
                    d_inputs["adflow_vol_coords"] += xVBar
                if "adflow_states" in d_inputs:
                    d_inputs["adflow_states"] += wBar

                for dv_name, dv_bar in xDVBar.items():
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += dv_bar.flatten()


FUNCS_UNITS = {
    "mdot": "kg/s",
    "mavgptot": "Pa",
    "mavgps": "Pa",
    "aavgptot": "Pa",
    "aavgps": "Pa",
    "mavgttot": "degK",
    "mavgvx": "m/s",
    "mavgvy": "m/s",
    "mavgvz": "m/s",
    "drag": "N",
    "lift": "N",
    "dragpressure": "N",
    "dragviscous": "N",
    "dragmomentum": "N",
    "fx": "N",
    "fy": "N",
    "fz": "N",
    "forcexpressure": "N",
    "forceypressure": "N",
    "forcezpressure": "N",
    "forcexviscous": "N",
    "forceyviscous": "N",
    "forcezviscous": "N",
    "forcexmomentum": "N",
    "forceymomentum": "N",
    "forcezmomentum": "N",
    "flowpower": "W",
    "area": "m**2",
}


class ADflowFunctions(ExplicitComponent):
    def initialize(self):
        self.options.declare("aero_solver", recordable=False)
        # flag to automatically add the AP functions as output
        self.options.declare("ap_funcs", default=True)
        self.options.declare("write_solution", default=True)

        # testing flag used for unit-testing to prevent the call to actually solve
        # NOT INTENDED FOR USERS!!! FOR TESTING ONLY
        self._do_solve = True

        self.extra_funcs = None

    def setup(self):
        self.solver = self.options["aero_solver"]
        self.ap_funcs = self.options["ap_funcs"]
        self.write_solution = self.options["write_solution"]
        # self.set_check_partial_options(wrt='*',directional=True)
        self.solution_counter = 0

        self.add_input("adflow_vol_coords", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])
        self.add_input("adflow_states", distributed=True, shape_by_conn=True, tags=["mphys_coupling"])

        # self.declare_partials(of=f_name, wrt='*')

    def _set_ap(self, inputs, print_dict=True):
        tmp = {}
        for args, _kwargs in self.ap_vars:
            name = args[0]
            tmp[name] = inputs[name]

        self.ap.setDesignVars(tmp)
        # self.options['solver'].setAeroProblem(self.options['ap'])
        if self.comm.rank == 0 and print_dict:
            pp(tmp)

    def mphys_set_ap(self, ap):
        # this is the external function to set the ap to this component
        self.ap = ap
        self.ap_vars, _ = get_dvs_and_cons(ap=ap)

        # parameter inputs
        # if self.comm.rank == 0:
        # print('adding ap var inputs:')
        for args, kwargs in self.ap_vars:
            name = args[0]
            size = args[1]
            self.add_input(name, distributed=False, shape=size, units=kwargs["units"], tags=["mphys_input"])
            # if self.comm.rank == 0:
            # print('%s with units %s'%(name, kwargs['units']))

        if self.ap_funcs:
            if self.comm.rank == 0:
                print("adding adflow funcs as output:")

            for f_name in sorted(list(self.ap.evalFuncs)):
                # get the function type. this is the first word before the first underscore
                f_type = f_name.split("_")[0]

                # check if we have a unit defined for this
                if f_type in FUNCS_UNITS:
                    units = FUNCS_UNITS[f_type]
                else:
                    units = None

                # print the function name and units
                if self.comm.rank == 0:
                    print("%s (%s)" % (f_name, units))

                self.add_output(f_name, distributed=False, shape=1, units=units, tags=["mphys_result"])

    # def mphys_add_prop_funcs(self, prop_funcs):
    #     save this list
    #     self.extra_funcs = prop_funcs

    #     if self.comm.rank == 0:
    #         print("adding adflow funcs as propulsion output:", prop_funcs)

    #     call the add_funcs function
    #     self.mphys_add_funcs(prop_funcs)

    def mphys_add_funcs(self, funcs):
        self.extra_funcs = funcs

        # loop over the functions here and create the output
        for f_name in funcs:
            # get the function type. this is the first word before the first underscore
            f_type = f_name.split("_")[0]

            # check if we have a unit defined for this
            if f_type in FUNCS_UNITS:
                units = FUNCS_UNITS[f_type]
            else:
                units = None

            # print the function name and units
            # if self.comm.rank == 0:
            #     print("%s (%s)" % (f_name, units))

            self.add_output(f_name, distributed=False, shape=1, units=units, tags=["mphys_result"])

    def _set_states(self, inputs):
        self.solver.setStates(inputs["adflow_states"])

    def _get_func_name(self, name):
        return "%s_%s" % (self.ap.name, name.lower())

    def nom_write_solution(self, **kwargs):
        # this writes the solution files and is callable from outside openmdao call routines
        solver = self.solver
        ap = self.ap

        # re-set the AP so that we are sure state is updated
        solver.setAeroProblem(ap)

        # write the solution files. Internally, this checks the
        # types of solution files specified in the options and
        # only outsputs these
        solver.writeSolution(number=self.solution_counter, **kwargs)
        self.solution_counter += 1

    def compute(self, inputs, outputs):
        solver = self.solver
        # print('funcs compute')
        # actually setting things here triggers some kind of reset, so we only do it if you're actually solving
        if self._do_solve:
            self._set_ap(inputs)
            # Set the warped mesh
            # solver.mesh.setSolverGrid(inputs['adflow_vol_coords'])
            # ^ This call does not exist. Assume the mesh hasn't changed since the last call to the warping comp for now
            self._set_states(inputs)

        if self.write_solution:
            # write the solution files. Internally, this checks the
            # types of solution files specified in the options and
            # only outsputs these
            solver.writeSolution(number=self.solution_counter)
            self.solution_counter += 1

        funcs = {}

        if self.ap_funcs:
            # without the sorted, each proc might get a different order...
            eval_funcs = sorted(list(self.ap.evalFuncs))
            solver.evalFunctions(self.ap, funcs, evalFuncs=eval_funcs)

            for name in self.ap.evalFuncs:
                f_name = self._get_func_name(name)
                if f_name in funcs:
                    outputs[name.lower()] = funcs[f_name]

        if self.extra_funcs is not None:
            # also do the prop
            solver.evalFunctions(self.ap, funcs, evalFuncs=self.extra_funcs)
            for name in self.extra_funcs:
                f_name = self._get_func_name(name)
                if f_name in funcs:
                    outputs[name.lower()] = funcs[f_name]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        solver = self.solver
        ap = self.ap
        self._set_ap(inputs, print_dict=False)

        # check if we changed APs, then we have to do a bunch of updates
        if ap != solver.curAP:
            # AP is changed, so we have to update the AP and
            # run a residual to make sure all intermediate vairables are up to date
            # we assume the AP has the last converged state information,
            # which is automatically set in the getResidual call
            solver.getResidual(ap)

        if mode == "fwd":
            xDvDot = {}
            for key in ap.DVs:
                if key in d_inputs:
                    mach_name = key.split("_")[0]
                    xDvDot[mach_name] = d_inputs[key]

            if "adflow_states" in d_inputs:
                wDot = d_inputs["adflow_states"]
            else:
                wDot = None

            if "adflow_vol_coords" in d_inputs:
                xVDot = d_inputs["adflow_vol_coords"]
            else:
                xVDot = None

            funcsdot = solver.computeJacobianVectorProductFwd(xDvDot=xDvDot, xVDot=xVDot, wDot=wDot, funcDeriv=True)

            for name in funcsdot:
                func_name = name.lower()
                if name in d_outputs:
                    d_outputs[name] += funcsdot[func_name]

        elif mode == "rev":
            funcsBar = {}

            if self.ap_funcs:
                for name in self.ap.evalFuncs:
                    func_name = name.lower()

                    # we have to check for 0 here, so we don't include any unnecessary variables in funcsBar
                    # becasue it causes ADflow to do extra work internally even if you give it extra variables, even if the seed is 0
                    if func_name in d_outputs and d_outputs[func_name] != 0.0:
                        funcsBar[func_name] = d_outputs[func_name][0]
                        # this stuff is fixed now. no need to divide
                        # funcsBar[func_name] = d_outputs[func_name][0] / self.comm.size
                        # print(self.comm.rank, func_name, funcsBar[func_name])

            # also do the same for prop functions
            if self.extra_funcs is not None:
                for name in self.extra_funcs:
                    func_name = name.lower()
                    if func_name in d_outputs and d_outputs[func_name] != 0.0:
                        funcsBar[func_name] = d_outputs[func_name][0]

            # print(funcsBar, flush=True)

            wBar = None
            xVBar = None
            xDVBar = None

            wBar, xVBar, xDVBar = solver.computeJacobianVectorProductBwd(
                funcsBar=funcsBar, wDeriv=True, xVDeriv=True, xDvDeriv=False, xDvDerivAero=True
            )
            if "adflow_states" in d_inputs:
                d_inputs["adflow_states"] += wBar
            if "adflow_vol_coords" in d_inputs:
                d_inputs["adflow_vol_coords"] += xVBar

            for dv_name, dv_bar in xDVBar.items():
                if dv_name in d_inputs:
                    d_inputs[dv_name] += dv_bar.flatten()


class ADflowGroup(Group):
    def initialize(self):
        self.options.declare("solver", recordable=False)
        self.options.declare("struct_coupling", default=False)
        self.options.declare("prop_coupling", default=False)
        self.options.declare("heat_transfer", default=False)
        self.options.declare("use_warper", default=True)
        self.options.declare("restart_failed_analysis", default=False)
        self.options.declare("err_on_convergence_fail", default=False)
        self.options.declare("balance_group", default=None, recordable=False)

    def setup(self):
        self.aero_solver = self.options["solver"]
        self.struct_coupling = self.options["struct_coupling"]
        self.prop_coupling = self.options["prop_coupling"]
        self.restart_failed_analysis = self.options["restart_failed_analysis"]
        self.err_on_convergence_fail = self.options["err_on_convergence_fail"]

        self.use_warper = self.options["use_warper"]

        balance_group = self.options["balance_group"]
        self.heat_transfer = self.options["heat_transfer"]

        if self.use_warper:
            # if we dont have geo_disp, we also need to promote the x_a as x_a0 from the deformer component
            self.add_subsystem(
                "deformer",
                ADflowWarper(
                    aero_solver=self.aero_solver,
                ),
                promotes_inputs=["x_aero"],
                promotes_outputs=["adflow_vol_coords"],
            )

        self.add_subsystem(
            "solver",
            ADflowSolver(
                aero_solver=self.aero_solver,
                restart_failed_analysis=self.restart_failed_analysis,
                err_on_convergence_fail=self.err_on_convergence_fail,
            ),
            promotes_inputs=["adflow_vol_coords"],
            promotes_outputs=["adflow_states"],
        )

        if self.struct_coupling:
            self.add_subsystem(
                "force",
                ADflowForces(aero_solver=self.aero_solver),
                promotes_inputs=["adflow_vol_coords", "adflow_states"],
                promotes_outputs=["f_aero"],
            )
        if self.prop_coupling:
            self.add_subsystem(
                "prop",
                ADflowFunctions(
                    aero_solver=self.aero_solver,
                    ap_funcs=False,
                    write_solution=False,
                ),
                promotes_inputs=["adflow_vol_coords", "adflow_states"],
            )

        if self.heat_transfer:
            self.add_subsystem(
                "heat_xfer", AdflowHeatTransfer(aero_solver=self.aero_solver), promotes_outputs=["q_convect"]
            )

        if balance_group is not None:
            self.add_subsystem("balance", balance_group)

    def mphys_set_ap(self, ap):
        # set the ap, add inputs and outputs, promote?
        self.solver.set_ap(ap)
        # self.funcs.set_ap(ap)
        if self.struct_coupling:
            self.force.set_ap(ap)
        if self.prop_coupling:
            self.prop.mphys_set_ap(ap)

        if self.heat_transfer:
            self.heat_xfer.set_ap(ap)

        # promote the DVs for this ap
        ap_vars, _ = get_dvs_and_cons(ap=ap)

        for args, _kwargs in ap_vars:
            name = args[0]
            self.promotes("solver", inputs=[name])
            # self.promotes('funcs', inputs=[name])
            if self.struct_coupling:
                self.promotes("force", inputs=[name])
            if self.prop_coupling:
                self.promotes("prop", inputs=[name])

    def mphys_add_prop_funcs(self, prop_funcs):
        # this is the main routine to enable outputs from the propulsion element

        # call the method of the prop element
        self.prop.mphys_add_funcs(prop_funcs)

        # promote these variables to the aero group level
        self.promotes("prop", outputs=prop_funcs)


class ADflowMeshGroup(Group):
    def initialize(self):
        self.options.declare("aero_solver", recordable=False)

    def setup(self):
        aero_solver = self.options["aero_solver"]

        self.add_subsystem("surface_mesh", ADflowMesh(aero_solver=aero_solver), promotes=["*"])
        self.add_subsystem(
            "volume_mesh",
            ADflowWarper(aero_solver=aero_solver),
            promotes_inputs=[("x_aero", "x_aero0")],
            promotes_outputs=["adflow_vol_coords"],
        )

    def mphys_add_coordinate_input(self):
        # just pass through the call
        return self.surface_mesh.mphys_add_coordinate_input()

    def mphys_get_triangulated_surface(self):
        # just pass through the call
        return self.surface_mesh.mphys_get_triangulated_surface()


class ADflowBuilder(Builder):
    def __init__(
        self,
        options,  # adflow options
        mesh_options=None,  # idwarp options
        scenario="aerodynamic",  # scenario type to configure the groups
        mesh_type="USMesh",  # mesh type option. USMesh or  MultiUSMesh
        restart_failed_analysis=False,  # retry after failed analysis
        err_on_convergence_fail=False,  # raise an analysis error if the solver stalls
        balance_group=None,
        user_family_groups=None,  # Dictonary of {group: surfs} to add
    ):
        # options dictionary for ADflow
        self.options = options

        # Check the user family groups to see if anything was set
        if user_family_groups is None:
            self.user_family_groups = None
        else:
            # Throw a type error if the user family group is not a dictionary
            if isinstance(user_family_groups, dict):
                self.user_family_groups = user_family_groups
            else:
                raise TypeError("ADflowBuilder parameter 'user_family_groups' must be a dictionary or NoneType")

        # MACH tools require separate option dictionaries for solver and mesh
        # if user did not provide a separate mesh_options dictionary, we just use
        # the grid file option from the aero options.
        if mesh_options is None:
            if "gridFile" in options:
                self.mesh_options = {
                    "gridFile": options["gridFile"],
                }
            elif "gridfile" in options:
                self.mesh_options = {
                    "gridFile": options["gridfile"],
                }
        else:
            self.mesh_options = mesh_options

        if mesh_type == "USMesh":
            self.multi_us_mesh = False

            if "multi_us_mesh_components" in self.mesh_options.keys():
                raise TypeError(
                    "'multi_us_mesh_components' is only for 'MultiUSMesh' mesh_type . Please don't provide any multi_us_mesh_components dictionary for 'USMesh' mesh_type."
                )

        elif mesh_type == "MultiUSMesh":
            if "multi_us_mesh_components" not in self.mesh_options.keys():
                raise TypeError(
                    "The 'multi_us_mesh_components' keyword argument is *NOT* "
                    "optional. A multi_us_mesh_components dictionary must be passed for 'MultiUSMesh' mesh_type. "
                )

            self.multi_us_mesh = True
            # self.multi_us_mesh_instances = multi_us_mesh_instancesdict

        else:
            raise TypeError(
                "Available options for mesh_type: 'USMesh' and 'MultiUSMesh'. By default, 'USMesh' is used. Choose 'MultiUSMesh' When several mesh components are considered for independent deformation of multiple overset component meshes."
            )

        # defaults:

        # flag to determine if the mesh warping component is added
        # in the nonlinear solver loop (e.g. for aerostructural)
        # or as a preprocessing step like the surface mesh coordinates
        # (e.g. for aeropropulsive). This will avoid doing extra work
        # for mesh deformation when the volume mesh does not change
        # during nonlinear iterations
        self.warp_in_solver = False
        # flag for aerostructural coupling variables
        self.struct_coupling = False
        # flag to enable propulsion coupling variables
        self.prop_coupling = False
        # flag to enable heat transfer coupling variables
        # TODO AY-JA: Can you rename heat_transfer to thermal_coupling to be consistent with other flags?
        self.heat_transfer = False

        # depending on the scenario we are building for, we adjust a few internal parameters:
        if scenario.lower() == "aerodynamic":
            # default
            pass

        elif scenario.lower() == "aerostructural":
            # volume mesh warping needs to be inside the coupling loop for aerostructural
            self.warp_in_solver = True
            self.struct_coupling = True

        elif scenario.lower() == "aeropropulsive":
            self.prop_coupling = True

        elif scenario.lower() == "aerothermal":
            self.heat_transfer = True

        # flag to determine if we want to restart a failed solution from free stream
        self.restart_failed_analysis = restart_failed_analysis

        # flag for raising an error on convergence stall
        self.err_on_convergence_fail = err_on_convergence_fail

        # balance group for propulsion
        self.balance_group = balance_group

        # TODO AY-JA: Check if we still need this
        # if self.heat_transfer:
        #     self.promotes('heat_xfer', inputs=[name])

    # api level method for all builders
    def initialize(self, comm):
        self.solver = ADFLOW(options=self.options, comm=comm)

        # We need to set the family groups here before we initialize the mesh.
        # The user can't do this in the setup method because we set the mesh
        # within the initialize method.  Allowing the user to set the
        # 'user_family_group' variable gives them a way to specificy custom
        # family groups within the setup when they create the builder.
        if self.user_family_groups is not None:
            # Loop over the surface family groups.  The keys are the
            # new group names and the vals are the subsurface to add
            # to the group.
            for key, val in self.user_family_groups.items():
                self.solver.addFamilyGroup(key, val)

        if self.multi_us_mesh:
            mesh = MultiUSMesh(self.mesh_options["gridFile"], self.mesh_options["multi_us_mesh_components"], comm=comm)
        else:
            mesh = USMesh(options=self.mesh_options, comm=comm)

        self.solver.setMesh(mesh)

    def get_solver(self):
        # this method is only used by the RLT transfer scheme
        return self.solver

    # api level method for all builders
    def get_coupling_group_subsystem(self, scenario_name=None):
        adflow_group = ADflowGroup(
            solver=self.solver,
            use_warper=self.warp_in_solver,
            struct_coupling=self.struct_coupling,
            prop_coupling=self.prop_coupling,
            restart_failed_analysis=self.restart_failed_analysis,
            err_on_convergence_fail=self.err_on_convergence_fail,
            balance_group=self.balance_group,
        )
        return adflow_group

    def get_mesh_coordinate_subsystem(self, scenario_name=None):
        # TODO modify this so that we can move the volume mesh warping to the top level
        # we need this to do mesh warping only once for all serial points.
        # for now, volume warping is duplicated on all scenarios, which is not efficient

        # use_warper = not self.warp_in_solver
        # # if we do warper in the mesh element, we will do a group thing
        # if use_warper:
        #     return ADflowMeshGroup(aero_solver=self.solver)
        # else:

        # just return the component that outputs the surface mesh.
        return ADflowMesh(aero_solver=self.solver)

    def get_pre_coupling_subsystem(self, scenario_name=None):
        if self.warp_in_solver:
            # if we warp in the solver, then we wont have any pre-coupling systems
            return None
        else:
            # we warp as a pre-processing step
            return ADflowWarper(aero_solver=self.solver)

    def get_post_coupling_subsystem(self, scenario_name=None):
        return ADflowFunctions(aero_solver=self.solver)

    # TODO the get_nnodes is deprecated. will remove
    def get_nnodes(self, groupName=None):
        return int(self.solver.getSurfaceCoordinates(groupName=groupName).size / 3)

    def get_number_of_nodes(self, groupName=None):
        return int(self.solver.getSurfaceCoordinates(groupName=groupName).size / 3)
