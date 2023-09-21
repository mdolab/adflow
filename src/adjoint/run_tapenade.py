#!/usr/bin/python
"""
Differentiate Fortran source code.

Differentiates fortran files with tapenade. It compares the timestamp of the
source file with the differentiated file and only selects the files that appear
to have changed.

There are 3 different modes: forward, reverse and reverse_fast. The routines
defined in `full_routines` are used for the forward and reverse mode. The
routines defined in `state_only_routines` are only used for the reverse_fast
mode.

Each mode is treated independendly as follows:
    1. Create temporary directories
    2. Add all files to the list of `to be differentiate files`
    3. loop through that list
        A. If the source file is older than the differentiated file, remove it
        from that list.
    4. preprocess all `to be differentiated files` witch cpp
    5. run tapenade on all `to be differentiated files`. This genereates
    differentiated fortran code
    6. postprocess all differentiated file with the scripts in the 'autoEdit'
    folder. This also copies the finished files to the correct folder
    7. Delete temporary directories

It is possible to force differentiation of all files, a single file or one mode
only. Please look at the options using `pyhton run_tapenade.py --help`

"""


import argparse
import os
import shutil
import os.path as osp
import subprocess

SRC_PATH = ".."
TAPENADE_PRECISION = "-i4 -dr8 -r8"


full_routines = {
    osp.join(SRC_PATH, "adjoint", "adjointExtra.F90"): [
        "adjointExtra%xhalo_block(x) >  (x)",
        "adjointExtra%volume_block(x) > (x, vol)",
        "adjointExtra%metric_block(x) > (x, si, sj, sk)",
        "adjointExtra%boundaryNormals(si, sj, sk) > (si, sj, sk, bcdata%norm)",
        "adjointExtra%sumDwandFw(fw, dw) > (dw)",
        "adjointExtra%resScale(dw) > (dw)",
    ],
    osp.join(SRC_PATH, "utils", "flowUtils.F90"): [
        "flowUtils%adjustInflowAngle(alpha, beta) > (veldirfreestream, dragdirection, liftdirection)",
        "flowUtils%computePressureSimple(w, pInfCorr) > (w, pInfCorr, p)",
        "flowUtils%computeLamViscosity(w, p, Tref, muRef, rGas) > (w, p, Tref, muRef, rGas, rlv)",
        "flowUtils%computeEtotBlock(w, p) > (w, p)",
        "flowUtils%computeSpeedOfSoundSquared(w, p) > (w, p, aa)",
        "flowUtils%allNodalGradients(w, vol, si, sj, sk, aa, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz) > (w, vol, si, sj, sk, aa, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "initFlow", "initializeFlow.F90"): [
        "initializeFlow%referenceState(mach, machcoef, pInfDim, TinfDim, rhoInfDim, velDirFreestream) > (mach, machcoef, pInfDim, TinfDim, rhoInfDim, velDirFreeStream, pRef, rhoRef, Tref, muRef, timeRef, pInf, pinfCorr, rhoInf, uInf, rGas, muInf,  winf) "  # noqa: E501
    ],
    osp.join(SRC_PATH, "bcdata", "BCData.F90"): [
        "bcData%BCDataIsoThermalwall(bcVarArray, Tref) > (BCData%TNS_wall, Tref)",
        "bcData%BCDataSubsonicInflow(bcVarArray, muRef, Pref, rhoRef, Tref, Href, wInf) > (bcData%PtInlet, bcData%ttInlet, bcData%htInlet, bcData%turbInlet, muRef, Pref, rhoRef, Tref, Href, wInf)",  # noqa: E501
        "bcData%BCDataSubsonicOutflow(bcVarArray, Pref) > (bcData%Ps,  Pref)",
        "bcData%BCDataSupersonicInflow(bcVarArray, muRef, rhoRef, Pref, uRef, wInf, pInfCorr) > (bcData%rho, bcData%velx, bcData%vely, bcData%velz, bcData%ps, bcData%turbInlet, muRef, rhoRef, Pref, uRef, wInf, pInfCorr)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "wallDistance", "wallDistance.F90"): [
        "wallDistance%updateWallDistancesQuickly(x, xSurf, d2wall) > (x, xSurf, d2wall)",
    ],
    osp.join(SRC_PATH, "turbulence", "turbUtils.F90"): [
        "turbUtils%computeEddyViscosity(w,rlv)>(w,rlv,rev)",
        "turbutils%turbAdvection(w, vol, si, sj, sk, scratch, sfacei, sfacej, sfacek) > (w, vol, si, sj, sk, scratch, sfacei, sfacej, sfacek)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "solver", "BCRoutines.F90"): [
        "BCRoutines%bcSymm1stHalo(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2, bcData%norm) > (ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2, bcData%norm)",  # noqa: E501
        "BCRoutines%bcSymm2ndHalo(ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3, bcData%norm) > (ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3, bcData%norm)",  # noqa: E501
        "BCRoutines%bcSymmPolar1stHalo(xx, ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2) > (xx, ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)",  # noqa: E501
        "BCRoutines%bcSymmPolar2ndHalo(xx, ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3) > (xx, ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3)",  # noqa: E501
        "BCRoutines%bcNSWallAdiabatic(ww0, ww1, ww2, pp0, pp1, pp2, pp3, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%uSlip) > (ww0, ww1, ww2, pp0, pp1, pp2, pp3, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%uSlip)",  # noqa: E501
        "BCRoutines%bcNSWallIsothermal(ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, rGas, BCData%TNS_Wall, bcData%uSlip) > (ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, rGas, BCData%TNS_Wall, bcData%uSlip)",  # noqa: E501
        "BCRoutines%bcFarField(ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, pInfCorr, wInf, bcData%norm, bcData%rface) > (ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, pInfCorr, wInf, bcData%norm, bcData%rface)",  # noqa: E501
        "BCRoutines%bcEulerWall(ww0, ww1, ww2, pp0, pp1, pp2, pp3, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, ssi, ssj, ssk, ss) > (ww0, ww1, ww2, pp0, pp1, pp2, pp3, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, ssi, ssj, ssk, ss)",  # noqa: E501
        "BCRoutines%bcSubsonicInflow(ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, bcData%PtInlet, bcData%ttInlet, bcData%htInlet) > (ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, bcData%PtInlet, bcData%ttInlet, bcData%htInlet)",  # noqa: E501
        "BCRoutines%bcSubsonicOutflow(ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, bcData%Ps) > (ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, rev0, rev1, rev2, bcData%norm, bcData%Ps)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "turbulence", "turbBCRoutines.F90"): [
        "turbBCRoutines%applyAllTurbBCThisBlock(rev, w, bvtj1, bvtj2, bvtk1, bvtk2, bvti1, bvti2) > (rev, w)",
        "turbBCRoutines%bcTurbTreatment(w, rlv, d2wall, winf) > (w, rlv, d2wall, winf, bvtj1, bvtj2, bvtk1, bvtk2, bvti1, bvti2)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "solver", "solverUtils.F90"): [
        "solverUtils%timeStep_block(w, pInfCorr, rhoInf, si, sj, sk, sFaceI, sFaceJ, sFaceK, p, radi, radj, radk, dtl, rlv, rev, vol) > (w, pInfCorr, rhoInf, si, sj, sk, sFaceI, sFaceJ, sFaceK, p, radi, radj, radk, dtl, rlv, rev, vol)",  # noqa: E501
        "solverutils%slipVelocitiesFineLevel_block(pinf, timeref, rhoinf, veldirfreestream, machgrid, bcdata%uslip, flowDoms%x) > (veldirfreestream, bcdata%uslip, flowDoms%x)",  # noqa: E501
        "solverutils%normalvelocities_block(sfacei, sfacej, sfacek, si, sj, sk, bcdata%rface) > (sfacei, sfacej, sfacek, si, sj, sk, bcdata%rface)",  # noqa: E501
        "solverutils%gridVelocitiesFineLevel_block(flowDoms%x, si, sj, sk, s, sfacei, sfacej, sfacek, velDirFreestream, pinf, timeref, rhoinf, machgrid) > (flowDoms%x, si, sj, sk, s, sfacei, sfacej, sfacek, velDirFreestream)",  # noqa: E501
        "solverutils%cellFaceVelocities(xc, rotCenter, rotRate, velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc) > (xc, velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "solver", "residuals.F90"): [
        "residuals%sourceTerms_block(w, pref, uref, plocal, dw, actuatorRegions%force, actuatorRegions%heat) > (w, pref, uref, plocal, dw, actuatorRegions%force, actuatorRegions%heat)",  # noqa: E501
        "residuals%initres_block(dw, fw, flowDoms%w, flowDoms%vol, dscalar, dvector) > (dw, fw, flowDoms%w, flowDoms%vol, dscalar, dvector)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "solver", "fluxes.F90"): [
        "fluxes%inviscidCentralFlux(w, timeRef, vol, si, sj, sk, p, dw, sfacei, sfacej, sfacek) > (w, timeRef, vol, si, sj, sk, p, dw, sfacei, sfacej, sfacek)",  # noqa: E501
        "fluxes%inviscidDissFluxScalar(w, pInfCorr, rhoInf, p, radi, radj, radk, fw) > (w, pInfCOrr, rhoInf, p, radi, radj, radk, fw)",  # noqa: E501
        "fluxes%inviscidDissFluxMatrix(w, pInfCorr, si, sj, sk, p, fw, sfacei, sfacej, sfacek) > (w, pInfCorr, si, sj, sk, p, fw, sfacei, sfacej, sfacek)",  # noqa: E501
        "fluxes%inviscidDissFluxScalarApprox(w, p, radi, radj, radk, fw, rhoinf, pinfcorr) > (w, p, fw)",
        "fluxes%inviscidDissFluxMatrixApprox(w, p, fw, pinfcorr, sfacei, sfacej, sfacek, si, sj, sk) > (w, p, fw, pinfcorr, sfacei, sfacej, sfacek, si, sj, sk)",  # noqa: E501
        "fluxes%inviscidUpwindFlux(w, si, sj, sk, p, fw, sfacei, sfacej, sfacek) > (w, si, sj, sk, p, fw, sfacei, sfacej, sfacek)",  # noqa: E501
        "fluxes%viscousFlux(w, x, si, sj, sk, rlv, rev, aa, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz, fw, viscsubface%tau, viscSubface%q) > (w, x, si, sj, sk, fw, viscsubface%tau, viscSubface%q)",  # noqa: E501
        "fluxes%viscousFluxApprox(w, x, rlv, rev, aa, fw, si, sj, sk) > (w, x, fw)",
    ],
    osp.join(SRC_PATH, "solver", "surfaceIntegrations.F90"): [
        "surfaceIntegrations%wallIntegrationFace(Pinf, Pref, ww2, pp1, pp2, xx, ssi, velDirFreeStream, pointRef, machCoef, viscSubface%tau, localValues, BCData%Fp, bcData%Fv, BCData%area, localValues) > (Pinf, Pref, ww2, pp1, pp2, xx, ssi, velDirFreeStream, pointRef, machCoef, viscSubface%tau, localValues, BCdata%Fp, bcData%Fv, BCData%area, localValues)",  # noqa: E501
        "surfaceIntegrations%flowIntegrationFace(Pref, rhoRef, timeRef, TRef, rGas, ww1, ww2, pp1, pp2, xx, ssi, pointRef, localValues) > (Pref, rhoRef, timeRef, TRef, rGas, ww1, ww2, pp1, pp2, xx, ssi, pointRef, localValues)",  # noqa: E501
        "surfaceIntegrations%getCostFunctions(globalVals, liftDirection, dragDirection, pRef, rhoRef, machCoef) > (liftDirection, dragDirection, pRef, rhoRef, machCoef, funcValues)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "solver", "zipperIntegrations.F90"): [
        "zipperIntegrations%wallIntegrationZipper(vars, pointRef, localValues) > (vars, pointRef, localValues)",
        "zipperIntegrations%flowIntegrationZipper(vars, Pref, rhoRef, timeRef, TRef, rGas, pointRef, localValues) > (vars, Pref, rhoRef, timeRef, TRef, rGas, pointRef, localValues)",  # noqa: E501
    ],
    osp.join(SRC_PATH, "turbulence", "sa.F90"): [
        "sa%saSource(w, rlv, vol, si, sj, sk, timeRef, d2wall) > (w, rlv, vol, si, sj, sk, timeREf, scratch)",
        "sa%saViscous(w, vol, si, sj, sk, rlv, scratch) > (w, vol, si, sj, sk, rlv, scratch)",
        "sa%saResScale(scratch, dw) > (dw)",
    ],
    osp.join(SRC_PATH, "overset", "oversetUtilities.F90"): [
        "oversetUtilities%newtonUpdate(xCen, blk, frac) > (blk, frac)",
        "oversetUtilities%fracToWeights(frac) > (weights)",
    ],
}

state_only_routines = {
    osp.join(SRC_PATH, "utils", "flowUtils.F90"): [
        "flowutils%allNodalGradients(aa, w) > (aa, w, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)",
    ],
    osp.join(SRC_PATH, "turbulence", "turbUtils.F90"): [
        "turbutils%saEddyViscosity(w,rlv) > (w,rlv,rev)",
        "turbutils%turbAdvection(w,scratch)>(w,scratch)",
    ],
    osp.join(SRC_PATH, "solver", "solverUtils.F90"): [
        "solverutils%timeStep_block(w, p) > (w, p, radi, radj, radk)",
    ],
    osp.join(SRC_PATH, "solver", "residuals.F90"): [
        "residuals%sourceTerms_block(w, dw) > (w, dw)",
        "residuals%initres_block(dw, fw, flowDoms%w) > (dw, fw, flowDoms%w)",
    ],
    osp.join(SRC_PATH, "solver", "fluxes.F90"): [
        "fluxes%inviscidCentralFlux(w, p, dw) > (w, p, dw)",
        "fluxes%inviscidDissFluxScalar(p, w, radi, radj, radk, fw) > (p, w, fw)",
        "fluxes%inviscidDissFluxMatrix(p, w, fw) > (p, w, fw)",
        "fluxes%inviscidUpwindFlux(w, p, fw) > (w, p, fw)",
        "fluxes%viscousFlux(w, aa, rlv, rev, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz, fw) > (w, aa, fw)",
    ],
    osp.join(SRC_PATH, "turbulence", "sa.F90"): [
        "sa%saSource(w, rlv) > (w, rlv, scratch)",
        "sa%saViscous(w, rlv, scratch) > (w, rlv, scratch)",
        "sa%saResScale(scratch) > (dw)",
    ],
}

dependencies = [
    osp.join(SRC_PATH, "modules", "block.F90"),
    osp.join(SRC_PATH, "modules", "surfaceFamilies.F90"),
    osp.join(SRC_PATH, "modules", "monitor.f90"),
    osp.join(SRC_PATH, "modules", "diffSizes.f90"),
    osp.join(SRC_PATH, "modules", "inputParam.F90"),
    osp.join(SRC_PATH, "modules", "constants.F90"),
    osp.join(SRC_PATH, "modules", "precision_tapenade.f90"),
    osp.join(SRC_PATH, "modules", "iteration.f90"),
    osp.join(SRC_PATH, "modules", "section.f90"),
    osp.join(SRC_PATH, "modules", "communication.F90"),
    osp.join(SRC_PATH, "modules", "paramTurb.F90"),
    osp.join(SRC_PATH, "modules", "cgnsGrid.F90"),
    osp.join(SRC_PATH, "modules", "CpCurveFits.f90"),
    osp.join(SRC_PATH, "modules", "blockPointers.F90"),
    osp.join(SRC_PATH, "modules", "BCPointers.F90"),
    osp.join(SRC_PATH, "modules", "flowVarRefState.F90"),
    osp.join(SRC_PATH, "modules", "wallDistanceData.F90"),
    osp.join(SRC_PATH, "modules", "actuatorRegionData.F90"),
    osp.join(SRC_PATH, "modules", "overset.F90"),
    osp.join(SRC_PATH, "modules", "cgnsNames.f90"),
    osp.join(SRC_PATH, "modules", "su_cgns.F90"),
    osp.join(SRC_PATH, "modules", "BCDataMod.F90"),
    osp.join(SRC_PATH, "solver", "ALEUtils.F90"),
    osp.join(SRC_PATH, "turbulence", "turbMod.F90"),
    osp.join(SRC_PATH, "utils", "utils.F90"),
    osp.join(SRC_PATH, "utils", "sorting.F90"),
]


class Main:
    """Main logic of differentiation script."""

    dirs = {
        "forward": "outputForward",
        "reverse": "outputReverse",
        "reverse_fast": "outputReverseFast",
    }

    tmp_dirs = {
        "forward": "tmp_forward",
        "reverse": "tmp_reverse",
        "reverse_fast": "tmp_reverse_fast",
    }

    preprocess_dirs = {
        "forward": "preprocess_forward",
        "reverse": "preprocess_reverse",
        "reverse_fast": "preprocess_reverse_fast",
    }

    identifiers = {
        "forward": "_d",
        "reverse": "_b",
        "reverse_fast": "_fast_b",
    }

    def __init__(self):
        """Runs Main logic of this script."""
        self.modes = ["forward", "reverse", "reverse_fast"]

        self.routines = {"forward": {}, "reverse": {}, "reverse_fast": {}}

        self.tmpdirs = []
        for d in self.tmp_dirs.values():
            self.tmpdirs.append(d)
        for d in self.preprocess_dirs.values():
            self.tmpdirs.append(d)

        # find files to process
        self.parse_args()
        self.filter_files()

        # if one mode is forced, set the other modes to zero
        if self.args.only_forward:
            self.routines["reverse"] = dict()
            self.routines["reverse_fast"] = dict()
        if self.args.only_reverse:
            self.routines["forward"] = dict()
            self.routines["reverse_fast"] = dict()
        if self.args.only_reverse_fast:
            self.routines["forward"] = dict()
            self.routines["reverse"] = dict()

        # report the files that are differentiated
        for mode in self.modes:
            print(f"Differentiating the following files in '{mode}' mode:")
            for routine in self.routines[mode].keys():
                print(routine)
            print(" ")

        self.create_tmp_dirs()

        # actually process files
        for mode in self.modes:
            if len(self.routines[mode]) == 0:
                continue

            self.preprocess_routines(mode)
            self.differentiate_routines(mode)

        self.postprocess_files()

        if not self.args.no_cleanup:
            self.delete_tmp_dirs()

    def parse_args(self):
        """Parses the command line arguments."""
        parser = argparse.ArgumentParser()
        parser.add_argument("-force_all", action="store_true", help="Forcefully differentiate all files when set.")
        parser.add_argument("-force_file", type=str, default=None, help="Forcefully differentaited specified file.")
        parser.add_argument("-no_cleanup", action="store_true", help="Do not remove temporary directories")

        for mode in self.modes:
            parser.add_argument(f"-only_{mode}", action="store_true", help=f"only run '{mode}' routines.")

        self.args = parser.parse_args()

    def filter_files(self):
        """Decides what files should be differentiated."""
        self.routines = {"forward": full_routines, "reverse": full_routines, "reverse_fast": state_only_routines}

        # if all files should be differentiated, return here
        if self.args.force_all:
            return

        # if only a specific file is specified
        if self.args.force_file is not None:
            self.routines = {"forward": {}, "reverse": {}, "reverse_fast": {}}
            f = self.args.force_file
            if f in list(full_routines.keys()):
                self.routines["forward"] = {f: full_routines[f]}
                self.routines["reverse"] = {f: full_routines[f]}
            if f in list(state_only_routines.keys()):
                self.routines["reverse_fast"] = {f: state_only_routines[f]}
            return

        # actually filter files
        self.routines["forward"] = self.drop_older_files("forward")
        self.routines["reverse"] = self.drop_older_files("reverse")
        self.routines["reverse_fast"] = self.drop_older_files("reverse_fast")

    def delete_tmp_dirs(self):
        """Deletes temporary directories."""
        for directory in self.tmpdirs:
            try:
                shutil.rmtree(directory)
            except FileNotFoundError:
                pass

    def create_tmp_dirs(self):
        """Creates temporary directories."""
        self.delete_tmp_dirs()
        for directory in self.tmpdirs:
            os.mkdir(directory)

    def postprocess_files(self):
        """Runs the post-processing scripts in 'autoEdit' on the differentiated files."""
        subprocess.run(
            ["python", osp.join("autoEdit", "autoEditForward.py"), self.tmp_dirs["forward"], self.dirs["forward"]],
            shell=False,
        )

        subprocess.run(
            ["python", osp.join("autoEdit", "autoEditReverse.py"), self.tmp_dirs["reverse"], self.dirs["reverse"]],
            shell=False,
        )

        subprocess.run(
            [
                "python",
                osp.join("autoEdit", "autoEditReverseFast.py"),
                self.tmp_dirs["reverse_fast"],
                self.dirs["reverse_fast"],
            ],
            shell=False,
        )

    def drop_older_files(self, mode):
        """
        Compares the timestamp of the differentiated file with its source. If
        the source is older, this file is dropped from the 'to be differentiated
        files' list

        Parameters
        ----------
        mode : str
            The mode to use, either 'forward', 'reverse' or 'reverse_fast'

        Returns
        -------
        return_routines : dict
            Holds the `files to be differentiated` in the same format as
            `state_only_routines` or `full_routines`
        """
        return_routines = dict()

        for source in self.routines[mode].keys():
            destination = osp.join(
                self.dirs[mode], osp.basename(osp.splitext(source)[0]) + self.identifiers[mode] + ".f90"
            )

            # get timestamps of files
            source_timestamp = osp.getmtime(source)
            destination_timestamp = osp.getmtime(destination)

            # only keep file if source is newer
            if source_timestamp > destination_timestamp:
                return_routines[source] = self.routines[mode][source]

        return return_routines

    def preprocess_routines(self, mode):
        """
        Runs the preprocessor on the `files to be differentiated` per mode.

        Parameters
        ----------
        mode : str
            The mode to use, either 'forward', 'reverse' or 'reverse_fast'
        """
        if mode == "forward":
            tapenade_strig = "-DTAPENADE_FORWARD"
        else:
            tapenade_strig = "-DTAPENADE_REVERSE"

        def run_preprocessor(file_name):
            subprocess.run(
                [
                    "cpp",
                    "-DUSE_TAPENADE",
                    tapenade_strig,
                    "-traditional",
                    "-P",
                    file_name,
                    osp.join(self.preprocess_dirs[mode], osp.basename(file_name)),
                ],
                shell=False,
            )

        # preprocess dependencies
        for dependency in dependencies:
            run_preprocessor(dependency)

        # preprocess routines
        for file_name in self.routines[mode].keys():
            run_preprocessor(file_name)

    def differentiate_routines(self, mode):
        """
        Runs tapenade to differentiate the preprocessed files.

        Parameters
        ----------
        mode : str
            The mode to use, either 'forward', 'reverse' or 'reverse_fast'
        """
        # assemble the list of files and the list of routines
        tapenade_files = ""
        tapenade_routines = ""
        # first, add dependencies
        for dependency in dependencies:
            tapenade_files += osp.join(self.preprocess_dirs[mode], osp.basename(dependency)) + " "

        for file_name, file_routines in self.routines[mode].items():
            tapenade_files += osp.join(self.preprocess_dirs[mode], osp.basename(file_name)) + " "

            for file_routine in file_routines:
                tapenade_routines += file_routine + " "

        # figure out the correct modes for tapenade
        if mode == "forward":
            tapenade_mode = "-forward"
            tapenade_isize = ""
            tapenade_targets = [
                "-tgtvarname",
                "%d",
                "-tgtfuncname",
                f"%{self.identifiers[mode]}",
                "-tgtmodulename",
                f"%{self.identifiers[mode]}",
            ]
        else:
            tapenade_mode = "-reverse"
            tapenade_isize = "-noisize"
            tapenade_targets = [
                "-adjvarname",
                "%d",
                "-adjfuncname",
                f"%{self.identifiers[mode]}",
                "-adjmodulename",
                f"%{self.identifiers[mode]}",
            ]

        tapenade_msglevel = [
            "-msglevel",
            "30",
        ]
        if mode != "reverse_fast":
            tapenade_msglevel = []

        subprocess.run(
            [
                osp.join(os.environ["TAPENADE_HOME"], "bin", "tapenade"),
                "-html",
                "-head",
                tapenade_routines,
                tapenade_mode,
                *tapenade_targets,
                *tapenade_msglevel,
                tapenade_isize,
                *TAPENADE_PRECISION.split(),
                "-outputdirectory",
                self.tmp_dirs[mode],
                *tapenade_files.split(),
            ],
            shell=False,
        )


if __name__ == "__main__":
    Main()
