#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 03-10-2003                                      *
#      * Last modified: 11-21-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to create the object files of this       *
#      * directory.                                                     *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Include the settings for the entire code.                      *
#      *                                                                *
#      ******************************************************************

SUMB_DIR = ../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FILES_TO_COMPLEXIFY = \
	       adjustInflowAngle.f90\
	       applyAllBC.F90\
	       bcBleedInflow.f90\
	       bcDomainInterface.f90\
	       bcEulerWall.f90\
	       bcExtrap.f90\
	       bcFarfield.f90\
	       bcMdot.f90\
	       bcNSWallAdiabatic.f90\
	       bcNSWallIsothermal.f90\
	       bcSubsonicInflow.f90\
	       bcSubsonicOutflow.f90\
	       bcSupersonicInflow.f90\
	       bcSymm.f90\
	       bcSymmPolar.f90\
	       bcThrust.f90\
	       bleedFlowParameters.f90\
	       calMassFlux.f90\
	       coeffConvergenceCheck.f90\
	       convergenceInfo.F90\
	       computeInletMassFlowFullWheel.f90\
	       computeLamViscosity.f90\
	       computeUtau.f90\
	       componentsBreakDownPrintout.f90\
	       derivativeRotMatrixRigid.f90\
	       eulerWallsPresent.f90\
	       executeMGCycle.f90\
	       extrapolate2ndHalo.f90\
               farFieldDrag.f90\
	       farFieldInducedDrag.f90\
	       forcesAndMoments.f90\
	       forcesAndMomentsCBD.f90\
	       farFieldDrag.f90\
	       getDirAngle.f90\
	       getDirVector.f90\
	       gridVelocitiesCoarseLevels.f90\
	       gridVelocitiesFineLevel.f90\
	       initres.F90\
	       inviscidCentralFlux.F90\
	       inviscidDissFluxCusp.f90\
	       inviscidDissFluxCuspCoarse.f90\
	       inviscidDissFluxMatrix.F90\
	       inviscidDissFluxMatrixCoarse.f90\
	       inviscidDissFluxScalar.F90\
	       inviscidDissFluxScalarCoarse.f90\
	       inviscidUpwindFlux.F90\
	       killFunctions.F90\
	       level0CoolingModel.f90\
	       maxEddyv.f90\
	       maxHdiffMach.f90\
	       normalVelocities.f90\
	       residual.F90\
	       residual2.f90\
	       residualAveraging.f90\
	       rungeKuttaSmoother.f90\
	       setCoefTimeIntegrator.f90\
	       setCornerRowHalos.f90\
	       setCorrectionsCoarseHalos.f90\
	       setCycleStrategy.f90\
	       shiftSolution.f90\
	       slipVelocities.f90\
	       solveState.F90\
	       solver.F90\
	       solverSteady.f90\
               solverUnsteadyExplicitRK.F90\
               solverUnsteadyImplicitRK.F90\
	       solverUnsteadyBDF.F90\
	       solverUnsteadyMD.F90\
	       timeStep.f90\
	       transferToCoarseGrid.f90\
	       transferToFineGrid.f90\
	       unsteadyHeader.f90\
	       utauWF.f90\
	       vectorRotation.f90\
	       viscousFlux.f90\
	       writeFamilyMassflow.f90

FILES_TO_COPY = connect_signals.c\
		convergenceHeader.f90\

ALL_FILES = $(FILES_TO_COMPLEXIFY) $(FILES_TO_COPY)
OFILES := $(ALL_FILES:%.f90=%.o) $(ALL_FILES:%.F90=%.o) $(ALL_FILES:%.c=%.o)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/solver/$$file -o ! -e ../../src_cs/solver/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/solver/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/solver/$$file -o ! -e ../../src_cs/solver/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/solver/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/solver

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(OFILES)