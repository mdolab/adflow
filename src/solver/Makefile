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

FILES_TO_COMPLEXIFY1 = \
	       BCRoutines.F90\

FILES_TO_COMPLEXIFY2 = \
	       applyAllBC.F90\
	       adjustInflowAngle.f90\
	       bleedFlowParameters.f90\
	       convergenceInfo.F90\
	       computedwD3ADI.f90\
	       computeInletMassFlowFullWheel.f90\
	       computeLamViscosity.F90\
	       computeUtau.f90\
               coorUtilsALE.F90\
	       DADISmoother.f90\
	       derivativeRotMatrixRigid.f90\
	       eulerWallsPresent.f90\
	       executeMGCycle.f90\
               farFieldDrag.f90\
	       farFieldInducedDrag.f90\
	       forcesAndMoments.F90\
	       farFieldDrag.f90\
               geomUtilsALE.F90\
	       getDirAngle.f90\
	       getDirVector.f90\
	       gridVelocitiesCoarseLevels.f90\
	       gridVelocitiesFineLevel.F90\
	       gridVelocitiesFineLevel_ALE.F90\
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
	       normalVelocities.F90\
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
	       slipVelocities_ALE.f90\
	       solveState.F90\
	       solver.F90\
	       solverSteady.f90\
               solverUnsteadyExplicitRK.F90\
               solverUnsteadyImplicitRK.F90\
	       solverUnsteadyBDF.F90\
	       solverUnsteady_ALE.F90\
	       solverUnsteadyWrap.F90\
	       timeStep.F90\
	       transferToCoarseGrid.f90\
	       transferToFineGrid.f90\
	       unsteadyHeader.f90\
	       utauWF.F90\
	       vectorRotation.f90\
	       viscousFlux.F90\
	       writeFamilyMassflow.f90

FILES_TO_COPY = connect_signals.c\
		convergenceHeader.f90\

ALL_FILES1 = $(FILES_TO_COMPLEXIFY1)
ALL_FILES2 = $(FILES_TO_COMPLEXIFY2) $(FILES_TO_COPY)
OFILES1 := $(ALL_FILES1:%.F90=%.o) $(ALL_FILES1:%.F90=%.o) $(ALL_FILES1:%.c=%.o)
OFILES2 := $(ALL_FILES2:%.f90=%.o) $(ALL_FILES2:%.F90=%.o) $(ALL_FILES2:%.c=%.o)
ALL_COMPLEXIFY = $(FILES_TO_COMPLEXIFY1) $(FILES_TO_COMPLEXIFY2)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src1
	cp -f *.mod $(SUMB_MODDIR)
	$(PMAKE) src2
clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(ALL_COMPLEXIFY); do \
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

src1:	$(OFILES1)
src2:	$(OFILES2)
