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
FF90_OBJECTS = applyAllBC.o\
	       bcBleedInflow.o\
	       bcDomainInterface.o\
	       bcEulerWall.o\
	       bcExtrap.o\
	       bcFarfield.o\
	       bcMdot.o\
	       bcNSWallAdiabatic.o\
	       bcNSWallIsothermal.o\
	       bcSubsonicInflow.o\
	       bcSubsonicOutflow.o\
	       bcSupersonicInflow.o\
	       bcSymm.o\
	       bcSymmPolar.o\
	       bcThrust.o\
	       bleedFlowParameters.o\
	       calMassFlux.o\
	       coeffConvergenceCheck.o\
	       computeInletMassFlowFullWheel.o\
	       computeLamViscosity.o\
	       computeUtau.o\
	       convergenceHeader.o\
	       convergenceInfo.o\
	       componentsBreakDownPrintout.o\
	       derivativeRotMatrixRigid.o\
	       eulerWallsPresent.o\
	       executeMGCycle.o\
	       extrapolate2ndHalo.o\
               farFieldDrag.o\
	       farFieldInducedDrag.o\
	       forcesAndMoments.o\
	       forcesAndMomentsCBD.o\
	       farFieldDrag.o\
	       gridVelocitiesCoarseLevels.o\
	       gridVelocitiesFineLevel.o\
	       initres.o\
	       inviscidCentralFlux.o\
	       inviscidDissFluxCusp.o\
	       inviscidDissFluxCuspCoarse.o\
	       inviscidDissFluxMatrix.o\
	       inviscidDissFluxMatrixCoarse.o\
	       inviscidDissFluxScalar.o\
	       inviscidDissFluxScalarCoarse.o\
	       inviscidUpwindFlux.o\
	       killFunctions.o\
	       level0CoolingModel.o\
	       maxEddyv.o\
	       maxHdiffMach.o\
	       normalVelocities.o\
	       residual.o\
	       residual2.o\
	       residualAveraging.o\
	       rungeKuttaSmoother.o\
	       setCoefTimeIntegrator.o\
	       setCornerRowHalos.o\
	       setCorrectionsCoarseHalos.o\
	       setCycleStrategy.o\
	       shiftSolution.o\
	       slipVelocities.o\
	       solveState.o\
	       solver.o\
	       solverSteady.o\
               solverUnsteadyExplicitRK.o\
               solverUnsteadyImplicitRK.o\
	       solverUnsteadyBDF.o\
	       solverUnsteadyMD.o\
	       timeStep.o\
	       transferToCoarseGrid.o\
	       transferToFineGrid.o\
	       unsteadyHeader.o\
	       utauWF.o\
	       viscousFlux.o\
	       writeFamilyMassflow.o

CC_OBJECTS   = connect_signals.o

FILES_TO_COMPLEXIFY = applyAllBC.f90\
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
	       computeInletMassFlowFullWheel.f90\
	       computeLamViscosity.f90\
	       computeUtau.f90\
	       convergenceInfo.f90\
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
	       gridVelocitiesCoarseLevels.f90\
	       gridVelocitiesFineLevel.f90\
	       initres.f90\
	       inviscidCentralFlux.f90\
	       inviscidDissFluxCusp.f90\
	       inviscidDissFluxCuspCoarse.f90\
	       inviscidDissFluxMatrix.f90\
	       inviscidDissFluxMatrixCoarse.f90\
	       inviscidDissFluxScalar.f90\
	       inviscidDissFluxScalarCoarse.f90\
	       inviscidUpwindFlux.f90\
	       killFunctions.F90\
	       level0CoolingModel.f90\
	       maxEddyv.f90\
	       maxHdiffMach.f90\
	       normalVelocities.f90\
	       residual.f90\
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
	       viscousFlux.f90\
	       writeFamilyMassflow.f90

FILES_TO_COPY = connect_signals.c\
		convergenceHeader.f90\


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

src:   $(FF90_OBJECTS) $(CC_OBJECTS)
