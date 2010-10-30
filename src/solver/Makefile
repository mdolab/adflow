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
	       computeInletMassFlowFullWheel.o\
	       computeLamViscosity.o\
	       computeUtau.o\
	       convergenceHeader.o\
	       convergenceInfo.o\
	       derivativeRotMatrixRigid.o\
	       eulerWallsPresent.o\
	       executeMGCycle.o\
	       extrapolate2ndHalo.o\
	       forcesAndMoments.o\
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
	       timeStep.o\
	       transferToCoarseGrid.o\
	       transferToFineGrid.o\
	       unsteadyHeader.o\
	       utauWF.o\
	       viscousFlux.o\
	       writeFamilyMassflow.o

CC_OBJECTS   = connect_signals.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS) $(CC_OBJECTS)
