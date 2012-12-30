#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
#      * Starting date: 01-14-2008                                      *
#      * Last modified: 01-14-2008                                      *
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

SUMB_DIR = ../../..
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
		adjustInflowAngleAdjTS.f90   \
		applyAllBCAdjTS.f90          \
		bcSymmAdjTS.f90              \
		bcEulerWallAdjTS.f90         \
		bcFarfieldAdjTS.f90          \
		checkInputParamAdjTS.f90     \
		checkOverlapAdjTS.f90        \
		checkXOverlapAdjTS.f90       \
		computeGammaAdjTS.f90        \
		computePressureAdjTS.f90     \
		computeRAdjTS.f90            \
		computeVolAdjTS.f90        \
		derivativeRigidRotAngleAdjTS.f90\
		derivativeRotMatrixRigidAdjTS.f90\
		etotArrayAdjTS.f90           \
		extractBCStatesAdjTS.f90     \
		extrapolate2ndHaloAdjTS.f90  \
		getDirAngleTS.f90            \
		getDirVectorTS.f90           \
		gridVelocitiesFineLevelAdjTS.f90 \
		initresAdjTS.f90             \
		inviscidCentralFluxAdjTS.f90 \
		inviscidDissFluxScalarAdjTS.f90 \
		inviscidUpwindFluxAdjTS.f90  \
		leftRightStateAdjTS.f90      \
		metricAdjTS.f90              \
		normalVelocitiesAdjTS.f90    \
		referenceStateAdjTS.f90      \
		replaceBCStatesAdjTS.f90     \
		residualAdjTS.f90            \
		riemannFluxAdjTS.f90         \
		rigidRotAngleAdjTS.f90       \
		rotMatrixRigidBodyAdjTS.f90     \
		setFlowInfinityStateAdjTS.f90\
		timeStepAdjTS.f90            \
		vectorRotationTS.f90         \
		xhaloAdjTS.f90

FILES_TO_COPY = 

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
	@echo "        Making clean in adjoint... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../../src_cs/adjoint/residualInputTS/$$file -o ! -e ../../../src_cs/adjoint/residualInputTS/$$file ]; then \
			python ../../complexify.py $$file; \
			mv c_$$file ../../../src_cs/adjoint/residualInputTS/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../../src_cs/adjoint/residualInputTS/$$file -o ! -e ../../../src_cs/adjoint/residualInputTS/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../../src_cs/adjoint/residualInputTS/$$file; \
		fi;\
		done

	cp Makefile ../../../src_cs/adjoint/residualInputTS


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(OFILES)
