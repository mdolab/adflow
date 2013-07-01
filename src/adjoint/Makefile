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
		adjointPythonInterfaceFunctions.F90\
		computeAeroCoef.f90        \
		computeObjPartials.F90\
		copyADjointStencil.f90     \
		createPETScVars.F90        \
		destroyPETScVars.F90       \
		preallocation.F90          \
		finalizePETSc.F90          \
		getSolution.f90            \
		initializePETSc.F90        \
		preprocessingADjoint.F90   \
		releaseMemADjoint.f90      \
		setupAllResidualMatrices.F90\
		setupRestrictionMatrices.F90\
		setupProlongationMatrix.F90\
		setupPETScKsp.F90         \
		solverADjoint.F90         \
		updateFlow.f90            \
		updateReferencePoint.f90   \
		updateRotationRate.f90     \
		setupCouplingMatrixStruct.F90\
		matrixInfo.F90\
		adjoint_stencils.f90\

FILES_TO_COPY = \
		solveADjointTransposePETSc.F90\
		initPETScWrap.c\
		saveAdjoint.F90\

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
		if [ $$file -nt ../../src_cs/adjoint/$$file -o ! -e ../../src_cs/adjoint/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/adjoint/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/adjoint/$$file -o ! -e ../../src_cs/adjoint/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/adjoint/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/adjoint


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(OFILES)
