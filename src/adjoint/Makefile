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

FF90_OBJECTS = computeAeroCoef.o        \
	       computeObjPartials.o\
	       copyADjointStencil.o     \
	       createPETScKsp.o         \
	       createPETScVars.o        \
	       createStatePETScVars.o   \
	       createSpatialPETScVars.o \
	       createCouplingPETScVars.o\
	       createExtraPETScVars.o   \
	       destroyPETScVars.o       \
	       preallocation.o          \
               finalizePETSc.o          \
	       getSolution.o            \
	       getADjoint.o             \
	       getdRdXvPsi.o            \
	       getdRdaPsi.o             \
	       getdRdwTVec.o            \
	       getdFdxVec.o             \
	       initializePETSc.o        \
	       preprocessingADjoint.o   \
	       releaseMemADjoint.o      \
	       setADjoint.o             \
	       setupAllResidualMatrices.o\
               setupPETScKsp.o         \
               solveADjointTransposePETSc.o\
               solverADjoint.o          \
	       updateFlow.o             \
	       updateReferencePoint.o   \
	       updateRotationRate.o     \
	       setupCouplingMatrixStruct.o\
	       agumentRHS.o\
	       matrixInfo.o\
               adjoint_stencils.o\
 #             verifydRdwFile.o\
 #             verifydRdxFile.o  \
 #             verifydRdExtraFile.o    
#	       verifyRAdj.o\
#	       verifyResiduals.o\
#	       verifydRdw.o             \
#              verifydRdwFile.o         \
#              verifydRdwFileFD.o       \
#	       verifydRdx.o             \
#	       verifydRdExtra.o         \
#	       verifydRdExtraFDAD.o      \
#	       verifydRdxFile.o             \
#              verifydRdxsFile.o             \
#	       verifydRdxFileFD.o             \

CC_OBJECTS   = initPETScWrap.o

FILES_TO_COMPLEXIFY = \
		computeAeroCoef.f90        \
		computeObjPartials.F90\
		copyADjointStencil.f90     \
		createPETScKsp.F90         \
		createPETScVars.F90        \
		destroyPETScVars.F90       \
		preallocation.F90          \
		finalizePETSc.F90          \
		getSolution.f90            \
		getADjoint.F90             \
		getdRdXvPsi.F90            \
		getdRdaPsi.F90             \
		getdRdwTVec.F90            \
		getdFdxVec.F90             \
		initializePETSc.F90        \
		preprocessingADjoint.F90   \
		releaseMemADjoint.f90      \
		setADjoint.F90             \
		setupAllResidualMatrices.F90\
		setupPETScKsp.F90         \
		solverADjoint.F90         \
		updateFlow.f90            \
		updateReferencePoint.f90   \
		updateRotationRate.f90     \
		setupCouplingMatrixStruct.F90\
		agumentRHS.F90\
		matrixInfo.F90\
		adjoint_stencils.f90\

FILES_TO_COPY = \
		solveADjointTransposePETSc.F90\
		initPETScWrap.c\
		createSpatialPETScVars.F90 \
		createStatePETScVars.F90   \
		createCouplingPETScVars.F90\
		createExtraPETScVars.F90   \

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

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
