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
	       setGlobalCellsAndNodes.o \
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

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
