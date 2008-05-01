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

FF90_MODULES = 

FF90_OBJECTS = applyAllBCAdj.o          \
	       bcSymmAdj.o              \
	       bcEulerWallAdj.o         \
	       bcFarfieldAdj.o          \
	       checkOverlapAdj.o        \
	       computePressureAdj.o     \
	       computeRAdj.o            \
	       etotArrayAdj.o           \
	       extractBCStatesAdj.o     \
	       extrapolate2ndHaloAdj.o  \
	       initresAdj.o             \
	       inviscidCentralFluxAdj.o \
	       inviscidUpwindFluxAdj.o  \
	       leftRightStateAdj.o      \
	       metricAdj.o              \
	       replaceBCStatesAdj.o     \
	       residualAdj.o            \
	       riemannFluxAdj.o         \

#FF90_OBJECTS = allocMemADjoint.o\
#	       assertionsADjoint.o\
#	       computeADjointGradientExtra.o\
#               computeADjointGradientSpatial.o\
#	       computeAeroCoef.o\
#	       convertDegToRad.o\
#	       convertRadToDeg.o\
#	       copyADjointStencil.o\
#	       copyLocalIndexHalo.o\
#	       createPETScKsp.o\
#	       createPETScMat.o\
#	       createPETScVars.o\
#	       createPETScVec.o\
#	       designExport.o\
#	       designImport.o\
#	       designInit.o\
#	       designWrite.o\
#	       destroyPETScVars.o\
#               destroyPETScKSP.o\
#	       finalizePETSc.o\
#	       getDirAngle.o\
#	       getDirVector.o\
#	       haloIndexCommunication.o\
#	       initializePETSc.o\
#	       preprocessingADjoint.o\
#	       getGlobalNodes.o\
#	       releaseMemADjoint.o\
#	       sendIndexHalo.o\
#	       setGlobalNodes.o\
#	       setPointersAdj.o\
#	       setupADjointMatrix.o\
#	       setupADjointRHSAeroCoeff.o\
#	       setupADjointRHS.o\
#	       setupGradientMatrixExtra.o\
#              setupGradientMatrixSpatial.o\
#	       setupGradientRHSExtra.o\
#               setupGradientRHSSpatial.o\
#	       solveADjointPETSc.o\
#	       solverADjoint.o\
#	       storeReceivedIndexHalo.o\
#	       vectorRotation.o\
#	       verifydRdW.o\
#	       verifydRdx.o\
#	       verifyRAdj.o\
#	       writeADjoint.o\
#	       writeADjointGradient.o\
#	       writeADjointHistory.o\
#	       writeADjointSol.o

CC_OBJECTS   = 

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) mod
#	cp -f *.mod $(SUMB_MODDIR)
	$(MAKE) src

clean:
	@echo "        Making clean in adjoint... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

mod:	$(FF90_MODULES)
src:	$(FF90_OBJECTS) $(CC_OBJECTS)
