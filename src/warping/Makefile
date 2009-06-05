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

FF90_MODULES = warpingPETSc.o

FF90_OBJECTS = warp_local.o        \
	       warp-mb.o      \
               integratedWarp.o\
               integratedWarpDeriv.o \
               integratedWarpDerivParallel.o \
               integratedWarpDerivFD.o \
               updateFaces.o\
               createWarpingPETScMat.o\
	       flagImplicitFacesAndEdges.o\
	       flagImplicitFacesAndEdgesDeriv.o\
	       blockRelations.o\
	       blockLogicalFunctions.o\
               storeGlobalSurfaceIndices.o\
               synchronizeFaces.o\
               synchronizeIndices.o\
               synchronizeSurfaceIndices.o\
               warpingInitializePetsc.o\
	       testwarping.o


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
