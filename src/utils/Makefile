#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 10-21-2007                                      *
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
#      * Names of the Fortran object files to be created.               *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)
FF90_OBJECTS = bsearchIntegers.o\
	       bsearchReals.o\
	       bsearchStrings.o\
	       byteswap.o\
	       computeEtot.o\
	       computeGamma.o\
	       computePtot.o\
	       computeTtot.o\
	       convertToLowerCase.o\
	       delta.o\
	       derivativeRigidRotAngle.o\
	       initExec.o\
	       myIsNAN.o\
	       nullifyCGNSDomPointers.o\
	       nullifyFlowDomPointers.o\
	       orphanAverage.o\
	       qsortIntegers.o\
	       qsortReals.o\
	       qsortStrings.o\
	       reallocateIntegers.o\
	       reallocateReals.o\
	       releaseMemory.o\
	       releaseOversetCGNSGrid.o\
	       removeTiesIntegers.o\
	       removeTiesStrings.o\
	       replaceTabsAndReturns.o\
	       rigidRotAngle.o\
	       rotMatrixRigidBody.o\
	       setBCPointers.o\
               setBufferSizes.o\
	       setCGNSRealType.o\
	       setPointers.o\
	       siAngle.o\
	       siDensity.o\
	       siLen.o\
	       siPressure.o\
	       siTemperature.o\
	       siTurb.o\
	       siVelocity.o\
	       spectralInterpolCoef.o\
	       tempMemory.o\
	       terminate.o\
	       whalo1to1.o\
	       whalo.o\
	       whaloMixing.o\
	       whaloSliding.o\
	       wOverset.o\
	       writeVarsToBuffer.o

CC_OBJECTS   = myIsNaNC.o\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in utils... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
