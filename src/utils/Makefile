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
	       getBlockData.o\
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
               secondDerivativeRigidRotAngle.o\
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
	       SUmb_finalize.o\
	       SUmb_init.o\
	       tempMemory.o\
	       terminate.o\
               TSalpha.o\
	       TSalphadot.o\
               TSBeta.o\
               TSBetadot.o\
               TSMach.o\
               TSMachdot.o\
	       unique.o\
	       whalo1to1.o\
	       whalo.o\
	       whaloMixing.o\
	       whaloSliding.o\
	       wOverset.o\
	       writeVarsToBuffer.o

CC_OBJECTS   = myIsNaNC.o\
	       f77flush.o

FILES_TO_COMPLEXIFY =  \
	       bsearchIntegers.f90\
	       bsearchReals.f90\
	       bsearchStrings.f90\
	       byteswap.f90\
	       computeEtot.f90\
	       computeGamma.f90\
	       computePtot.f90\
	       computeTtot.f90\
	       convertToLowerCase.f90\
	       delta.f90\
	       derivativeRigidRotAngle.f90\
	       getBlockData.f90\
	       initExec.F90\
	       myIsNAN.f90\
	       nullifyCGNSDomPointers.f90\
	       nullifyFlowDomPointers.f90\
	       orphanAverage.f90\
	       qsortIntegers.f90\
	       qsortReals.f90\
	       qsortStrings.f90\
	       reallocateIntegers.f90\
	       reallocateReals.f90\
	       releaseMemory.f90\
	       releaseOversetCGNSGrid.f90\
	       removeTiesIntegers.f90\
	       removeTiesStrings.f90\
	       replaceTabsAndReturns.f90\
	       rigidRotAngle.f90\
	       rotMatrixRigidBody.f90\
               secondDerivativeRigidRotAngle.f90\
	       setBCPointers.f90\
               setBufferSizes.f90\
	       setCGNSRealType.F90\
	       setPointers.f90\
	       siAngle.f90\
	       siDensity.f90\
	       siLen.f90\
	       siPressure.f90\
	       siTemperature.f90\
	       siTurb.f90\
	       siVelocity.f90\
	       spectralInterpolCoef.f90\
	       SUmb_finalize.f90\
	       SUmb_init.f90\
	       tempMemory.f90\
	       terminate.f90\
               TSalpha.f90\
	       TSalphadot.f90\
               TSBeta.f90\
               TSBetadot.f90\
               TSMach.f90\
               TSMachdot.f90\
	       unique.f90\
	       whalo1to1.f90\
	       whalo.f90\
	       whaloMixing.f90\
	       whaloSliding.f90\
	       wOverset.f90\
	       writeVarsToBuffer.f90


FILES_TO_COPY = myIsNaNC.c\
		f77flush.c\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src

clean:
	@echo "        Making clean in utils... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/utils/$$file -o ! -e ../../src_cs/utils/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/utils/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/utils/$$file -o ! -e ../../src_cs/utils/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/utils/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/utils

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
