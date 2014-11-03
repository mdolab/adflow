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
	       initExec.F90\
	       myIsNAN.f90\
	       nullifyCGNSDomPointers.f90\
	       nullifyFlowDomPointers.f90\
	       orphanAverage.f90\
	       pointReduce.F90\
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
	       setBCPointers.F90\
	       resetBCPointers.F90\
	       setss.f90\
	       resetss.f90\
	       setallx.f90\
	       resetallx.f90\
	       setssMetric.f90\
	       resetssMetric.f90\
	       setpp3pp4.f90\
	       resetpp3pp4.f90\
	       setgamma.f90\
	       resetgamma.f90\
	       setxxssrhodd2Wall.f90\
	       resetxxssrhodd2Wall.f90\
	       setww0pp0rlv0rev0.f90\
	       resetww0pp0rlv0rev0.f90\
	       setBCPointersBwd.F90\
	       resetBCPointersBwd.F90\
	       setssBwd.f90\
	       resetssBwd.f90\
	       setallxBwd.f90\
	       resetallxBwd.f90\
	       setssMetricBwd.f90\
	       resetssMetricBwd.f90\
	       setpp3pp4Bwd.f90\
	       resetpp3pp4Bwd.f90\
	       setgammaBwd.f90\
	       resetgammaBwd.f90\
	       setxxssrhodd2WallBwd.f90\
	       resetxxssrhodd2WallBwd.f90\
	       setww0pp0rlv0rev0Bwd.f90\
	       resetww0pp0rlv0rev0Bwd.f90\
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
	       terminate.F90\
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
	       writeVarsToBuffer.f90\
	       getLiftDirFromSymmetry.F90

FILES_TO_COPY = myIsNaNC.c\
		f77flush.c\


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

src:	$(OFILES)
