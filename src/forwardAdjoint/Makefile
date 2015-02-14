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
FILES_TO_COMPLEXIFY = \
		alloc_derivative_values.F90\
		dealloc_derivative_values.F90\
		determineStencil.F90\
		setupAllResidualMatricesFwd.F90\
		setupColoring.F90\
		setPointers_d.f90\
		setDiffSizes.F90\
		zeroADSeeds.F90\
		setFDReference.F90\
		referenceShockSensor.F90\

FILES_TO_COPY = setupStateResidualMatrix.F90\
		verifyAD.F90\
		test.F90\
#		verifyOneValue.F90
#		verifyMatProd.F90\


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
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/forwardAdjoint/$$file -o ! -e ../../src_cs/forwardAdjoint/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/forwardAdjoint/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/forwardAdjoint/$$file -o ! -e ../../src_cs/forwardAdjoint/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/forwardAdjoint/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/forwardAdjoint


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(OFILES)
