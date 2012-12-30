#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        C.A.(Sandy) Mader                               *
#      * Starting date: 01-15-2007                                      *
#      * Last modified: 02-19-2007                                      *
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

FF90_OBJECTS = computeTSStabilityDerivAdj.o\
	       verifyTSStabilityDerivAdj.o\
	       verifydIdwFile.o\
               verifydIdxFile.o

CC_OBJECTS   = 

FILES_TO_COMPLEXIFY = \
		computeTSStabilityDerivAdj.f90\
		verifyTSStabilityDerivAdj.f90\
		verifydIdwFile.F90\
		verifydIdxFile.F90

FILES_TO_COPY = 
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
		if [ $$file -nt ../../../src_cs/adjoint/stabilityInput/$$file -o ! -e ../../../src_cs/adjoint/stabilityInput/$$file ]; then \
			python ../../complexify.py $$file; \
			mv c_$$file ../../../src_cs/adjoint/stabilityInput/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../../src_cs/adjoint/stabilityInput/$$file -o ! -e ../../../src_cs/adjoint/stabilityInput/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../../src_cs/adjoint/stabilityInput/$$file; \
		fi;\
		done

	cp Makefile ../../../src_cs/adjoint/stabilityInput


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
