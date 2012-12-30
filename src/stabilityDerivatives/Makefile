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

FF90_MODULES = 

FF90_OBJECTS = computeLeastSquaresRegression.o\
	       computeTSDerivatives.o   \
               stabilityDerivativeDriver.o

CC_OBJECTS   = 

FILES_TO_COMPLEXIFY = computeLeastSquaresRegression.f90\
	       computeTSDerivatives.f90 \
               stabilityDerivativeDriver.f90
FILES_TO_COPY = 

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) mod
	$(PMAKE) src

clean:
	@echo "        Making clean in adjoint... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/stabilityDerivatives/$$file -o ! -e ../../src_cs/stabilityDerivatives/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/stabilityDerivatives/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/stabilityDerivatives/$$file -o ! -e ../../src_cs/stabilityDerivatives/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/stabilityDerivatives/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/stabilityDerivatives


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

mod:	$(FF90_MODULES)
src:	$(FF90_OBJECTS) $(CC_OBJECTS)
