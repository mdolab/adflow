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

FF90_OBJECTS = 	computePressureForcesAdj.o\
	        computeForcesPressureAdj.o\
	        computeForcesAdj.o\
	        bcEulerWallForcesAdj.o\
	        getSurfaceNormalsForcesAdj.o\
	        forcesAdj.o\
	        computeForceAndMomentAdj.o\
		adjustInflowAngleForcesAdj.o\
		getDirVectorForces.o\
		vectorRotationForces.o\
                verifydCfdx.o

FILES_TO_COMPLEXIFY  = \
		computePressureForcesAdj.f90\
	        computeForcesPressureAdj.f90\
	        computeForcesAdj.f90\
	        bcEulerWallForcesAdj.f90\
	        getSurfaceNormalsForcesAdj.f90\
	        forcesAdj.f90\
	        verifyForces.F90\
	        computeForceAndMomentAdj.f90\
		adjustInflowAngleForcesAdj.f90\
		getDirVectorForces.f90\
		vectorRotationForces.f90\
                verifydCfdx.F90

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
		if [ $$file -nt ../../../src_cs/adjoint/forcesInput/$$file -o ! -e ../../../src_cs/adjoint/forcesInput/$$file ]; then \
			python ../../complexify.py $$file; \
			mv c_$$file ../../../src_cs/adjoint/forcesInput/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../../src_cs/adjoint/forcesInput/$$file -o ! -e ../../../src_cs/adjoint/forcesInput/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../../src_cs/adjoint/forcesInput/$$file; \
		fi;\
		done

	cp Makefile ../../../src_cs/adjoint/forcesInput

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
