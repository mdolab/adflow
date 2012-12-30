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

FF90_OBJECTS = adjustInflowAngle.o\
               applyAllBC_block.o\
	       block_res.o\
               block_res_extra.o\
               block_res_spatial.o\
               checkInputParam_mod.o\
	       computeUtau_block.o\
               gridVelocitiesFineLevel_block.o\
	       initres_block.o\
	       initres_block_mod.o\
               metric_block.o\
               normalVelocities_block.o\
               referenceState_mod.o\
	       residual_block.o\
	       sa_block.o\
	       saSolve_mod.o\
	       timeStep_block.o\
               turbResidual_block.o\
	       xhalo_block.o\

FILES_TO_COMPLEXIFY = \
		adjustInflowAngle.f90\
		applyAllBC_block.f90\
		block_res.f90\
		block_res_extra.f90\
		block_res_spatial.f90\
		checkInputParam_mod.f90\
		computeUtau_block.f90\
		gridVelocitiesFineLevel_block.f90\
		initres_block.f90\
		initres_block_mod.f90\
		metric_block.f90\
		normalVelocities_block.f90\
		referenceState_mod.f90\
		residual_block.f90\
		sa_block.f90\
		saSolve_mod.f90\
		timeStep_block.f90\
		turbResidual_block.f90\
		xhalo_block.f90\

FILES_TO_COPY = 
#      ******************************************************************
#      *                                                        *
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
		if [ $$file -nt ../../../src_cs/forwardAdjoint/residualInput/$$file -o ! -e ../../../src_cs/forwardAdjoint/residualInput/$$file ]; then \
			python ../../complexify.py $$file; \
			mv c_$$file ../../../src_cs/forwardAdjoint/residualInput/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../../src_cs/forwardAdjoint/residualInput/$$file -o ! -e ../../../src_cs/forwardAdjoint/residualInput/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../../src_cs/forwardAdjoint/residualInput/$$file; \
		fi;\
		done

	cp Makefile ../../../src_cs/forwardAdjoint/residualInput


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS) 

