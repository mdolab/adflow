#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Steve Repsher                                   *
#      * Starting date: 04-06-2005                                      *
#      * Last modified: 02-24-2006                                      *
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
#      * These are split in sets, such that no dependencies are present *
#      * within a set. This guarantees that parallel make works         *
#      * correctly.                                                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)
FF90_OBJECTS_1 = boundaryList.o\
		 searchMod.o

FILES_TO_COMPLEXIFY =  boundaryList.f90\
		 searchMod.f90\
		 addHalosToBoundary.f90\
		 badDonorQuality.f90\
		 buildCoarseBoundaryList.f90\
		 buildFineBoundaryList.f90\
		 cellCentroid.f90\
		 changeIblanks.f90\
		 checkOverset.f90\
		 checkSizeBoundaryList.f90\
		 createCoarseBoundary.f90\
		 determineOversetComm.f90\
		 donorSearchCycle.f90\
		 donorWeights.f90\
		 estimateCoarseDonor.f90\
		 exchangeIblanks.f90\
		 fringeToField.f90\
		 mergeComm.f90\
		 numAdjacentField.f90\
		 oversetComm.f90\
		 releaseMemOverset.f90\
		 scan1to1HalosForFringe.f90\
		 setBlockOversetData.f90\
		 stencilSearch.f90\
		 transformCellIndices.f90\
		 updateOversetAllLevels.f90

ALL_FILES = $(FILES_TO_COMPLEXIFY) $(FILES_TO_COPY)
OFILES := $(ALL_FILES:%.f90=%.o) $(ALL_FILES:%.F90=%.o) $(ALL_FILES:%.c=%.o)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src_1
	$(PMAKE) src_2

clean:
	@echo "        Making clean in overset... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/overset/$$file -o ! -e ../../src_cs/overset/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/overset/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/overset/$$file -o ! -e ../../src_cs/overset/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/overset/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/overset


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(OFILES)
