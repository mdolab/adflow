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

FF90_OBJECTS_2 = addHalosToBoundary.o\
		 badDonorQuality.o\
		 buildCoarseBoundaryList.o\
		 buildFineBoundaryList.o\
		 cellCentroid.o\
		 changeIblanks.o\
		 checkOverset.o\
		 checkSizeBoundaryList.o\
		 createCoarseBoundary.o\
		 determineOversetComm.o\
		 donorSearchCycle.o\
		 donorWeights.o\
		 estimateCoarseDonor.o\
		 exchangeIblanks.o\
		 fringeToField.o\
		 mergeComm.o\
		 numAdjacentField.o\
		 oversetComm.o\
		 releaseMemOverset.o\
		 scan1to1HalosForFringe.o\
		 setBlockOversetData.o\
		 stencilSearch.o\
		 transformCellIndices.o\
		 updateOversetAllLevels.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src_1
	$(MAKE) src_2

clean:
	@echo "        Making clean in overset... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(FF90_OBJECTS_2)
