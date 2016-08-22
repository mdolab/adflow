#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-18-2002                                      *
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
#      * Set the additional compiler flags for this library.            *
#      *                                                                *
#      ******************************************************************

ADDITIONAL_CC_FLAGS = -I../metis-4.0

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
#      * These are split in sets, such that no dependencies are present *
#      * within a set. This guarantees that parallel make works         *
#      * correctly.                                                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)
FF90_OBJECTS_1 = partitionMod.o

FILES_TO_COMPLEXIFY = \
		allocCoorFineGrid.f90\
		 blockDistribution.f90\
		 check1to1Subfaces.f90\
		 checkFaces.f90\
		 checkLoadBalance.f90\
		 checkPartitioning.f90\
		 checkTransform.f90\
		 countConnectivities.F90\
		 determineBleedFlowRegions.f90\
		 determineComputeBlocks.f90\
		 determineGridFileNames.f90\
		 determineInterfaceIDs.f90\
		 determineNeighborIDs.f90\
		 fineGridSpectralCoor.f90\
		 initFlowDoms.f90\
		 initFineGridIblank.f90\
		 internalBC.f90\
		 loadBalance.f90\
		 overwriteFamilyData.f90\
		 overwriteWallBCData.f90\
		 partitionAndReadGrid.f90\
		 qsortSubblocksOfCGNSType.f90\
		 qsortSubfaceNonMatchType.f90\
		 read1to1Conn.F90\
		 readBlockSizes.F90\
		 readBocos.F90\
		 readFamilyInfo.F90\
		 readGeneralConn.F90\
		 readLevel0CoolingParameters.f90\
		 readPeriodicSubface.F90\
		 readZoneInfo.F90\
		 reallocSubfaceMemory.f90\
		 sortSubfaces.f90\
		 splitBlock.f90\
		 timePeriodSpectral.f90\
		 timeRotMatricesSpectral.f90\


FILES_TO_COPY = partitionMod.f90\
		sortRangesSplitInfo.f90\
		sortBadEntities.F90\
		graphPartitioning.f90\
		metisInterface.c\
		readGrid.F90\

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
	@echo "        Making clean in partitioning... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/partitioning/$$file -o ! -e ../../src_cs/partitioning/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/partitioning/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/partitioning/$$file -o ! -e ../../src_cs/partitioning/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/partitioning/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/partitioning


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1) 
src_2:  $(OFILES)
