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

FF90_OBJECTS_2 = allocCoorFineGrid.o\
		 blockDistribution.o\
		 cgnsBCForPlot3D.o\
		 check1to1Subfaces.o\
		 checkFaces.o\
		 checkLoadBalance.o\
		 checkTransform.o\
		 countConnectivities.o\
		 determineBleedFlowRegions.o\
		 determineComputeBlocks.o\
		 determineGridFileNames.o\
		 determineInterfaceIDs.o\
		 determineNeighborIDs.o\
		 distributeOversetCells.o\
		 fineGridSpectralCoor.o\
		 graphPartitioning.o\
		 initFlowDoms.o\
		 initFineGridIblank.o\
		 internalBC.o\
		 loadBalance.o\
		 overwriteFamilyData.o\
		 partitionAndReadGrid.o\
		 qsortSubblocksOfCGNSType.o\
		 qsortSubfaceNonMatchType.o\
		 read1to1Conn.o\
		 readBlockSizes.o\
		 readBlockSizesPlot3D.o\
		 readBocos.o\
		 readConnectivityPlot3D.o\
		 readFamilyInfo.o\
		 readGeneralConn.o\
		 readGrid.o\
		 readGridPlot3D.o\
		 readLevel0CoolingParameters.o\
		 readOversetHoles.o\
		 readOversetPlot3D.o\
		 readPeriodicSubface.o\
		 readZoneInfo.o\
		 reallocSubfaceMemory.o\
		 setChangingOverset.o\
		 sortSubfaces.o\
		 sortBadEntities.o\
		 sortRangesSplitInfo.o\
		 splitBlock.o\
		 timePeriodSpectral.o\
		 timeRotMatricesSpectral.o

CC_OBJECTS     = metisInterface.o

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
	@echo "        Making clean in partitioning... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1) $(CC_OBJECTS)
src_2:  $(FF90_OBJECTS_2)
