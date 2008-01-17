#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-21-2003                                      *
#      * Last modified: 11-30-2007                                      *
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
FF90_OBJECTS_1 = bcHalo.o\
		 coarse1to1Subface.o\
		 haloList.o\
		 indirectHalo.o\
		 periodicInfo.o

FF90_OBJECTS_2 = allocMemHaloList.o\
		 bsearchCGNSPeriodicType.o\
		 cellRangeSubface.o\
		 checkCoarse1to1.o\
		 checkSymmetry.o\
		 closestDirectHalos.o\
		 coarseDonorInfo.o\
		 coarseLevel0CoolingParameters.o\
		 createCoarseBlocks.o\
		 determineAreaLevel0Cooling.o\
		 determineCommPattern.o\
		 determineFaceHalos.o\
		 determineIndirectHalos.o\
		 determineNumberOfHalos.o\
		 determineNcellGlobal.o\
		 determinePeriodicData.o\
		 determinePeriodicFaces.o\
		 exchangeCoor.o\
		 faceRotationMatrices.o\
		 finalCommStructures.o\
		 indirectHalosPerLevel.o\
		 init2ndLevelCellHalos.o\
		 mdUpdateRoutines.o\
		 metric.o\
		 preprocessing.o\
		 qsortHaloListType.o\
		 qsortIndHaloType.o\
		 qsortPeriodicSubfacesHaloType.o\
		 setFamilyInfoFaces.o\
		 setPorosities.o\
		 shiftCoorAndVolumes.o\
		 unitVectorsInAxialPlane.o\
		 update1to1Coarse.o\
		 updateCoorFineMesh.o\
		 viscSubfaceInfo.o\
		 xhalo.o

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
	@echo "        Making clean in preprocessing... "
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
