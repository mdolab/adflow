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

FILES_TO_COMPLEXIFY = \
		 coarse1to1Subface.f90\
		 allocMemHaloList.f90\
		 cellRangeSubface.f90\
		 checkCoarse1to1.f90\
		 checkSymmetry.f90\
		 closestDirectHalos.f90\
		 coarseDonorInfo.f90\
		 createCoarseBlocks.f90\
		 determineCommPattern.f90\
		 determineFaceHalos.f90\
		 determineIndirectHalos.f90\
		 determineNumberOfHalos.f90\
		 determineNcellGlobal.f90\
		 determinePeriodicFaces.f90\
		 exchangeCoor.f90\
		 exchangeCoor_d.f90\
		 exchangeCoor_b.f90\
		 faceRotationMatrices.f90\
		 finalCommStructures.f90\
		 indirectHalosPerLevel.f90\
		 init2ndLevelCellHalos.f90\
		 setSurfaceFamilyInfo.F90\
		 mdUpdateRoutines.f90\
		 metric.f90\
		 preprocessing.f90\
		 preprocessingCustomOverset.F90\
		 qsortHaloListType.f90\
		 qsortIndHaloType.f90\
		 qsortPeriodicSubfacesHaloType.f90\
		 setFamilyInfoFaces.f90\
		 setGlobalCellsAndNodes.f90\
		 setPorosities.f90\
		 shiftCoorAndVolumes.f90\
		 unitVectorsInAxialPlane.f90\
		 update1to1Coarse.f90\
		 updateCoorFineMesh.f90\
		 viscSubfaceInfo.f90\
		 xhalo.f90\
		 setReferenceVolume.F90\

FILES_TO_COPY = bcHalo.f90\
		indirectHalo.f90\
		haloList.F90\
		periodicInfo.f90\
		bsearchCGNSPeriodicType.f90\
		determinePeriodicData.F90\

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
	@echo "        Making clean in preprocessing... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/preprocessing/$$file -o ! -e ../../src_cs/preprocessing/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/preprocessing/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/preprocessing/$$file -o ! -e ../../src_cs/preprocessing/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/preprocessing/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/preprocessing


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(OFILES)
