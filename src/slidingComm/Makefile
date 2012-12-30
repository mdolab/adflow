#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 09-26-2003                                      *
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
FF90_OBJECTS_1 = localSubfaces.o\
		 mixingData.o\
		 thisSlide.o\
		 tmpSliding.o\
		 updateComm.o

FILES_TO_COMPLEXIFY =  \
		 localSubfaces.f90\
		 thisSlide.f90\
		 tmpSliding.f90\
		 cumulativeNSendReceives.f90\
		 cylCoorNodesOnBlockFace.f90\
		 determineNSlices.f90\
		 determineRotationMatrices.f90\
		 determineSections.f90\
		 haloType.f90\
		 initMemMixingPlane.f90\
		 initMemSliding.f90\
		 initThisSlide.f90\
		 interpolateSlide.f90\
		 localNodesAndQuads.f90\
		 localSurfaceGrids.f90\
		 mixingDonorInterpol.f90\
		 mixingHaloInterpol.f90\
		 mixingIntervals.f90\
		 mixingPlaneComm.f90\
		 mixingPlaneInterpol.f90\
		 mySubfacesSlide.f90\
		 nLocalNodesAndQuads.f90\
		 procSliding.f90\
		 qsortMixingIntervalType.f90\
		 qsortSortedDonorCommType.f90\
		 releaseMemSliding.f90\
		 slidingComm.f90\
		 slidingIntersections.f90\
		 slidingMesh.f90\
		 slidingMeshGroups.f90\
		 updateExternalCommSlide.f90\
		 updateInterpolSendBuf.f90\
		 updateLocalCommSlide.f90\
		 updateRotationInfo.f90\
		 updateSlidingAllLevels.f90\
		 updateSlidingCommPattern.f90\
		 volumeInterpol.f90

FILES_TO_COPY =  mixingData.F90\
		 updateComm.F90\
		 bsearchSortedDonorCommType.f90

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
	@echo "        Making clean in slidingComm... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/slidingComm/$$file -o ! -e ../../src_cs/slidingComm/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/slidingComm/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/slidingComm/$$file -o ! -e ../../src_cs/slidingComm/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/slidingComm/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/slidingComm


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(OFILES)
