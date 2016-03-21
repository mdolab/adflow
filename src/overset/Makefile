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

FILES_TO_COMPLEXIFY =	computeInterpolationParallel.F90\
			createZipperMesh.F90\
			bowtieAndIsolationElimination.F90\
			determineClusters.F90 \
			determineClusterAreas.F90 \
			checkOverset.F90 \
			exchangeIblanks.f90\
			exchangeIsCompute.f90\
			floodInteriorCells.F90\
			oversetLoadBalance.F90\
			qsortFringeType.F90\
			initializeOBlock.F90\
			initializeOWall.F90\
			initializeFringes.F90\
			determineDonors.F90\
			computeHolesInsideBody.F90\
			packingRoutines.F90\
			fringeSearch.F90\
			oversetUtilities.F90\
			irregularCellCorrection.F90\
			fringeReduction.F90\
			finalOversetCommStructures.F90\
			setIBlankArray.F90\
		        exchangeDonorStatus.F90\
			exchangeFringes.f90\
			surfaceCorrection.F90\
			oversetPresent.F90\
			emptyOversetComm.F90\
			wallSearch.F90\
			initBCDataiBlank.F90\
			exchangeSurfaceIBlanks.F90\
			makeBoundaryStrings.F90\
			writeDualGrid.F90\
		       # computeQuadraticWeights.F90 \
		       matrixinv3by3.F90 \



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
