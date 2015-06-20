#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 03-07-2003                                      *
#      * Last modified: 09-13-2007                                      *
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
FF90_OBJECTS_1 = restartMod.o

FILES_TO_COMPLEXIFY= restartMod.f90\
		 allocConvArrays.f90\
		 allocMemFlovar.f90\
		 allocTimeArrays.f90\
		 computePressure.f90\
		 copySpectralSolution.f90\
		 determineSolFileNames.f90\
		 getSortedVarNumbers.F90\
		 getSortedZoneNumbers.F90\
		 initBleedsFromHalos.f90\
		 initDepvarAndHalos.f90\
		 initFlow.f90\
		 initFlowfield.f90\
		 initializeHalos.f90\
		 interpolateSpectralSolution.f90\
		 readConvHistory.F90\
		 readConvHistoryPlot3D.f90\
		 readDensity.f90\
		 readDensityPlot3D.f90\
		 readEnergy.f90\
		 readEnergyPlot3D.f90\
		 readPressure.f90\
		 readPressurePlot3D.f90\
		 readRestartFile.F90\
		 readRestartFilePlot3D.f90\
		 readRestartVariable.F90\
		 readTimeHistory.F90\
		 readTimeHistoryPlot3D.f90\
		 readTurbEddyVis.f90\
		 readTurbEddyVisPlot3D.f90\
		 readTurbKwType.f90\
		 readTurbKwTypePlot3D.f90\
		 readTurbSA.f90\
		 readTurbSAPlot3D.f90\
		 readTurbV2f.f90\
		 readTurbV2fPlot3D.f90\
		 readTurbvar.f90\
		 readTurbvarPlot3D.f90\
		 readXmomentum.f90\
		 readXmomentumPlot3D.f90\
		 readXvelocity.f90\
		 readXvelocityPlot3D.f90\
		 readYmomentum.f90\
		 readYmomentumPlot3D.f90\
		 readYvelocity.f90\
		 readYvelocityPlot3D.f90\
		 readZmomentum.f90\
		 readZmomentumPlot3D.f90\
		 readZvelocity.f90\
		 readZvelocityPlot3D.f90\
		 referenceState.F90\
		 releaseExtraMemBCs.f90\
		 scaleFactors.F90\
		 setFlowInfinityState.F90\
		 setIOVar.f90\
		 setGlobalCellsAndNodes.f90\
		 setPressureAndComputeEnergy.f90\
		 setUniformFlow.f90\
		 timeSpectralCoef.f90\
		 timeSpectralMatrices.f90\
		 updateGamma.F90\
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
	@echo "        Making clean in initFlow... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/initFlow/$$file -o ! -e ../../src_cs/initFlow/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/initFlow/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/initFlow/$$file -o ! -e ../../src_cs/initFlow/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/initFlow/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/initFlow

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(OFILES)
