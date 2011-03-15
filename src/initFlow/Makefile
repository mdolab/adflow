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

FF90_OBJECTS_2 = allocConvArrays.o\
		 allocMemFlovar.o\
		 allocTimeArrays.o\
		 computePressure.o\
		 copySpectralSolution.o\
		 determineSolFileNames.o\
		 getSortedVarNumbers.o\
		 getSortedZoneNumbers.o\
		 initBleedsFromHalos.o\
		 initDepvarAndHalos.o\
		 initFlow.o\
		 initFlowfield.o\
		 initializeHalos.o\
		 interpolateSpectralSolution.o\
		 readConvHistory.o\
		 readConvHistoryPlot3D.o\
		 readDensity.o\
		 readDensityPlot3D.o\
		 readEnergy.o\
		 readEnergyPlot3D.o\
		 readPressure.o\
		 readPressurePlot3D.o\
		 readRestartFile.o\
		 readRestartFilePlot3D.o\
		 readRestartVariable.o\
		 readTimeHistory.o\
		 readTimeHistoryPlot3D.o\
		 readTurbEddyVis.o\
		 readTurbEddyVisPlot3D.o\
		 readTurbKwType.o\
		 readTurbKwTypePlot3D.o\
		 readTurbSA.o\
		 readTurbSAPlot3D.o\
		 readTurbV2f.o\
		 readTurbV2fPlot3D.o\
		 readTurbvar.o\
		 readTurbvarPlot3D.o\
		 readXmomentum.o\
		 readXmomentumPlot3D.o\
		 readXvelocity.o\
		 readXvelocityPlot3D.o\
		 readYmomentum.o\
		 readYmomentumPlot3D.o\
		 readYvelocity.o\
		 readYvelocityPlot3D.o\
		 readZmomentum.o\
		 readZmomentumPlot3D.o\
		 readZvelocity.o\
		 readZvelocityPlot3D.o\
		 referenceState.o\
		 releaseExtraMemBCs.o\
		 scaleFactors.o\
		 setFlowInfinityState.o\
		 setIOVar.o\
		 setPressureAndComputeEnergy.o\
		 setUniformFlow.o\
		 timeSpectralCoef.o\
		 timeSpectralMatrices.o

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
	@echo "        Making clean in initFlow... "
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
