#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 04-11-2003                                      *
#      * Last modified: 03-29-2006                                      *
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
FF90_OBJECTS_1 = outputMod.o


FILES_TO_COMPLEXIFY = outputMod.f90\
		 copyDataBuf.f90\
		 computeIsoVariable.F90\
		 describeScheme.f90\
		 gridFileNamesWrite.f90\
		 isoSurfNames.f90\
		 marchingCubesData.F90\
		 numberOfSurfSolVariables.f90\
		 numberOfVolSolVariables.f90\
		 numberOfIsoVariables.f90\
		 openCGNSVolumeSol.F90\
		 releaseHelpVariablesWriting.f90\
		 resHalo1.f90\
		 setHelpVariablesWriting.f90\
		 storeOldSolInBuffer.f90\
		 storeSolInBuffer.f90\
		 storeSurfsolInBuffer.f90\
		 surfSolFileNamesWrite.f90\
		 surfSolNames.f90\
		 updateParamfile.f90\
		 volSolFileNamesWrite.f90\
		 volSolNames.f90\
		 writeCGNSBaldwinLomaxInfo.F90\
		 writeCGNSConvInfo.F90\
		 writeCGNSGridFile.F90\
		 writeCGNSGridFrame.F90\
		 writeCGNSHeader.F90\
		 writeCGNSKomegaModifiedInfo.F90\
		 writeCGNSKomegaWilcoxInfo.F90\
		 writeCGNSKtauInfo.F90\
		 writeCGNSMenterSSTInfo.F90\
		 writeCGNSReferenceState.F90\
		 writeCGNSSaInfo.F90\
		 writeCGNSSaeInfo.F90\
		 writeCGNSSurfaceSol.F90\
		 writeCGNSTimeHistory.F90\
		 writeCGNSV2fInfo.F90\
		 writeCGNSVolumeSol.F90\
		 writeCoorCGNSZone.F90\
		 writeInterpolFromSpectral.f90\
		 writeIsoSurface.F90\
		 writeOversetCGNSZone.F90\
		 writePlot3DConnFile.f90\
		 writePlot3DConvInfo.f90\
		 writePlot3DGridFile.f90\
		 writePlot3DSurfaceSol.f90\
		 writePlot3DTimeHistory.f90\
		 writePlot3DVolumeSol.f90\
		 writeSol.F90\
		 writeSolCGNSZone.F90\
		 writeSurfsolCGNSZone.F90\
		 writeUnsteadyFromSpectral.f90

FILES_TO_COPY = 

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
	@echo "        Making clean in output... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/output/$$file -o ! -e ../../src_cs/output/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/output/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/output/$$file -o ! -e ../../src_cs/output/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/output/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/output


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:	$(FF90_OBJECTS_1)
src_2:	$(OFILES)
