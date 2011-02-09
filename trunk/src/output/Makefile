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

FF90_OBJECTS_2 = copyDataBuf.o\
		 describeScheme.o\
		 gridFileNamesWrite.o\
		 numberOfSurfSolVariables.o\
		 numberOfVolSolVariables.o\
		 openCGNSVolumeSol.o\
		 releaseHelpVariablesWriting.o\
		 resHalo1.o\
		 setHelpVariablesWriting.o\
		 storeOldSolInBuffer.o\
		 storeSolInBuffer.o\
		 storeSurfsolInBuffer.o\
		 surfSolFileNamesWrite.o\
		 surfSolNames.o\
		 updateParamfile.o\
		 volSolFileNamesWrite.o\
		 volSolNames.o\
		 writeCGNSBaldwinLomaxInfo.o\
		 writeCGNSConvInfo.o\
		 writeCGNSGridFile.o\
		 writeCGNSGridFrame.o\
		 writeCGNSHeader.o\
		 writeCGNSKomegaModifiedInfo.o\
		 writeCGNSKomegaWilcoxInfo.o\
		 writeCGNSKtauInfo.o\
		 writeCGNSMenterSSTInfo.o\
		 writeCGNSReferenceState.o\
		 writeCGNSSaInfo.o\
		 writeCGNSSaeInfo.o\
		 writeCGNSSurfaceSol.o\
		 writeCGNSTimeHistory.o\
		 writeCGNSV2fInfo.o\
		 writeCGNSVolumeSol.o\
		 writeCoorCGNSZone.o\
		 writeInterpolFromSpectral.o\
		 writeOversetCGNSZone.o\
		 writePlot3DConnFile.o\
		 writePlot3DConvInfo.o\
		 writePlot3DGridFile.o\
		 writePlot3DSurfaceSol.o\
		 writePlot3DTimeHistory.o\
		 writePlot3DVolumeSol.o\
		 writeSol.o\
		 writeSolCGNSZone.o\
		 writeSurfsolCGNSZone.o\
		 writeUnsteadyFromSpectral.o

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
	@echo "        Making clean in output... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:	$(FF90_OBJECTS_1)
src_2:	$(FF90_OBJECTS_2)
