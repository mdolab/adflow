#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-13-2002                                      *
#      * Last modified: 11-27-2007                                      *
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
FF90_OBJECTS_1 = localModules.o

FF90_OBJECTS_2 = analyzeString.o\
		 checkInputParam.o\
		 checkMonitor.o\
		 checkOutput.o\
		 createTemplate.o\
		 defaultMonitor.o\
		 defaultSurfaceOut.o\
		 defaultVolumeOut.o\
		 extractMGInfo.o\
		 monitorVariables.o\
		 readCpTempCurveFits.o\
		 readParamFile.o\
		 dummyReadParamFile.o\
		 setDefaultValues.o\
		 setEquationParameters.o\
		 setIterationParam.o\
		 setStageCoeffExplicitRK.o\
		 setStageCoeffImplicitRK.o\
		 surfaceVariables.o\
		 volumeVariables.o\
		 writeInputParam.o

FILES_TO_COMPLEXIFY =localModules.f90\
		 analyzeString.f90\
		 checkInputParam.F90\
		 checkMonitor.f90\
		 checkOutput.f90\
		 createTemplate.f90\
		 defaultMonitor.f90\
		 defaultSurfaceOut.f90\
		 defaultVolumeOut.f90\
		 extractMGInfo.f90\
		 monitorVariables.f90\
		 readCpTempCurveFits.f90\
		 readParamFile.f90\
		 dummyReadParamFile.f90\
		 setDefaultValues.F90\
		 setEquationParameters.f90\
		 setIterationParam.f90\
		 setStageCoeffExplicitRK.f90\
		 setStageCoeffImplicitRK.f90\
		 surfaceVariables.f90\
		 volumeVariables.f90\
		 writeInputParam.f90


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
	@echo "        Making clean in inputParam... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/inputParam/$$file -o ! -e ../../src_cs/inputParam/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/inputParam/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/inputParam/$$file -o ! -e ../../src_cs/inputParam/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/inputParam/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/inputParam
#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(FF90_OBJECTS_2)
