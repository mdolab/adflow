#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 10-27-2005                                      *
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
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FILES_TO_COMPLEXIFY = determineBytesCommChunk.f90\
	       determineCommSizesIO.f90\
	       determineMyChunk.f90\
	       determineVarLocalIO.f90\
	       getIndices.f90\
	       prepareReadPlot3D.f90\
	       prepareWritePlot3D.f90\
	       readPlot3DVar.f90\
	       releaseMemIOPlot3D.f90\
	       sendBufNonLocalIO.f90\
	       writePlot3DVar.f90

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
	$(PMAKE) src

clean:
	@echo "        Making clean in parallelIO... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/parallelIO/$$file -o ! -e ../../src_cs/parallelIO/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/parallelIO/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/parallelIO/$$file -o ! -e ../../src_cs/parallelIO/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/parallelIO/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/parallelIO


#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(OFILES)
