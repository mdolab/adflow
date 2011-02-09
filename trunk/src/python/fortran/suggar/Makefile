#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn,           *
#      *                Steve Repsher                                   *
#      * Starting date: 01-17-2004                                      *
#      * Last modified: 03-02-2006                                      *
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
#      * Include the settings for the directory containing the routines *
#      * for running SUmb for multi-disciplinary applications.          *
#      *                                                                *
#      ******************************************************************

SUMB_MD_DIR = ..
SUMB_MD_COMMON_FILE = $(SUMB_MD_DIR)/SUmb_MD_Common.mk
include ${SUMB_MD_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran object files to be created.               *
#      * These are split in sets, such that no dependencies are present *
#      * within a set. This guarantees that parallel make works         *
#      * correctly.                                                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_MD_OBJDIR)
FF90_OBJECTS_1 = suggarData.o

FF90_OBJECTS_2 = initSuggarInterface.o\
                 loadSuggarDCIFile.o\
                 writePlot3DZoneFile.o

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
	@echo "        Making clean in python/fortran/suggar... "
	rm -f *~ *.o *.mod *.il *.stb */*~ */*.o */*.mod */*.il */*.stb

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(FF90_OBJECTS_2)
