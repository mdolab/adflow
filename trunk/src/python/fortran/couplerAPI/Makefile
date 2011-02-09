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
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_MD_OBJDIR)
FF90_OBJECTS = sumb_coupler_m.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src
	cp -f *.mod $(SUMB_MODDIR)

clean:
	@echo "        Making clean in python/fortran/couplerAPI.... "
	rm -f *~ *.o *.mod *.il *.stb

#      ******************************************************************
#      *                                                                *
#      * Rules to make src.                                             *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS)
