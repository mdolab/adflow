#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 04-26-2006                                      *
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
#      * Include the common file.                                       *
#      *                                                                *
#      ******************************************************************

ADT_DIR = ..
ADT_COMMON_FILE = $(ADT_DIR)/ADT_Common.mk
include ${ADT_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran object files to be created.               *
#      * These are split in sets, such that no dependencies are present *
#      * within a set. This guarantees that parallel make works         *
#      * correctly.                                                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(ADT_OBJDIR)
FF90_OBJECTS_1 = adtPrecision.o
FF90_OBJECTS_2 = adtData.o\
		 gradD2Hexa.o\
		 hessD2Hexa.o\
		 minD2Hexa.o\
		 newtonStep.o
FF90_OBJECTS_3 = adtUtils.o
FF90_OBJECTS_4 = adtBuild.o\
		 adtLocalSearch.o
FF90_OBJECTS_5 = adtSearch.o
FF90_OBJECTS_6 = adtAPI.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all: 
	$(MAKE) src_1
	$(MAKE) src_2
	$(MAKE) src_3
	$(MAKE) src_4
	$(MAKE) src_5
	$(MAKE) src_6
	cp -f *.mod $(ADT_MODDIR)

clean:
	@echo "        Making clean in ADT/src ... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1, etc.       *
#      *                                                                *
#      ******************************************************************

src_1:	$(FF90_OBJECTS_1)
src_2:	$(FF90_OBJECTS_2)
src_3:	$(FF90_OBJECTS_3)
src_4:	$(FF90_OBJECTS_4)
src_5:	$(FF90_OBJECTS_5)
src_6:	$(FF90_OBJECTS_6)
