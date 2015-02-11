#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 03-10-2003                                      *
#      * Last modified: 11-21-2007                                      *
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

SUMB_DIR = ../../../
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FF90_OBJECTS1 = bcroutines_fast_b.o\
		samodule_fast_b.o\

FF90_OBJECTS2 =  \
		allnodalgradients_fast_b.o\
		computelamviscosity_fast_b.o\
		computespeedofsoundsquared_fast_b.o\
		timestep_block_fast_b.o\
		invisciddissfluxscalar_fast_b.o\
		invisciddissfluxmatrix_fast_b.o\
		dim_fast_b.o\
		viscousflux_fast_b.o \
	        inviscidcentralflux_fast_b.o\
		computepressuresimple_fast_b.o\
		BCExtra.o\
		prodsmag2_fast_b.o\
		prodwmag2_fast_b.o\
		prodkatolaunder_fast_b.o\
		turbadvection_fast_b.o\
		saeddyviscosity_fast_b.o\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src1
	cp -f *.mod $(SUMB_MODDIR)
	$(PMAKE) src2

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src1:   $(FF90_OBJECTS1) 
src2:   $(FF90_OBJECTS2) 
