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

FF90_OBJECTS1 =\
		utils_d.o\
		flowutils_d.o\
		bcroutines_d.o\
		samodule_d.o\
		walldistance_d.o\

FF90_OBJECTS2 = allnodalgradients_d.o\
		viscousflux_d.o \
		viscousfluxapprox_d.o\
		invisciddissfluxmatrix_d.o \
		diffSizes.o\
		computelamviscosity_d.o\
		computepressuresimple_d.o\
		timestep_block_d.o\
		invisciddissfluxscalar_d.o\
		invisciddissfluxscalarapprox_d.o\
		invisciddissfluxmatrixapprox_d.o\
		inviscidupwindflux_d.o  \
		block_res_d.o  \
		residual_block_d.o \
	        inviscidcentralflux_d.o\
                saeddyviscosity_d.o\
		xhalo_block_d.o\
		metric_block_d.o\
		gridvelocitiesfinelevel_block_d.o\
		normalvelocities_block_d.o\
		adjustinflowangle_d.o\
		referencestate_d.o\
		forcesandmoments_d.o\
		getdirangle_d.o\
		getdirvector_d.o\
		vectorrotation_d.o\
		computeeddyviscosity_d.o\
		sanuknowneddyratio_d.o\
		prodkatolaunder_d.o\
		prodwmag2_d.o\
		prodsmag2_d.o\
		turbadvection_d.o\
		unsteadyturbterm_d.o\
		slipvelocitiesfinelevel_block_d.o\
		getcostfunction2_d.o\
		bcturbwall_d.o\
                bcturbtreatment_d.o\
                applyallturbbcthisblock_d.o\
                bceddynowall_d.o\
                bceddywall_d.o\
		bcturbfarfield_d.o\
		bcturbinterface_d.o\
		turb2ndhalo_d.o\
		boundarynormals_d.o\
		volume_block_d.o\
		derivativerotmatrixrigid_d.o\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src1
	$(PMAKE) src2
	cp -f *.mod $(SUMB_MODDIR)
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
