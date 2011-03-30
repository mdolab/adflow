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
FF90_OBJECTS =  blockpointers_d.o\
		diffSizes.o\
		applyallbc_block_d.o\
		bcsymm_d.o\
                computeetot_d.o\
		extrapolate2ndhalo_d.o\
	        invisciddissfluxscalarcoarse_d.o\
		timestep_block_d.o\
		bceulerwall_d.o\
	        bcsymmpolar_d.o\
		curveupre_d.o\
		initres_block_d.o\
		invisciddissfluxscalar_d.o\
		turbbcnswall_d.o\
		bcfarfield_d.o\
		bcturbwall_d.o\
		dim_d.o\
		inviscidcentralflux_d.o \
		inviscidupwindflux_d.o  \
		utauwf_d.o \
		bcnswalladiabatic_d.o \
		block_res_d.o  \
		invisciddissfluxmatrixcoarse_d.o \
		residual_block_d.o \
		viscousflux_d.o \
		bcnswallisothermal_d.o \
		computeetotcellcpfit_d.o \
		etotarray_d.o \
		invisciddissfluxmatrix_d.o \
		setbcpointers_d.o\
		alloc_derivative_values.o\
		dealloc_derivative_values.o\
		setPointersd.o\
		eintarray_d.o  \

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS) 
