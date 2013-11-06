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
		computegamma_d.o\
		computelamviscosity_d.o\
		extrapolate2ndhalo_d.o\
	        invisciddissfluxscalarcoarse_d.o\
		timestep_block_d.o\
		bceulerwall_d.o\
		curveupre_d.o\
		invisciddissfluxscalar_d.o\
		turbbcnswall_d.o\
		bcfarfield_d.o\
		bcturbwall_d.o\
		dim_d.o\
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
	        sa_block_d.o\
	        inviscidcentralflux_d.o\
                saeddyviscosity_d.o\
                bcturbtreatment_d.o\
                applyallturbbcthisblock_d.o\
                bceddynowall_d.o\
                bceddywall_d.o\
		bcturbfarfield_d.o\
		bcturbsymm_d.o\
		eintarray_d.o  \
		xhalo_block_d.o\
		metric_block_d.o\
		gridvelocitiesfinelevel_block_d.o\
		normalvelocities_block_d.o\
		adjustinflowangle_d.o\
		referencestate_d.o\
		setflowinfinitystate_d.o\
		forcesandmoments_d.o\
		derivativerotmatrixrigid_d.o\
		getdirangle_d.o\
		getdirvector_d.o\
		vectorrotation_d.o\
		derivativerigidrotangle_d.o\
		computeeddyviscosity_d.o\
		sanuknowneddyratio_d.o\
		unsteadyturbspectral_block_d.o\
		turb2ndhalo_d.o\
		sasolve_d.o\
		prodkatolaunder_d.o\
		prodwmag2_d.o\
		prodsmag2_d.o\
		turbadvection_d.o\
		unsteadyturbterm_d.o\
		curvetupyp_d.o\
		slipvelocitiesfinelevel_block_d.o\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS) 
