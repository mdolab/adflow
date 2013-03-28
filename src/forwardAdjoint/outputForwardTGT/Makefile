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
FF90_OBJECTS =  blockpointers_t.o\
		applyallbc_block_t.o\
		bcsymm_t.o\
                computeetot_t.o\
		computelamviscosity_t.o\
		extrapolate2ndhalo_t.o\
	        invisciddissfluxscalarcoarse_t.o\
		timestep_block_t.o\
		bceulerwall_t.o\
	  	getcostfuncmat_t.o\
		curveupre_t.o\
		invisciddissfluxscalar_t.o\
		turbbcnswall_t.o\
		bcfarfield_t.o\
		bcturbwall_t.o\
		dim_t.o\
		inviscidupwindflux_t.o  \
		utauwf_t.o \
		bcnswalladiabatic_t.o \
		block_res_t.o  \
		invisciddissfluxmatrixcoarse_t.o \
		residual_block_t.o \
		viscousflux_t.o \
		bcnswallisothermal_t.o \
		computeetotcellcpfit_t.o \
		etotarray_t.o \
		invisciddissfluxmatrix_t.o \
                setbcpointers_t.o\
	        sa_block_t.o\
	        inviscidcentralflux_t.o\
                saeddyviscosity_t.o\
                bcturbtreatment_t.o\
                applyallturbbcthisblock_t.o\
                bceddynowall_t.o\
                bceddywall_t.o\
                bcturbfarfield_t.o\
               	bcturbsymm_t.o\
		eintarray_t.o  \
		xhalo_block_t.o\
		metric_block_t.o\
		gridvelocitiesfinelevel_block_t.o\
		normalvelocities_block_t.o\
		adjustinflowangle_t.o\
		referencestate_t.o\
		setflowinfinitystate_t.o\
		forcesandmoments_t.o\
		derivativerotmatrixrigid_t.o\
		getdirangle_t.o\
		getdirvector_t.o\
		vectorrotation_t.o\
		derivativerigidrotangle_t.o\
		computeeddyviscosity_t.o\
		sanuknowneddyratio_t.o\
		unsteadyturbspectral_block_t.o\
		turb2ndhalo_t.o\
		sasolve_t.o\
		prodkatolaunder_t.o\
		prodwmag2_t.o\
		prodsmag2_t.o\
		turbadvection_t.o\
		unsteadyturbterm_t.o\
		curvetupyp_t.o\


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
