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

FF90_OBJECTS =  computeetot_b.o\
		computelamviscosity_b.o\
		extrapolate2ndhalo_b.o\
		timestep_block_b.o\
		metric_block_b.o\
		xhalo_block_b.o\
		applyallbc_block_b.o\
		residual_block_b.o\
		bceulerwall_b.o\
		bcfarfield_b.o\
		bcsymm_b.o\
		bcnswallisothermal_b.o\
		bcnswalladiabatic_b.o\
		invisciddissfluxscalar_b.o\
		invisciddissfluxscalarcoarse_b.o\
		dim_b.o\
		block_res_b.o  \
		viscousflux_b.o \
		computeetotcellcpfit_b.o \
                setbcpointersbwd_b.o\
		setssbwd_b.o\
		resetssbwd_b.o\
		setpp3pp4bwd_b.o\
		resetpp3pp4bwd_b.o\
		setallxbwd_b.o\
		resetallxbwd_b.o\
		setssmetricbwd_b.o\
		resetssmetricbwd_b.o\
		setgammabwd_b.o\
		resetgammabwd_b.o\
                setww0pp0rlv0rev0bwd_b.o\
		resetww0pp0rlv0rev0bwd_b.o\
                setxxssrhodd2wallbwd_b.o\
                resetxxssrhodd2wallbwd_b.o\
                resetbcpointersbwd_b.o\
	        sa_block_b.o\
	        inviscidcentralflux_b.o\
                saeddyviscosity_b.o\
		computeeddyviscosity_b.o\
		sasolve_b.o\
		unsteadyturbterm_b.o\
		prodkatolaunder_b.o\
		prodwmag2_b.o\
		prodsmag2_b.o\
		turbadvection_b.o\
		forcesandmoments_b.o\
		setflowinfinitystate_b.o\
		sanuknowneddyratio_b.o\
		etotarray_b.o\
		eintarray_b.o\
		getcostfunction_b.o\
		getcostfunction2_b.o\

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
