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
FF90_OBJECTS =  adjustinflowangleadj_extra_d.o\
                adjustinflowangle_extra_d.o\
		applyallbc_block_extra_d.o\
                bceulerwall_extra_d.o\
                bcfarfield_extra_d.o\
		bcsymm_extra_d.o\
                bcsymmpolar_extra_d.o\
                block_res_extra_extra_d.o\
                checkinputparam_mod_extra_d.o\
                computeetotcellcpfit_extra_d.o\
                computeetot_extra_d.o\
		curveupre_extra_d.o\
                derivativerigidrotangle_extra_d.o\
                derivativerotmatrixrigid_extra_d.o\
                dim_extra_d.o\
                eintarray_extra_d.o\
                etotarray_extra_d.o\
                extrapolate2ndhalo_extra_d.o\
                getdirangle_extra_d.o\
                getdirvector_extra_d.o\
                gridvelocitiesfinelevel_block_extra_d.o\
                inviscidcentralflux_extra_d.o\
	        invisciddissfluxmatrixcoarse_extra_d.o\
	        invisciddissfluxmatrix_extra_d.o\
	        invisciddissfluxscalarcoarse_extra_d.o\
	        invisciddissfluxscalar_extra_d.o\
                inviscidupwindflux_extra_d.o\
                normalvelocities_block_extra_d.o\
                referencestate_mod_extra_d.o\
                residual_block_extra_d.o\
                rotmatrixrigidbody_extra_d.o\
                sanuknowneddyratio_extra_d.o\
                setbcpointers_extra_d.o\
                setflowinfinitystate_extra_d.o\
                timestep_block_extra_d.o\
                utauwf_extra_d.o \
                vectorrotation_extra_d.o\
                viscousflux_extra_d.o


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
