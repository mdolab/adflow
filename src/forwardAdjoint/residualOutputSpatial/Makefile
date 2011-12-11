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
#SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
SUMB_RULES_FILE  = adjointRules.mk

include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)
FF90_OBJECTS =  accuracy_spatial_d.o\
		applyallbc_block_spatial_d.o\
		bceulerwall_spatial_d.o\
		bcfarfield_spatial_d.o\
		bcsymmpolar_spatial_d.o\
		bcsymm_spatial_d.o\
		bctypes_spatial_d.o\
		block_res_spatial_spatial_d.o\
		block_spatial_d.o\
		cgnsgrid_spatial_d.o\
		communication_spatial_d.o\
		computeetotcellcpfit_spatial_d.o\
		computeetot_spatial_d.o\
		constants_spatial_d.o\
		cpcurvefits_spatial_d.o\
		curveupre_spatial_d.o\
		dim_spatial_d.o\
		eintarray_spatial_d.o\
		etotarray_spatial_d.o\
		extrapolate2ndhalo_spatial_d.o\
		flowvarrefstate_spatial_d.o\
		initres_block_ts_spatial_d.o\
		inputadjoint_spatial_d.o\
		inputdiscretization_spatial_d.o\
		inputio_spatial_d.o\
		inputiteration_spatial_d.o\
		inputmotion_spatial_d.o\
		inputoverset_spatial_d.o\
		inputparallel_spatial_d.o\
		inputphysics_spatial_d.o\
		inputtimespectral_spatial_d.o\
		inputtsstabderiv_spatial_d.o\
		inputunsteady_spatial_d.o\
		inputvisualization_spatial_d.o\
		inviscidcentralflux_spatial_d.o\
		invisciddissfluxmatrixcoarse_spatial_d.o\
		invisciddissfluxmatrix_spatial_d.o\
		invisciddissfluxscalarcoarse_spatial_d.o\
		invisciddissfluxscalar_spatial_d.o\
		inviscidupwindflux_spatial_d.o\
		iteration_spatial_d.o\
		metric_block_spatial_d.o\
		paramturb_spatial_d.o\
		precision_spatial_d.o\
		residual_block_spatial_d.o\
		section_spatial_d.o\
		setbcpointers_spatial_d.o\
		su_mpi_spatial_d.o\
		timestep_block_spatial_d.o\
		turbmod_spatial_d.o\
		utauwf_spatial_d.o\
		viscousflux_spatial_d.o\
		xhalo_block_spatial_d.o\
                normalvelocities_block_spatial_d.o\
                gridvelocitiesfinelevel_block_spatial_d.o


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
