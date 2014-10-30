#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 10-17-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to build the SUmb library and the        *
#      * executable.                                                    *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * The subdirectories where the sources are located.              *
#      * Note that PV3_INT_SRC_DIR is defined in the config and is      *
#      * empty if no pV3 support is desired.                            *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src/modules       \
		src/ADT           \
		src/bcdata        \
	        src/initFlow      \
	        src/inputParam    \
	        src/metis-4.0     \
	        src/output        \
	        src/overset       \
		src/parallelIO    \
	        src/partitioning  \
	        src/preprocessing \
	        src/slidingComm   \
	        src/solver        \
		src/NKsolver      \
	        src/stabilityDerivatives \
	        src/turbulence    \
	        src/utils         \
	        src/wallDistance  \
		src/warping       \
		src/adjoint               \
                src/adjoint/ADFirstAidKit \
                src/forwardAdjoint \
                src/forwardAdjoint/residualInput \
                src/forwardAdjoint/outputForwardTGT \
                src/forwardAdjoint/outputForward \
		src/forwardAdjoint/costInput \
		src/forwardAdjoint/outputCost\
		src/bendingMomentAnalysis \
                src/adjoint/forcesOutput\
                src/adjoint/forcesInput\
#                src/forwardAdjoint/outputReverse \
#                src/forwardAdjoint/residualOutputReverse\

SUBDIR_EXEC   = src/exec
SUBDIR_PV3    = src/pv3Interface
CONFIG_DEFAULT_DIR = config/defaults
CONFIG_DIR         = config
SUMB_SUBDIRS       = $(SUBDIR_SRC) $(PV3_INT_SRC_DIR)
SUMB_CLEAN_SUBDIRS = $(SUBDIR_SRC)  $(SUBDIR_PV3) $(SUBDIR_EXEC)

#      ******************************************************************
#      *                                                                *
#      * General targets.                                               *
#      *                                                                *
#      ******************************************************************

default:
	@echo "Usage: make <arch> or make <arch>_PYTHON for python wrapper"
	@echo "Currently Supported Architectures:"
	@echo " NOTE: There has been a change in the config files. Please "
	@echo "       use a default config file or modifiy an existing one "
	@echo "       to match the defaults"
	@echo "                         LINUX_INTEL_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI_SCINET"
	@echo "                         LINUX_INTEL_INTELMPI_SCINET"
	@echo "                         LINUX_GFORTRAN_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI_NYX"
	@echo " Previously Supported Architectures"
	@echo "                         LINUX_INTEL_MPICH"
	@echo "                         LINUX_INTEL"
	@echo "                         ABLATION_INTEL_IB"
	@echo "                         ABLATION_PG_IB"
	@echo "                         ALC_INTEL"
	@echo "                         ALTIX"
	@echo "                         ALTIX_MPI"
	@echo "                         ALTIX_MPICH2"
	@echo "                         APPLE_MAC_NAG"
	@echo "                         APPLE_MAC_NAG_MPICH"
	@echo "                         APPLE_MAC_XLF"
	@echo "                         APPLE_MAC_XLF_MPICH"
	@echo "                         ASCI_QSC"
	@echo "                         FLASH_INTEL"
	@echo "                         FLASH_PG"
	@echo "                         IBM_BLUEGENE"
	@echo "                         IBM_DATASTAR"
	@echo "                         LINUX_ABSOFT"
	@echo "                         LINUX_G95"
	@echo "                         LINUX_G95_MPICH"
	@echo "                         LINUX_PG"
	@echo "                         LINUX_PG_MPICH"
	@echo "                         REDHOT_IFC_ETHERNET"
	@echo "                         REDHOT_IFC_MYRINET"
	@echo "                         REDHOT_PG_ETHERNET"
	@echo "                         REDHOT_PG_MYRINET"
	@echo "                         REDSTORM"
	@echo "                         SGI"
	@echo "                         SGI_MPI_ORIGIN"
	@echo "                         SGI_N32"
	@echo "                         SGI_N32_MPICH"
	@echo "                         SGI_N32_MPI_ORIGIN"

all:	 default

dirs:	
	mkdir -p bin
	mkdir -p obj
	mkdir -p mod

clean:
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	@echo " Making clean ... "
	@for subdir in $(SUMB_CLEAN_SUBDIRS) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*

#      ******************************************************************
#      *                                                                *
#      * The actual make. This is not a direct target, but is called    *
#      * from the architectures.                                        *
#      *                                                                *
#      ******************************************************************

sumb:
	@for subdir in $(SUMB_SUBDIRS) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)
	(cd $(SUBDIR_EXEC) && make)

#      ******************************************************************
#      *                                                                *
#      * Currently Supported Platforms                                  *
#      *                                                                *
#      ******************************************************************

LINUX_INTEL_OPENMPI:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	make sumb
	(cd src/python/f2py && make)

LINUX_GFORTRAN_OPENMPI:
	make dirs
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	make sumb
	(cd src/python/f2py &&  make)

LINUX_INTEL_OPENMPI_SCINET:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	make sumb
	(cd src/python/f2py && make)

LINUX_INTEL_INTELMPI_SCINET:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_INTELMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_INTELMPI_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_INTELMPI_SCINET.mk config.mk
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	make sumb
	(cd src/python/f2py && make)

LINUX_INTEL_OPENMPI_NYX:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_NYX.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_NYX.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI_NYX.mk config.mk
	ln -sf SUmb_Common_real.mk SUmb_Common.mk
	make sumb
	(cd src/python/f2py && make)

LINUX_INTEL_OPENMPI_PYTHON:
	@echo "Calling with _PYTHON is no longer necessary."

LINUX_GFORTRAN_OPENMPI_PYTHON:
	@echo "Calling with _PYTHON is no longer necessary."

LINUX_INTEL_OPENMPI_NYX_PYTHON:
	@echo "Calling with _PYTHON is no longer necessary."

LINUX_INTEL_OPENMPI_SCINET_PYTHON:
	@echo "Calling with _PYTHON is no longer necessary."

LINUX_INTEL_INTELMPI_SCINET_PYTHON:
	@echo "Calling with _PYTHON is no longer necessary."

#      ******************************************************************
#      *                                                                *
#      * Previously Supported Platforms                                 *
#      *                                                                *
#      ******************************************************************

ABLATION_INTEL_IB:
ABLATION_PG_IB:
ALC_INTEL:
ALTIX:
ALTIX_MPI:
ALTIX_MPICH2:
APPLE_MAC_NAG:
APPLE_MAC_NAG_MPICH:
APPLE_MAC_XLF:
APPLE_MAC_XLF_MPICH:
ASCI_QSC:
FLASH_INTEL:
FLASH_PG:
IBM_BLUEGENE:
IBM_DATASTAR:
LINUX_ABSOFT:
LINUX_G95:
LINUX_G95_MPICH:
LINUX_G95_OPENMPI:
LINUX_G95_OPENMPI_PYTHON:
LINUX_INTEL:
LINUX_INTEL_MPICH:
LINUX_PG:
LINUX_PG_MPICH:
REDHOT_IFC_ETHERNET:
REDHOT_IFC_MYRINET:
REDHOT_PG_ETHERNET:
REDHOT_PG_MYRINET:
REDSTORM:
SGI:
SGI_MPI_ORIGIN:
SGI_N32:
SGI_N32_MPICH:
SGI_N32_MPI_ORIGIN:
SICORTEX:
SICORTEX_MPI:
