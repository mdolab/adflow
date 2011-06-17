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
	        src/wallDistance

SUBDIR_ADJOINT = src/adjoint               \
		 src/adjoint/ADFirstAidKit \
		 src/adjoint/residualInput \
		 src/adjoint/residualOutput\
                 src/adjoint/stabilityInput\
	         src/adjoint/stabilityOutput\
		 src/adjoint/forcesInput\
	         src/adjoint/forcesOutput\
		 src/adjoint/residualNKPCInput \
	         src/adjoint/residualNKPCOutput \
		 src/forwardAdjoint \
		 src/forwardAdjoint/residualInput \
		 src/forwardAdjoint/residualOutput


SUBDIR_WARPING = src/warping \
		 src/warping/warpingOutput

SUBDIR_EXEC   = src/exec
SUBDIR_PYTHON = src/python/fortran
SUBDIR_PV3    = src/pv3Interface
CONFIG_DEFAULT_DIR = config/defaults
CONFIG_DIR         = config
#SUMB_SUBDIRS       = $(SUBDIR_SRC) $(PV3_INT_SRC_DIR)
SUMB_SUBDIRS       = $(SUBDIR_SRC) $(PV3_INT_SRC_DIR)\
		     $(SUBDIR_PYTHON) $(SUBDIR_WARPING)\
                     $(SUBDIR_ADJOINT)

SUMB_CLEAN_SUBDIRS = $(SUBDIR_SRC) $(SUBDIR_PYTHON) $(SUBDIR_PV3) \
		     $(SUBDIR_EXEC) $(SUBDIR_WARPING) $(SUBDIR_ADJOINT)

#      ******************************************************************
#      *                                                                *
#      * General targets.                                               *
#      *                                                                *
#      ******************************************************************

default:
	@echo "Usage: make <arch>"
	@echo "Supported architectures: ABLATION_INTEL_IB"
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
	@echo "                         LINUX_INTEL"
	@echo "                         LINUX_INTEL_MPICH"
	@echo "                         LINUX_INTEL_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI_SCINET"
	@echo "                         LINUX_INTEL_OPENMPI_SCINET_PYTHON"
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
	mkdir -p externals/ADT/obj
	mkdir -p externals/ADT/mod
	mkdir -p externals/SU_MPI/obj
	mkdir -p externals/SU_MPI/mod
	mkdir -p src/python/fortran/obj

clean:
	@echo " Making clean ... "

	(cd externals/SU_MPI && gmake clean)
	(cd externals/ADT && gmake clean)

	@for subdir in $(SUMB_CLEAN_SUBDIRS) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && gmake $@) || exit 1; \
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
			(cd $$subdir && gmake) || exit 1; \
		done
	(cd lib && gmake)
	(cd $(SUBDIR_EXEC) && gmake)

#      ******************************************************************
#      *                                                                *
#      * Platform specific targets.                                     *
#      *                                                                *
#      ******************************************************************

ABLATION_INTEL_IB:
	make dirs
	if [ ! -f "config/config.ABLATION_INTEL_IB.mk" ]; then cp "config/defaults/config.ABLATION_INTEL_IB.mk" ./config; fi
	(cd externals/SU_MPI && gmake ABLATION_INTEL_IB)
	(cd externals/ADT && gmake ABLATION_INTEL_IB)
	ln -sf config/config.ABLATION_INTEL_IB.mk config.mk
	gmake sumb

ABLATION_PG_IB:
	make dirs
	if [ ! -f "config/config.ABLATION_PG_IB.mk" ]; then cp "config/defaults/config.ABLATION_PG_IB.mk" ./config; fi
	(cd externals/SU_MPI && gmake ABLATION_PG_IB)
	(cd externals/ADT && gmake ABLATION_PG_IB)
	ln -sf config/config.ABLATION_PG_IB.mk config.mk
	gmake sumb

ALC_INTEL:
	make dirs
	if [ ! -f "config/config.ALC_INTEL.mk" ]; then cp "config/defaults/config.ALC_INTEL.mk" ./config; fi
	(cd externals/SU_MPI && gmake ALC_INTEL)
	(cd externals/ADT && gmake ALC_INTEL)
	ln -sf config/config.ALC_INTEL.mk config.mk
	gmake sumb

ALTIX:
	make dirs
	if [ ! -f "config/config.ALTIX.mk" ]; then cp "config/defaults/config.ALTIX.mk" ./config; fi
	(cd externals/SU_MPI && gmake ALTIX)
	(cd externals/ADT && gmake ALTIX)
	ln -sf config/config.ALTIX.mk config.mk
	gmake sumb

ALTIX_MPI:
	make dirs
	if [ ! -f "config/config.ALTIX_MPI.mk" ]; then cp "config/defaults/config.ALTIX_MPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake ALTIX_MPI)
	(cd externals/ADT && gmake ALTIX_MPI)
	ln -sf config/config.ALTIX_MPI.mk config.mk
	gmake sumb

ALTIX_MPICH2:
	make dirs
	if [ ! -f "config/config.ALTIX_MPICH2.mk" ]; then cp "config/defaults/config.ALTIX_MPICH2.mk" ./config; fi
	(cd externals/SU_MPI && gmake ALTIX_MPICH2)
	(cd externals/ADT && gmake ALTIX_MPICH2)
	ln -sf config/config.ALTIX_MPICH2.mk config.mk
	gmake sumb

APPLE_MAC_NAG:
	make dirs
	if [ ! -f "config/config.APPLE_MAC_NAG.mk" ]; then cp "config/defaults/config.APPLE_MAC_NAG.mk" ./config; fi
	(cd externals/SU_MPI && gmake APPLE_MAC_NAG)
	(cd externals/ADT && gmake APPLE_MAC_NAG)
	ln -sf config/config.APPLE_MAC_NAG.mk config.mk
	gmake sumb

APPLE_MAC_NAG_MPICH:
	make dirs
	if [ ! -f "config/config.APPLE_MAC_NAG_MPICH.mk" ]; then cp "config/defaults/config.APPLE_MAC_NAG_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake APPLE_MAC_NAG_MPICH)
	(cd externals/ADT && gmake APPLE_MAC_NAG_MPICH)
	ln -sf config/config.APPLE_MAC_NAG_MPICH.mk config.mk
	gmake sumb

APPLE_MAC_XLF:
	make dirs
	if [ ! -f "config/config.APPLE_MAC_XLF.mk" ]; then cp "config/defaults/config.APPLE_MAC_XLF.mk" ./config; fi
	(cd externals/SU_MPI && gmake APPLE_MAC_XLF)
	(cd externals/ADT && gmake APPLE_MAC_XLF)
	ln -sf config/config.APPLE_MAC_XLF.mk config.mk
	gmake sumb

APPLE_MAC_XLF_MPICH:
	make dirs
	if [ ! -f "config/config.APPLE_MAC_XLF_MPICH.mk" ]; then cp "config/defaults/config.APPLE_MAC_XLF_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake APPLE_MAC_XLF_MPICH)
	(cd externals/ADT && gmake APPLE_MAC_XLF_MPICH)
	ln -sf config/config.APPLE_MAC_XLF_MPICH.mk config.mk
	gmake sumb

ASCI_QSC:
	make dirs
	if [ ! -f "config/config.ASCI_QSC.mk" ]; then cp "config/defaults/config.ASCI_QSC.mk" ./config; fi
	(cd externals/SU_MPI && gmake ASCI_QSC)
	(cd externals/ADT && gmake ASCI_QSC)
	ln -sf config/config.ASCI_QSC.mk config.mk
	gmake sumb

FLASH_INTEL:
	make dirs
	if [ ! -f "config/config.FLASH_INTEL.mk" ]; then cp "config/defaults/config.FLASH_INTEL.mk" ./config; fi
	(cd externals/SU_MPI && gmake FLASH_INTEL)
	(cd externals/ADT && gmake FLASH_INTEL)
	ln -sf config/config.FLASH_INTEL.mk config.mk
	gmake sumb

FLASH_PG:
	make dirs
	if [ ! -f "config/config.FLASH_PG.mk" ]; then cp "config/defaults/config.FLASH_PG.mk" ./config; fi
	(cd externals/SU_MPI && gmake FLASH_PG)
	(cd externals/ADT && gmake FLASH_PG)
	ln -sf config/config.FLASH_PG.mk config.mk
	gmake sumb

IBM_BLUEGENE:
	make dirs
	if [ ! -f "config/config.IBM_BLUEGENE.mk" ]; then cp "config/defaults/config.IBM_BLUEGENE.mk" ./config; fi
	(cd externals/SU_MPI && gmake IBM_BLUEGENE)
	(cd externals/ADT && gmake IBM_BLUEGENE)
	ln -sf config/config.IBM_BLUEGENE.mk config.mk
	gmake sumb

IBM_DATASTAR:
	make dirs
	if [ ! -f "config/config.IBM_DATASTAR.mk" ]; then cp "config/defaults/config.IBM_DATASTAR.mk" ./config; fi
	(cd externals/SU_MPI && gmake IBM_DATASTAR)
	(cd externals/ADT && gmake IBM_DATASTAR)
	ln -sf config/config.IBM_DATASTAR.mk config.mk
	gmake sumb

LINUX_ABSOFT:
	make dirs
	if [ ! -f "config/config.LINUX_ABSOFT.mk" ]; then cp "config/defaults/config.LINUX_ABSOFT.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_ABSOFT)
	(cd externals/ADT && gmake LINUX_ABSOFT)
	ln -sf config/config.LINUX_ABSOFT.mk config.mk
	gmake sumb

LINUX_G95:
	make dirs
	if [ ! -f "config/config.LINUX_G95.mk" ]; then cp "config/defaults/config.LINUX_G95.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_G95)
	(cd externals/ADT && gmake LINUX_G95)
	ln -sf config/config.LINUX_G95.mk config.mk
	gmake sumb

LINUX_G95_MPICH:
	make dirs
	if [ ! -f "config/config.LINUX_G95_MPICH.mk" ]; then cp "config/defaults/config.LINUX_G95_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_G95_MPICH)
	(cd externals/ADT && gmake LINUX_G95_MPICH)
	ln -sf config/config.LINUX_G95_MPICH.mk config.mk
	gmake sumb

LINUX_G95_OPENMPI:
	make dirs
	if [ ! -f "config/config.LINUX_G95_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_G95_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_G95_OPENMPI)
	(cd externals/ADT && gmake LINUX_G95_OPENMPI)
	ln -sf config/config.LINUX_G95_OPENMPI.mk config.mk
	gmake sumb

LINUX_G95_OPENMPI_PYTHON:
	make dirs
	if [ ! -f "config/config.LINUX_G95_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_G95_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_G95_OPENMPI)
	(cd externals/ADT && gmake LINUX_G95_OPENMPI)
	ln -sf config/config.LINUX_G95_OPENMPI.mk config.mk
	gmake sumb
	(cd src/python/f2py && gmake)

LINUX_INTEL:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL.mk" ]; then cp "config/defaults/config.LINUX_INTEL.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL)
	(cd externals/ADT && gmake LINUX_INTEL)
	ln -sf config/config.LINUX_INTEL.mk config.mk
	gmake sumb

LINUX_INTEL_MPICH:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_MPICH.mk" ]; then cp "config/defaults/config.LINUX_INTEL_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_MPICH)
	(cd externals/ADT && gmake LINUX_INTEL_MPICH)
	ln -sf config/config.LINUX_INTEL_MPICH.mk config.mk
	gmake sumb

LINUX_INTEL_OPENMPI:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI)
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	gmake sumb

LINUX_INTEL_OPENMPI_PYTHON:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI)
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	gmake sumb
	(cd src/python/f2py && gmake)

LINUX_GFORTRAN_OPENMPI:
	make dirs
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_GFORTRAN)
	(cd externals/ADT && gmake LINUX_GFORTRAN)
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	gmake sumb

LINUX_GFORTRAN_OPENMPI_PYTHON:
	make dirs
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_GFORTRAN)
	(cd externals/ADT && gmake LINUX_GFORTRAN)
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	gmake sumb
	(cd src/python/f2py && gmake)


LINUX_INTEL_OPENMPI_SCINET:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_SCINET.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI_SCINET)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI_SCINET)
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	gmake sumb

LINUX_INTEL_OPENMPI_SCINET_PYTHON:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_SCINET.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI_SCINET)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI_SCINET)
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	gmake sumb
	(cd src/python/f2py && gmake)

LINUX_INTEL_OPENMPI_NYX:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_NYX.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_NYX.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI_SCINET)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI_SCINET)
	ln -sf config/config.LINUX_INTEL_OPENMPI_NYX.mk config.mk
	gmake sumb

LINUX_INTEL_OPENMPI_NYX_PYTHON:
	make dirs
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_NYX.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_NYX.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_INTEL_OPENMPI_SCINET)
	(cd externals/ADT && gmake LINUX_INTEL_OPENMPI_SCINET)
	ln -sf config/config.LINUX_INTEL_OPENMPI_NYX.mk config.mk
	gmake sumb
	(cd src/python/f2py && gmake)

LINUX_PG:
	make dirs
	if [ ! -f "config/config.LINUX_PG.mk" ]; then cp "config/defaults/config.LINUX_PG.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_PG)
	(cd externals/ADT && gmake LINUX_PG)
	ln -sf config/config.LINUX_PG.mk config.mk
	gmake sumb

LINUX_PG_MPICH:
	make dirs
	if [ ! -f "config/config.LINUX_PG_MPICH.mk" ]; then cp "config/defaults/config.LINUX_PG_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake LINUX_PG_MPICH)
	(cd externals/ADT && gmake LINUX_PG_MPICH)
	ln -sf config/config.LINUX_PG_MPICH.mk config.mk
	gmake sumb

REDHOT_IFC_ETHERNET:
	make dirs
	if [ ! -f "config/config.REDHOT_IFC_ETHERNET.mk" ]; then cp "config/defaults/config.REDHOT_IFC_ETHERNET.mk" ./config; fi
	(cd externals/SU_MPI && gmake REDHOT_IFC_ETHERNET)
	(cd externals/ADT && gmake REDHOT_IFC_ETHERNET)
	ln -sf config/config.REDHOT_IFC_ETHERNET.mk config.mk
	gmake sumb

REDHOT_IFC_MYRINET:
	make dirs
	if [ ! -f "config/config.REDHOT_IFC_MYRINET.mk" ]; then cp "config/defaults/config.REDHOT_IFC_MYRINET.mk" ./config; fi
	(cd externals/SU_MPI && gmake REDHOT_IFC_MYRINET)
	(cd externals/ADT && gmake REDHOT_IFC_MYRINET)
	ln -sf config/config.REDHOT_IFC_MYRINET.mk config.mk
	gmake sumb

REDHOT_PG_ETHERNET:
	make dirs
	if [ ! -f "config/config.REDHOT_PG_ETHERNET.mk" ]; then cp "config/defaults/config.REDHOT_PG_ETHERNET.mk" ./config; fi
	(cd externals/SU_MPI && gmake REDHOT_PG_ETHERNET)
	(cd externals/ADT && gmake REDHOT_PG_ETHERNET)
	ln -sf config/config.REDHOT_PG_ETHERNET.mk config.mk
	gmake sumb

REDHOT_PG_MYRINET:
	make dirs
	if [ ! -f "config/config.REDHOT_PG_MYRINET.mk" ]; then cp "config/defaults/config.REDHOT_PG_MYRINET.mk" ./config; fi
	(cd externals/SU_MPI && gmake REDHOT_PG_MYRINET)
	(cd externals/ADT && gmake REDHOT_PG_MYRINET)
	ln -sf config/config.REDHOT_PG_MYRINET.mk config.mk
	gmake sumb

REDSTORM:
	make dirs
	if [ ! -f "config/config.REDSTORM.mk" ]; then cp "config/defaults/config.REDSTORM.mk" ./config; fi
	(cd externals/SU_MPI && gmake REDSTORM)
	(cd externals/ADT && gmake REDSTORM)
	ln -sf config/config.REDSTORM.mk config.mk
	gmake sumb

SGI:
	make dirs
	if [ ! -f "config/config.SGI.mk" ]; then cp "config/defaults/config.SGI.mk" ./config; fi
	(cd externals/SU_MPI && gmake SGI)
	(cd externals/ADT && gmake SGI)
	ln -sf config/config.SGI.mk config.mk
	gmake sumb

SGI_MPI_ORIGIN:
	make dirs
	if [ ! -f "config/config.SGI_MPI_ORIGIN.mk" ]; then cp "config/defaults/config.SGI_MPI_ORIGIN.mk" ./config; fi
	(cd externals/SU_MPI && gmake SGI_MPI_ORIGIN)
	(cd externals/ADT && gmake SGI_MPI_ORIGIN)
	ln -sf config/config.SGI_MPI_ORIGIN.mk config.mk
	gmake sumb

SGI_N32:
	make dirs
	if [ ! -f "config/config.SGI_N32.mk" ]; then cp "config/defaults/config.SGI_N32.mk" ./config; fi
	(cd externals/SU_MPI && gmake SGI_N32)
	(cd externals/ADT && gmake SGI_N32)
	ln -sf config/config.SGI_N32.mk config.mk
	gmake sumb

SGI_N32_MPICH:
	make dirs
	if [ ! -f "config/config.SGI_N32_MPICH.mk" ]; then cp "config/defaults/config.SGI_N32_MPICH.mk" ./config; fi
	(cd externals/SU_MPI && gmake SGI_N32_MPICH)
	(cd externals/ADT && gmake SGI_N32_MPICH)
	ln -sf config/config.SGI_N32_MPICH.mk config.mk
	gmake sumb

SGI_N32_MPI_ORIGIN:
	make dirs
	if [ ! -f "config/config.SGI_N32_MPI_ORIGIN.mk" ]; then cp "config/defaults/config.SGI_N32_MPI_ORIGIN.mk" ./config; fi
	(cd externals/SU_MPI && gmake SGI_N32_MPI_ORIGIN)
	(cd externals/ADT && gmake SGI_N32_MPI_ORIGIN)
	ln -sf config/config.SGI_N32_MPI_ORIGIN.mk config.mk
	gmake sumb

SICORTEX:
	make dirs
	if [ ! -f "config/config.SICORTEX.mk" ]; then cp "config/defaults/config.SICORTEX.mk" ./config; fi
	(cd externals/SU_MPI && gmake SICORTEX)
	(cd externals/ADT && gmake SICORTEX)
	ln -sf config/config.SICORTEX.mk config.mk
	gmake sumb

SICORTEX_MPI:
	make dirs
	if [ ! -f "config/config.SICORTEX_MPI.mk" ]; then cp "config/defaults/config.SICORTEX_MPI.mk" ./config; fi
	(cd externals/SU_MPI && gmake SICORTEX_MPI)
	(cd externals/ADT && gmake SICORTEX_MPI)
	ln -sf config/config.SICORTEX_MPI.mk config.mk
	gmake sumb
