#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 02-21-2006                                      *
#      * Last modified: 10-17-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to build the su_mpi library.             *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * The subdirectory where the sources are located.                *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC = src/modules
SUBDIR_SRC_SEQ = $(SUBDIR_SRC) src/mpi

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
	@echo "                         APPLE_MAC_NAG"
	@echo "                         APPLE_MAC_NAG_MPICH"
	@echo "                         APPLE_MAC_XLF"
	@echo "                         APPLE_MAC_XLF_MPICH"
	@echo "                         ASCI_QSC"
	@echo "                         CELL_G95"
	@echo "                         FLASH_INTEL"
	@echo "                         FLASH_PG"
	@echo "                         GOUGING"
	@echo "                         MELTING_INTEL"
	@echo "                         IBM_BLUEGENE"
	@echo "                         IBM_DATASTAR"
	@echo "                         LINUX_ABSOFT"
	@echo "                         LINUX_G95"
	@echo "                         LINUX_G95_MPICH"
	@echo "                         LINUX_INTEL"
	@echo "                         LINUX_INTEL_MPICH"
	@echo "                         LINUX_PG"
	@echo "                         LINUX_PG_MPICH"
	@echo "                         LINUX_XLF"
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

clean:
	@echo " Making clean ... "
	@for subdir in $(SUBDIR_SRC_SEQ) ; \
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
#      * Parallel and sequential makes. These are not direct targets,   *
#      * but are called from the architectures.                         *
#      *                                                                *
#      ******************************************************************

parallel:
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)

sequential:
	@for subdir in $(SUBDIR_SRC_SEQ) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)

#      ******************************************************************
#      *                                                                *
#      * Platform specific targets.                                     *
#      *                                                                *
#      ******************************************************************

ABLATION_INTEL_IB:
	ln -sf config/config.ABLATION_INTEL_IB.mk config.mk
	make parallel

ABLATION_PG_IB:
	ln -sf config/config.ABLATION_PG_IB.mk config.mk
	make parallel

ALC_INTEL:
	ln -sf config/config.ALC_INTEL.mk config.mk
	make parallel

ALTIX:
	ln -sf config/config.ALTIX.mk config.mk
	make sequential

ALTIX_MPI:
	ln -sf config/config.ALTIX_MPI.mk config.mk
	make parallel

ALTIX_MPICH2:
	ln -sf config/config.ALTIX_MPICH2.mk config.mk
	make parallel

APPLE_MAC_NAG:
	ln -sf config/config.APPLE_MAC_NAG.mk config.mk
	make sequential

APPLE_MAC_NAG_MPICH:
	ln -sf config/config.APPLE_MAC_NAG_MPICH.mk config.mk
	make parallel

APPLE_MAC_XLF:
	ln -sf config/config.APPLE_MAC_XLF.mk config.mk
	make sequential

APPLE_MAC_XLF_MPICH:
	ln -sf config/config.APPLE_MAC_XLF_MPICH.mk config.mk
	make parallel

ASCI_QSC:
	ln -sf config/config.ASCI_QSC.mk config.mk
	make parallel

CELL_G95:
	ln -sf config/config.CELL_G95.mk config.mk
	make sequential

FLASH_INTEL:
	ln -sf config/config.FLASH_INTEL.mk config.mk
	make parallel

FLASH_PG:
	ln -sf config/config.FLASH_PG.mk config.mk
	make parallel

GOUGING_INTEL_IB:
	ln -sf config/config.GOUGING_INTEL_IB.mk config.mk
	make parallel

MELTING_INTEL:
	ln -sf config/config.MELTING_INTEL.mk config.mk
	make parallel

IBM_BLUEGENE:
	ln -sf config/config.IBM_BLUEGENE.mk config.mk
	make parallel

IBM_DATASTAR:
	ln -sf config/config.IBM_DATASTAR.mk config.mk
	make parallel

LINUX_ABSOFT:
	ln -sf config/config.LINUX_ABSOFT.mk config.mk
	make sequential

LINUX_G95:
	ln -sf config/config.LINUX_G95.mk config.mk
	make sequential

LINUX_G95_MPICH:
	ln -sf config/config.LINUX_G95_MPICH.mk config.mk
	make parallel

LINUX_G95_OPENMPI:
	ln -sf config/config.LINUX_G95_OPENMPI.mk config.mk
	make parallel

LINUX_INTEL:
	ln -sf config/config.LINUX_INTEL.mk config.mk
	make sequential

LINUX_INTEL_MPICH:
	ln -sf config/config.LINUX_INTEL_MPICH.mk config.mk
	make parallel

LINUX_INTEL_OPENMPI:
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	make parallel

LINUX_INTEL_OPENMPI_SCINET:
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	make parallel

LINUX_GFORTRAN:
	ln -sf config/config.LINUX_GFORTRAN.mk config.mk
	make parallel

LINUX_PG:
	ln -sf config/config.LINUX_PG.mk config.mk
	make sequential

LINUX_PG_MPICH:
	ln -sf config/config.LINUX_PG_MPICH.mk config.mk
	make parallel

LINUX_XLF:
	ln -sf config/config.LINUX_XLF.mk config.mk
	make sequential

REDHOT_IFC_ETHERNET:
	ln -sf config/config.REDHOT_IFC_ETHERNET.mk config.mk
	make parallel

REDHOT_IFC_MYRINET:
	ln -sf config/config.REDHOT_IFC_MYRINET.mk config.mk
	make parallel

REDHOT_PG_ETHERNET:
	ln -sf config/config.REDHOT_PG_ETHERNET.mk config.mk
	make parallel

REDHOT_PG_MYRINET:
	ln -sf config/config.REDHOT_PG_MYRINET.mk config.mk
	make parallel

REDSTORM:
	ln -sf config/config.REDSTORM.mk config.mk
	make parallel

SGI:
	ln -sf config/config.SGI.mk config.mk
	make sequential

SGI_MPI_ORIGIN:
	ln -sf config/config.SGI_MPI_ORIGIN.mk config.mk
	make parallel

SGI_N32:
	ln -sf config/config.SGI_N32.mk config.mk
	make sequential

SGI_N32_MPICH:
	ln -sf config/config.SGI_N32_MPICH.mk config.mk
	make parallel

SGI_N32_MPI_ORIGIN:
	ln -sf config/config.SGI_N32_MPI_ORIGIN.mk config.mk
	make parallel

SICORTEX:
	ln -sf config/config.SICORTEX.mk config.mk
	make sequential

SICORTEX_MPI:
	ln -sf config/config.SICORTEX_MPI.mk config.mk
	make parallel