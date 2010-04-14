#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 10-17-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to build the ADT library.                *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * The subdirectory where the sources are located.                *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC = src

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
	@echo "                         CYGWIN_INTEL"
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
	@for subdir in $(SUBDIR_SRC) ; \
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
#      * The actual make. This is not a direct target, but it is called *
#      * from the architectures.                                        *
#      *                                                                *
#      ******************************************************************

adt:
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && gmake) || exit 1; \
		done
	(cd lib && gmake)

#      ******************************************************************
#      *                                                                *
#      * Platform specific targets. Note that for the compilation of    *
#      * the ADT library it is not needed to know whether this is a     *
#      * sequential or parallel build. The distinction is present for   *
#      * the build of the entire suite of source codes.                 *
#      *                                                                *
#      ******************************************************************

ABLATION_INTEL_IB:
	ln -sf config/config.ABLATION_INTEL_IB.mk config.mk
	gmake adt

ABLATION_PG_IB:
	ln -sf config/config.ABLATION_PG_IB.mk config.mk
	gmake adt

ALC_INTEL:
	ln -sf config/config.ALC_INTEL.mk config.mk
	gmake adt

ALTIX:
	ln -sf config/config.ALTIX.mk config.mk
	gmake adt

ALTIX_MPI:
	ln -sf config/config.ALTIX.mk config.mk
	gmake adt

ALTIX_MPICH2:
	ln -sf config/config.ALTIX.mk config.mk
	gmake adt

APPLE_MAC_NAG:
	ln -sf config/config.APPLE_MAC_NAG.mk config.mk
	gmake adt

APPLE_MAC_NAG_MPICH:
	ln -sf config/config.APPLE_MAC_NAG.mk config.mk
	gmake adt

APPLE_MAC_XLF:
	ln -sf config/config.APPLE_MAC_XLF.mk config.mk
	gmake adt

APPLE_MAC_XLF_MPICH:
	ln -sf config/config.APPLE_MAC_XLF.mk config.mk
	gmake adt

ASCI_QSC:
	ln -sf config/config.ASCI_QSC.mk config.mk
	gmake adt

CELL_G95:
	ln -sf config/config.CELL_G95.mk config.mk
	gmake adt

FLASH_INTEL:
	ln -sf config/config.FLASH_INTEL.mk config.mk
	gmake adt

FLASH_PG:
	ln -sf config/config.FLASH_PG.mk config.mk
	gmake adt

CYGWIN_INTEL:
	ln -sf config/config.CYGWIN_INTEL.mk config.mk
	gmake adt

GOUGING_INTEL_IB:
	ln -sf config/config.GOUGING_INTEL_IB.mk config.mk
	make adt

MELTING_INTEL:
	ln -sf config/config.MELTING_INTEL.mk config.mk
	make adt

IBM_BLUEGENE:
	ln -sf config/config.IBM_BLUEGENE.mk config.mk
	gmake adt

IBM_DATASTAR:
	ln -sf config/config.IBM_DATASTAR.mk config.mk
	gmake adt

LINUX_ABSOFT:
	ln -sf config/config.LINUX_ABSOFT.mk config.mk
	gmake adt

LINUX_G95:
	ln -sf config/config.LINUX_G95.mk config.mk
	gmake adt

LINUX_G95_MPICH:
	ln -sf config/config.LINUX_G95.mk config.mk
	gmake adt

LINUX_G95_OPENMPI:
	ln -sf config/config.LINUX_G95.mk config.mk
	gmake adt

LINUX_INTEL:
	ln -sf config/config.LINUX_INTEL.mk config.mk
	gmake adt

LINUX_INTEL_MPICH:
	ln -sf config/config.LINUX_INTEL.mk config.mk
	gmake adt

LINUX_INTEL_OPENMPI:
	ln -sf config/config.LINUX_INTEL.mk config.mk
	gmake adt

LINUX_INTEL_OPENMPI_SCINET:
	ln -sf config/config.LINUX_INTEL.mk config.mk
	gmake adt

LINUX_PG:
	ln -sf config/config.LINUX_PG.mk config.mk
	gmake adt

LINUX_PG_MPICH:
	ln -sf config/config.LINUX_PG.mk config.mk
	gmake adt

LINUX_XLF:
	ln -sf config/config.LINUX_XLF.mk config.mk
	gmake adt

REDHOT_IFC_ETHERNET:
	ln -sf config/config.REDHOT_INTEL.mk config.mk
	gmake adt

REDHOT_IFC_MYRINET:
	ln -sf config/config.REDHOT_INTEL.mk config.mk
	gmake adt

REDHOT_PG_ETHERNET:
	ln -sf config/config.REDHOT_PG.mk config.mk
	gmake adt

REDHOT_PG_MYRINET:
	ln -sf config/config.REDHOT_PG.mk config.mk
	gmake adt

REDSTORM:
	ln -sf config/config.REDSTORM.mk config.mk
	gmake adt


SGI:
	ln -sf config/config.SGI.mk config.mk
	gmake adt

SGI_MPI_ORIGIN:
	ln -sf config/config.SGI.mk config.mk
	gmake adt

SGI_N32:
	ln -sf config/config.SGI_N32.mk config.mk
	gmake adt

SGI_N32_MPICH:
	ln -sf config/config.SGI_N32.mk config.mk
	gmake adt

SGI_N32_MPI_ORIGIN:
	ln -sf config/config.SGI_N32.mk config.mk
	gmake adt

SICORTEX:
	ln -sf config/config.SICORTEX.mk config.mk
	gmake adt

SICORTEX_MPI:
	ln -sf config/config.SICORTEX.mk config.mk
	gmake adt