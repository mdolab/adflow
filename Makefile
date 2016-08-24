#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Gaetan Kenway              *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 11-16-2014                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src/modules       \
	        src/solver        \
		src/ADT           \
		src/bcdata        \
	        src/initFlow      \
	        src/inputParam    \
	        src/metis-4.0     \
	        src/output        \
	        src/overset       \
	        src/partitioning  \
	        src/preprocessing \
	        src/slidingComm   \
	        src/stabilityDerivatives \
	        src/turbulence    \
	        src/wallDistance  \
		src/warping       \
		src/bendingMomentAnalysis \
                src/adjoint/ADFirstAidKit \
                src/adjoint/residualInput \
		src/NKSolver      \
		src/ANKSolver     \
		src/adjoint/outputForward \
		src/adjoint/outputReverse \
                src/adjoint/outputReverseFast \
		src/adjoint \

SUMB_SUBDIRS       = $(SUBDIR_SRC)
SUMB_CLEAN_SUBDIRS = $(SUBDIR_SRC)

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_INTEL_OPENMPI.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. Typically the CGNS directory "; \
	echo "will have to be modified. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make sumb;\
	fi;

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

sumb:
	mkdir -p bin;
	mkdir -p obj;
	mkdir -p mod;
	ln -sf config/config.mk config.mk;
	ln -sf SUmb_Common_real.mk SUmb_Common.mk;

	@for subdir in $(SUMB_SUBDIRS) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)
	(cd src/python/f2py && make)
