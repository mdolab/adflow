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
		src/parallelIO    \
	        src/partitioning  \
	        src/preprocessing \
	        src/slidingComm   \
		src/NKsolver      \
	        src/stabilityDerivatives \
	        src/turbulence    \
	        src/utils         \
	        src/wallDistance  \
		src/warping       \
		src/bendingMomentAnalysis \
		src/adjoint               \
                src/adjoint/ADFirstAidKit \
                src/adjoint/forcesOutput\
                src/adjoint/forcesInput\
                src/forwardAdjoint/outputReverse \
                src/forwardAdjoint/outputReverseFast \
                src/forwardAdjoint \
                src/forwardAdjoint/residualInput \
                src/forwardAdjoint/outputForward \
		src/forwardAdjoint/costInput \
		src/forwardAdjoint/outputCost\

SUBDIR_EXEC   = src/exec
SUBDIR_PV3    = src/pv3Interface
SUMB_SUBDIRS       = $(SUBDIR_SRC) $(PV3_INT_SRC_DIR)
SUMB_CLEAN_SUBDIRS = $(SUBDIR_SRC)  $(SUBDIR_PV3) $(SUBDIR_EXEC)

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
	(cd $(SUBDIR_EXEC) && make)
	(cd src/python/f2py && make)
