#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Gaetan Kenway              *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 11-16-2014                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src \
		src/utils \
		src/ADT \
		src/bcdata \
		src/partitioning \
		src/modules       \
	        src/solver        \
	        src/metis-4.0     \
		src/NKSolver \
		src/initFlow \
		src/inputParam\
	        src/output        \
	        src/overset       \
	        src/preprocessing \
		src/wallDistance \
	        src/slidingComm   \
	        src/turbulence    \
		src/warping       \
                src/adjoint/ADFirstAidKit \
                src/adjoint/residualInput \
		src/adjoint/outputForward \
		src/adjoint/outputReverse \
                src/adjoint/outputReverseFast \
		src/adjoint 

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
	@for subdir in $(SUBDIR_SRC)  ; \
		do \
			echo "making $@ in $$subdir"; \
			rm -fr $$subdir/*.o; \
		done
	rm -fr src/*.mod
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*

sumb:
	mkdir -p bin;
	mkdir -p obj;
	mkdir -p mod;
	ln -sf config/config.mk config.mk;
	ln -sf SUmb_Common_real.mk SUmb_Common.mk;
	(cd src && make)
