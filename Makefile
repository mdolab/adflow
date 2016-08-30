#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Gaetan Kenway              *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 11-16-2014                                      *
#      *                                                                *
#      ******************************************************************

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
	rm -fr src/build/*.mod
	rm -fr src/build/*.o
	rm -f *~ config.mk;


sumb:
	ln -sf config/config.mk config.mk;
	(cd src/build/ && make)

