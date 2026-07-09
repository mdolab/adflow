#!/bin/bash

# Need one argument at least
if [ -z $1 ]; then
    echo "$0 needs to be called with the compiler name as an argument"; exit 1
fi

# Initialize variables
FF90=$1
FFLAGS=""


# Allow argument mismatch for gfortran >= 10

# NOTE: The primary reason for this is that the CGNS library does not provide explicit
# Fortran interface at this time for some of the functions as mentioned
# in the docs https://cgns.github.io/CGNS_docs_current/midlevel/general.html
# and source https://github.com/CGNS/CGNS/blob/develop/src/cgns_f.F90
# Once CGNS lib supports explicit interface this script should be removed
# and any issues with the code addressed.

fc=$("$FF90" --version 2>&1 | grep -i 'gnu')
if [ ! -z "$fc" ]; then
    # Get the GNU compiler version
    version=$("$FF90" -v 2>&1 | grep 'gcc version' | cut -d' ' -f3)
    if [ ! -z "$version" ]; then
      # Get the major version
      version_major=`echo $version | cut -f1 -d.`
    fi

    if [ $version_major -ge 10 ]; then
        FFLAGS="-fallow-argument-mismatch"
    fi
fi

# Print at end to add to the makefile
echo "$FFLAGS"
