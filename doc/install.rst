.. _adflow_install:

Installation
============
All the core computations in ``ADflow`` are coded in Fortran.
It is therefore necessary to build this library before using ``ADflow``.

Requirements
------------
ADflow requires the following dependencies:
- CGNS Library
- PETSc
- MPI

See the MDO Lab installation guide `here <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/installInstructions/install3rdPartyPackages.html>`_ for the supported versions and installation instructions.

Building
--------
ADflow follows the standard MDO Lab build procedure.
To start, find a configuration file close to your current setup in::

    $ config/defaults

and copy it to ''config/config.mk''. For example::

    $ cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

If you are a beginner user installing the packages on a linux desktop, 
you should use the ``config.LINUX_GFORTRAN`` versions of the configuration 
files. The ``config.LINUX_INTEL`` versions are usually used on clusters.
ADflow has been successfully compiled on LINUX with either
``ifort`` or ``gfortran``.

Once you have copied the config file, compile ADflow by running::

    $ make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module adflow can be imported...
   Module adflow was successfully imported.

If you don't see this, it will be necessary to configure the build manually.
To configure manually, open ``config/config.mk`` and modify options as necessary.

It is most likely that you need to modify the ``CGNS_INCLUDE_FLAGS`` and the ``CGNS_LINKER_FLAGS`` variables.
After changes to the configuration file, run ``make clean`` before attempting a new build.

.. NOTE::

    Compiling ADflow on HPC clusters requires additional care, as some systems do not have a homogeneous CPU architecture across all nodes.
    For example, the architecture of the login nodes may differ from the architecture of compute nodes available on the same cluster. 

    Compiling the code on/for a specific login node type may result in unexpected crashes if the compute nodes have an incompatible (newer) architecture.
    You can append the ``-march=<HPC-ARCH>`` flag to the ``config.mk`` file to specify the architecture of the compute node and avoid such issues. 
    To optimize the compiled code for a specific architecture, one can add the ``-mtune=<HPC-ARCH>`` flag. However, this is rarely needed. 
    An example of the updated flags in the config file is:: 

        FF90_FLAGS = <normal-flags> -march=<HPC-ARCH> -mtune=<HPC-ARCH>

    Note that the ``<HPC-ARCH>`` should be replaced with the name of the correct compute node architecture. This name is compiler-specific and also depends on the compiler version.

    For example, if the login nodes use newer Cascade Lake CPUs while the compute nodes are based on older Sandy Bridge CPUs it may be necessary to set ``-march=sandybridge`` in your ``config.mk`` file::

        FF90_FLAGS = <normal-flags> -march=sandybridge

    We recommend to contact your local HPC team to get more information about hardware specific issues. 

Lastly, to build and install the Python interface, type::

    pip install .


Verification
------------
ADflow contains a set of simple tests that can be run automatically to ensure ADflow reproduces the expected reference results.
First, install `idwarp <https://github.com/mdolab/idwarp/>`__ following the instructions there.
Next, install the testing dependencies by going to the root ADflow directory and typing::

    $ pip install .[testing]

With all of these packages installed, you can fully verify your ADflow installation.
First, run the script::

    $ inputFiles/get-input-files.sh

to download and extract the necessary files.
Then in the root directory run::

    $ testflo .


Complex Build
-------------
ADflow contains scripts to automatically build a "complexified"
version of ADflow directly from the real version.

ADflow_CS REQUIRES a complex build of petsc to build and run. The
petsc configuration script must be re-run with the following
options::

    $ ./configure --with-shared-libraries --download-superlu_dist=yes --download-parmetis=yes --download-metis=yes --with-fortran-interfaces=1 --with-debugging=yes --with-scalar-type=complex --PETSC_ARCH=complex-debug

Follow instructions as before to complete complex build.

Now, to build complex ADflow do::

    $ export PETSC_ARCH=complex-debug
    $ make -f Makefile_CS

Note that the correct, complex PETSC_ARCH MUST be set before the code is
compiled and also must be set when running in complex mode.

To run the complex tests, first set the ``$PETSC_ARCH`` to the complex architecture.
Then run::

    $ testflo . -m "cmplx_test*"
