!
!     ******************************************************************
!     *                                                                *
!     * File:          adtPrecision.F90                                *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 11-27-2004                                      *
!     * Last modified: 06-12-2005                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtPrecision
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the kinds used for the integer and real types.   *
!     * Due to mpi, it is a bit messy to use the compiler options -r8  *
!     * and -r4 (and it is depreciated as well). Therefore the kind    *
!     * construction is used here, where the precision is set using    *
!     * compiler flags of -D type.                                     *
!     *                                                                *
!     * This is the only source file that should be changed when a     *
!     * user wants single precision instead of double precision. All   *
!     * other routines use the definitions in this file whenever       *
!     * possible. If other definitions are used, there is a good       *
!     * reason to do so, e.g. when calling MPI functions.              *
!     *                                                                *
!     * The actual types used are determined by compiler flags like    *
!     * -DUSE_LONG_INT and -DUSE_SINGLE_PRECISION. If these are        *
!     * omitted the default integer and double precision are used.     *
!     *                                                                *
!     ******************************************************************
!
!     ******************************************************************
!     *                                                                *
!     * Include the su_mpi module; inside this module it is            *
!     * controlled whether a sequential or a parallel executable is    *
!     * built. For the ADT sources this is completely transparent,     *
!     * although a completely new build must be performed when a       *
!     * change is made from sequential to parallel and vice versa.     *
!     *                                                                *
!     ******************************************************************
!
      use su_mpi
      implicit none
      save
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the integer type. There might be a more elegant  *
!     * solution to do this, but be sure that compatability with MPI   *
!     * must be guaranteed. Note that adtDummyInt is a private         *
!     * variable, only used for the definition of the integer type.    *
!     * Note furthermore that the parameters defining the MPI types    *
!     * are integers. This is because of the definition in MPI.        *
!     *                                                                *
!     ******************************************************************
!

#ifdef USE_LONG_INT

      ! Long, i.e. 8 byte, integers are used as default integers.

      integer(kind=8), private :: adtDummyInt
      integer, parameter       :: adt_integer = mpi_integer8
#else

      ! Standard 4 byte integer types are used as default integers.

      integer(kind=4), private :: adtDummyInt
      integer, parameter       :: adt_integer = mpi_integer4

#endif

!
!     ******************************************************************
!     *                                                                *
!     * Definition of the float type used in the entire code. The      *
!     * remarks mentioned before the integer type definition also      *
!     * apply here.                                                    *
!     *                                                                *
!     ******************************************************************
!

#ifdef USE_SINGLE_PRECISION

      ! Single precision reals are used as default real types.

      real(kind=4), private :: adtDummyReal
      integer, parameter    :: adt_real = mpi_real4

#elif USE_QUADRUPOLE_PRECISION

      ! Quadrupole precision reals are used as default real types.
      ! This may not be supported on all platforms.

      real(kind=16), private :: adtDummyReal
      integer, parameter     :: adt_real = mpi_real16

#else

      ! Double precision reals are used as default real types.

      real(kind=8), private :: adtDummyReal
      integer, parameter    :: adt_real = mpi_real8

#endif

!
!     ******************************************************************
!     *                                                                *
!     * Definition of the integer type for the element types. As only  *
!     * a limited number element types are present, a 1 byte integer   *
!     * is enough.                                                     *
!     *                                                                *
!     ******************************************************************
!
      integer(kind=1), private :: adtDummyElementInt
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the kind parameters for the integer and real     *
!     * types.                                                         *
!     *                                                                *
!     ******************************************************************
!
      integer, parameter :: adtIntType     = kind(adtDummyInt)
      integer, parameter :: adtElementType = kind(adtDummyElementInt)
      integer, parameter :: adtRealType    = kind(adtDummyReal)
!
!     ******************************************************************
!     *                                                                *
!     * Set the parameter adtDebug, depending on the compiler option.  *
!     *                                                                *
!     ******************************************************************
!

#ifdef DEBUG_MODE
      logical, parameter :: adtDebug = .true.
#else
      logical, parameter :: adtDebug = .false.
#endif

      end module adtPrecision
