module precision
  !
  !       Definition of the kinds used for the integer and real types.
  !       Due to MPI, it is a bit messy to use the compiler options -r8
  !       and -r4 and therefore the kind construction is used here,
  !       where the precision is set using compiler flags of -d type.
  !       This is the only file of the code that should be changed when
  !       a user wants single precision instead of double precision. All
  !       other routines use the definitions in this file whenever
  !       possible. If other definitions are used, there is a good
  !       reason to do so, e.g. when calling the cgns or MPI functions.
  !       The actual types used are determined by compiler flags like
  !       -DUSE_LONG_INT and -DUSE_SINGLE_PRECISION. If these are
  !       omitted the default integer and double precision are used.
  !
  !

  use complexify
  use mpi
  implicit none
  save

  !
  !       Definition of the integer type used in the entire code. There
  !       might be a more elegant solution to do this, but be sure that
  !       compatability with MPI must be guaranteed. Note that dummyInt
  !       is a private variable, only used for the definition of the
  !       integer type. Note furthermore that the parameters defining
  !       the MPI types are integers. This is because of the definition
  !       in MPI.
  !

#ifdef USE_LONG_INT

  ! Long, i.e. 8 byte, integers are used as default integers

  integer(kind=8), private :: dummyInt
  integer, parameter       :: adflow_integer  = mpi_integer8
  integer, parameter       :: sizeOfInteger = 8
#else

  ! Standard 4 byte integer types are used as default integers.

  integer(kind=4), private :: dummyInt
  integer, parameter       :: adflow_integer  = mpi_integer4
  integer, parameter       :: sizeOfInteger = 4
#endif

  !
  !       Definition of the float type used in the entire code. The
  !       remarks mentioned before the integer type definition also
  !       apply here.
  !

#ifdef USE_SINGLE_PRECISION

  ! Single precision reals are used as default real types.

  complex(kind=4), private :: dummyReal
  integer, parameter    :: adflow_real  = MPI_DOUBLE_COMPLEX
  integer, parameter    :: sizeOfReal = 4
  real(kind=4), private :: dummyCGNSReal

#elif USE_QUADRUPLE_PRECISION

  ! Quadrupole precision reals are used as default real types.
  ! This may not be supported on all platforms.
  ! As cgns does not support quadrupole precision, double
  ! precision is used instead.

  complex(kind=16), private :: dummyReal
  integer, parameter     :: adflow_real  = mpi_DOUBLE_COMPLE16
  integer, parameter     :: sizeOfReal = 16
  real(kind=8), private :: dummyCGNSReal

#else

  ! Double precision reals are used as default real types.

  complex(kind=8), private :: dummyReal
  integer, parameter    :: adflow_real   = MPI_DOUBLE_COMPLEX
  integer, parameter    :: sizeOfReal = 8
  real(kind=8), private :: dummyCGNSReal

#endif

 ! Dummy single and double types
  complex(kind=4) :: dummySingle
  complex(kind=8) :: dummyDouble

  !
  !       Definition of the porosity type. As this is only a flag to
  !       indicate whether or not fluxes must be computed, an integer1
  !       is perfectly okay.
  !
  integer(kind=1), private :: dummyPor

  !      Definition of the integer type for the element types. As only
  !      a limited number element types are present, a 1 byte integer
  !      is enough.
  !
  integer(kind=1), private :: adtDummyElementInt

  !       Definition of the cgns periodic type.
  !
  real(kind=4), private :: dummyCGNSPer
  !
  !       Definition of the kind parameters for the integer and real
  !       types.
  !
  integer, parameter :: intType      = kind(dummyInt)
  integer, parameter :: porType      = kind(dummyPor)
  integer, parameter :: realType     = kind(dummyReal)
  integer, parameter :: adtElementType = kind(adtDummyElementInt)
  integer, parameter :: cgnsRealType = kind(dummyCGNSReal)
  integer, parameter :: cgnsPerType  = kind(dummyCGNSPer)
  integer, parameter :: alwaysRealType = kind(dummyReal)
  integer, parameter :: singleType   = kind(dummySingle)
  integer, parameter :: doubleType   = kind(dummyDouble)

  !
  !       Set the parameter debug, depending on the compiler option.
  !
#ifdef DEBUG_MODE
  logical, parameter :: debug = .true.
#else
  logical, parameter :: debug = .false.
#endif

end module precision
