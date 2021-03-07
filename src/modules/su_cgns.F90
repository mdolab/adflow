       module su_cgns
!
!       Module that contains the definition of the cgns parameters.
!       Depending on the compiler flags either the cgns module is
!       included or the functionality is faked by just defining the
!       parameters.
!

#ifdef USE_TAPENADE
!      ******************************************************************
!      *                                                                *
!      * Dimensional units                                              *
!      *                                                                *
!      ******************************************************************
!

       implicit none
       save

       integer, parameter :: Null = 0
       integer, parameter :: UserDefined = 1

       integer, parameter :: Kilogram  = 2
       integer, parameter :: Gram      = 3
       integer, parameter :: Slug      = 4
       integer, parameter :: PoundMass = 5

       integer, parameter :: Meter      = 2
       integer, parameter :: Centimeter = 3
       integer, parameter :: Millimeter = 4
       integer, parameter :: Foot       = 5
       integer, parameter :: Inch       = 6
       integer, parameter :: Second = 2

       integer, parameter :: Kelvin     = 2
       integer, parameter :: Celcius    = 3
       integer, parameter :: Celsius    = 3
       integer, parameter :: Rankine    = 4
       integer, parameter :: Fahrenheit = 5

       integer, parameter :: Degree = 2
       integer, parameter :: Radian = 3
#else

       use cgns
       implicit none

#endif


       end module su_cgns
