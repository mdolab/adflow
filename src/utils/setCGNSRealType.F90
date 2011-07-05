!
!      ******************************************************************
!      *                                                                *
!      * File:          setCGNSRealType.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-12-2003                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       integer function setCGNSRealType()
!
!      ******************************************************************
!      *                                                                *
!      * setCGNSRealType sets the cgns real type, depending on the      *
!      * compiler options. Note that quadrupole precision is not        *
!      * supported by CGNS; double precision is used instead for the    *
!      * CGNS IO.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use su_cgns
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("setCGNSRealType", &
                      "Function should not be called if no cgns support &
                      &is selected.")

#else

# ifdef USE_SINGLE_PRECISION
       setCGNSRealType = RealSingle
# else
       setCGNSRealType = RealDouble
# endif

#endif

       end function setCGNSRealType
