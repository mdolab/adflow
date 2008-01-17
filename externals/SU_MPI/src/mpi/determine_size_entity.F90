!
!      ******************************************************************
!      *                                                                *
!      * File:          determine_size_entity.F90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 09-27-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       integer function determine_size_entity(datatype)
!
!      ******************************************************************
!      *                                                                *
!      * The function determine_size_entity determines the size of an   *
!      * entity of the type corresponding to datatype.                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Function arguments
!
       integer :: datatype
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the size of each individual entity. Take the situation
       ! of a Cray into account in case of an integer; an integer is 4
       ! bytes on all the machines I know, except cray where it is 8.

       select case (datatype)
         case (mpi_character, mpi_byte)
           determine_size_entity = 1
         case (mpi_integer)
#ifdef _CRAY
           determine_size_entity = 8
#else
           determine_size_entity = 4
#endif
         case (mpi_integer2)
           determine_size_entity = 2
         case (mpi_integer4)
           determine_size_entity = 4
         case (mpi_integer8)
           determine_size_entity = 8
         case (mpi_real)
           determine_size_entity = 4
         case (mpi_double_precision)
           determine_size_entity = 8
         case (mpi_complex)
           determine_size_entity = 8
         case (mpi_double_complex)
           determine_size_entity = 16
         case (mpi_real4)
           determine_size_entity = 4
         case (mpi_real8)
           determine_size_entity = 8
         case (mpi_real16)
           determine_size_entity = 16
         case default
           call su_mpi_terminate("determine_size_entity", &
                                 "Unknown dataype")
       end select

       end function determine_size_entity
