!
!      ******************************************************************
!      *                                                                *
!      * File:          myIsNAN.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-28-2006                                      *
!      * Last modified: 10-20-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       logical function myIsNAN(val)
!
!      ******************************************************************
!      *                                                                *
!      * myIsNAN determines whether or not the given value is a NAN and *
!      * returns the according logical.                                 *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Function arguments.
!
       real(kind=realType), intent(in) :: val
!
!      Local variable.
!
       integer(kind=intType) :: res
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call myIsNaNC(val, res)
       if(res == 1) then
         myIsNAN = .true.
       else
         myIsNAN = .false.
       endif

       end function myIsNAN
