!
!      ******************************************************************
!      *                                                                *
!      * File:          delta.f90                                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-28-2003                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function delta(val1,val2)
!
!      ******************************************************************
!      *                                                                *
!      * delta is a function used to determine the contents of the full *
!      * transformation matrix from the shorthand form. It returns 1    *
!      * if the absolute value of the two arguments are identical.      *
!      * Otherwise it returns 0.                                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Function type.
!
       integer(kind=intType) :: delta
!
!      Function arguments.
!
       integer(kind=intType) :: val1, val2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if(abs(val1) == abs(val2)) then
         delta = 1_intType
       else
         delta = 0_intType
       endif

       end function delta
