!
!      ******************************************************************
!      *                                                                *
!      * File:          bcThrust.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcThrust(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcThrust applies the duct outflow boundary condition to a      *
!      * block. It is assumed that the pointers in blockPointers are    *
!      * already set to the correct block on the correct grid level.    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo, correctForK
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for duct outflow boundary conditions.

         ductOutflow: if(BCType(nn) == thrust) then

           call terminate("bcThrust", "Not implemented yet")

         endif ductOutflow
       enddo bocos

       end subroutine bcThrust
