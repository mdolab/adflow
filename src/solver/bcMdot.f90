!
!      ******************************************************************
!      *                                                                *
!      * File:          bcMdot.f90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcMdot(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcMdot applies the duct inflow boundary condition to a block.  *
!      * It is assumed that the pointers in blockPointers are already   *
!      * set to the correct block on the correct grid level.            *
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

         ! Check for duct inflow boundary conditions.

         ductInflow: if(BCType(nn) == mdot) then

           call terminate("bcMdot", "Not implemented yet")

         endif ductInflow
       enddo bocos

       end subroutine bcMdot
