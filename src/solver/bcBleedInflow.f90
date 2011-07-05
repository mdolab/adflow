!
!      ******************************************************************
!      *                                                                *
!      * File:          bcBleedInflow.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-16-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcBleedInflow(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcBleedInflow applies the boundary conditions to an inflow     *
!      * bleed region of a block.                                       *
!      * It is assumed that the pointers in blockPointers are already   *
!      * set to the correct block on the correct grid level.            *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use bleedFlows
       use blockPointers
       use flowVarRefState
       use iteration
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

         ! Check for an inflow bleed region.

         inflowBleed: if(BCType(nn) == MassBleedInflow) then

           call terminate("bcBleedInflow", "Not implemented yet")

         endif inflowBleed
       enddo bocos

       end subroutine bcBleedInflow
