!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBC.f90                                  *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine applyAllBC(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllBC applies all boundary conditions for the all         *
!      * blocks on the grid level currentLevel.                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputDiscretization
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo
!
!      Local variables.
!
       integer(kind=intType) :: nn, sps

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not the total energy must be corrected
       ! for the presence of the turbulent kinetic energy.

       if( kPresent ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
           correctForK = .true.
         else
           correctForK = .false.
         endif
       else
         correctForK = .false.
       endif

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Loop over the number of blocks.

         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, currentLevel, sps)

           ! Apply all the boundary conditions. The order is important.

           ! The symmetry boundary conditions.

           call bcSymm(secondHalo)
           call bcSymmPolar(secondHalo)

           ! call bcEulerWall(secondHalo, correctForK)

           ! The viscous wall boundary conditions.

           call bcNSWallAdiabatic( secondHalo, correctForK)
           call bcNSWallIsothermal(secondHalo, correctForK)

           ! The farfield is a special case, because the treatment
           ! differs when preconditioning is used. Make that distinction
           ! and call the appropriate routine.

           select case (precond)

             case (noPrecond)
               call bcFarfield(secondHalo, correctForK)

             case (Turkel)
               call terminate("applyAllBC", &
                              "Farfield boundary conditions for Turkel &
                              &preconditioner not implemented")

             case (ChoiMerkle)
               call terminate("applyAllBC", &
                              "Farfield boundary conditions for Choi and &
                              &Merkle preconditioner not implemented")

           end select

           ! Subsonic outflow and bleed outflow boundaries.

           call bcSubsonicOutflow(secondHalo, correctForK)

           ! Subsonic inflow boundary.

           call bcSubsonicInflow(secondHalo, correctForK)

           ! Bleed inflow regions.

           call bcBleedInflow( secondHalo, correctForK)

           ! Engine boundary conditions. Not implemented yet.

           call bcMdot(secondHalo, correctForK)
           call bcThrust(secondHalo, correctForK)

           ! Extrapolation boundary conditions; this also includes
           ! the supersonic outflow boundary conditions. The difference
           ! between the two is that the extrap boundary conditions
           ! correspond to singular lines and supersonic outflow
           ! boundaries to physical boundaries. The treatment however
           ! is identical.

           call bcExtrap(secondHalo, correctForK)

           ! Inviscid wall boundary conditions.

           call bcEulerWall(secondHalo, correctForK)

           ! Domain-interface boundary conditions,
           ! when coupled with other solvers.

           call bcDomainInterface(secondHalo, correctForK)

           ! Supersonic inflow boundary conditions.

           call bcSupersonicInflow(secondHalo, correctForK)

         enddo domains
       enddo spectralLoop

       end subroutine applyAllBC
