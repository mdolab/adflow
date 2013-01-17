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
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  logical, intent(in) :: secondHalo
  !
  !      Local Variables
  integer(kind=intType) :: sps, nn
  
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
 
  ! Loop over the number of spectral solutions.

  spectralLoop: do sps=1,nTimeIntervalsSpectral

     ! Loop over the number of blocks.

     domains: do nn=1,nDom

        ! Set the pointers for this block.

        call setPointers(nn, currentLevel, sps)

        call applyAllBC_block(secondHalo)
  
     enddo domains
  enddo spectralLoop

end subroutine applyAllBC

subroutine applyAllBC_block(secondHalo)

  ! Apply BC's for a single block

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
  
  ! Apply all the boundary conditions. The order is important.

  ! The symmetry boundary conditions.
  
  call bcSymm(secondHalo)
#ifndef USE_TAPENADE
  call bcSymmPolar(secondHalo)
#endif
  !call bcEulerWall(secondHalo, correctForK)

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
          "Farfield Turkel preconditioner not implemented")
     
  case (ChoiMerkle)
     call terminate("applyAllBC", &
          "Farfield Choi and Merkle preconditioner not implemented")
     
  end select
  
  ! Subsonic outflow and bleed outflow boundaries.
#ifndef USE_TAPENADE  
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
#endif  
  
  ! Inviscid wall boundary conditions.

  call bcEulerWall(secondHalo, correctForK)
  
  ! Domain-interface boundary conditions,
  ! when coupled with other solvers.
#ifndef USE_TAPENADE  
  call bcDomainInterface(secondHalo, correctForK)
  
  ! Supersonic inflow boundary conditions.
  
  call bcSupersonicInflow(secondHalo, correctForK)
#endif  
end subroutine applyAllBC_block
