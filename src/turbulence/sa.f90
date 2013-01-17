subroutine sa(resOnly)
  !
  ! Shell function to call sa_block on al blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  logical, intent(in) :: resOnly
  !      Local variables.
  !
  integer(kind=intType) :: nn, sps
  !
  ! Compute the time derivative for the time spectral mode.
  !
  call unsteadyTurbSpectral(itu1,itu1)

  ! Loop over the number of spectral solutions.

  spectralLoop: do sps=1,nTimeIntervalsSpectral
     
     ! Loop over the number of blocks.

     domains: do nn=1,nDom

        ! Set the pointers for this block.

        call setPointers(nn, currentLevel, sps)

        call sa_block(resOnly)

     end do domains

  end do spectralLoop

end subroutine sa

!
!      ******************************************************************
!      *                                                                *
!      * File:          sa.f90                                          *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine sa_block(resOnly)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * sa solves the transport equation for the Spalart-Allmaras      *
  !      * turbulence model in a segregated manner using a diagonal       *
  !      * dominant ADI-scheme.                                           *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Subroutine argument.
  !
  logical, intent(in) :: resOnly
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, sps
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
 
  ! Set the arrays for the boundary condition treatment.
  
  call bcTurbTreatment
  
  ! Solve the transport equation for nuTilde.
  
  call saSolve(resOnly)
  
  ! The eddy viscosity and the boundary conditions are only
  ! applied if an actual update has been computed in saSolve.
  
  if(.not. resOnly ) then
     
     ! Compute the corresponding eddy viscosity.
     
     call saEddyViscosity
     
     ! Set the halo values for the turbulent variables.
     ! We are on the finest mesh, so the second layer of halo
     ! cells must be computed as well.
     
     call applyAllTurbBCThisBlock(.true.)
     
     ! Write the loglaw for a flat plate boundary layer.
     
     ! call writeLoglaw
     
  endif

end subroutine sa_block
