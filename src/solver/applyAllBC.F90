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
  use constants
  use blockPointers, only : nDom
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use iteration, only :currentLevel
  use BCRoutines, only : applyAllBC_block
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

        call interpLevelALEBC_block
        call applyAllBC_block(secondHalo)
        call recoverLevelALEBC_block
  
     enddo domains
  enddo spectralLoop

end subroutine applyAllBC
