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
