subroutine getCurrentResidual(rhoRes,totalRRes)

  use constants
  use communication, only : sumb_comm_world 
  use blockPointers, only : nDom, nCellGlobal
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use iteration, only : currentLevel
  use monitor, only: monLoc, monGlob, nMonSum
  use utils, only : setPointers
  implicit none

  ! Compute rhoRes and totalR. The actual residual must have already
  ! been evaluated

  real(kind=realType), intent(out) :: rhoRes,totalRRes
  integer(kind=intType) :: sps,nn,ierr

  monLoc = zero
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, currentLevel, sps)
        call sumResiduals(1, 1) ! Sum 1st state res into first mon location
        call sumAllResiduals(2) ! Sum into second mon location
     end do
  end do
  
  ! This is the same calc as in convergence info, just for rehoRes and
  ! totalR only. 
  call mpi_allreduce(monLoc, monGlob, nMonSum, sumb_real, &
       mpi_sum, SUmb_comm_world, ierr)

  rhoRes = sqrt(monGlob(1)/nCellGlobal(currentLevel))
  totalRRes = sqrt(monGlob(2))

end subroutine getCurrentResidual
