subroutine getCurrentResidual(rhoRes,totalRRes)

  use communication
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use iteration
  use inputPhysics
  use inputIteration
  use monitor
  implicit none

  ! Compute rhoRes and totalR. The actual residual must have already
  ! been evaluated

  real(kind=realType), intent(out) :: rhoRes,totalRRes
  integer(kind=intType) :: sps,nn,ierr

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
