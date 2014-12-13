subroutine setupAllResidualMatricesfwd
  use ADjointPETSc
  use ADjointVars 
  use communication    
  use inputADjoint       

  implicit none
 
  logical :: useAD, useTranspose, usePC, useObjective
  real(kind=realType) :: timeAdjLocal, timeAdj, time(2)
  integer(kind=intType) :: ierr
  if( myid ==0 ) &
       write(*, 10) "Assembling All Residual Matrices in Forward mode..."

  call cpu_time(time(1))

  ! If we are assembling matrices...we ned to assemble the
  ! 'transpose', with 'AD', we want the exact matrix not the 'PC',
  ! and will compute objective RHS
  useAD = .True.
  usePC = .False.
  useTranspose = .True.
  useObjective = .True.

  if (.not. useMatrixFreedRdw) then 
     call setupStateResidualMatrix(drdwT, useAD, usePC, useTranspose, &
          useObjective, 1_intType)
  end if

  if (.not. useMatrixFreedRdx) then 
     call setupSpatialResidualMatrix(drdx, useAD, useObjective)
     call setupExtraResidualMatrix(drda, useAD)
  end if

  call cpu_time(time(2))
  timeAdjLocal = time(2)-time(1)

  call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
       mpi_max, 0, SUMB_COMM_WORLD, ierr)
  call EChk(ierr,  __FILE__, __LINE__)
  if(myid ==0) &
       write(*, 20) "Assembling All Residaul Matrices Fwd time (s) = ", timeAdj

  ! Output formats.
10 format(a)
20 format(a, 1x, f8.2)

end subroutine setupAllResidualMatricesfwd
