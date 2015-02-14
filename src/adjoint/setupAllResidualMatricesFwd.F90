subroutine setupAllResidualMatricesfwd
  use ADjointPETSc
  use ADjointVars 
  use communication    
  use inputADjoint       

  implicit none
 
  logical :: useAD, useTranspose, usePC, useObjective
  real(kind=realType) :: timeAdjLocal, timeAdj, time(2)
  integer(kind=intType) :: ierr
#ifndef USE_COMPLEX

  ! If we are assembling matrices...we ned to assemble the
  ! 'transpose', with 'AD', we want the exact matrix not the 'PC',
  ! and will compute objective RHS
  useAD = .True.
  usePC = .False.
  useTranspose = .True.
  useObjective = .True.

  if (.not. useMatrixFreedRdw) then 
     if( myid ==0 ) then
        write(*, 10) "Assembling State Residual Matrix in Forward mode..."
     end if
     time(1) = mpi_wtime()
     call setupStateResidualMatrix(drdwT, useAD, usePC, useTranspose, &
          useObjective, frozenTurbulence, 1_intType)
     time(2) = mpi_wtime()
     timeAdjLocal = time(2)-time(1)

     call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
          mpi_max, 0, SUMB_COMM_WORLD, ierr)
     call EChk(ierr,  __FILE__, __LINE__)
     
     if(myid ==0)  then 
       write(*, 20) "Assembling State Residaul Matrices Fwd time (s) = ", timeAdj
    end if
  end if

  ! Output formats.
10 format(a)
20 format(a, 1x, f8.2)
#endif
end subroutine setupAllResidualMatricesfwd
