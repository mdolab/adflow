subroutine setupSpatialMatrix(useAD)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Assembles the dRdx Matrix                                      *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none
  !
  !     Local variables.

  logical, intent(in) :: useAD

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj

  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  if( myid ==0 ) &
       write(*,10) "Assembling dRdx Matrix..."

  call cpu_time(time(1))

  call setupSpatialResidualMatrix(drdx,useAD)

  call cpu_time(time(2))
  timeAdjLocal = time(2)-time(1)

  ! MPI Reduce time
  call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
       mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  if(myid ==0) &
       write(*,20) "Assembling dRdx Matrix Time (s) = ", timeAdj
#endif

  ! Output formats.
 
10 format(a)
20 format(a,1x,f8.2)

end subroutine setupSpatialMatrix
