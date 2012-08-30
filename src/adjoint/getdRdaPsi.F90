!
!     ******************************************************************
!     *                                                                *
!     * File:          getdRdaPsi.F90                                 *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-08-2010                                      *
!     * Last modified: 10-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getdRdaPsi(output,ndv,adjoint,nstate)
#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRda                    *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc
  use ADjointVars
  use blockPointers
  use inputADjoint
  use section
  use inputTimeSpectral 
  use monitor 

  implicit none
  integer(kind=intType), intent(in) :: ndv, nstate
  real(kind=realType), intent(in) :: adjoint(nstate)
  real(kind=realType),intent(out) :: output(ndv)
  integer(kind=intType) :: ierr, i

  ! put adjoint arry in w_like1
  call VecPlaceArray(w_like1, adjoint, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the result vector for dRda^T * psi
  call VecCreateMPI(SUMB_PETSC_COMM_WORLD, PETSC_DECIDE, ndv, dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do the Multiplication
  call MatMultTranspose(dRda, w_like1, dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! This is a little wonkly, since dRdaTPsi only contains a handful of
  ! variables.
  ! Use a proper vec scatter here:
  call VecScatterCreateToAll(dRdaTPsi, dRdaTpsi_scatter, dRdaTPsi_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterBegin(dRdaTpsi_scatter, dRdaTPsi, dRdaTpsi_local, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd (dRdaTpsi_scatter, dRdaTPsi, dRdaTpsi_local, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now just pluck off the local values
  do i=1, nDv
     call VecGetValues(dRdaTPsi_local, 1, i-1, output(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

  ! No longer need vectors/scatter
  call VecScatterDestroy(dRdaTpsi_scatter, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(dRdaTpsi_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(w_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif
end subroutine getdRdaPsi
