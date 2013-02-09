Subroutine computeObjectivePartialsFwd(costFunction, usedJdw, usedJdx)

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the spatial derivative matrix using a forward mode calc*
  !     * There is one flag to determine how this routine is run:        *
  !     *                                                                *
  !     * useAD: if True, AD is used for derivative calculation, if      *
  !     *        False, FD is used.                                      *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointVars
  use ADjointPETSc
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use costFunctions
  implicit none

  ! Input Variables
  logical, intent(in) :: usedJdw, usedJdx
  integer(kind=intType), intent(in) :: costFunction

  ! Working Variables
  integer(kind=intType) :: i
  !******************************************!!
  !**        dIdw CALCULATIONS            * *!!
  !******************************************!!
  if (usedJdw) then
     ! Zero Entries and multiply through by costFuncVechen
     
     call VecZeroEntries(dJdw, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Write force and moment to dJdw
     do i = 1,6
        call VecAXPY(dJdw, FMw(i), costfuncvec(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do

     ! Assemble dJdw
     call VecAssemblyBegin(dJdw, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call VecAssemblyEnd(dJdw, ierr)
     call EChk(ierr, __FILE__, __LINE__)

  end if

  !******************************************!!
  !**        dIdx CALCULATIONS            * *!!
  !******************************************!!

  if (usedJdx) then
     ! Zero Entries and multiply through by costFuncVec
     call VecZeroEntriex(dJdx, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     do i = 1,6
        call VecAXPY(dJdx, FMx(i), costfuncvec(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do

     ! Assemble dJdx
     call VecAssemblyBegin(dJdx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call VecAssemblyEnd(dJdx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  
  end if

  !******************************************!!
  !**        dIda CALCULATIONS            * *!!
  !******************************************!!

  dIda = zero
  do i = 1,6
     dIda(:) = dFMdExtra(i,:) * costFuncMat(i, costFunction)
  end do

end subroutine computeObjectivePartialsFwd
