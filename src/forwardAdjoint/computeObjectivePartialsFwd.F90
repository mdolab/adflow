Subroutine computeObjectivePartialsFwd(costFunction)

  !     ******************************************************************
  !     *                                                                *
  !     * Assemble the objective partials from the data precomputed      *
  !     * from the forward mode computations. 
  !     *                                                                *
  !     * costFunction: The index of the cost function to use.           *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointVars
  use ADjointPETSc
  use blockPointers   
  use communication  
  use costFunctions
  use inputTimeSpectral
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: costFunction

  ! Working Variables
  integer(kind=intType) :: i, ierr, sps
  real(kind=realtype) :: val

  !******************************************! 
  !               dIdw                       ! 
  !******************************************! 

  ! Zero Entries and multiply through by costFuncMat
     
  call VecZeroEntries(dJdw, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Write force and moment to dJdw
  do sps=1,nTimeIntervalsSpectral
     do i = 1,6
        call VecAXPY(dJdw, costfuncmat(i, costFunction, sps), FMw(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Assemble dJdw
  call VecAssemblyBegin(dJdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(dJdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !******************************************!
  !          dIdx CALCULATIONS               !
  !******************************************!

  ! Zero Entries and multiply through by costFuncMat
  call VecZeroEntries(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do sps=1,nTimeIntervalsSpectral
     do i = 1,6
        call VecAXPY(dJdx, costfuncmat(i, costFunction, sps), FMx(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do
  ! Assemble dJdx
  call VecAssemblyBegin(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !******************************************!
  !          dIda CALCULATIONS               !
  !******************************************!
  if (nDesignExtra > 0) then
     dIda = zero
     do sps=1,nTimeIntervalsSpectral
        do i = 1,6
           ! This is actually doing a product rule of :
           ! I = sum(i=1,6): FMextra(i) * costFuncVec(i)
           dIda(:) = dIda(:) + dFMdExtra(i,:) * costFuncMat(i, costFunction, sps) + &
                FMextra(i) * dCostFuncMatdExtra(i, costFunction, 1, :)
        end do
     end do
  end if
     
end subroutine computeObjectivePartialsFwd
