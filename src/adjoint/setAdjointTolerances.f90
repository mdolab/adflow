!
!     ******************************************************************
!     *                                                                *
!     * File:          setAdjointTolerances.f90                        *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-13-2010                                      *
!     * Last modified: 10-13-2010                                      *
!     *                                                                *
!     ******************************************************************

subroutine setAdjointTolerances(L2,L2Rel,L2Abs)

  use ADjointPETSc
  use inputADjoint
  implicit none

  ! Subroutine arguments:
  real(kind=realType) :: L2,L2Rel,L2Abs

  ! Local Arguments
  real(kind=realType) :: curRes

  ! Get Current Residual
  call MatMult(dRdWT,psi,pvr,PETScIerr)
  call VecAXPY(pvr,PETScNegOne,dJdW,PETScIerr)
  call VecNorm(pvr,NORM_2,curRes,PETScIerr)  

  ! We are only going to overwrite adjRelTol and adjAbsTol

  adjAbsTol = curRes * L2Rel

  if (L2Abs > adjAbsTol) then
     adjAbsTol = L2Abs
  end if
 
  adjRelTol = L2
 
  call KSPSetTolerances(ksp,adjRelTol,adjAbsTol,adjDivTol,adjMaxIter,PETScIerr)

end subroutine setAdjointTolerances
