! there are a number of differenet functions in destroyPETScVars. They
! coorespond to the different creation functions:

subroutine destroyPETScVars

  use constants
  use ADjointPETSc, only : dRdWT, dRdwPreT, adjointKSP, adjointPETScVarsAllocated
  use inputAdjoint, only : approxPC
  use utils, only : EChk
  implicit none

  integer(kind=intType) ::  ierr
  
  if (adjointPETScVarsAllocated) then 

     ! Matrices
     call MatDestroy(dRdWT, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     if (ApproxPC) then
        call MatDestroy(dRdWPreT, ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
     
     call KSPDestroy(adjointKSP, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     adjointPETScVarsAllocated = .False.
  end if

end subroutine destroyPETScVars

