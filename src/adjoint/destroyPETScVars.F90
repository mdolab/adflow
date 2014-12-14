! there are a number of differenet functions in destroyPETScVars. They
! coorespond to the different creation functions:

subroutine destroyPETScVars

  use ADjointPETSc
  use inputADjoint    
  use blockPointers
  use costFunctions
  use inputTimeSpectral
  implicit none

  integer(kind=intType) :: i, nLevels, sps
  
#ifndef USE_NO_PETSC
  if (adjointPETScVarsAllocated) then 

     ! Matrices
     call MatDestroy(dRdWT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     if (ApproxPC) then
        call MatDestroy(dRdWPreT, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end if
     
     call MatDestroy(dFcdw, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call MatDestroy(dFcdx, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call MatDestroy(doAdx, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call MatDestroy(dFndFc, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call MatDestroy(dFdx, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call MatDestroy(dFdw, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
        
     if (.not. useMatrixFreedRdx) then
        call MatDestroy(dRdx, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call MatDestroy(dRda, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end if
     
     ! Vectors
     do sps=1,nTimeIntervalsSpectral
        do i=1,nFM
           call VecDestroy(FMw(i, sps), PETScIerr)
           call EChk(PETScIerr,__FILE__,__LINE__)
        end do
     end do
     deallocate(FMw)
       
     do sps=1,ntimeIntervalsSpectral
        do i=1,nFM
           call VecDestroy(FMx(i, sps), PETScIerr)
           call EChk(PETScIerr,__FILE__,__LINE__)
        end do
     end do
     deallocate(FMx)
     
     call vecDestroy(overArea, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call vecDestroy(fCell, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call vecDestroy(fNode, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     if (allocated(dRda_data)) then
        deallocate(dRda_data)
     end if
     
     call KSPDestroy(adjointKSP, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     adjointPETScVarsAllocated = .False.
  end if
#endif
end subroutine destroyPETScVars

