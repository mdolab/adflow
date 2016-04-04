subroutine FormJacobianANK
#ifndef USE_NO_PETSC

  use ANKSolverVars, only : dRdWPre, ANK_asmOverlap, &
       ANK_innerPreConIts,  MAT_FINAL_ASSEMBLY, ANK_iluFill, &
       ANK_subSpace, KSP_GMRES_CGS_REFINE_NEVER, ANK_KSP, ANK_CFL, deltaW
  use inputADjoint, only : viscPC
  use inputiteration
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
    implicit none
#define PETSC_AVOID_MPIF_H


#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Local Variables
  character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
  integer(kind=intType) ::ierr
  logical :: useAD, usePC, useTranspose, useObjective, tmp
  real(kind=realType) :: newCFL
  integer(kind=intType) :: i, j, k, l, ii, nn, sps, outerPreConIts
  real(kind=realType), pointer :: diag(:)

  ! Assemble the approximate PC (fine leve, level 1)
  useAD = .False.
  usePC = .True.
  useTranspose = .False.
  useObjective = .False.
  tmp = viscPC ! Save what is in viscPC and set to the NKvarible
  viscPC = .False.
  
  call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
       useObjective, .True., 1_intType)
  ! Reset saved value
  viscPC = tmp

  call VecGetArrayF90(deltaW, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now put the lumping on the digonal
  ii = 0 
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1, sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1, nwf
                    ii = ii + 1
                    diag(ii) = one/(ANK_CFL*dtl(i,j,k))
                    
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(deltaW, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatDiagonalSet(dRdwPre, deltaW, ADD_VALUES, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Setup KSP Options
  preConSide = 'right'
  localPCType = 'ilu'
  kspObjectType = 'gmres'
  globalPCType = 'asm'
  localOrdering = 'rcm'
  outerPreConIts = 1 
  ! Setup the KSP using the same code as used for the adjoint
  call setupStandardKSP(ANK_KSP, kspObjectType,  ANK_subSpace, &
       preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
       localOrdering, ANK_iluFill, ANK_innerPreConIts)

  ! Don't do iterative refinement for the NKSolver.
  call KSPGMRESSetCGSRefinementType(ANK_KSP, &
       KSP_GMRES_CGS_REFINE_NEVER, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
#endif

end subroutine FormJacobianANK
