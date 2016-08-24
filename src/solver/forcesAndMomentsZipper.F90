subroutine forcesAndMomentsZipper(cFp, cFv, cMp, cMv, sps)

  use communication
  use blockPointers
  use BCTypes
  use flowVarRefState
  use inputPhysics
  use costFunctions
  use overset, only : nodeZipperScatter, globalNodalVec, localZipperNodes, localZipperTp, localZipperTv
  use inputTimeSpectral
  use inputIteration
  use utils, only : EChk, setPointers
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

  ! Input/Output
  real(kind=realType), dimension(3), intent(out) :: cFp, cFv
  real(kind=realType), dimension(3), intent(out) :: cMp, cMv
  integer(kind=intType), intent(in) :: sps

  ! Working
  integer(kind=intType) :: i, j, k, l, ierr, nn, mm, gind, ind, iVar
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, rowStart, rowEnd
  real(kind=realType), dimension(:, :, :), pointer :: xx
  real(kind=realType), dimension(:), pointer :: xPtr, pPtr, vPtr, localPtr
  integer(kind=intType), dimension(:, :), pointer :: gnp

  real(kind=realType), dimension(3) :: x1, x2, x3, xc, ss, norm, refPoint
  real(kind=realType), dimension(3) :: pp1, pp2, pp3, presForce
  real(kind=realType), dimension(3) :: vv1, vv2, vv3, viscForce
  real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv, MTmp
  real(kind=realType) :: fact, scaleDim, triArea

  ! PETsc
  Vec, pointer :: localVec

  ! Set the actual scaling factor such that ACTUAL forces are computed
  scaleDim = pRef/pInf

  ! Determine the reference point for the moment computation in
  ! meters.
  refPoint(1) = LRef*pointRef(1)
  refPoint(2) = LRef*pointRef(2)
  refPoint(3) = LRef*pointRef(3)

  ! Communicate the nodes, pressure tractions and viscous tractions to
  ! the root proc for the triangles. We do a generic loop

  call VecGetOwnershipRange(globalNodalVec, rowStart, rowEnd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do iVar=1, 3

     call VecGetArrayF90(globalNodalVec, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     domainLoop: do nn=1, nDom
        call setPointers(nn, 1, sps)
        
        bocoLoop: do mm=1, nBocos

           if(BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              
              select case (BCFaceID(mm))
              case (iMin)
                 xx => x(1, :, :, :)
                 gnp => globalNode(1, :, :)
              case (iMax)
                 xx => x(il, :, :, :)
                 gnp => globalNode(il, :, :)
              case (jMin)
                 xx => x(:, 1, :, :)
                 gnp => globalNode(:, 1, :)
              case (jMax)
                 xx => x(:, jl, :, :)
                 gnp => globalNode(:, jl, :)
              case (kMin)
                 xx => x(:, :, 1, :)          
                 gnp => globalNode(:, :, 1)
              case (kMax)
                 xx => x(:, :, kl, :)
                 gnp => globalNode(:, :, kl)
              end select

              jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd  
              
              ! Loop over the nodes of the subface:
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    gInd = gnp(i+1, j+1)
                    ind = (gInd-rowStart/3)*3+1
                    if (ind < 0) then 
                       print *,'something wrong:', myid, gind, rowstart, ind
                       stop
                    end if

                    select case(iVar)
                    case (1) ! Nodes
                       localPtr(ind:ind+2) = xx(i+1, j+1, :)
                    case (2)
                       localPtr(ind:ind+2) = bcData(mm)%Tp(i, j, :)
                    case (3)
                       localPtr(ind:ind+2) = bcData(mm)%Tv(i, j, :)
                    end select
                 end do
              end do
           end if
        end do bocoLoop
     end do domainLoop
     
     call vecRestoreArrayF90(globalNodalVec, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     select case(iVar)
     case(1)
        localVec => localZipperNodes
     case(2)
        localVec => localZipperTp
     case(3)
        localVec => localZipperTv
     end select
     
     ! Perform the scatter from the global x vector to zipperNodes
     call VecScatterBegin(nodeZipperScatter, globalNodalVec, &
          localVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecScatterEnd(nodeZipperScatter, globalNodalVec, &
          localVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do

  ! Puck out pointers for the nodes and tractions
  cFp = zero
  cFv = zero
  cMp = zero
  cMv = zero
  if (myid == 0) then 

     Fp = zero
     Fv = zero
     Mp = zero
     Mv = zero
     
     call VecGetArrayF90(localZipperTp, pPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(localZipperTv, vPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(localZipperNodes, xPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Number of triangles is the length of localZipperNodes /9

     do i=1, size(xPtr)/9

        ! Nodes for the triangles
        x1 = xPtr((i-1)*9+1:(i-1)*9+3)
        x2 = xPtr((i-1)*9+4:(i-1)*9+6)
        x3 = xPtr((i-1)*9+7:(i-1)*9+9)

        ! Nodal pressure tractions
        pp1 = pPtr((i-1)*9+1:(i-1)*9+3)
        pp2 = pPtr((i-1)*9+4:(i-1)*9+6)
        pp3 = pPtr((i-1)*9+7:(i-1)*9+9)

        ! Nodal viscous tractions
        vv1 = vPtr((i-1)*9+1:(i-1)*9+3)
        vv2 = vPtr((i-1)*9+4:(i-1)*9+6)
        vv3 = vPtr((i-1)*9+7:(i-1)*9+9)

        ! Compute area
        call cross_prod(x2-x1, x3-x1, norm)
        ss = half * norm
        triArea = norm2(ss)
        
        ! This is the actual integration
        presForce = third*(pp1 + pp2 + pp3) * triArea
        viscForce = third*(vv1 + vv2 + vv3) * triArea

        ! Add to Fp and Fv
        Fp = Fp + presForce
        Fv = Fv + viscForce

        ! Compute cell center
        xc(:) = third*(x1 + x2 + x3) - refPoint(:)
        
        ! Add to Mp and Mv
        call cross_prod(xc, presForce, Mtmp)
        Mp = Mp + Mtmp

        call cross_prod(xc, viscForce, Mtmp)
        Mv = Mv + Mtmp

     end do

     ! Return the array pointers
     call VecRestoreArrayF90(localZipperNodes, xPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecRestoreArrayF90(localZipperTp, pPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecRestoreArrayF90(localZipperTv, vPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     fact = two/(gammaInf*pInf*MachCoef*MachCoef &
          *surfaceRef*LRef*LRef*scaleDim)

     ! Convert the values to coefficients
     cFp(:) =  fact*Fp
     cFv(:) =  fact*Fv

     fact = fact/(lengthRef*LRef)
     cMp(:) = Mp(:)*fact
     cMv(:) = Mv(:)*fact

  end if

end subroutine forcesAndMomentsZipper
