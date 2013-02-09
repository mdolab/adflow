subroutine setupExtraResidualMatrix(matrix, useAD)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the extra derivative matrix using a forward mode calc  *
  !     * There is one different flags that determine how this           *
  !     * routine is run:                                                *
  !     * useAD: if True, AD is used for derivative calculation, if      *
  !     *        False, FD is used.                                      *
  !     ******************************************************************

  use ADjointVars
  use ADjointPETSc, only : dFMdExtra
  use blockPointers      
  use communication      
  use inputDiscretization
  USE inputTimeSpectral  
  use iteration         
  use flowVarRefState   
  use inputTimeSpectral 
  use inputDiscretization
  use inputPhysics 
  use inputMotion     
  use stencils
  use cgnsGrid

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical :: useAD

  !     Local variables.
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ll
  integer(kind=intType) :: irow, icol
  real(kind=realType) :: delta_x, one_over_dx

  integer(kind=intType) :: n_stencil,i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: nColor, iColor, idxblk, level
  logical :: secondHalo
  integer(kind=intType) :: FMDim

  ! Values for block_res
  real(kind=realType) :: alpha, beta, Lift, Drag, CL, CD
  real(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
  real(kind=realType) :: alphad, betad, Liftd, Dragd, CLd, CDd
  real(kind=realType), dimension(3) :: Forced, Momentd, cForced, cMomentd
  integer(kind=intType) :: liftIndex
  
  !Reference values for FD
  real(kind=realType) :: alpharef, betaref, machref, machGridRef
  real(kind=realType), dimension(3) :: rotRateRef,rotcenterRef
  real(kind=realType), dimension(3) :: rotPointRef,pointRefRef

  rkStage = 0
  currentLevel =1 
  groundLevel = 1

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Call the residual to make sure its up to date withe current w
  call whalo2(1_intType, 1_intType, nw, .True.,.True.,.True.)
  call computeResidualNK

  ! call getDirAngle to get the baseline values for alpha and beta
  call getDirAngle(velDirFreestream,LiftDirection,liftIndex,alpha,beta)

  ! Set delta_x
  delta_x = 1e-8
  one_over_dx = 1.0/delta_x
  rkStage = 0
  secondHalo = .True. 
  ! Master Domain Loop
  level = 1_intType
  nColor = nDesignExtra
  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn,level,1)
     idxblk = nbkGlobal
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn, level)

     ! Save the reference values in case we are doing finite differencing
     alpharef = alpha
     betaref = beta
     machref = mach
     machGridRef = machGrid
     rotRateRef = cgnsDoms(idxblk)%rotRate
     rotcenterRef = cgnsDoms(idxblk)%rotCenter
     rotPointRef = rotPoint
     pointRefRef = pointRef

     ! Do 'Coloring' and extra varibales
     do iColor = 1,nColor !set colors based on extra vars....
        !zero derivatives
        do sps = 1,nTimeIntervalsSpectral
           flowDomsd(nn,1,sps)%dw_deriv(:,:,:,:,:) = zero
        end do
       
        !reset all of the seeds
        alphad = 0.0
        betad = 0.0
        machd = 0.0
        machGridd = 0.0
        cgnsDomsd(idxblk)%rotrate(:) = 0.0
        cgnsDomsd(idxblk)%rotcenter(:) = 0.0
        rotpointd(:) = 0.0
        pointrefd(:) = 0.0

        if (useAD) then
           if (nDesignAoA ==icolor-1) then
              alphad = 1.0
           elseif (nDesignSSA==icolor-1) then
              betad = 1.0
           elseif (nDesignMach ==icolor-1) then
              machd = 1.0
           elseif (nDesignMachGrid==icolor-1) then
              machGridd = 1.0
           elseif(nDesignRotX==icolor-1) then
              cgnsDomsd(idxblk)%rotrate(1) = 1.0
           elseif(nDesignRotY==icolor-1) then
              cgnsDomsd(idxblk)%rotrate(2) = 1.0
           elseif(nDesignRotZ==icolor-1) then
              cgnsDomsd(idxblk)%rotrate(3) = 1.0 
           elseif(nDesignRotCenX==icolor-1) then
              cgnsDomsd(idxblk)%rotcenter(1) = 1.0
              rotpointd(1) = 1.0
              !consider this!
              !+rotpointxcorrection
           elseif(nDesignRotCenY==icolor-1)then
              cgnsDomsd(idxblk)%rotcenter(2) = 1.0
              rotpointd(2) = 1.0
              !consider this!+rotpointxcorrection
           elseif(nDesignRotCenZ==icolor-1)then      
              cgnsDomsd(idxblk)%rotcenter(3)=1.0
              rotpointd(3) = 1.0
              !+rotpointzcorrection
           end if
        else
           if (nDesignAoA ==icolor-1) then
              alpha = alphaRef+delta_x
           elseif (nDesignSSA==icolor-1) then
              beta = betaRef+delta_x
           elseif (nDesignMach ==icolor-1) then
              mach = machRef +delta_x
           elseif (nDesignMachGrid==icolor-1) then
              machGrid = machGridRef+delta_x
           elseif(nDesignRotX==icolor-1) then
              cgnsDoms(idxblk)%rotrate(1) = rotrateref(1)+delta_x
           elseif(nDesignRotY==icolor-1) then
              cgnsDoms(idxblk)%rotrate(2) = rotrateref(2)+delta_x
           elseif(nDesignRotZ==icolor-1) then
              cgnsDoms(idxblk)%rotrate(3) = rotrateref(3)+delta_x 
           elseif(nDesignRotCenX==icolor-1) then
              !consider this!
              cgnsDoms(idxblk)%rotcenter(1) = rotcenterRef(1)+delta_x
              rotpoint(1)=rotPointRef(1)+delta_x
              !rotcenteradjb(1)+rotpointadjb(1)+rotpointxcorrection
           elseif(nDesignRotCenY==icolor-1)then
              !consider this!
              cgnsDoms(idxblk)%rotcenter(2) = rotcenterRef(2)+delta_x
              rotpoint(2)=rotPointRef(2)+delta_x
              !rotcenteradjb(2)+rotpointadjb(2)+rotpointxcorrection
           elseif(nDesignRotCenZ==icolor-1)then      
              cgnsDoms(idxblk)%rotcenter(3) = rotcenterRef(3)+delta_x
              rotpoint(3)=rotPointRef(3)+delta_x
              !+rotpointzcorrection
           end if
        end if

        ! Take all Derivatives
        do sps = 1,nTimeIntervalsSpectral
           call setPointers_d(nn, level, sps)

           ! Block-based residual
           if (useAD) then
#ifndef USE_COMPLEX         
              call block_res_d(nn, sps, .True., .False., &
                   alpha, alphad, beta, betad, liftIndex, Force, Forced, &
                   Moment, Momentd, lift, liftd, drag, dragd, cForce, &
                   cForced, cMoment, cMomentd, CL, CLd, CD, CDd)

#else
              print *,'Forward AD routines are not complexified'
              stop
#endif
           else
              call block_res(nn, sps, .True., .False., &
                   alpha, beta, liftIndex, Force, Moment, Lift, Drag, &
                   cForce, cMoment, CL, CD)
           end if
           
           do FMDim=1,3
              dFMdExtra(FMDim, iColor) = Forced(FMDim)
              dFMdExtra(FMDim+3, iColor) = Momentd(FMDim)
           end do

           ! Set the computed residual in dw_deriv. If using FD,
           ! actually do the FD calculation if AD, just copy out dw
           ! in flowdomsd
        
           do ll=1,nw
              do k=2,kl 
                 do j=2,jl
                    do i=2,il
                       if (useAD) then
                          flowDomsd(nn,1,sps)%dw_deriv(i,j,k,ll,1) = &
                               flowdomsd(nn,1,sps)%dw(i,j,k,ll)
                       else
                          flowDomsd(nn,1,sps)%dw_deriv(i,j,k,ll,1) = &
                               one_over_dx*(flowDoms(nn,1,sps)%dw(i,j,k,ll) - &
                               flowDomsd(nn,1,sps)%dwtmp(i,j,k,ll))
                       end if
                    end do
                 end do
              end do
           end do
        end do
        
        ! Set derivatives by block in "matrix" after we've peturbed
        ! all states in "color"
        do sps = 1,nTimeIntervalsSpectral
           call setPointers(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    !Assume the entire vector is present and set accordingly
                    irow = flowDoms(nn,1,sps)%globalCell(i,j,k)
                    if ( irow >= 0) then
                       icol = icolor-1
                       call setBlock(flowDomsd(nn, 1, sps)%dw_deriv(i, j, k, :, 1))
                    end if ! Color If check
                 end do ! i loop
              end do ! j loop
           end do ! k loop
        end do ! spectral Loop
     end do !color loop

     ! Deallocate and reset Values
     call dealloc_derivative_values(nn, level)

     ! Reset values to reference 
     alpha = alpharef
     beta = betaref
     mach = machref
     machGrid = machGridRef
     cgnsDoms(idxblk)%rotRate = rotRateRef
     cgnsDoms(idxblk)%rotcenter = rotCenterRef
     rotPoint = rotPointRef
     pointRef = pointRefRef
     
  end do domainLoopAD
  
  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol, if useTranspose is False
    ! Sets a block at icol,irow with transpose of blk if useTranspose is True

    implicit none
    real(kind=realType), dimension(nw,1) :: blk
    integer(kind=intType) :: iii
       
    do iii=1, nw
       call MatSetValues(matrix, 1, irow*nw+iii-1, 1, icol,blk(iii,1), &
            ADD_VALUES, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end do
  end subroutine setBlock

#endif
end subroutine setupExtraResidualMatrix
