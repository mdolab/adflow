subroutine setupExtraResidualMatrix(matrix,useAD)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the extra derivative matrix using a forward mode calc  *
  !     * There is one different flags that determine how this           *
  !     * routine is run:                                                *
  !     * useAD: if True, AD is used for derivative calculation, if      *
  !     *        False, FD is used.                                      *
  !     ******************************************************************
  !

  use ADjointVars
  use ADjointPETSc , only:localInfo
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use inputDiscretization ! spaceDiscr
  USE inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputTimeSpectral   ! spaceDiscr
  use inputDiscretization
  use inputPhysics 
  use inputMotion         ! rotpoint,rotpointd
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
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ll,ii,jj,kk
  integer(kind=intType) :: irow,icol,ilow,ihigh
  real(kind=realType) :: delta_x,one_over_dx

  real(kind=realType)::alpha,beta
  integer(kind=intType)::liftIndex

  integer(kind=intType) :: n_stencil,i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: nColor,iColor,idxblk
  logical :: secondHalo
  
  !Reference values for FD
  real(kind=realType)::alpharef,betaref,machref, machGridRef
  real(kind=realType), dimension(3) :: rotRateRef,rotcenterRef
  real(kind=realType), dimension(3) :: rotPointRef,pointRefRef

  !Seed values for AD
  real(kind=realType)::alphad,betad!,machd, machGridd
  !real(kind=realType), dimension(3) :: rotRated,rotcenterd
  !real(kind=realType), dimension(3) :: rotPointd!,pointRefd
  

  rkStage = 0
  currentLevel =1 
  groundLevel = 1

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  !The Extra variables are sometimes dense, therfore no stencil will be used

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

  nColor = nDesignExtra
  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn,1,1)
     idxblk = nbkGlobal
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn)

     ! Setup the coloring for this block depending on if its
     ! drdw or a PC

     !let the DV number be the color
    
     alpharef = alpha
     betaref = beta
     machref = mach
     machGridRef = machGrid
     rotRateRef = cgnsDoms(idxblk)%rotRate
     rotcenterRef = cgnsDoms(idxblk)%rotCenter
     rotPointRef = rotPoint
     pointRefRef = pointRef


     ! Do Coloring and perturb states
     do iColor = 1,nColor !set colors based on extra vars....
        !zero derivatives
        do sps = 1,nTimeIntervalsSpectral
           flowDomsd(nn,1,sps)%dw_deriv(:,:,:,:,:) = 0.0
        end do
       
        !reset all of the seeds
        alphad = 0.0
        betad = 0.0
        machd = 0.0
        machGridd = 0.0
        cgnsDomsd(idxblk)%rotrate(:) = 0.0
        cgnsDomsd(idxblk)%rotcenter(:) = 0.0
        rotpointd(:) = 0.0
        !pointrefd(:) = 0.0

        if (useAD) then
           !Set the seeds by color
           if (nDesignAoA ==icolor-1) then
              !Angle of Attack
              alphad = 1.0
           elseif (nDesignSSA==icolor-1) then
              ! Side slip angle
              betad = 1.0
           elseif (nDesignMach ==icolor-1) then
              !Mach Number
              machd = 1.0
           elseif (nDesignMachGrid==icolor-1) then
              !Mach NumberGrid
              machGridd = 1.0
           elseif(nDesignRotX==icolor-1) then
              !X Rotation
              cgnsDomsd(idxblk)%rotrate(1) = 1.0
           elseif(nDesignRotY==icolor-1) then
              !Y Rotation
              cgnsDomsd(idxblk)%rotrate(2) = 1.0
           elseif(nDesignRotZ==icolor-1) then
              !Z Rotation
              cgnsDomsd(idxblk)%rotrate(3) = 1.0 
           elseif(nDesignRotCenX==icolor-1) then
              !X Rotation Center
              cgnsDomsd(idxblk)%rotcenter(1) = 1.0
              rotpointd(1) = 1.0
              !consider this!
              !+rotpointxcorrection
           elseif(nDesignRotCenY==icolor-1)then
              !Y Rotation Center
              cgnsDomsd(idxblk)%rotcenter(2) = 1.0
              rotpointd(2) = 1.0
              !consider this!+rotpointxcorrection
           elseif(nDesignRotCenZ==icolor-1)then      
              !Z Rotation Center
              cgnsDomsd(idxblk)%rotcenter(3)=1.0
              rotpointd(3) = 1.0
              !+rotpointzcorrection
           end if

        else
           alpha = alpharef
           beta = betaref
           mach = machref
           machGrid = machGridRef
           cgnsDoms(idxblk)%rotRate = rotRateRef
           cgnsDoms(idxblk)%rotcenter = rotCenterRef
           rotPoint = rotPointRef
           pointRef = pointRefRef

           if (nDesignAoA ==icolor-1) then
              !Angle of Attack
              alpha = alphaRef+delta_x
           elseif (nDesignSSA==icolor-1) then
              ! Side slip angle
              beta = betaRef+delta_x
           elseif (nDesignMach ==icolor-1) then
              !Mach Number
              mach = machRef +delta_x
           elseif (nDesignMachGrid==icolor-1) then
              !Mach NumberGrid
              machGrid = machGridRef+delta_x
           elseif(nDesignRotX==icolor-1) then
              !X Rotation
              cgnsDoms(idxblk)%rotrate(1) = rotrateref(1)+delta_x
           elseif(nDesignRotY==icolor-1) then
              !Y Rotation
              cgnsDoms(idxblk)%rotrate(2) = rotrateref(2)+delta_x
           elseif(nDesignRotZ==icolor-1) then
              !Z Rotation
              cgnsDoms(idxblk)%rotrate(3) = rotrateref(3)+delta_x 
           elseif(nDesignRotCenX==icolor-1) then
              !X Rotation Center
              !consider this!
              cgnsDoms(idxblk)%rotcenter(1) = rotcenterRef(1)+delta_x
              rotpoint(1)=rotPointRef(1)+delta_x
              !rotcenteradjb(1)+rotpointadjb(1)+rotpointxcorrection
           elseif(nDesignRotCenY==icolor-1)then
              !Y Rotation Center
              !consider this!
              cgnsDoms(idxblk)%rotcenter(2) = rotcenterRef(2)+delta_x
              rotpoint(2)=rotPointRef(2)+delta_x
              !rotcenteradjb(2)+rotpointadjb(2)+rotpointxcorrection
           elseif(nDesignRotCenZ==icolor-1)then      
              !Z Rotation Center
              cgnsDoms(idxblk)%rotcenter(3) = rotcenterRef(3)+delta_x
              rotpoint(3)=rotPointRef(3)+delta_x
              !+rotpointzcorrection
           end if
        end if

        ! Take all Derivatives
        do sps = 1,nTimeIntervalsSpectral
           ! Block-based residual
           if (useAD) then
              call block_res_extra_extra_d(nn,sps,alpha,alphad,beta,&
                   betad,liftIndex)
           else
              call block_res_extra(nn,sps,alpha,beta,liftIndex)
           end if


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
                       call setBlock(flowDomsd(nn,1,sps)%dw_deriv(i,j,k,:,1))
                    end if ! Color If check
                 end do ! i loop
              end do ! j loop
           end do ! k loop
        end do ! spectral Loop
     end do !color loop

     ! Deallocate and reset Values
     call dealloc_derivative_values(nn)
     
  end do domainLoopAD
  
  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatSetOption(matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol, if useTranspose is False
    ! Sets a block at icol,irow with transpose of blk if useTranspose is True

    implicit none
    real(kind=realType), dimension(nw,1) :: blk
    integer(kind=intType) :: iii
       
    do iii=1,nw
       call MatSetValues(matrix,1,irow*nw+iii-1,1,icol,blk(iii,1),ADD_VALUES,ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do
  end subroutine setBlock

#endif
end subroutine setupExtraResidualMatrix
