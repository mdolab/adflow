subroutine setupExtraResidualMatrix(matrix,useAD)

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
  Mat mat_copy

  ! Input Variables
  logical :: useAD

  !     Local variables.
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ll,ii,jj,kk
  integer(kind=intType) :: irow,icol,ilow,ihigh
  real(kind=realType) :: delta_x,one_over_dx

  real(kind=realType)::alpha,beta
  integer(kind=intType)::liftIndex

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               ::setupTime,trace
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
  

  if (myid==0) print *,'Setting up dRda in forward Mode...'

  rkStage = 0
  currentLevel =1 
  groundLevel = 1
  ! Start Timer
  time(1) = mpi_wtime()

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  !The Extra variables are sometimes dense, therfore no stencil will be used

  ! Call the residual to make sure its up to date withe current w
  call computeResidualNK ! This is the easiest way to do this

  ! call getDirAngle to get the baseline values for alpha and beta
  call getDirAngle(velDirFreestream,LiftDirection,liftIndex,alpha,beta)

  ! Set delta_x
  delta_x = 1e-8
  one_over_dx = 1.0/delta_x
  rkStage = 0
  secondHalo = .True. 
  ! Master Domain Loop
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
     nColor = nDesignExtra
     !print *,'nColor',nColor
     !print *,'Alpha',alpha
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
        !print *,'icolor',icolor
        !zero derivatives
        do sps = 1,nTimeIntervalsSpectral
           flowDomsd(sps)%dw_deriv(:,:,:,:,:) = 0.0
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
           !print *,'icolor',icolor-1
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
           !print *,'icolor',icolor-1,nDesignAoA
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
              !print *,'alpha',alpha,alphad,beta,liftIndex,mach,machd
              call block_res_extra_extra_d(nn,sps,alpha,alphad,beta,&
                                          &betad,liftIndex)
              !print *,'liftdir',liftDirection
              !print *,'AD Not Implmented Yet'
              !stop
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
                          flowDomsd(sps)%dw_deriv(i,j,k,ll,1) = &
                               flowdomsd(sps)%dw(i,j,k,ll)
                       else
                          
                          flowDomsd(sps)%dw_deriv(i,j,k,ll,1) = &
                               one_over_dx*(flowDoms(nn,1,sps)%dw(i,j,k,ll) - &
                               flowDomsd(sps)%dwtmp(i,j,k,ll))
                          !print *,'deriv',flowDomsd(sps)%dw_deriv(i,j,k,ll,1)
                       end if
                       
                    end do
                 end do
              end do
           end do
        end do
     
        ! Set derivatives by block in "matrix" after we've peturbed
        ! all states in "color"
        do sps = 1,nTimeIntervalsSpectral
           call setPointersAdj(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    !Assume the entire vector is present and set accordingly
                    irow = flowDoms(nn,1,sps)%globalCell(i,j,k)
                    if ( irow >= 0) then
                       icol = icolor-1
                       call setBlock(flowDomsd(sps)%dw_deriv(i,j,k,:,1))
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

#ifdef USE_PETSC_3
  call MatSetOption(matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call MatSetOption(matrix,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  if (myid == 0) then
     print *,'Assembly time:',setupTime
  end if
  ! Debugging ONLY!
  !call writeOutMatrix()

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol, if useTranspose is False
    ! Sets a block at icol,irow with transpose of blk if useTranspose is True

    implicit none
    real(kind=realType), dimension(nw,1) :: blk
    integer(kind=intType) :: iii
       
    do iii=1,nw
       !print *,'output',abs(blk(iii,1)),irow,icol
       if (abs(blk(iii,1)).ne. 0.0)then
          call MatSetValues(matrix,1,irow*nw+iii-1,1,icol,blk(iii,1),ADD_VALUES,ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end if
    end do
    
    
  end subroutine setBlock



end subroutine setupExtraResidualMatrix
