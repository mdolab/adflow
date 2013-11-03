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
  use inputADjoint
  use diffSizes
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
  logical :: secondHalo, resetToRANS
  integer(kind=intType) :: FMDim, nState

  ! Values for block_res
  real(kind=realType) :: alpha, beta, force(3), moment(3)
  real(kind=realType) :: alphad, betad, forced(3), momentd(3)
  integer(kind=intType) :: liftIndex
  
  !Reference values for FD
  real(kind=realType) :: alpharef, betaref, machref, machGridRef, machCoefRef
  real(kind=realType), dimension(3) :: rotRateRef,rotcenterRef
  real(kind=realType), dimension(3) :: rotPointRef,pointRefRef
  real(kind=realType) :: lengthrefref, reynoldslengthref
  if (ndesignextra < 1) then
     ! No need to do anything here
     return
  end if

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  rkStage = 0
  currentLevel =1 
  groundLevel = 1

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! call getDirAngle to get the baseline values for alpha and beta
  call getDirAngle(velDirFreestream, LiftDirection, liftIndex, alpha, beta)

  ! Set delta_x
  delta_x = 1e-8
  one_over_dx = 1.0/delta_x
  rkStage = 0
  secondHalo = .True. 
  ! Master Domain Loop
  level = 1_intType
  nColor = nDesignExtra
  dFMdextra = zero

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use forzenTurbulent Adjiont
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn, level, 1)
     idxblk = nbkGlobal
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn, level)
     ISIZE1OFDrfbcdata = nBocos
     ISIZE1OFDrfviscsubface = nViscBocos

     ! Save the reference values in case we are doing finite differencing
     alpharef = alpha
     betaref = beta
     machref = mach
     machGridRef = machGrid
     machCoefRef = machCoef
     rotRateRef = cgnsDoms(idxblk)%rotRate
     rotcenterRef = cgnsDoms(idxblk)%rotCenter
     rotPointRef = rotPoint
     pointRefRef = pointRef
     lengthRefRef = lengthref
     reynoldslengthref = reynoldslength
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
        machCoefd = 0.0
        cgnsDomsd(idxblk)%rotrate(:) = 0.0
        cgnsDomsd(idxblk)%rotcenter(:) = 0.0
        rotpointd(:) = 0.0
        pointrefd(:) = 0.0
        lengthrefd = 0.0
        reynoldslengthd = 0.0
        reynoldsd = 0.0
        if (useAD) then
           if (nDesignAoA ==icolor-1) then
              alphad = one
           elseif (nDesignSSA==icolor-1) then
              betad = one
           elseif (nDesignMach ==icolor-1) then
              machd = one
              machCoefd = one
           elseif (nDesignMachGrid==icolor-1) then
              machGridd = one
              machCoefd = one
           elseif (nDesignPointRefX==icolor-1) then
              pointrefd(1) = one
           elseif (nDesignPointRefY==icolor-1) then
              pointrefd(2) = one
           elseif (nDesignPointRefZ==icolor-1) then
              pointrefd(3) = one
           elseif (nDesignLengthRef==icolor-1) then
              lengthrefd = one
              reynoldslengthd = one
           end if
        else
           if (nDesignAoA ==icolor-1) then
              alpha = alphaRef+delta_x
           elseif (nDesignSSA==icolor-1) then
              beta = betaRef+delta_x
           elseif (nDesignMach ==icolor-1) then
              mach = machRef +delta_x
              machCoef = machCoefRef + delta_x
           elseif (nDesignMachGrid==icolor-1) then
              machGrid = machGridRef+delta_x
              machCoef = machCoefRef + delta_x
           elseif (nDesignPointRefX==icolor-1) then
              pointref(1) = pointref(1) + delta_x
           elseif (nDesignPointRefY==icolor-1) then
              pointref(2) = pointref(2) + delta_x
           elseif (nDesignPointRefZ==icolor-1) then
              pointref(3) = pointref(3) + delta_x
           elseif(nDesignLengthRef==icolor-1) then
              lengthref = lengthref + delta_x
              reynoldslength = reynoldslength + delta_x
           end if
        end if

        ! Take all Derivatives
        do sps = 1,nTimeIntervalsSpectral
           call setPointers_d(nn, level, sps)

           ! Block-based residual
           if (useAD) then
#ifndef USE_COMPLEX       
              call block_res_d(nn, sps, .True., &
                   alpha, alphad, beta, betad, liftIndex, force, forced, &
                   moment, momentd)
#else
              print *,'Forward AD routines are not complexified'
              stop
#endif
           else
              call block_res(nn, sps, .True., &
                   alpha, beta, liftIndex, force, moment)
           end if

           ! Save the values of FMExtra and the derivatives

           do FMDim=1,3
              dFMdExtra(FMDim, iColor, sps) = dFMdExtra(FMDim, iColor, sps) + Forced(FMDim)
              dFMdExtra(FMDim+3, iColor, sps) = dFMdExtra(FMDim+3, iColor, sps) + Momentd(FMDim)
           end do

           ! Set the computed residual in dw_deriv. If using FD,
           ! actually do the FD calculation if AD, just copy out dw
           ! in flowdomsd
        
           do ll=1,nState
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
     machCoef = machCoefRef
     cgnsDoms(idxblk)%rotRate = rotRateRef
     cgnsDoms(idxblk)%rotcenter = rotCenterRef
     rotPoint = rotPointRef
     pointRef = pointRefRef
     LengthRef = LengthRefRef
     ReynoldsLength = ReynoldsLengthRef
  end do domainLoopAD

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the parameters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if
contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol, if useTranspose is False
    ! Sets a block at icol,irow with transpose of blk if useTranspose is True

    implicit none
    real(kind=realType), dimension(nState,1) :: blk
    integer(kind=intType) :: iii
       
    do iii=1, nState
       call MatSetValues(matrix, 1, irow*nState+iii-1, 1, icol,blk(iii,1), &
            ADD_VALUES, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end do
  end subroutine setBlock

#endif
end subroutine setupExtraResidualMatrix
