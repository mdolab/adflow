subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
     useObjective, frozenTurb, level)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the state derivative matrix using a forward mode calc  *
  !     * There are three different flags that determine how this        *
  !     * routine is run:                                                *
  !     * useAD: if True, AD is used for derivative calculation, if      *
  !     *        False, FD is used.                                      *
  !     * usePC: if True, the reduced 1st order stencil with dissipation *
  !     *        lumping is assembled instead of the actual exact        *
  !     *        full stencil jacobian                                   *
  !     * useTranspose: If true, the transpose of dRdw is assembled.     *
  !     *               For use with the adjoint this must be true.      *
  !     * useObjective: If true, the force matrix is assembled           *
  !     *                                                                *
  !     * level : What level to use to form the matrix. Level 1 is       *
  !     *         always the finest level                                *         
  !     ******************************************************************
  !
  use BCTypes
  use blockPointers
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use stencils
  use diffSizes
  use NKSolverVars, only : diag
  use communication
  use adjointVars
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
  integer(kind=intType), intent(in) :: level

  ! Local variables.
  integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
  integer(kind=intType) :: nColor, iColor, jColor, irow, icol, fmDim, frow
  integer(kind=intType) :: nTransfer, nState, tmp, icount
  integer(kind=intType) :: n_stencil, i_stencil
  integer(kind=intType), dimension(:, :), pointer :: stencil
  real(kind=realType) :: delta_x, one_over_dx

#ifdef USE_COMPLEX
  complex(kind=realType) :: alpha, beta, sepSensor, Cavitation
  complex(kind=realType) :: alphad, betad, sepSensord, Cavitationd
  complex(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment, forced, momentd
  complex(kind=realType), dimension(:,:), allocatable :: blk
#else
  real(kind=realType) :: alpha, beta, sepSensor, Cavitation
  real(kind=realType) :: alphad, betad,  sepSensord, Cavitationd
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment, forced, momentd
  real(kind=realType), dimension(:,:), allocatable :: blk
#endif
  integer(kind=intType) :: liftIndex
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, mm, colInd
  logical :: resetToRANS
  real :: val

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurb ) then
     nState = nwf
  else
     nState = nw
  endif

  allocate(blk(nState, nState))
  ! Exchange data and call the residual to make sure its up to date
  ! withe current w
  call whalo2(1_intType, 1_intType, nw, .True., .True., .True.)

  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here
  alphad = zero
  betad  = zero
  machd  = zero
  machGridd = zero
  machcoefd = zero
  pointRefd  = zero
  lengthRefd = zero
  prefd = zero
  tempfreestreamd = zero
  reynoldsd = zero
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

  rkStage = 0

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set a pointer to the correct set of stencil depending on if we are
  ! using the first order stencil or the full jacobian

  if (usePC) then
     if (viscous .and. viscPC) then
        stencil => visc_pc_stencil
        n_stencil = N_visc_pc
     else
        stencil => euler_pc_stencil
        n_stencil = N_euler_pc
     end if

     ! Very important to use only Second-Order dissipation for PC 
     lumpedDiss=.True.
  else
     if (viscous) then
        stencil => visc_drdw_stencil
        n_stencil = N_visc_drdw
     else
        stencil => euler_drdw_stencil
        n_stencil = N_euler_drdw
     end if
  end if

  ! Need to trick the residual evalution to use coupled (mean flow and
  ! turbulent) together.

  ! If we want to do the matrix on a coarser level, we must first
  ! restrict the fine grid solutions, since it is possible the
  ! NKsolver was used an the coarse grid solutions are (very!) out of
  ! date. 
  
  ! Assembling matrix on coarser levels is not entirely implemented yet. 
  currentLevel = level
  groundLevel = level

  ! Set delta_x
  delta_x = 1e-9_realType
  one_over_dx = one/delta_x
  rkStage = 0
  
  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurb .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Allocate the additional memory we need for doing forward mode AD
  !  derivatives and copy any required reference values:
  if (.not. derivVarsAllocated) then 
     call alloc_derivative_values(level)
  end if
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        call zeroADSeeds(nn,level, sps)
     end do
  end do

  if (usePC) then 
     call referenceShockSensor
  end if

  if (.not. useAD) then 
     call setFDReference(level)
  end if
  ! Master Domain Loop
  domainLoopAD: do nn=1, nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn, level, 1)
     ! Set unknown sizes in diffSizes for AD routine
     ISIZE1OFDrfbcdata = nBocos
     ISIZE1OFDrfviscsubface = nViscBocos

     ! Setup the coloring for this block depending on if its
     ! drdw or a PC
     
     ! List of all Coloring Routines:
     !   Debugging Colorings Below:
     !       call setup_3x3x3_coloring(nn, level,  nColor)
     !       call setup_5x5x5_coloring(nn, level,  nColor)
     !       call setup_BF_coloring(nn, level,  nColor)
     !   Regular:
     !       call setup_PC_coloring(nn, level,  nColor)
     !       call setup_dRdw_euler_coloring(nn, level,  nColor)
     !       call setup_dRdw_visc_coloring(nn, level,  nColor)

     if (usePC) then
      
        if (viscous .and. viscPC) then
           call setup_3x3x3_coloring(nn, level,  nColor) ! dense 3x3x3 coloring
        else
           call setup_PC_coloring(nn, level,  nColor) ! Euler Colorings
        end if
     else
        if (viscous) then
           !call setup_5x5x5_coloring(nn, level,  nColor)
           call setup_dRdw_visc_coloring(nn, level,  nColor)! Viscous/RANS
        else 
           call setup_dRdw_euler_coloring(nn, level,  nColor) ! Euler Colorings
        end if
     end if

     spectralLoop: do sps=1, nTimeIntervalsSpectral
        ! Set pointers and derivative pointers
        call setPointers_d(nn, level, sps)

        ! Do Coloring and perturb states
        colorLoop: do iColor = 1, nColor
           do sps2 = 1, nTimeIntervalsSpectral
              flowDomsd(nn, 1, sps2)%dw_deriv(:, :, :, :, :) = zero
           end do

           ! Master State Loop
           stateLoop: do l=1, nState

              ! Reset All States and possibe AD seeds
              do sps2 = 1, nTimeIntervalsSpectral
                 if (.not. useAD) then 
                    do ll=1,nw
                       do k=0,kb
                          do j=0,jb
                             do i=0,ib
                                flowDoms(nn, level, sps2)%w(i,j,k,ll) =  flowDomsd(nn, 1, sps2)%wtmp(i,j,k,ll)
                             end do
                          end do
                       end do
                    end do
                 end if
                 
                 if (useAD) then
                    flowdomsd(nn, 1, sps2)%w = zero ! This is actually w seed
                 end if
              end do
              
              ! Peturb w or set AD Seed according to iColor
              do k=0, kb
                 do j=0, jb
                    do i=0, ib
                       if (flowdomsd(nn, 1, 1)%color(i, j, k) == icolor) then
                          if (useAD) then
                             flowDomsd(nn, 1, sps)%w(i, j, k, l) = one
                          else

                             w(i, j, k, l) = w(i, j, k, l) + delta_x
                          end if
                       end if
                    end do
                 end do
              end do
             
              ! Run Block-based residual 
              if (useAD) then
#ifndef USE_COMPLEX
                 call block_res_d(nn, sps, .False., &
                      alpha, alphad, beta, betad, liftIndex, force, forced, moment, momentd,&
                      sepSensor, sepSensord, Cavitation, Cavitationd, frozenTurb)
#else
                 print *, 'Forward AD routines are not complexified'
                 stop
#endif
              else
                 call block_res(nn, sps, .False., alpha, beta, liftIndex, force, moment, &
                 sepSensor, Cavitation, frozenTurb)
              end if

              ! Set the computed residual in dw_deriv. If using FD, 
              ! actually do the FD calculation if AD, just copy out dw
              ! in flowdomsd

              ! Compute/Copy all derivatives
              do sps2 = 1, nTimeIntervalsSpectral
                 do ll=1, nState
                    do k=2, kl 
                       do j=2, jl
                          do i=2, il
                             if (useAD) then
                                flowDomsd(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                     flowdomsd(nn, 1, sps2)%dw(i, j, k, ll)
                             else
                                if (sps2 == sps) then
                                   ! If the peturbation is on this
                                   ! instance, we've computed the spatial
                                   ! contribution so subtrace dwtmp

                                   flowDomsd(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                        one_over_dx*&
                                        (flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                        flowDomsd(nn, 1, sps2)%dwtmp(i, j, k, ll))
                                else
                                   ! If the peturbation is on an off
                                   ! instance, only subtract dwtmp2
                                   ! which is the reference result
                                   ! after initres

                                   flowDomsd(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                        one_over_dx*(&
                                        flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                        flowDomsd(nn, 1, sps2)%dwtmp2(i, j, k, ll))
                                end if
                             end if
                          end do
                       end do
                    end do
                 end do
              end do
           end do stateLoop

           ! Set derivatives by block in "matrix" after we've peturbed
           ! all states in "color"

           kLoop: do k=0, kb
              jLoop: do j=0, jb
                 iLoop: do i=0, ib
                    iCol = flowDoms(nn, level, sps)%globalCell(i, j, k)
                    ! if (iCol < 0) then
                    !    iCol = -(iCol + 1)
                    ! end if
                    colorCheck: if (flowdomsd(nn, 1, 1)%color(i, j, k) == icolor) then! &
                       !.and. icol >= 0) then
                       !colorCheck: if (flowdomsd(nn, 1, 1)%color(i, j, k) == icolor) then

                       ! i, j, k are now the "Center" cell that we
                       ! actually petrubed. From knowledge of the
                       ! stencil, we can simply take this cell and
                       ! using the stencil, set the values around it
                       ! in PETSc

                       stencilLoop: do i_stencil=1, n_stencil
                          ii = stencil(i_stencil, 1)
                          jj = stencil(i_stencil, 2)
                          kk = stencil(i_stencil, 3)

                          ! Check to see if the cell in this
                          ! sentcil is on a physical cell, not a
                          ! halo/BC halo
                          onBlock: if ( i+ii >= 2 .and. i+ii <= il .and. &
                                        j+jj >= 2 .and. j+jj <= jl .and. &
                                        k+kk >= 2 .and. k+kk <= kl) then 

                             irow = flowDoms(nn, level, sps)%globalCell(&
                                 i+ii, j+jj, k+kk)

                             centerCell: if ( ii == 0 .and. jj == 0 &
                                  .and. kk == 0) then
                                useDiagPC: if (usePC .and. useDiagTSPC) then
                                   ! If we're doing the PC and we want
                                   ! to use TS diagonal form, only set
                                   ! values for on-time insintance
                                   blk = flowDomsd(nn, 1, sps)%dw_deriv(i+ii, j+jj, k+kk, &
                                        1:nstate, 1:nstate)
                                   call setBlock(blk)
                                else
                                   ! Otherwise loop over spectral
                                   ! instances and set all.
                                   do sps2=1, nTimeIntervalsSpectral
                                      irow = flowDoms(nn, level, sps2)%&
                                           globalCell(i+ii, j+jj, k+kk)
                                      blk = flowDomsd(nn, 1, sps2)%dw_deriv(i+ii, j+jj, k+kk, &
                                           1:nstate, 1:nstate)
                                      call setBlock(blk)
                                   end do
                                end if useDiagPC
                             else
                                ! ALl other cells just set.
                                blk = flowDomsd(nn, 1, sps)%dw_deriv(i+ii, j+jj, k+kk, &
                                     1:nstate, 1:nstate)
                                call setBlock(blk)
                             end if centerCell
                          end if onBlock
                       end do stencilLoop
                    end if colorCheck
                 end do iLoop
              end do jLoop
           end do kLoop
        end do colorLoop
     end do spectralLoop
  end do domainLoopAD

  ! Deallocate and reset values 
  if (.not. useAD) then 
     call resetFDReference(level)
  end if
  
  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd  (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !Return dissipation Parameters to normal -> VERY VERY IMPORTANT
  if (usePC) then
     lumpedDiss = .False.
  end if

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  deallocate(blk)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow, icol, if useTranspose is False
    ! Sets a block at icol, irow with transpose of blk if useTranspose is True

    implicit none
#ifndef USE_COMPLEX
    real(kind=realType), dimension(nState, nState) :: blk
#else
    complex(kind=realType), dimension(nState, nState) :: blk
#endif
    ! local variables
    integer(kind=intType) :: i, j
    logical :: zeroFlag
    zeroFlag = .False.
#ifndef USE_COMPLEX
    ! Check if the blk is all zeros

    zeroFlag = .True.
    do i = 1, nState
       do j = 1, nState
          if ( .not. blk(i,j) == zero) then
             zeroFlag = .False.
          end if
       end do
    end do
    
    ! Check if the blk has nan
    if (isnan(sum(blk))) then
       print *,'Bad Block:',blk
       print *,'irow:',irow
       print *,'icol',icol
       print *,'nn:',nn
       print *,'ijk:',i,j,k
       call EChk(1, __FILE__, __LINE__)
    end if
#endif

    if (.not. zeroFlag) then
       if (useTranspose) then
          blk = transpose(blk)
          call MatSetValuesBlocked(matrix, 1, icol, 1, irow, blk, &
               ADD_VALUES, ierr)
          call EChk(ierr, __FILE__, __LINE__)
       else
          call MatSetValuesBlocked(matrix, 1, irow, 1, icol, blk, &
               ADD_VALUES, ierr)
          call EChk(ierr, __FILE__, __LINE__)
       end if
    end if

  end subroutine setBlock
#endif
end subroutine setupStateResidualMatrix
