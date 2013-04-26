subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
     useObjective, level)
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
  !     * useObjective: If true, the derivative of Fx,Fy,Fz and Mx, My,Mz*
  !     *               are assembled into the size FMw vectors          *
  !     * level : What level to use to form the matrix. Level 1 is       *
  !     *         always the finest level                                *         
  !     ******************************************************************
  !
  use ADjointPetsc, only : FMw, dFcdW
  use BCTypes
  use blockPointers_d      
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use stencils
  use diffSizes
  use NKSolverVars, only : diag
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical, intent(in) :: useAD, usePC, useTranspose, useObjective
  integer(kind=intType), intent(in) :: level

  ! Local variables.
  integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
  integer(kind=intType) :: nColor, iColor, jColor, irow, icol, fmDim, frow
  integer(kind=intType) :: nTransfer, nState
  integer(kind=intType) :: n_stencil, i_stencil, n_force_stencil
  integer(kind=intType), dimension(:, :), pointer :: stencil, force_stencil
  real(kind=realType) :: delta_x, one_over_dx

#ifdef USE_COMPLEX
  complex(kind=realType) :: alpha, beta, Lift, Drag, CL, CD
  complex(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
  complex(kind=realType) :: alphad, betad, Liftd, Dragd, CLd, CDd
  complex(kind=realType), dimension(3) :: Forced, Momentd, cForced, cMomentd
#else
  real(kind=realType) :: alpha, beta, Lift, Drag, CL, CD
  real(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
  real(kind=realType) :: alphad, betad, Liftd, Dragd, CLd, CDd
  real(kind=realType), dimension(3) :: Forced, Momentd, cForced, cMomentd
#endif
  integer(kind=intType) :: liftIndex
  integer(kind=intType), dimension(:,:), pointer ::  colorPtr1, colorPtr2
  integer(kind=intType), dimension(:,:), pointer ::  globalCellPtr1, globalCellPtr2
  integer(kind=intType), dimension(:,:), pointer ::  colorPtr, globalCellPtr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, mm, colInd
  logical :: resetToRANS

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here
  alphad = zero
  betad  = zero
  machd  = zero
  machGridd = zero
  lengthRefd = zero
  pointRefd  = zero
  surfaceRefd = zero
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
        force_stencil => visc_force_w_stencil
        n_force_stencil = N_visc_force_w
     else
        stencil => euler_drdw_stencil
        n_stencil = N_euler_drdw
        force_stencil => euler_force_w_stencil
        n_force_stencil = N_euler_force_w
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

  ! Exchange data and call the residual to make sure its up to date
  ! withe current w
  call whalo2(1_intType, 1_intType, nw, .True., .True., .True.)
  call computeResidualNK ! This is the easiest way to do this

  ! Set delta_x
  delta_x = 1e-6_realType
  one_over_dx = one/delta_x
  rkStage = 0
  
  if (useObjective .and. useAD) then
     do fmDim=1,6
        call VecZeroEntries(FMw(fmDim), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Master Domain Loop
  domainLoopAD: do nn=1, nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn, level, 1)
     
     ! Set unknown sizes in diffSizes for AD routine
     ISIZE1OFDrfbcdata = nBocos
     ISIZE1OFDrfviscsubface = nViscBocos
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn, level)

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
        ! Note: The lumped dissipation doesn't quite result in a
        !3-cell stencil in each direction but we will still use PC
        !coloring. Not really a big deal.
     
        if (viscous) then
           call setup_3x3x3_coloring(nn, level,  nColor) ! dense 3x3x3 coloring
        else
           call setup_PC_coloring(nn, level,  nColor) ! Euler Colorings
        end if
     else
        if (viscous) then
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
                 flowDoms(nn, level, sps2)%w(:, :, :, :) =  flowDomsd(nn, 1, sps2)%wtmp
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
                 call block_res_d(nn, sps, .False., useObjective, &
                      alpha, alphad, beta, betad, liftIndex, Force, Forced, &
                      Moment, Momentd, lift, liftd, drag, dragd, cForce, &
                      cForced, cMoment, cMomentd, CL, CLD, CD, CDd)
#else
                 print *, 'Forward AD routines are not complexified'
                 stop
#endif
              else
                 call block_res(nn, sps, .False., .False., &
                      alpha, beta, liftIndex, Force, Moment, Lift, Drag, &
                      cForce, cMoment, CL, CD)
              end if

              ! If required, set values in the 6 vectors defined in
              ! FMw. We have to be a little carful actually. What we
              ! have to do is loop over subfaces where the cell
              ! centered forces are defined. Then for each force, we
              ! loop over its stencil. There should be at most 1
              ! peturbed cell/state in its stencil. Then, we can take
              ! the derivatives out F and M defined on the face. The
              ! derivatives are correct, since for objective we are
              ! using a simple sum. 
              sotreObjectivePartials: if (useObjective .and. useAD .and. .not. usePC) then
                 bocos: do mm=1,nBocos
                    if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic &
                         .or. BCType(mm) == NSWallIsothermal) then

                       ! Set the globalNodePtr depending on what face
                       ! we are on:

                       select case (BCFaceID(mm))
                       case (iMin)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(2, :, :)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(3, :, :)
                          globalCellPtr1 => globalCell(2, :, :)
                          globalCellPtr2 => globalCell(3, :, :)
                       case (iMax)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(il, :, :)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(il-1, :, :)
                          globalCellPtr1 => globalCell(il, :, :)
                          globalCellPtr2 => globalCell(il-1, :, :)
                       case (jMin)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(:, 2, :)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(:, 3, :)
                          globalCellPtr1 => globalCell(:, 2, :)
                          globalCellPtr2 => globalCell(:, 3, :)
                       case (jMax)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(:, jl, :)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(:, jl-1, :)
                          globalCellPtr1 => globalCell(:, jl, :)
                          globalCellPtr2 => globalCell(:, jl-1, :)
                       case (kMin)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(:, :, 2)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(:, :, 3)
                          globalCellPtr1 => globalCell(:, :, 2)
                          globalCellPtr2 => globalCell(:, :, 3)
                       case (kMax)
                          colorPtr1 => flowDomsd(nn, 1, 1)%color(:, :, kl)
                          colorPtr2 => flowDomsd(nn, 1, 1)%color(:, :, kl-1)
                          globalCellPtr1 => globalCell(:, :, kl)
                          globalCellPtr2 => globalCell(:, :, kl-1)
                       end select

                       ! These are the indices for the INTERNAL CELLS!
                       jBeg = BCData(mm)%jnBeg+1; jEnd = BCData(mm)%jnEnd
                       iBeg = BCData(mm)%inBeg+1; iEnd = BCData(mm)%inEnd
                       
                       do j=jBeg, jEnd ! This is a cell loop
                          do i=iBeg, iEnd ! This is a cell loop
                             ! Basically what we are doing is
                             ! searching the force stencil for this
                             ! face to see if one of them has been
                             ! petrubed. If one has, then we set the
                             ! values appropriately
                             forceStencilLoop: do i_stencil=1, n_force_stencil
                                ii = force_stencil(i_stencil, 1)
                                jj = force_stencil(i_stencil, 2)
                                kk = force_stencil(i_stencil, 3)
                     
                                if (kk == 0) then
                                   colorPtr => colorPtr1
                                   globalCellPtr => globalCellPtr1
                                else if (kk==1) then
                                   colorPtr => colorPtr2
                                   globalCellPtr => globalCellPtr2
                                end if

                                colInd = globalCellPtr(i+1+ii, j+1+jj)*nState + l -1
                                ! The extra + 1 is due to the pointer offset       
                                if (colorPtr(i+1+ii, j+1+jj) == iColor .and. &
                                     colInd >= 0) then
                                   ! This real cell has been peturbed!
                                   do fmDim = 1,3
                                      call VecSetValues(FMw(fmDim), 1, colInd, &
                                           bcDatad(mm)%F(i,j,fmDim), &
                                           ADD_VALUES, ierr) 
                                      call EChk(ierr, __FILE__, __LINE__)

                                      call VecSetValues(FMw(fmDim+3), 1, Colind, &
                                           bcDatad(mm)%M(i,j,fmDim), &
                                           ADD_VALUES, ierr) 
                                      call EChk(ierr, __FILE__, __LINE__)
                                      
                                      ! While we are at it, we have
                                      ! all the info we need for dFcdw
                                      fRow = BCData(mm)%FMCellIndex(i,j)*3 + fmDim - 1
                                      !call MatSetValues(dFcdw, 1, &
                                      !     fRow, 1, colInd, &
                                      !     bcDatad(mm)%F(i,j,fmDim), &
                                      !     ADD_VALUES, ierr)
                                      !call EChk(ierr, __FILE__, __LINE__)
                                   end do
                                end if
                             end do forceStencilLoop
                          end do
                       end do
                    end if
                 end do bocos
              end if sotreObjectivePartials
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
                    icol = flowDoms(nn, level, sps)%globalCell(i, j, k)
                    colorCheck: if (flowdomsd(nn, 1, 1)%color(i, j, k) == icolor&
                         .and. icol >= 0) then

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
                                   call setBlock(&
                                        flowDomsd(nn, 1, sps)%&
                                        dw_deriv(i+ii, j+jj, k+kk, :, :))
                                else
                                   ! Otherwise loop over spectral
                                   ! instances and set all.
                                   do sps2=1, nTimeIntervalsSpectral
                                      irow = flowDoms(nn, level, sps2)%&
                                           globalCell(i+ii, j+jj, k+kk)
                                      call setBlock(&
                                           flowDomsd(nn, 1, sps2)%&
                                           dw_deriv(i+ii, j+jj, k+kk, :, :))
                                   end do
                                end if useDiagPC
                             else
                                ! ALl other cells just set.
                                call setBlock(&
                                     flowDomsd(nn, 1, sps)%&
                                     dw_deriv(i+ii, j+jj, k+kk, :, :))
                             end if centerCell
                          end if onBlock
                       end do stencilLoop
                    end if colorCheck
                 end do iLoop
              end do jLoop
           end do kLoop
        end do colorLoop
     end do spectralLoop

     ! Deallocate and reset values for block nn
     call dealloc_derivative_values(nn, level)
  end do domainLoopAD

  if (useObjective .and. useAD) then
     do fmDim=1,6
        call VecAssemblyBegin(FMw(fmDim), ierr) 
        call EChk(ierr, __FILE__, __LINE__)
        call VecAssemblyEnd(FMw(fmDim), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do

     call MatAssemblyBegin(dFcdw, MAT_FINAL_ASSEMBLY, ierr)
     call MatAssemblyEnd(dFcdw, MAT_FINAL_ASSEMBLY, ierr)
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

  ! Reset the paraters to use segrated turbulence solve. 
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
       call EChk(1, __FILE__, __LINE__)
    end if
#endif
    
    if (.not. zeroFlag) then
       if (useTranspose) then
          call MatSetValuesBlocked(matrix, 1, icol, 1, irow, transpose(blk), &
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
