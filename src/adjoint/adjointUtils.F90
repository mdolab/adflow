! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.
module adjointUtils

contains

  subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
       useObjective, frozenTurb, level, useTurbOnly, useCoarseMats)

    !      Compute the state derivative matrix using a forward mode calc
    !      There are three different flags that determine how this
    !      routine is run:
    !      useAD: if True, AD is used for derivative calculation, if
    !             False, FD is used.
    !      usePC: if True, the reduced 1st order stencil with dissipation
    !             lumping is assembled instead of the actual exact
    !             full stencil jacobian
    !      useTranspose: If true, the transpose of dRdw is assembled.
    !                    For use with the adjoint this must be true.
    !      useObjective: If true, the force matrix is assembled
    !      level : What level to use to form the matrix. Level 1 is
    !              always the finest level
    !
    use block, only : flowDomsd
    use blockPointers
    use inputDiscretization
    use inputTimeSpectral
    use inputPhysics
    use iteration
    use flowVarRefState
    use inputAdjoint
    use stencils
    use diffSizes
    use communication
    use adjointVars
    use turbMod
    use surfaceFamilies, only : fullFamList
    use oversetUtilities, only : fracToWeights
    use utils, only : EChk, setPointers, getDirAngle, setPointers_d
    use haloExchange, only : whalo2
    use masterRoutines, only : block_res_state, master
    use agmg, only : agmgLevels, A, coarseIndices, coarseOversetIndices
#ifndef USE_COMPLEX
    use masterRoutines, only : block_res_state_d
#endif
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Matrix Variable
    Mat :: matrix

    ! Input Variables
    logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
    logical, intent(in), optional :: useTurbOnly, useCoarseMats
    integer(kind=intType), intent(in) :: level

    ! Local variables.
    integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
    integer(kind=intType) :: nColor, iColor, jColor, irow, icol, fmDim, frow
    integer(kind=intType) :: nTransfer, nState, lStart, lEnd, tmp, icount, cols(8), nCol
    integer(kind=intType) :: n_stencil, i_stencil, m, iFringe, fInd, lvl, orderturbsave
    integer(kind=intType), dimension(:, :), pointer :: stencil
    real(kind=alwaysRealType) :: delta_x, one_over_dx
    real(kind=realType) :: weights(8)
    real(kind=realType), dimension(:,:), allocatable :: blk
    integer(kind=intType), dimension(2:10) :: coarseRows
    integer(kind=intType), dimension(8, 2:10) :: coarseCols
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, mm, colInd
    logical :: resetToRANS, turbOnly, flowRes, turbRes, buildCoarseMats

    ! Determine if we are assembling a turb only PC
    turbOnly = .False.
    if (present(useTurbOnly)) then
       turbOnly = useTurbOnly
    end if

    buildCoarseMats = .False.
    if (present(useCoarseMats)) then
       buildCoarseMats = useCoarseMats
    end if

    if (turbOnly) then
       ! we are making a PC for turbulence only KSP
       flowRes = .False.
       turbRes = .True.
       lStart = nt1
       lEnd = nt2
       nState = nt2 - nt1 + 1
    else
       ! We are making a "matrix" for either NK or adjoint
       ! Setup number of state variable based on turbulence assumption
       if ( frozenTurb ) then
          flowRes = .True.
          turbRes = .False.
          lStart = 1
          lEnd =nwf
          nState = nwf
       else
          flowRes = .True.
          turbRes = .True.
          lStart = 1
          lEnd =nw
          nState = nw
       endif
    end if

    ! Generic block to use while setting values
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
    pinfdimd = zero
    rhoinfdimd = zero
    tinfdimd = zero

    rkStage = 0

    ! Zero out the matrix before we start
    call MatZeroEntries(matrix, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the diagonal to 1 if the cell is not a compute cell:

    ! Make an identity block
    blk = zero
    do i=1, nState
       blk(i,i) = one
    end do

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   if (iblank(i, j, k) /= 1) then
                      iRow = flowDoms(nn, level, sps)%globalCell(i, j, k)
                      cols(1) = irow
                      nCol = 1

                      if (buildCoarseMats) then
                         do lvl=1, agmgLevels-1
                            coarseRows(lvl+1) = coarseIndices(nn, lvl)%arr(i, j, k)
                            coarseCols(1, lvl+1) = coarseRows(lvl+1)
                         end do
                      end if

                      call setBlock(blk)
                   end if
                end do
             end do
          end do
       end do
    end do

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
       ! also use first order advection terms for turbulence
       orderturbsave = orderturb
       orderturb = firstOrder
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
    if (.not. derivVarsAllocated .and. useAD) then
       call allocDerivativeValues(level)
    end if

    ! For AD the initial seeds must all be zeroed.
    if (useAD) then
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn, level, sps)
             call zeroADSeeds(nn, level, sps)
          end do
       end do
    end if

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)

          ! Allocate some extra routines used only for assembly
          allocate(&
               flowDoms(nn, level, sps)%dw_deriv(2:il, 2:jl, 2:kl, 1:nw, 1:nw), &
               flowDoms(nn, level, sps)%wtmp(0:ib,0:jb,0:kb,1:nw), &
               flowDoms(nn, level, sps)%dwtmp(0:ib,0:jb,0:kb,1:nw), &
               flowDoms(nn, level, sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw), &
          stat=ierr)
          call EChk(ierr, __FILE__, __LINE__)

          if (sps == 1) then
             allocate(flowDoms(nn, level, sps)%color(0:ib, 0:jb, 0:kb), stat=ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if
       end do
    end do

    ! For the PC we don't linearize the shock sensor so it must be
    ! computed here.

    if (usePC) then
       call referenceShockSensor
    end if

    ! For FD, the initial reference values must be computed and stored.
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
          ! Set pointers and (possibly derivative pointers)
          if (useAD) then
             call setPointers_d(nn, level, sps)
          else
             call setPointers(nn, level, sps)
          end if

          ! Do Coloring and perturb states
          colorLoop: do iColor = 1, nColor
             do sps2 = 1, nTimeIntervalsSpectral
                flowDoms(nn, 1, sps2)%dw_deriv(:, :, :, :, :) = zero
             end do

             ! Master State Loop
             stateLoop: do l=lStart, lEnd

                ! Reset All States and possibe AD seeds
                do sps2 = 1, nTimeIntervalsSpectral
                   if (.not. useAD) then
                      do ll=1,nw
                         do k=0,kb
                            do j=0,jb
                               do i=0,ib
                                  flowDoms(nn, level, sps2)%w(i,j,k,ll) =  flowDoms(nn, 1, sps2)%wtmp(i,j,k,ll)
                               end do
                            end do
                         end do
                      end do
                   end if

                   if (useAD) then
                      flowdomsd(nn, 1, sps2)%w = zero ! This is actually w seed
                   end if
                end do

                ! Peturb w or set AD Seed according to iColor. Note:
                ! Do NOT try to putt he useAD if check inside the
                ! color if check. ifort barfs hard-core on that and it
                ! segfaults with AVX2
                if (useAD) then
                   do k=0, kb
                      do j=0, jb
                         do i=0, ib
                            if (flowdoms(nn, 1, 1)%color(i, j, k) == icolor) then
                               flowDomsd(nn, 1, sps)%w(i, j, k, l) = one
                            end if
                         end do
                      end do
                   end do
                else
                   do k=0, kb
                      do j=0, jb
                         do i=0, ib
                            if (flowdoms(nn, 1, 1)%color(i, j, k) == icolor) then
                               w(i, j, k, l) = w(i, j, k, l) + delta_x
                            end if
                         end do
                      end do
                   end do
                end if

                ! Run Block-based residual
                if (useAD) then
#ifndef USE_COMPLEX
                   call block_res_state_d(nn, sps)
#else
                   print *, 'Forward AD routines are not complexified'
                   stop
#endif
                else
                   call block_res_state(nn, sps, useFlowRes=flowRes, useTurbRes=turbRes)
                end if

                ! Set the computed residual in dw_deriv. If using FD,
                ! actually do the FD calculation if AD, just copy out dw
                ! in flowdomsd

                ! Compute/Copy all derivatives
                do sps2 = 1, nTimeIntervalsSpectral
                   do ll=lStart, lEnd
                      do k=2, kl
                         do j=2, jl
                            do i=2, il
                               if (useAD) then
                                  flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                       flowdomsd(nn, 1, sps2)%dw(i, j, k, ll)
                               else
                                  if (sps2 == sps) then
                                     ! If the peturbation is on this
                                     ! instance, we've computed the spatial
                                     ! contribution so subtrace dwtmp

                                     flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                          one_over_dx * &
                                          (flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                          flowDoms(nn, 1, sps2)%dwtmp(i, j, k, ll))
                                  else
                                     ! If the peturbation is on an off
                                     ! instance, only subtract dwtmp2
                                     ! which is the reference result
                                     ! after initres

                                     flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                          one_over_dx*(&
                                          flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                          flowDoms(nn, 1, sps2)%dwtmp2(i, j, k, ll))
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
                      colBlank: if (flowDoms(nn, level, sps)%iblank(i, j, k) == 1 .or. &
                           flowDoms(nn, level, sps)%iBlank(i, j, k) == -1) then

                         ! If the cell we perturned ('iCol') is an
                         ! interpolated cell, we don't actually use
                         ! iCol, rather we use the 8 real donors that
                         ! comprise the cell's value.
                         if (flowDoms(nn, level, sps)%iblank(i, j, k) == 1) then
                            cols(1) = flowDoms(nn, level, sps)%globalCell(i, j, k)
                            nCol = 1

                            if (buildCoarseMats) then
                               do lvl=1, agmgLevels-1
                                  coarseCols(1, lvl+1) = coarseIndices(nn, lvl)%arr(i, j, k)
                               end do
                            end if

                         else
                            do m=1,8
                               cols(m) = flowDoms(nn, level, sps)%gInd(m, i, j, k)

                               if (buildCoarseMats) then
                                  do lvl=1, agmgLevels-1
                                     coarseCols(m, lvl+1) = coarseOversetIndices(nn, lvl)%arr(m, i, j, k)
                                  end do
                               end if
                            end do

                            fInd = fringePtr(1,i,j,k)
                            call fracToWeights(flowDoms(nn, level, sps)%fringes(fInd)%donorFrac, &
                                 weights)
                            nCol = 8
                         end if

                         colorCheck: if (flowdoms(nn, 1, 1)%color(i, j, k) == icolor) then

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

                                  if (buildCoarseMats) then
                                     do lvl=1, agmgLevels-1
                                        coarseRows(lvl+1) = coarseIndices(nn, lvl)%arr(i+ii, j+jj, k+kk)
                                     end do
                                  end if

                                  rowBlank: if (flowDoms(nn, level, sps)%iBlank(i+ii, j+jj, k+kk) == 1) then

                                     centerCell: if ( ii == 0 .and. jj == 0  .and. kk == 0) then
                                        useDiagPC: if (usePC .and. useDiagTSPC) then
                                           ! If we're doing the PC and we want
                                           ! to use TS diagonal form, only set
                                           ! values for on-time insintance
                                           blk = flowDoms(nn, 1, sps)%dw_deriv(i+ii, j+jj, k+kk, &
                                                lStart:lEnd, lStart:lEnd)
                                           call setBlock(blk)
                                        else
                                           ! Otherwise loop over spectral
                                           ! instances and set all.
                                           do sps2=1, nTimeIntervalsSpectral
                                              irow = flowDoms(nn, level, sps2)%&
                                                   globalCell(i+ii, j+jj, k+kk)
                                              blk = flowDoms(nn, 1, sps2)%dw_deriv(i+ii, j+jj, k+kk, &
                                                   lStart:lEnd, lStart:lEnd)
                                              call setBlock(blk)
                                           end do
                                        end if useDiagPC
                                     else
                                        ! ALl other cells just set.
                                        blk = flowDoms(nn, 1, sps)%dw_deriv(i+ii, j+jj, k+kk, &
                                             lStart:lEnd, lStart:lEnd)
                                        call setBlock(blk)
                                     end if centerCell
                                  end if rowBlank
                               end if onBlock
                            end do stencilLoop
                         end if colorCheck
                      end if colBlank
                   end do iLoop
                end do jLoop
             end do kLoop
          end do colorLoop
       end do spectralLoop
    end do domainLoopAD

    ! PETSc Matrix Assembly begin
    call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (buildCoarseMats) then
       do lvl=2, agmgLevels
          call MatAssemblyBegin(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
          call EChk(ierr, __FILE__, __LINE__)
       end do
    end if

    ! Maybe we can do something useful while the communication happens?
    ! Deallocate the temporary memory used in this routine.

    ! Deallocate and reset values
    if (.not. useAD) then
       call resetFDReference(level)
    end if

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          deallocate(&
               flowDoms(nn, 1, sps)%dw_deriv, &
               flowDoms(nn, 1, sps)%wTmp, &
               flowDoms(nn, 1, sps)%dwTmp, &
               flowDoms(nn, 1, sps)%dwTmp2)
          if (sps==1) then
             deallocate(flowDoms(nn, 1, 1)%color)
          end if
       end do
    end do

    ! Return dissipation Parameters to normal -> VERY VERY IMPORTANT
    if (usePC) then
       lumpedDiss = .False.
       ! also recover the turbulence advection order
       orderturb = orderturbsave
    end if

    ! Reset the correct equation parameters if we were useing the frozen
    ! Turbulent
    if (resetToRANS) then
       equations = RANSEquations
    end if

    deallocate(blk)

    ! Complete the matrix assembly.
    call MatAssemblyEnd  (matrix, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (buildCoarseMats) then
       do lvl=2, agmgLevels
          call MatAssemblyEnd(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
          call EChk(ierr, __FILE__, __LINE__)
       end do
    end if

    call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! ================= Important ===================

    ! We must run the residual computation to make sure that all
    ! intermediate variables are up to date. We can just call master
    ! for this. No need to recompute spatial terms.
    call master(.false.)

  contains

    subroutine setBlock(blk)
      ! Sets a block at irow, icol, if useTranspose is False
      ! Sets a block at icol, irow with transpose of blk if useTranspose is True

      use utils, only : myisnan
      implicit none
      real(kind=realType), dimension(nState, nState) :: blk

      ! local variables
      integer(kind=intType) :: i, j, tmp, iRowSet, iColSet
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
      if (myisnan(sum(blk))) then
         print *,'Bad Block:',blk
         print *,'irow:',irow
         print *,'icol',cols(1:ncol)
         print *,'nn:',nn
         print *,'ijk:',i,j,k
         call EChk(1, __FILE__, __LINE__)
      end if
#endif

      if (.not. zeroFlag) then
         if (nCol == 1) then
            if (useTranspose) then
               blk = transpose(blk)
               call MatSetValuesBlocked(matrix, 1, cols(1), 1, irow, blk, &
                    ADD_VALUES, ierr)
               call EChk(ierr, __FILE__, __LINE__)
            else
               call MatSetValuesBlocked(matrix, 1, irow, 1, cols(1), blk, &
                    ADD_VALUES, ierr)
               call EChk(ierr, __FILE__, __LINE__)
            end if
         else
            if (useTranspose) then
               blk = transpose(blk)
               do m=1, ncol
                  if (cols(m) >= 0) then
                     call MatSetValuesBlocked(matrix, 1, cols(m), 1, irow, blk*weights(m), &
                          ADD_VALUES, ierr)
                     call EChk(ierr, __FILE__, __LINE__)
                  end if
               end do
            else
               do m=1, ncol
                  if (cols(m) >= 0) then
                     call MatSetValuesBlocked(matrix, 1, irow, 1, cols(m), blk*weights(m), &
                          ADD_VALUES, ierr)
                     call EChk(ierr, __FILE__, __LINE__)
                  end if
               end do
            end if
         end if

         ! Extension for setting coarse grids:
         if (buildCoarseMats) then
            if (nCol == 1) then
               do lvl=2, agmgLevels
                  if (useTranspose) then
                     ! Loop over the coarser levels
                     call MatSetValuesBlocked(A(lvl), 1, coarseCols(1, lvl), 1, coarseRows(lvl), &
                          blk, ADD_VALUES, ierr)
                  else
                     call MatSetValuesBlocked(A(lvl), 1, coarseRows(lvl), 1, coarseCols(1, lvl), &
                          blk, ADD_VALUES, ierr)
                  end if
               end do
            else
               do m=1, nCol
                  do lvl=2, agmgLevels
                     if (coarseCols(m, lvl) >= 0) then
                        if (useTranspose) then
                           ! Loop over the coarser levels
                           call MatSetValuesBlocked(A(lvl), 1, coarseCols(m, lvl), 1, coarseRows(lvl), &
                                blk*weights(m), ADD_VALUES, ierr)
                        else
                           call MatSetValuesBlocked(A(lvl), 1, coarseRows(lvl), 1, coarseCols(m, lvl), &
                                blk*weights(m), ADD_VALUES, ierr)
                        end if
                     end if
                  end do
               end do
            end if
         end if
      end if
    end subroutine setBlock
  end subroutine setupStateResidualMatrix

  subroutine allocDerivativeValues(level)

    use constants
    use block, only : flowDoms, flowDomsd
    use blockPointers, only : nDom, il, jl, kl, ie, je, ke, ib, jb, kb, BCData, &
         nBOcos, nViscBocos
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : winf, winfd, nw, nt1, nt2
    use inputDiscretization, only : useApproxWallDistance
    use inputPhysics, only : wallDistanceNeeded
    use communication, only : adflow_comm_world
    use wallDistanceData, only : xSurfVec, xSurfVecd!, PETSC_DETERMINE
    use BCPointers_b
    use adjointVars, only : derivVarsAllocated
    use utils, only : EChk, setPointers, getDirAngle
    use cgnsGrid, only : cgnsDoms, cgnsDomsd, cgnsNDom
    implicit none

    ! Input parameters
    integer(kind=intType) :: level

    ! Local variables
    integer(kind=intType) :: sps,ierr,i,j,k,l, mm, nn
    integer(kind=intType) :: iBeg, jBeg, iStop, jStop, isizemax, jsizemax
    integer(kind=intType) :: inBeg, jnBeg, inStop, jnStop
    integer(kind=intType) :: massShape(2), max_face_size
    integer(kind=intType) :: iBoco, nDataset, iData, nDirichlet, iDirichlet, nArray
    ! First create the derivative flowdoms structure flowDomsd. Note we
    ! only allocate information for the finest grid.

    allocate(flowDomsd(nDom, 1, nTimeIntervalsSpectral), stat=ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! winfd hasn't be allocated so we'll do it here
    allocate(winfd(size(winf)))

    ! If we are not using RANS with walDistance create a dummy xSurfVec
    ! since one does not yet exist
    if (.not. wallDistanceNeeded .or. .not. useApproxWallDistance) then
       do sps=1, nTimeIntervalsSpectral
          call VecCreateMPI(ADFLOW_COMM_WORLD, 1, PETSC_DETERMINE, xSurfVec(1, sps), ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end do
    end if

    ! Duplicate the PETSc Xsurf Vec, but only on the first level:
    allocate(xSurfVecd(nTimeIntervalsSpectral))
    do sps=1, nTimeIntervalsSpectral
       call VecDuplicate(xSurfVec(1, sps), xSurfVecd(sps), ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)

          ! Allocate d2wall if not already done so
          if (.not. associated(flowDoms(nn, 1, sps)%d2wall)) then
             allocate(flowDoms(nn, 1, sps)%d2wall(2:il, 2:jl, 2:kl))
             call EChk(ierr,__FILE__,__LINE__)
          end if

          ! Now allocate all valus that have a differentiable
          ! dependence.
          allocate(&
               flowDomsd(nn, level, sps)%x(0:ie,0:je,0:ke,3), &
               flowDomsd(nn, level, sps)%vol(0:ib,0:jb,0:kb), &
               flowDomsd(nn, level, sps)%si(0:ie,1:je,1:ke,3), &
               flowDomsd(nn, level, sps)%sj(1:ie,0:je,1:ke,3), &
               flowDomsd(nn, level, sps)%sk(1:ie,1:je,0:ke,3), &
               flowDomsd(nn, level, sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
               flowDomsd(nn, level, sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
               flowDomsd(nn, level, sps)%rotMatrixK(2:il,2:jl,kl,3,3), &
               flowDomsd(nn, level, sps)%s(ie,je,ke,3),      &
               flowDomsd(nn, level, sps)%sFaceI(0:ie,je,ke), &
               flowDomsd(nn, level, sps)%sFaceJ(ie,0:je,ke), &
               flowDomsd(nn, level, sps)%sFaceK(ie,je,0:ke), &
               flowDomsd(nn, level, sps)%w (0:ib,0:jb,0:kb,1:nw), &
               flowDomsd(nn, level, sps)%dw(0:ib,0:jb,0:kb,1:nw), &
               flowDomsd(nn, level, sps)%fw(0:ib,0:jb,0:kb,1:nw), &
               flowDomsd(nn, level, sps)%scratch(0:ib,0:jb,0:kb,5), &
               flowDomsd(nn, level, sps)%p(0:ib,0:jb,0:kb), &
               flowDomsd(nn, level, sps)%gamma(0:ib,0:jb,0:kb),  &
               flowDomsd(nn, level, sps)%aa(0:ib,0:jb,0:kb), &
               flowDomsd(nn, level, sps)%ux(il,jl,kl), &
               flowDomsd(nn, level, sps)%uy(il,jl,kl), &
               flowDomsd(nn, level, sps)%uz(il,jl,kl), &
               flowDomsd(nn, level, sps)%vx(il,jl,kl), &
               flowDomsd(nn, level, sps)%vy(il,jl,kl), &
               flowDomsd(nn, level, sps)%vz(il,jl,kl), &
               flowDomsd(nn, level, sps)%wx(il,jl,kl), &
               flowDomsd(nn, level, sps)%wy(il,jl,kl), &
               flowDomsd(nn, level, sps)%wz(il,jl,kl), &
               flowDomsd(nn, level, sps)%qx(il,jl,kl), &
               flowDomsd(nn, level, sps)%qy(il,jl,kl), &
               flowDomsd(nn, level, sps)%qz(il,jl,kl), &
               flowDomsd(nn, level, sps)%rlv(0:ib,0:jb,0:kb), &
               flowDomsd(nn, level, sps)%rev(0:ib,0:jb,0:kb), &
               flowDomsd(nn, level, sps)%radI(1:ie,1:je,1:ke), &
               flowDomsd(nn, level, sps)%radJ(1:ie,1:je,1:ke), &
               flowDomsd(nn, level, sps)%radK(1:ie,1:je,1:ke), &
               flowDomsd(nn, level, sps)%BCData(nBocos), &
               flowDomsd(nn, level, sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvti1(je,ke,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvti2(je,ke,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvtj1(ie,ke,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvtj2(ie,ke,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvtk1(ie,je,nt1:nt2), &
               flowDomsd(nn, level, sps)%bvtk2(ie,je,nt1:nt2), &
               flowDomsd(nn, level, sps)%d2Wall(2:il,2:jl,2:kl), &
               stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          allocate(flowDomsd(nn, level, sps)%viscSubface(nviscBocos), &
               stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Set the number of bocos/viscbocs
          flowdomsd(nn, level, sps)%nBocos = flowdoms(nn, level, sps)%nbocos
          flowDomsd(nn, level, sps)%nViscBocos = flowDoms(nn, level, sps)%nViscBocos
          bocoLoop: do mm=1,nBocos

             ! Store the cell range of the boundary subface
             ! a bit easier.

             iBeg = BCData(mm)%icbeg; iStop = BCData(mm)%icend
             jBeg = BCData(mm)%jcbeg; jStop = BCData(mm)%jcend

             inBeg = BCData(mm)%inBeg; inStop = BCData(mm)%inEnd
             jnBeg = BCdata(mm)%jnBeg; jnStop = BCData(mm)%jnEnd
             allocate(&
                  flowDomsd(nn, level, sps)%BCData(mm)%norm(iBeg:iStop,jBeg:jStop,3), &
                  flowDomsd(nn, level, sps)%BCData(mm)%rface(iBeg:iStop,jBeg:jStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%Fp(inBeg+1:inStop, jnBeg+1:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%Fv(inBeg+1:inStop, jnBeg+1:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%Tp(inBeg:inStop, jnBeg:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%Tv(inBeg:inStop, jnBeg:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%F(inBeg:inStop, jnBeg:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%T(inBeg:inStop, jnBeg:jnStop, 3),&
                  flowDomsd(nn, level, sps)%BCData(mm)%area(inBeg+1:inStop, jnBeg+1:jnStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%uSlip(iBeg:iStop,jBeg:jStop,3), &
                  flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall(iBeg:iStop,jBeg:jStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%ptInlet(iBeg:iStop,jBeg:jStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%htInlet(iBeg:iStop,jBeg:jStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%ttInlet(iBeg:iStop,jBeg:jStop), &
                  flowDomsd(nn, level, sps)%BCData(mm)%turbInlet(iBeg:iStop,jBeg:jStop,nt1:nt2), &
                  flowDomsd(nn, level, sps)%BCData(mm)%ps(iBeg:iStop,jBeg:jStop), stat=ierr)

             call EChk(ierr,__FILE__,__LINE__)
          end do bocoLoop

          viscbocoLoop: do mm=1,nviscBocos

             iBeg = BCData(mm)%inBeg + 1; iStop = BCData(mm)%inEnd
             jBeg = BCData(mm)%jnBeg + 1; jStop = BCData(mm)%jnEnd

             allocate(&
                  flowDomsd(nn, level, sps)%viscSubface(mm)%tau(iBeg:iStop,jBeg:jStop,6), &
                  flowDomsd(nn, level, sps)%viscSubface(mm)%q(iBeg:iStop,jBeg:jStop,6), &
                  stat=ierr)
             call EChk(ierr,__FILE__,__LINE__)
          enddo viscbocoLoop
       end do
    end do

    ! Allocate the derivatives values for the CGNS data structure used
    ! to store boundary condition values
    do nn=1, cgnsNDom
       nBocos = cgnsDoms(nn)%nBocos
       allocate(cgnsDomsd(nn)%bocoInfo(nBocos))
       do iBoco = 1,nBocos
          if (associated(cgnsDoms(nn)%bocoInfo(iBoco)%dataSet)) then
             nDataSet = size(cgnsDoms(nn)%bocoInfo(iBoco)%dataSet)
             allocate(cgnsDomsd(nn)%bocoInfo(iBoco)%dataSet(nDataSet))

             do iData=1, nDataSet
                if (associated(cgnsDoms(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)) then
                   nDirichlet = size(cgnsDoms(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)
                   allocate(cgnsDomsd(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(nDirichlet))

                   do iDirichlet = 1, nDirichlet
                      nArray = size(cgnsDoms(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr)
                      allocate(cgnsDomsd(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(nArray))
                      cgnsDomsd(nn)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(nArray) = zero
                   end do
                end if
             end do
          end if
       end do
    end do

    derivVarsAllocated = .True.
  end subroutine allocDerivativeValues

  subroutine zeroADSeeds(nn, level, sps)

    use constants
    use block, only : flowDomsd, flowDoms
    use blockPointers
    use inputTimeSpectral
    use flowVarRefState
    use inputPhysics
    use BCPointers_b
    use communication
    use oversetData, only : oversetPresent
    use cgnsGrid, only : cgnsDoms, cgnsDomsd, cgnsNDom
    use actuatorRegionData, only : nActuatorRegions, actuatorRegionsd
    implicit none

    ! Input parameters
    integer(kind=intType) :: nn, level, sps

    ! Working parameters
    integer(kind=intType) :: mm, i, iDom
    integer(kind=intType) :: iBoco, iData, iDirichlet
    flowDomsd(nn, level, sps)%d2wall = zero
    flowDomsd(nn, level, sps)%x = zero
    flowDomsd(nn, level, sps)%si = zero
    flowDomsd(nn, level, sps)%sj = zero
    flowDomsd(nn, level, sps)%sk = zero
    flowDomsd(nn, level, sps)%vol = zero

    flowDomsd(nn, level, sps)%s = zero
    flowDomsd(nn, level, sps)%sFaceI = zero
    flowDomsd(nn, level, sps)%sFaceJ = zero
    flowDomsd(nn, level, sps)%sFaceK = zero

    flowDomsd(nn, level, sps)%w = zero
    flowDomsd(nn, level, sps)%dw = zero
    flowDomsd(nn, level, sps)%fw = zero
    flowDomsd(nn, level, sps)%scratch = zero

    flowDomsd(nn, level, sps)%p = zero
    flowDomsd(nn, level, sps)%gamma = zero
    flowDomsd(nn, level, sps)%aa = zero

    flowDomsd(nn, level, sps)%rlv = zero
    flowDomsd(nn, level, sps)%rev = zero

    flowDomsd(nn, level, sps)%radI = zero
    flowDomsd(nn, level, sps)%radJ = zero
    flowDomsd(nn, level, sps)%radK = zero

    flowDomsd(nn, level, sps)%ux = zero
    flowDomsd(nn, level, sps)%uy = zero
    flowDomsd(nn, level, sps)%uz = zero
    flowDomsd(nn, level, sps)%vx = zero
    flowDomsd(nn, level, sps)%vy = zero
    flowDomsd(nn, level, sps)%vz = zero
    flowDomsd(nn, level, sps)%wx = zero
    flowDomsd(nn, level, sps)%wy = zero
    flowDomsd(nn, level, sps)%wz = zero
    flowDomsd(nn, level, sps)%qx = zero
    flowDomsd(nn, level, sps)%qy = zero
    flowDomsd(nn, level, sps)%qz = zero

    flowDomsd(nn, level, sps)%bmti1 = zero
    flowDomsd(nn, level, sps)%bmti2 = zero
    flowDomsd(nn, level, sps)%bmtj1 = zero
    flowDomsd(nn, level, sps)%bmtj2 = zero
    flowDomsd(nn, level, sps)%bmtk1 = zero
    flowDomsd(nn, level, sps)%bmtk2 = zero
    flowDomsd(nn, level, sps)%bvti1 = zero
    flowDomsd(nn, level, sps)%bvti2 = zero
    flowDomsd(nn, level, sps)%bvtj1 = zero
    flowDomsd(nn, level, sps)%bvtj2 = zero
    flowDomsd(nn, level, sps)%bvtk1 = zero
    flowDomsd(nn, level, sps)%bvtk2 = zero

    bocoLoop: do mm=1, flowDoms(nn, level, sps)%nBocos
       flowDomsd(nn, level, sps)%BCData(mm)%norm= zero
       flowDomsd(nn, level, sps)%bcData(mm)%rface = zero
       flowDomsd(nn, level, sps)%bcData(mm)%Fv = zero
       flowDomsd(nn, level, sps)%bcData(mm)%Fp = zero
       flowDomsd(nn, level, sps)%bcData(mm)%Tv = zero
       flowDomsd(nn, level, sps)%bcData(mm)%Tp = zero
       flowDomsd(nn, level, sps)%bcData(mm)%area = zero
       flowDomsd(nn, level, sps)%BCData(mm)%uSlip = zero
       flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall = zero
       flowDomsd(nn, level, sps)%BCData(mm)%ptInlet = zero
       flowDomsd(nn, level, sps)%BCData(mm)%htInlet = zero
       flowDomsd(nn, level, sps)%BCData(mm)%ttInlet = zero
       flowDomsd(nn, level, sps)%BCData(mm)%turbInlet = zero
       flowDomsd(nn, level, sps)%BCData(mm)%ps = zero
    end do bocoLoop


    viscbocoLoop: do mm=1,flowDoms(nn, level, sps)%nViscBocos
       flowDomsd(nn, level, sps)%viscSubface(mm)%tau = zero
       flowDomsd(nn, level, sps)%viscSubface(mm)%q = zero
    end do viscbocoLoop

    ! For overset, the weights may be active in the comm structure. We
    ! need to zero them before we can accumulate.
    if (oversetPresent) then
       ! Pointers to the overset comms to make it easier to read
       sends: do i=1,commPatternOverset(level, sps)%nProcSend
          commPatternOverset(level, sps)%sendList(i)%interpd = zero
       end do sends
       internalOverset(level, sps)%donorInterpd = zero
    end if

    alphad = zero
    betad = zero
    machd = zero
    machGridd = zero
    machCoefd = zero
    pinfdimd = zero
    tinfdimd = zero
    rhoinfdimd = zero
    rgasdimd = zero
    pointrefd = zero
    prefd = zero
    rhoRefd = zero
    Trefd = zero
    murefd = zero
    urefd = zero
    hrefd = zero
    timerefd = zero
    pinfd = zero
    pinfCorrd = zero
    rhoinfd = zero
    uinfd = zero
    rgasd = zero
    muinfd = zero
    gammainfd = zero
    winfd = zero
    veldirfreestreamd = zero
    liftdirectiond = zero
    dragdirectiond = zero

    ! Zero all the reverse seeds in the dirichlet input arrays
    do iDom=1, cgnsNDom
       do iBoco=1, cgnsDoms(iDom)%nBocos
          if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)) then
             do iData=1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)
                if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)) then
                   do iDirichlet = 1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)
                      cgnsDomsd(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(:) = zero
                   end do
                end if
             end do
          end if
       end do
    end do

    ! And the reverse seeds in the actuator zones
    do i=1, nActuatorRegions
       actuatorRegionsd(i)%F = zero
       actuatorRegionsd(i)%T = zero
    end do

  end subroutine zeroADSeeds
  ! This is a special function that is sued to dealloc derivative values
  ! in blockpointers_d for use with the AD code.

  ! This routine setups "coloring" stencils for various sizes

  subroutine setup_PC_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k

    call setPointers(nn, level, 1)
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)
             flowDoms(nn, level, 1)%color(i, j, k) = &
                  mod(i + 5*j + 4*k, 7) + 1
          end do
       end do
    end do

    nColor = 7

  end subroutine setup_PC_coloring

  subroutine setup_dRdw_euler_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k

    call setPointers(nn, level, 1)
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)
             flowDoms(nn, level, 1)%color(i, j, k) = &
                  mod( i + 3*j + 4*k , 13) + 1
          end do
       end do
    end do

    nColor = 13

  end subroutine setup_dRdw_euler_coloring

  subroutine setup_dRdw_visc_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k

    call setPointers(nn, level, 1) ! Just to get the correct sizes
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)
             flowDoms(nn, level, 1)%color(i,j,k) = &
                  mod( i + 19*j + 11*k ,35) + 1
          end do
       end do
    end do

    nColor = 35

  end subroutine setup_dRdw_visc_coloring

  ! -------------------------------------------------------------
  !                   Debugging Color Colorings
  ! -------------------------------------------------------------

  subroutine setup_3x3x3_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! This is a dense 3x3x3 cube for debugging only
    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k, modi, modj, modk

    call setPointers(nn, level, 1)
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)
             modi = mod(i, 3)
             modj = mod(j, 3)
             modk = mod(k, 3)

             flowDoms(nn, level, 1)%color(i, j, k) = modi + 3*modj + 9*modk + 1

          end do
       end do
    end do

    nColor = 27
  end subroutine setup_3x3x3_coloring

  subroutine setup_5x5x5_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! This is a dense 5x5x5 cube for debugging only
    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k, modi, modj, modk

    call setPointers(nn, level, 1)
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)
             modi = mod(i, 5)
             modj = mod(j, 5)
             modk = mod(k, 5)

             flowDoms(nn, level, 1)%color(i, j, k) = modi + 5*modj + 25*modk + 1

          end do
       end do
    end do

    nColor = 125
  end subroutine setup_5x5x5_coloring

  subroutine setup_BF_coloring(nn, level, nColor)

    use constants
    use blockPointers, only : flowDoms, ib, jb, kb
    use utils, only : setPointers
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn, level

    ! Output parameters
    integer(kind=intTYpe), intent(out) :: nColor

    ! Working
    integer(kind=intType) :: i, j, k

    ! This is a REALLY brute force coloring for debugging

    call setPointers(nn, level, 1)
    !DIR$ NOVECTOR
    do k=0, kb
       do j=0, jb
          do i=0, ib
             ! Add the extra one for 1-based numbering (as opposed to zero-based)

             flowDoms(nn, level, 1)%color(i, j, k) = i + j*(ib+1) + k*((ib+1)*(jb+1)) + 1
          end do
       end do
    end do

    nColor = (ib+1)*(jb+1)*(kb+1)
  end subroutine setup_BF_coloring


  subroutine myMatCreate(matrix, blockSize, m, n, nnzDiagonal, nnzOffDiag, &
       file, line)
    ! Function to create petsc matrix to make stuff a little cleaner in
    ! the code above. Also, PETSc always thinks is a good idea to
    ! RANDOMLY change syntax between versions so this way there is only
    ! one place to make a change based on petsc version.

    use constants
    use communication, only : adflow_comm_world
    use utils, only : EChk, setPointers
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    Mat matrix
    integer(kind=intType), intent(in) :: blockSize, m, n
    integer(kind=intType), intent(in), dimension(*) :: nnzDiagonal, nnzOffDiag
    character*(*) :: file
    integer(kind=intType) :: ierr, line
    ! if (blockSize > 1) then
       call MatCreateBAIJ(ADFLOW_COMM_WORLD, blockSize, &
            m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
            0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
    ! else
       ! call MatCreateAIJ(ADFLOW_COMM_WORLD,&
            ! m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
            ! 0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
       call EChk(ierr, file, line)
    ! end if

    ! Warning: The array values is logically two-dimensional,
    ! containing the values that are to be inserted. By default the
    ! values are given in row major order, which is the opposite of
    ! the Fortran convention, meaning that the value to be put in row
    ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
    ! the insertion of values in column major order, one can call the
    ! command MatSetOption(Mat A, MAT COLUMN ORIENTED);

    call MatSetOption(matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine myMatCreate

  subroutine MyKSPMonitor(myKsp, n, rnorm, dummy, ierr)
    !
    !      This is a user-defined routine for monitoring the KSP
    !      iterative solvers. Instead of outputing the L2-norm at every
    !      iteration (default PETSc monitor), it only does it every
    !      'adjMonStep' iterations.
    !
    use ADjointPETSc
    use inputADjoint
    use communication
    implicit none
    !
    !     Subroutine arguments.
    !
    ! myKsp - Iterative context
    ! n     - Iteration number
    ! rnorm - 2-norm (preconditioned) residual value
    ! dummy - Optional user-defined monitor context (unused here)
    ! ierr  - Return error code

    real(kind=realType), pointer, dimension(:, :) :: myKsp
    integer(kind=intType) :: n, dummy, ierr
    real(kind=realType)   :: rnorm

    ! Write the residual norm to stdout every adjMonStep iterations.

    if(mod(n, adjMonStep) ==0 ) then
       if( myid==0 ) write(*, 10) n, rnorm
    end if

    ierr = 0

    ! Output format.

10  format(i4, 1x, 'KSP Residual norm', 1x, e16.10)

  end subroutine MyKSPMonitor

  subroutine setupStandardKSP(kspObject, kspObjectType, gmresRestart, preConSide, &
       globalPCType, ASMOverlap, globalPreConIts, localPCType, &
       localMatrixOrdering, localFillLevel, localPreConIts)

    ! This function sets up the supplied kspObject in the followin
    ! specific fashion. The reason this setup is in
    ! its own function is that it is used in the following places:
    ! 1. Setting up the preconditioner to use for the NKsolver
    ! 2. Setting up the preconditioner to use for the adjoint solver
    ! 3. Setting up the smoothers on the coarse multigrid levels.
    !
    ! The hierarchy of the setup is:
    !  kspObject --> Supplied KSP object
    !  |
    !  --> master_PC --> Preconditioner type set to KSP
    !      |
    !      --> master_PC_KSP --> KSP type set to Richardson with 'globalPreConIts'
    !          |
    !           --> globalPC --> PC type set to 'globalPCType'
    !               |            Usually Additive Schwartz and overlap is set
    !               |            with 'ASMOverlap'. Use 0 to get BlockJacobi
    !               |
    !               --> subKSP --> KSP type set to Richardon with 'LocalPreConIts'
    !                   |
    !                   --> subPC -->  PC type set to 'loclaPCType'.
    !                                  Usually ILU. 'localFillLevel' is
    !                                  set and 'localMatrixOrder' is used.
    !
    ! Note that if globalPreConIts=1 then maser_PC_KSP is NOT created and master_PC=globalPC
    ! and if localPreConIts=1 then subKSP is set to preOnly.
    use constants
    use utils, only : ECHk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Params
    KSP kspObject
    character(len=*), intent(in) :: kspObjectType, preConSide
    character(len=*), intent(in) :: globalPCType, localPCType
    character(len=*), intent(in) :: localMatrixOrdering
    integer(kind=intType), intent(in) :: ASMOverlap, localFillLevel, gmresRestart
    integer(kind=intType), intent(in) :: globalPreConIts, localPreConIts

    ! Working Variables
    PC  master_PC, globalPC, subpc
    KSP master_PC_KSP, subksp
    integer(kind=intType) :: nlocal, first, ierr


    ! First, KSPSetFromOptions MUST be called
    call KSPSetFromOptions(kspObject, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the type of solver to use:
    call KSPSetType(kspObject, kspObjectType, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! If we're using GMRES set the possible gmres restart
    call KSPGMRESSetRestart(kspObject, gmresRestart, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! If you're using GMRES, set refinement type
    call KSPGMRESSetCGSRefinementType(kspObject, &
         KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the preconditioner side from option:
    if (trim(preConSide) == 'right') then
       call KSPSetPCSide(kspObject, PC_RIGHT, ierr)
    else
       call KSPSetPCSide(kspObject, PC_LEFT, ierr)
    end if
    call EChk(ierr, __FILE__, __LINE__)

    if (trim(kspObjectType) == 'richardson') then
       call KSPSetPCSide(kspObject, PC_LEFT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Since there is an extraneous matMult required when using the
    ! richardson precondtiter with only 1 iteration, only use it we need
    ! to do more than 1 iteration.
    if (globalPreConIts > 1) then
       ! Extract preconditioning context for main KSP solver: (master_PC)
       call KSPGetPC(kspObject, master_PC, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the type of master_PC to ksp. This lets us do multiple
       ! iterations of preconditioner application
       call PCSetType(master_PC, 'ksp', ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Get the ksp context from master_PC which is the actual preconditioner:
       call PCKSPGetKSP(master_PC, master_PC_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! master_PC_KSP type will always be of type richardson. If the
       ! number  of iterations is set to 1, this ksp object is transparent.

       call KSPSetType(master_PC_KSP, 'richardson', ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call KSPMonitorSet(master_PC_KSP, MyKSPMonitor, PETSC_NULL_FUNCTION, &
            PETSC_NULL_FUNCTION, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Important to set the norm-type to None for efficiency.
       call kspsetnormtype(master_PC_KSP, KSP_NORM_NONE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Do one iteration of the outer ksp preconditioners. Note the
       ! tolerances are unsued since we have set KSP_NORM_NON
       call KSPSetTolerances(master_PC_KSP, PETSC_DEFAULT_REAL, &
            PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
            globalPreConIts, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Get the 'preconditioner for master_PC_KSP, called 'globalPC'. This
       ! preconditioner is potentially run multiple times.
       call KSPgetPC(master_PC_KSP, globalPC, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    else
       ! Just pull out the pc-object if we are not using kspRichardson
       call KSPGetPC(kspObject, globalPC, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Set the type of 'globalPC'. This will almost always be additive schwartz
    call PCSetType(globalPC, 'asm', ierr)!globalPCType, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the overlap required
    call PCASMSetOverlap(globalPC, ASMOverlap, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    !Setup the main ksp context before extracting the subdomains
    call KSPSetUp(kspObject, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Extract the ksp objects for each subdomain
    call PCASMGetSubKSP(globalPC, nlocal, first, subksp, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Since there is an extraneous matMult required when using the
    ! richardson precondtiter with only 1 iteration, only use it we need
    ! to do more than 1 iteration.
    if (localPreConIts > 1) then
       ! This 'subksp' object will ALSO be of type richardson so we can do
       ! multiple iterations on the sub-domains
       call KSPSetType(subksp, 'richardson', ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the number of iterations to do on local blocks. Tolerances are ignored.

       call KSPSetTolerances(subksp, PETSC_DEFAULT_REAL, &
            PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
            localPreConIts, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Again, norm_type is NONE since we don't want to check error
       call kspsetnormtype(subksp, KSP_NORM_NONE, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    else
       call KSPSetType(subksp, 'preonly', ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Extract the preconditioner for subksp object.
    call KSPGetPC(subksp, subpc, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! The subpc type will almost always be ILU
    call PCSetType(subpc, localPCType, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Setup the matrix ordering for the subpc object:
    call PCFactorSetMatOrderingtype(subpc, localMatrixOrdering, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the ILU parameters
    call PCFactorSetLevels(subpc, localFillLevel , ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine setupStandardKSP

  subroutine setupStandardMultigrid(kspObject, kspObjectType, gmresRestart, &
       preConSide, ASMoverlap, outerPreconIts, localMatrixOrdering, fillLevel)

    ! and if localPreConIts=1 then subKSP is set to preOnly.
    use constants
    use utils, only : ECHk
    use agmg, only : agmgOuterIts, agmgASMOverlap, agmgFillLevel, agmgMatrixOrdering, &
         setupShellPC, destroyShellPC, applyShellPC
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Params
    KSP kspObject
    character(len=*), intent(in) :: kspObjectType, preConSide
    character(len=*), intent(in) :: localMatrixOrdering
    integer(kind=intType), intent(in) :: ASMOverlap, fillLevel, gmresRestart
    integer(kind=intType), intent(in) :: outerPreconIts

    ! Working Variables
    PC  shellPC
    integer(kind=intType) :: ierr

    call KSPSetType(kspObject, kspObjectType, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the preconditioner side from option:
    if (trim(preConSide) == 'right') then
       call KSPSetPCSide(kspObject, PC_RIGHT, ierr)
    else
       call KSPSetPCSide(kspObject, PC_LEFT, ierr)
    end if
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGMRESSetRestart(kspObject, gmresRestart, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetPC(kspObject, shellPC, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PCSetType(shellPC, PCSHELL, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PCShellSetSetup(shellPC, setupShellPC, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PCShellSetDestroy(shellPC, destroyShellPC, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PCShellSetApply(shellPC, applyShellPC, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Just save the remaining pieces ofinformation in the agmg module.
    agmgOuterIts = outerPreConIts
    agmgASMOverlap = asmOverlap
    agmgFillLevel = fillLevel
    agmgMatrixOrdering = localMatrixOrdering
  end subroutine setupStandardMultigrid

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

  subroutine initializePETSc

    ! Call the C-version of the petsc initialize routine

    use ADjointPETSc, only : petsc_comm_world
    use communication, only : adflow_comm_world
    implicit none

    PETSC_COMM_WORLD= ADFLOW_COMM_WORLD
    call initPETScWrap()

  end subroutine initializePETSc

  subroutine finalizePETSc
    !
    !      Finalize PETSc by calling the appropriate routine
    !      PetscFinalize provided in the PETSc library. This
    !      automatically calls MPI_Finalize().
    !
    use ADjointPETSc, only : PETScIerr
    implicit none
    call PetscFinalize(PETScIerr)
  end subroutine finalizePETSc

subroutine statePreAllocation(onProc, offProc, wSize, stencil, N_stencil, &
     level, transposed)

  ! This is a generic function that determines the correct
  ! pre-allocation for on and off processor parts of the TRANSPOSED
  ! matrix. With overset, it is quite tricky to determine the
  ! transpose sparsity structure exactly, so we use an alternative
  ! approach. We proceed to determine the non-zeros of the untranposed
  ! matrix, but instead of assigning a non-zero the row we're looping
  ! over, we assign it to the column, which will become a row in the
  ! tranposed matrix. Since this requires communication we use a petsc
  ! vector for doing off processor values. This is not strictly
  ! correct since we will be using the real values as floats, but
  ! since the number of non-zeros per row is always going to be
  ! bounded, we don't have to worry about the integer/floating point
  ! conversions.

  use constants
  use blockPointers, only : nDom, il, jl, kl, fringes, flowDoms, globalCell, &
       iBlank, gInd
  use communication, only : adflow_comm_world
  use inputTimeSpectral , only : nTimeIntervalsSpectral
  use utils, only : setPointers, EChk
  use sorting, only : unique
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

  ! Subroutine Arguments
  integer(kind=intType), intent(in)  :: wSize
  integer(kind=intType), intent(in)  :: N_stencil
  integer(kind=intType), intent(in)  :: stencil(N_stencil, 3)
  integer(kind=intType), intent(out) :: onProc(wSize), offProc(wSize)
  integer(kind=intType), intent(in)  :: level
  logical, intent(in) :: transposed

  ! Local Variables
  integer(kind=intType) :: nn, i, j, k, sps, ii, jj, kk, iii, jjj, kkk, n, m, gc
  integer(kind=intType) :: iRowStart, iRowEnd, ierr, fInd
  integer(kind=intType), dimension((N_stencil-1)*8+1) :: cellBuffer, dummy
  Vec offProcVec
  logical :: overset
  real(kind=realType), pointer :: tmpPointer(:)


  call vecCreateMPI(adflow_comm_world, wSize, PETSC_DETERMINE, offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  ! Zero the cell movement counter
  ii = 0

  ! Set the onProc values for each cell to the number of "OFF" time
  ! spectral instances. The "on" spectral instances are accounted for
  ! in the stencil
  onProc(:) = nTimeIntervalsSpectral-1
  offProc(:) = 0
  ! Determine the range of onProc in dRdwT
  iRowStart = flowDoms(1, 1, 1)%globalCell(2,2,2)
  call setPointers(nDom, 1, nTimeIntervalsSpectral)
  iRowEnd   = flowDoms(nDom, 1, nTimeIntervalsSpectral)%globalCell(il, jl, kl)

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Loop over each real cell
        do k=2, kl
           do j=2, jl
              do i=2, il

                 ! Increment the running ii counter ONLY for each each
                 ! movement of center cell
                 ii = ii + 1

                 ! Reset the running tally of the number of neighbours
                 n = 0

                 blankedTest: if (iblank(i, j, k) == 1) then

                    ! Short-cut flag for cells without interpolated
                    ! cells in it's stencil
                    overset = .False.

                    ! Loop over the cells in the provided stencil:
                    do jj=1, N_stencil

                       ! Determine the cell we are dealing with
                       iii = stencil(jj, 1) + i
                       jjj = stencil(jj, 2) + j
                       kkk = stencil(jj, 3) + k

                       ! Index of the cell we are dealing with. Make
                       ! code easier to read
                       gc = globalCell(iii, jjj, kkk)

                       ! Check if the cell in question is a fringe or not:
                       if (iblank(iii, jjj, kkk) == 1) then
                          ! regular cell, add to our list, if it is
                          ! not a boundary
                          if (gc >= 0) then
                             n = n + 1
                             cellBuffer(n) = gc
                          end if

                       else if (iblank(iii, jjj, kkk) == -1) then
                          ! Fringe cell. What we do here is loop over
                          ! the donors for this cell and add any
                          ! entries that are real cells
                          overset = .True.
                          do kk=1,8
                             gc = gInd(kk, iii, jjj, kkk)
                             if (gc >= 0) then
                                n = n + 1
                                cellBuffer(n) = gc
                             end if
                          end do
                       end if
                    end do

                    ! We have now added 'n' cells to our buffer. For
                    ! the overset interpolation case, it is possible
                    ! (actually highly likely) that the same donor
                    ! cells are used in multiple fringes. To avoid
                    ! allocating more space than necessary, we
                    ! unique-ify the values, producing 'm' unique
                    ! values. If overset wasn't present, we can be
                    ! sure that m=n and we simply don't do the unique
                    ! operation.

                    if (overset) then
                       call unique(cellBuffer, n, m, dummy)
                    else
                       m = n
                    end if

                    ! -------------------- Non-transposed code ----------------
                    if (.not. transposed) then
                       ! Now we loop over the total number of
                       ! (unique) neighbours we have and assign them
                       ! to either an on-proc or an off-proc entry:
                       do jj=1, m
                          gc = cellBuffer(jj)

                          if (gc >= irowStart .and. gc <= iRowEnd) then
                             onProc(ii) = onProc(ii) + 1
                          else
                             offProc(ii) = offProc(ii) + 1
                          end if
                       end do
                    else
                       ! -------------------- Ttransposed code ----------------

                       ! Now we ALSO loop over the total number of
                       ! (unique) neighbours. However, instead of
                       ! adding to the non-zeros to the on/offproc for
                       ! row 'ii', we add them to the column index
                       ! which will be the row index for the
                       ! transposed matrix.
                       do jj=1, m
                          gc = cellBuffer(jj)

                          if (gc >= irowStart .and. gc <= iRowEnd) then
                             ! On processor values can be dealt with
                             ! directly since the diagonal part is square.
                             onProc(gc-iRowStart + 1) = onProc(gc-iRowStart+1)  +1
                          else
                             ! The offproc values need to be sent to
                             ! the other processors and summed.
                             call VecSetValue(offProcVec, gc, real(1), ADD_VALUES, ierr)
                             call EChk(ierr, __FILE__, __LINE__)
                          end if
                       end do
                    end if
                 else
                    ! Blanked and interpolated cells only need a single
                    ! non-zero per row for the identity on the diagonal.
                    onProc(ii) = onProc(ii) + 1
                 end if blankedTest
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

  ! Assemble the offproc vector. This doesn't take any time for the
  ! non-transposed operation.
  call VecAssemblyBegin(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (transposed) then
     ! Pull the local vector out and convert it back to integers.
     call VecGetArrayF90(offProcVec, tmpPointer, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     do i=1,wSize
        offProc(i) = int(tmpPointer(i) + half) ! Make sure, say 14.99999 is 15.
     end do

     call VecRestoreArrayF90(offProcVec, tmpPointer, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  ! Done with the temporary offProcVec
  call vecDestroy(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine statePreAllocation
  subroutine referenceShockSensor

    ! Compute the reference shock sensor for PC computations
    use constants
    use blockPointers, only : ib, jb, kb, il, jl, kl, ie, je, ke, shockSensor, &
         w, gamma, p, nDom, flowDoms
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputPhysics, only : equations
    use inputDiscretization, only : spaceDiscr
    use utils, only : setPointers, EChk
    implicit none

    ! Working variables
    integer(kind=intType) :: nn, level, sps, i, j, k

    level = 1

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)

          if (equations == EulerEquations .or. spaceDiscr == dissMatrix) then
             !shockSensor is Pressure
             do k=0, kb
                do j=0, jb
                   do i=0, ib
                      shockSensor(i,j,k) = P(i,j,k)
                   end do
                end do
             end do
          else
             ! Enthalpy is used instead
             do k=0, kb
                do j=2, jl
                   do i=2, il
                      shockSensor(i, j, k) = p(i, j, k)/(w(i, j, k, irho)**gamma(i, j, k))
                   enddo
                enddo
             enddo

             do k=2, kl
                do j=2, jl
                   shockSensor(0,  j, k) = p(0,  j, k)/(w(0,  j, k, irho)**gamma(0,  j, k))
                   shockSensor(1,  j, k) = p(1,  j, k)/(w(1,  j, k, irho)**gamma(1,  j, k))
                   shockSensor(ie, j, k) = p(ie, j, k)/(w(ie, j, k, irho)**gamma(ie, j, k))
                   shockSensor(ib, j, k) = p(ib, j, k)/(w(ib, j, k, irho)**gamma(ib, j, k))
                enddo
             enddo

             do k=2, kl
                do i=2, il
                   shockSensor(i, 0,  k) = p(i, 0,  k)/(w(i, 0,  k, irho)**gamma(i, 0,  k))
                   shockSensor(i, 1,  k) = p(i, 1,  k)/(w(i, 1,  k, irho)**gamma(i, 1,  k))
                   shockSensor(i, je, k) = p(i, je, k)/(w(i, je, k, irho)**gamma(i, je, k))
                   shockSensor(i, jb, k) = p(i, jb, k)/(w(i, jb, k, irho)**gamma(i, jb, k))
                enddo
             enddo
          end if
       end do
    end do
  end subroutine referenceShockSensor

  subroutine setFDReference(level)
    use constants
    use blockPointers, only : nDom, flowDoms, ib, jb, kb, il, jl, kl, &
         shockSensor, w, volRef, dw
    use inputPhysics, only : liftDirection, velDirFreeStream
    use flowVarRefState, only : nw, nwf
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : EChk, setPointers, getDirAngle
    use residuals, only : initRes_block
    use masterRoutines, only : block_res_state

    implicit none

    ! Input Parameters
    integer(kind=intType) :: level
    ! Working Parameters
    integer(kind=intType) :: i, j, k, l, nn, sps

    ! Compute the reference values for doing jacobian with FD
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral

          call setPointers(nn, level, sps)
          call block_res_state(nn, sps)
          ! Set the values
          do l=1, nw
             do k=0, kb
                do j=0, jb
                   do i=0, ib
                      flowdoms(nn, 1, sps)%wtmp(i,j,k,l)  = w(i, j, k, l)
                      flowdoms(nn, 1, sps)%dwtmp(i, j, k, l) = dw(i, j, k, l)
                   end do
                end do
             end do
          end do

          call initRes_block(1, nwf, nn, sps)

          ! Note: we have to divide by the volume for dwtmp2 since
          ! normally, dw would have been mulitpiled by 1/Vol in block_res_state

          do l=1, nw
             do k=0, kb
                do j=0, jb
                   do i=0, ib
                      flowdoms(nn, 1, sps)%dwtmp2(i, j, k, l) = &
                           dw(i, j, k, l)/volRef(i, j, k)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine setFDReference

  subroutine resetFDReference(level)

    use constants
    use blockPointers, only : nDom, flowDoms, ib, jb, kb, w, dw
    use flowVarRefState, only : nw, nwf
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    implicit none

    ! Input Parameters
    integer(kind=intType) :: level

    ! Working Parameters
    integer(kind=intType) :: i, j, k, l, nn, sps
    real(kind=realType) :: sepSensor, Cavitation, axisMoment

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          ! Reset w and dw
          do l=1, nw
             do k=0, kb
                do j=0, jb
                   do i=0, ib
                      w(i, j, k, l) = flowdoms(nn, 1, sps)%wtmp(i, j, k, l)
                      dw(i, j, k, l) = flowdoms(nn, 1, sps)%dwtmp(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine resetFDReference

  subroutine setDiffSizes
    !
    !       This routine set the sizes for the pointers that will be
    !       used in the forward debug mode and reverse mode AD.
    !
    use constants
    use blockPointers, only : flowDoms, ib, jb, kb, ie, je, ke, ib, jb, ke, &
         nBocos, nViscBocos, nDom
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputPhysics, only : equations
    use diffSizes
    implicit none

    ! local variables
    integer(kind=intType) :: nLevels

    ! Compute nlevels
    nLevels = ubound(flowDoms, 2)


    ! Set the size for dynamic pointers to zero for debug purpose
    ! bcdata%norm
    ISIZE1OFDrfDrfbcdata_norm = 0
    ISIZE2OFDrfDrfbcdata_norm = 0
    ISIZE3OFDrfDrfbcdata_norm = 0

    ! bcdata%rface
    ISIZE1OFDrfDrfbcdata_rface = 0
    ISIZE2OFDrfDrfbcdata_rface = 0

    ! bcdata%m
    ISIZE1OFDrfDrfbcdata_m = 0
    ISIZE2OFDrfDrfbcdata_m = 0
    ISIZE3OFDrfDrfbcdata_m = 0

    ! bcdata%fp
    ISIZE1OFDrfDrfbcdata_fp = 0
    ISIZE2OFDrfDrfbcdata_fp = 0
    ISIZE3OFDrfDrfbcdata_fp = 0

    ! bcdata%fv
    ISIZE1OFDrfDrfbcdata_fv = 0
    ISIZE2OFDrfDrfbcdata_fv = 0
    ISIZE3OFDrfDrfbcdata_fv = 0

    ! bcdata%fp
    ISIZE1OFDrfDrfbcdata_oarea = 0
    ISIZE2OFDrfDrfbcdata_oarea = 0
    ISIZE3OFDrfDrfbcdata_oarea = 0

    ! sepSensor
    ISIZE1OFDrfDrfbcdata_sepSensor = 0
    ISIZE2OFDrfDrfbcdata_sepSensor = 0

    ! Cavitation
    ISIZE1OFDrfDrfbcdata_Cavitation = 0
    ISIZE2OFDrfDrfbcdata_Cavitation = 0

    ! AxisMoment
    ISIZE1OFDrfDrfbcdata_axisMoment = 0
    ISIZE2OFDrfDrfbcdata_axisMoment = 0

    ! viscsubface%tau
    ISIZE1OFDrfDrfviscsubface_tau = 0
    ISIZE2OFDrfDrfviscsubface_tau = 0
    ISIZE3OFDrfDrfviscsubface_tau = 0

    ! prod
    ISIZE1OFDrfprod = 0
    ISIZE2OFDrfprod = 0
    ISIZE3OFDrfprod = 0

    ! vort
    ISIZE1OFDrfvort = 0
    ISIZE2OFDrfvort = 0
    ISIZE3OFDrfvort = 0

    ! dvt
    ISIZE1OFDrfdvt = 0
    ISIZE2OFDrfdvt = 0
    ISIZE3OFDrfdvt = 0
    ISIZE4OFDrfdvt = 0

    ! vol
    ISIZE1OFDrfvol = 0
    ISIZE2OFDrfvol = 0
    ISIZE3OFDrfvol = 0

    ! rho,etot,u,v,w,p,k
    ISIZE1OFrho = 0
    ISIZE1OFetot = 0
    ISIZE1OFu = 0
    ISIZE1OFv = 0
    ISIZE1OFw = 0
    ISIZE1OFp = 0
    ISIZE1OFk = 0

    ! Du1, Du2, Du3
    ISIZE1OFDu1 = 0
    ISIZE1OFDu2 = 0
    ISIZE1OFDu3 = 0

    ! Left, Right, Flux
    ISIZE1OFLeft = 0
    ISIZE1OFRight = 0
    ISIZE1OFFlux = 0

    ! bcdata
    ISIZE1OFDrfbcdata = 0!nbocos

    ! s
    ISIZE1OFDrfs = 0!ie
    ISIZE2OFDrfs = je
    ISIZE3OFDrfs = ke
    ISIZE4OFDrfs = 3

    ! sfacei
    ISIZE3OFDrfsfaceI = 0!ie + 1
    ISIZE2OFDrfsfaceI = je
    ISIZE1OFDrfsfaceI = ke

    ! sfacej
    ISIZE3OFDrfsfaceJ = 0!ie
    ISIZE2OFDrfsfaceJ = je + 1
    ISIZE1OFDrfsfaceJ = ke

    ! sfacek
    ISIZE3OFDrfsfaceK = 0!ie
    ISIZE2OFDrfsfaceK = je
    ISIZE1OFDrfsfaceK = ke + 1

    ! Define size for the pointers
    ! flowdoms
    ISIZE1OFDrfflowdoms = nDom
    ISIZE2OFDrfflowdoms = nLevels
    ISIZE3OFDrfflowdoms = nTimeIntervalsSpectral

    !viscSubface
    ISIZE1OFDrfviscsubface = nViscBocos
    ISIZE1OFDrfflowdoms_bcdata = nBocos

    ! x
    ISIZE4OFDrfx = 3
    ISIZE1OFDrfx = ie + 1
    ISIZE2OFDrfx = je + 1
    ISIZE3OFDrfx = ke + 1

    ! flowdoms_x
    ISIZE4OFDRFFLOWDOMS_X = 3
    ISIZE1OFDRFFLOWDOMS_X = ie + 1
    ISIZE2OFDRFFLOWDOMS_X = je + 1
    ISIZE3OFDRFFLOWDOMS_X = ke + 1

    if ( equations == RANSEquations ) then
       ! rev
       ISIZE1OFDrfrev = ib + 1
       ISIZE2OFDrfrev = jb + 1
       ISIZE3OFDrfrev = kb + 1
    else
       ! rev
       ISIZE1OFDrfrev = 0
       ISIZE2OFDrfrev = 0
       ISIZE3OFDrfrev = 0
    end if

    ! rlv
    ISIZE1OFDrfrlv = ib + 1
    ISIZE2OFDrfrlv = jb + 1
    ISIZE3OFDrfrlv = kb + 1

    ! w
    ISIZE4OFDrfw = nw
    ISIZE1OFDrfw = ib + 1
    ISIZE2OFDrfw = jb + 1
    ISIZE3OFDrfw = kb + 1

    ! flowdoms_x
    ISIZE4OFDRFFLOWDOMS_W = nw
    ISIZE1OFDRFFLOWDOMS_W = ib + 1
    ISIZE2OFDRFFLOWDOMS_W = jb + 1
    ISIZE3OFDRFFLOWDOMS_W = kb + 1

    ! flowdoms_dw
    ISIZE4OFDRFFLOWDOMS_dw = nw
    ISIZE1OFDRFFLOWDOMS_dw = ib + 1
    ISIZE2OFDRFFLOWDOMS_dw = jb + 1
    ISIZE3OFDRFFLOWDOMS_dw = kb + 1

    ! flowdoms_vol
    ISIZE1OFDRFFLOWDOMS_vol = ib + 1
    ISIZE2OFDRFFLOWDOMS_vol = jb + 1
    ISIZE3OFDRFFLOWDOMS_vol = kb + 1

    ! fw
    ISIZE4OFDrffw = nwf
    ISIZE1OFDrffw = ib + 1
    ISIZE2OFDrffw = jb + 1
    ISIZE3OFDrffw = kb + 1

    ! dw
    ISIZE4OFDrfdw = nw
    ISIZE1OFDrfdw = ib + 1
    ISIZE2OFDrfdw = jb + 1
    ISIZE3OFDrfdw = kb + 1

    ! p
    ISIZE1OFDrfp = ib + 1
    ISIZE2OFDrfp = jb + 1
    ISIZE3OFDrfp = kb + 1

    ! gamma
    ISIZE1OFDrfgamma = ib + 1
    ISIZE2OFDrfgamma = jb + 1
    ISIZE3OFDrfgamma = kb + 1

    ! radI
    ISIZE1OFDrfradI = ie
    ISIZE2OFDrfradI = je
    ISIZE3OFDrfradI = ke

    ! radJ
    ISIZE1OFDrfradJ = ie
    ISIZE2OFDrfradJ = je
    ISIZE3OFDrfradJ = ke

    ! radK
    ISIZE1OFDrfradK = ie
    ISIZE2OFDrfradK = je
    ISIZE3OFDrfradK = ke

    ! sI
    ISIZE1OFDRFsI = ie + 1
    ISIZE2OFDRFsI = je
    ISIZE3OFDRFsI = ke
    ISIZE4OFDRFsI = 3

    ! sJ
    ISIZE1OFDRFsJ = ie
    ISIZE2OFDRFsJ = je + 1
    ISIZE3OFDRFsJ = ke
    ISIZE4OFDRFsJ = 3

    ! sK
    ISIZE1OFDRFsK = ie
    ISIZE2OFDRFsK = je
    ISIZE3OFDRFsK = ke + 1
    ISIZE4OFDRFsK = 3

    !bmti1
    ISIZE1OFDrfbmti1 = je
    ISIZE2OFDrfbmti1 = ke
    ISIZE3OFDrfbmti1 = nt2 - nt1 + 1
    ISIZE4OFDrfbmti1 = nt2 - nt1 + 1

    !bmti2
    ISIZE1OFDrfbmti2 = je
    ISIZE2OFDrfbmti2 = ke
    ISIZE3OFDrfbmti2 = nt2 - nt1 + 1
    ISIZE4OFDrfbmti2 = nt2 - nt1 + 1

    !bmtj1
    ISIZE1OFDrfbmtj1 = ie
    ISIZE2OFDrfbmtj1 = ke
    ISIZE3OFDrfbmtj1 = nt2 - nt1 + 1
    ISIZE4OFDrfbmtj1 = nt2 - nt1 + 1

    !bmtj2
    ISIZE1OFDrfbmtj2 = ie
    ISIZE2OFDrfbmtj2 = ke
    ISIZE3OFDrfbmtj2 = nt2 - nt1 + 1
    ISIZE4OFDrfbmtj2 = nt2 - nt1 + 1

    !bmtk1
    ISIZE1OFDrfbmtk1 = ie
    ISIZE2OFDrfbmtk1 = je
    ISIZE3OFDrfbmtk1 = nt2 - nt1 + 1
    ISIZE4OFDrfbmtk1 = nt2 - nt1 + 1

    !bmtk2
    ISIZE1OFDrfbmtk2 = ie
    ISIZE2OFDrfbmtk2 = je
    ISIZE3OFDrfbmtk2 = nt2 - nt1 + 1
    ISIZE4OFDrfbmtk2 = nt2 - nt1 + 1

    !bvti1
    ISIZE1OFDrfbvti1 = je
    ISIZE2OFDrfbvti1 = ke
    ISIZE3OFDrfbvti1 = nt2 - nt1 + 1

    !bvti2
    ISIZE1OFDrfbvti2 = je
    ISIZE2OFDrfbvti2 = ke
    ISIZE3OFDrfbvti2 = nt2 - nt1 + 1
    !bvti1
    ISIZE1OFDrfbvti1 = je
    ISIZE2OFDrfbvti1 = ke
    ISIZE3OFDrfbvti1 = nt2 - nt1 + 1

    !bvti2
    ISIZE1OFDrfbvti2 = je
    ISIZE2OFDrfbvti2 = ke
    ISIZE3OFDrfbvti2 = nt2 - nt1 + 1

    !bvtk1
    ISIZE1OFDrfbvtk1 = ie
    ISIZE2OFDrfbvtk1 = je
    ISIZE3OFDrfbvtk1 = nt2 - nt1 + 1

    !bvtk2
    ISIZE1OFDrfbvtk2 = ie
    ISIZE2OFDrfbvtk2 = je
    ISIZE3OFDrfbvtk2 = nt2 - nt1 + 1

  end subroutine setDiffSizes

end module adjointUtils
