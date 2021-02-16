! This module contains the majoriy of the implemention of a
! preconditioner based on our own memory managment and matrix vector
! operation. The idea is that this could be used to form a DADI based
! preconditioner for the NK/adjoint systems.  It is experiemental and
! not production ready at all.

module fortranPC
  contains
  subroutine allocPCMem(level)

    ! This routine allocates memory for the fortran-based PC. It is
    ! currently not used anywhere, but it become useful in the future.

    use constants
    use blockPointers, only : nDom, nx, ny, nz, il, jl, kl, ie, je, ke, flowDoms
    use inputTimeSpectral, onlY : nTimeIntervalsSpectral
    use flowVarRefState, only : nw
    use utils, only : EChk, setPointers
    implicit none

    integer(kind=intType), intent(in) :: level
    integer(kind=intType) :: nn, sps


    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)

          if (.not. associated(flowDoms(nn, level, sps)%PCMat)) then

             allocate(flowDoms(nn, level, sps)%PCMat(nw, 7*nw, 2:il, 2:jl, 2:kl))

             ! I direciton data
             allocate(flowDoms(nn, level, sps)%i_D_Fact(nw, nx*nw, 2:jl, 2:kl))
             allocate(flowDoms(nn, level, sps)%i_L_Fact(nw, (nx-1)*nw, 2:jl, 2:kl))
             allocate(flowDoms(nn, level, sps)%i_U_Fact(nw, (nx-1)*nw, 2:jl, 2:kl))
             allocate(flowDoms(nn, level, sps)%i_U2_Fact(nw, (nx-2)*nw, 2:jl, 2:kl))

             ! J direction data
             allocate(flowDoms(nn, level, sps)%j_D_Fact(nw, ny*nw, 2:il, 2:kl))
             allocate(flowDoms(nn, level, sps)%j_L_Fact(nw, (ny-1)*nw, 2:il, 2:kl))
             allocate(flowDoms(nn, level, sps)%j_U_Fact(nw, (ny-1)*nw, 2:il, 2:kl))
             allocate(flowDoms(nn, level, sps)%j_U2_Fact(nw, (ny-2)*nw, 2:il, 2:kl))

             ! K direciton data
             allocate(flowDoms(nn, level, sps)%k_D_Fact(nw, nz*nw, 2:il, 2:jl))
             allocate(flowDoms(nn, level, sps)%k_L_Fact(nw, (nz-1)*nw, 2:il, 2:jl))
             allocate(flowDoms(nn, level, sps)%k_U_Fact(nw, (nz-1)*nw, 2:il, 2:jl))
             allocate(flowDoms(nn, level, sps)%k_U2_Fact(nw, (nz-2)*nw, 2:il, 2:jl))

             ! iPIV arrays
             allocate(flowDoms(nn, level, sps)%i_ipiv(nw, nx, 2:jl, 2:kl))
             allocate(flowDoms(nn, level, sps)%j_ipiv(nw, ny, 2:il, 2:kl))
             allocate(flowDoms(nn, level, sps)%k_ipiv(nw, nz, 2:il, 2:jl))

             ! Vectors
             allocate(flowDoms(nn, level, sps)%pcVec1(nw, 1:ie, 1:je, 1:ke))
             allocate(flowDoms(nn, level, sps)%pcVec2(nw, 1:ie, 1:je, 1:ke))

          end if
       end do
    end do
  end subroutine allocPCMem


  subroutine setupPCMatrix(useAD,  useTranspose, frozenTurb, level)
#ifndef USE_NO_PETSC

    ! This routine generates a fortran form of the PCmatrix. It is
    ! currently not used anywhere, but it become useful in the future.
    use block, only : flowDomsd, flowDoms
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
    use utils, only : setPointers, EChk, getDirAngle, setPointers_d
    use haloExchange, only : whalo2
    use masterRoutines, only : block_res_state
    use adjointUtils, only : zeroADSeeds
#ifndef USE_COMPLEX
    use masterRoutines, only : block_res_state_d
#endif
    implicit none

    ! Input Variables
    logical, intent(in) :: useAD,  useTranspose, frozenTurb
    integer(kind=intType), intent(in) :: level

    ! Local variables.
    integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
    integer(kind=intType) :: nColor, iColor, jColor, irow, icol, fmDim, frow
    integer(kind=intType) :: nTransfer, nState, tmp, icount
    integer(kind=intType) :: n_stencil, i_stencil, ind1, orderturbsave
    integer(kind=intType), dimension(:, :), pointer :: stencil
    real(kind=realType) :: delta_x, one_over_dx
    real(kind=realType) :: delta_x_turb, one_over_dx_turb

#ifdef USE_COMPLEX
    complex(kind=realType), dimension(:,:), allocatable :: blk
#else
    real(kind=realType), dimension(:,:), allocatable :: blk
#endif
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, mm, colInd
    logical :: resetToRANS,  splitMat
    real :: val

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurb ) then
       nState = nwf
    else
       nState = nw
    endif

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
    tinfdimd = zero
    rhoinfdimd = zero


    ! Set a pointer to the correct set of stencil depending on if we are
    ! using the first order stencil or the full jacobian

    stencil => euler_pc_stencil
    n_stencil = N_euler_pc

    ! Very important to use only Second-Order dissipation for PC
    lumpedDiss=.True.
    ! also use first order advection terms for turbulence
    orderturbsave = orderturb
    orderturb = firstOrder

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


    delta_x_turb = 1e-14
    one_over_dx_turb = one/delta_x_turb

    rkStage = 0

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False.
    if (frozenTurb .and. equations == RANSEquations) then
       equations = NSEquations
       resetToRANS = .True.
    end if

    ! Allocate the additional memory we need for doing forward mode AD
    !  derivatives and copy any required reference values:

    call allocPCMem(level)
    if (.not. derivVarsAllocated .and. useAD) then
       call alloc_derivative_values(level)
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

          ! Allocate temporary space only needed while assembling.
          allocate(flowDoms(nn, 1, sps)%dw_deriv(2:il, 2:jl, 2:kl, 1:nw, 1:nw), stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          allocate(flowDoms(nn, 1, sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          allocate(flowDoms(nn, 1, sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          allocate(flowDoms(nn, 1, sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw),stat=ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Only need 1 set of colors on the first sps instance.
          if (sps == 1) then
             allocate(flowDoms(nn, 1, 1)%color(0:ib, 0:jb, 0:kb), stat=ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end if

          ! Zero the matrix
          pcMat = zero

       end do
    end do

    ! For the PC we don't linearize the shock sensor so it must be
    ! computed here.
    call referenceShockSensor

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

       call setup_PC_coloring(nn, level,  nColor) ! Euler Colorings

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
             stateLoop: do l=1, nState

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

                ! Peturb w or set AD Seed according to iColor
                do k=0, kb
                   do j=0, jb
                      do i=0, ib
                         if (flowdoms(nn, 1, 1)%color(i, j, k) == icolor) then
                            if (useAD) then
                               flowDomsd(nn, 1, sps)%w(i, j, k, l) = one
                            else
                               if (l <= nwf) then
                                  w(i, j, k, l) = w(i, j, k, l) + delta_x
                               else
                                  w(i, j, k, l) = w(i, j, k, l) + delta_x_turb
                               end if
                            end if
                         end if
                      end do
                   end do
                end do

                ! Run Block-based residual
                if (useAD) then
#ifndef USE_COMPLEX
                   call block_res_state_d(nn, sps)
#else
                   print *, 'Forward AD routines are not complexified'
                   stop
#endif
                else
                   call block_res(nn, sps)
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
                                  flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                       flowdomsd(nn, 1, sps2)%dw(i, j, k, ll)
                               else
                                  if (sps2 == sps) then
                                     ! If the peturbation is on this
                                     ! instance, we've computed the spatial
                                     ! contribution so subtrace dwtmp

                                     if (l <= nwf) then
                                        flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                             one_over_dx * &
                                             (flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                             flowDoms(nn, 1, sps2)%dwtmp(i, j, k, ll))
                                     else
                                        flowDoms(nn, 1, sps2)%dw_deriv(i, j, k, ll, l) = &
                                             one_over_dx_turb * &
                                             (flowDoms(nn, 1, sps2)%dw(i, j, k, ll) - &
                                             flowDoms(nn, 1, sps2)%dwtmp(i, j, k, ll))
                                     end if
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
                      if (flowdoms(nn, 1, 1)%color(i, j, k) == icolor .and. globalCell(i,j,k)>=0) then

                         ! Diagonal block is easy.
                         if (onBlock(i, j, k)) then
                            PCMat(:, 1:nw, i, j, k) = flowDoms(nn, 1, sps)%dw_deriv(i, j, k, 1:nstate, 1:nstate)
                         end if


                         if (onBlock(i-1, j, k)) then
                            PCMat(:, 2*nw+1:3*nw, i-1, j, k) = flowDoms(nn, 1, sps)%dw_deriv(i-1, j, k, 1:nstate, 1:nstate)
                         end if

                         if (onBlock(i+1, j, k)) then
                            PCMat(:, 1*nw+1:2*nw, i+1, j, k) = flowDoms(nn, 1, sps)%dw_deriv(i+1, j, k, 1:nstate, 1:nstate)
                         end if

                         if (onBlock(i, j-1, k)) then
                            PCMat(:, 4*nw+1:5*nw, i, j-1, k) = flowDoms(nn, 1, sps)%dw_deriv(i, j-1, k, 1:nstate, 1:nstate)
                         end if

                         if (onBlock(i, j+1, k)) then
                            PCMat(:, 3*nw+1:4*nw, i, j+1, k) = flowDoms(nn, 1, sps)%dw_deriv(i, j+1, k, 1:nstate, 1:nstate)
                         end if

                         if (onBlock(i, j, k-1)) then
                            PCMat(:, 6*nw+1:7*nw, i, j, k-1) = flowDoms(nn, 1, sps)%dw_deriv(i, j, k-1, 1:nstate, 1:nstate)
                         end if

                         if (onBlock(i, j, k+1)) then
                            PCMat(:, 5*nw+1:6*nw, i, j, k+1) = flowDoms(nn, 1, sps)%dw_deriv(i, j, k+1, 1:nstate, 1:nstate)
                         end if
                      end if
                   end do iLoop
                end do jLoop
             end do kLoop
          end do colorLoop
       end do spectralLoop
    end do domainLoopAD

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
          if (sps == 1) then
             deallocate(flowDoms(nn, 1, sps)%color)
          end if
       end do
    end do

    ! Return dissipation Parameters to normal -> VERY VERY IMPORTANT
    lumpedDiss = .False.
    orderturb = orderturbsave

    ! Reset the correct equation parameters if we were useing the frozen
    ! Turbulent
    if (resetToRANS) then
       equations = RANSEquations
    end if

    deallocate(blk)
#endif
  contains
    function onBlock(i, j, k)

      use precision
      implicit none

      integer(kind=intType), intent(in) :: i, j, k
      logical :: onBlock

      if (i >= 2 .and. i <= il .and. j >= 2 .and. j<= jl .and. k >= 2 .and. k <= kl) then
         onBlock = .True.
      else
         onBlock = .False.
      end if

    end function onBlock

    function getIndex(i, j, k)

      use precision
      implicit none

      integer(kind=intType), intent(in) :: i, j, k
      integer(kind=intType) :: getIndex

      getIndex = ((k-2)*nx*ny + (j-2)*nx + (i-2))*nw  + 1

    end function getIndex

  end subroutine setupPCMatrix
  subroutine testpc()

    ! use constants
    ! !use communication
    ! !use adjointPETSc
    ! use adjointVars
    ! use flowvarrefstate
    ! use inputtimespectral
    ! use blockPointers
    ! use utils, only : Echk
    ! use adjointUtils, only : setupStateResidualMatrix
    ! implicit none

    ! real(kind=alwaysRealType) , dimension(:), allocatable :: v1, v2
    ! real(kind=realType) :: timeA, timeB, timeC, val
    ! integer(kind=intType) :: nDimw, ierr, i, nn, sps, j, k, ii, jj, kk, info
    ! logical :: useAD

    ! useAD = .True.
    ! call setupPCMatrix(useAD, .False., .False., 1)

    ! ! Now create all the factors
    ! call createPETscVars
    ! call setupStateResidualMatrix(drdwpret, usead, .True. , .False., &
    !      .False., .False., 1)

    ! nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
    ! allocate(v1(nDimw), v2(nDimW))

    ! ! Set random into V1
    ! call random_number(v1)

    ! ! Dump psi into psi_like1 and RHS into psi_like2
    ! call VecPlaceArray(psi_like1, v1, ierr)
    ! call EChk(ierr,__FILE__,__LINE__)

    ! call VecPlaceArray(psi_like2, v2, ierr)
    ! call EChk(ierr,__FILE__,__LINE__)

    ! call mpi_barrier(adflow_comm_world, ierr)
    ! timeA = mpi_wtime()
    ! call matMult(drdwpret, psi_like1, psi_like2, ierr)
    ! call EChk(ierr,__FILE__,__LINE__)
    ! timeB =mpi_wtime()

    ! print *,'Petsc Time:', myid, timeB-timeA
    ! call vecNorm(psi_like2, NORM_2, val, ierr)
    ! print *,'PETsc Norm:', val


    ! call mpi_barrier(adflow_comm_world, ierr)
    ! timeA = mpi_wtime()
    ! call PCMatMult(drdwpret, psi_like1, psi_like2, ierr)
    ! timeB = mpi_wtime()

    ! print *,'My Time:', myid, timeB-timeA
    ! call vecNorm(psi_like2, NORM_2, val, ierr)
    ! print *,'My Norm:', val
    ! call VecResetArray(psi_like1, ierr)
    ! call VecResetArray(psi_like2, ierr)
  end subroutine testpc

  subroutine factorPCMatrix()

    use communication
    use adjointPETSc
    use adjointVars
    use flowvarrefstate
    use inputtimespectral
    use blockPointers
    use utils, only : Echk, setPointers
    implicit none

    ! integer(kind=intType) :: ierr, i, nn, sps, j, k, ii, jj, kk, info, timeA, timeB

    ! timeA = mpi_wtime()
    ! do nn=1,nDom
    !    do sps=1,nTimeIntervalsSpectral
    !       call setPointers(nn, 1, sps)

    !       ! ========= I-Lines ============
    !       do k=2, kl
    !          do j=2, jl

    !             ! Copy data from PCMat
    !             ii = 0
    !             jj = 0
    !             kk = 0
    !             do i=2, il
    !                i_D_fact(:, nw*ii+1:nw*(ii+1), j, k) = PCMat(:,      1:1*nw, i, j, k)
    !                ii = ii + 1

    !                if (i > 2) then
    !                   i_L_fact(:, nw*jj+1:nw*(jj+1), j, k) = PCMat(:, 1*nw+1:2*nw, i, j, k)
    !                   jj = jj + 1
    !                end if

    !                if (i < il) then
    !                   i_U_fact(:, nw*kk+1:nw*(kk+1), j, k) = PCMat(:, 2*nw+1:3*nw, i, j, k)
    !                   kk = kk + 1
    !                end if

    !             end do

    !             ! ! Perform factorization
    !             ! call dgeblttrf(nx, nw, i_D_fact(:, :, j, k), i_L_fact(:, :, j, k), &
    !             !      i_U_fact(:, :, j, k), i_U2_fact(:, :, j, k), i_ipiv(:, :, j, k), info)

    !          end do
    !       end do

    !       ! ========= J-Lines ============
    !       do k=2, kl
    !          do i=2, il

    !             ! Copy data from PCMat
    !             ii = 0
    !             jj = 0
    !             kk = 0
    !             do j=2, jl
    !                j_D_fact(:, nw*ii+1:nw*(ii+1), i, k) = PCMat(:,      1:1*nw, i, j, k)
    !                ii = ii + 1

    !                if (j > 2) then
    !                   j_L_fact(:, nw*jj+1:nw*(jj+1), i, k) = PCMat(:, 3*nw+1:4*nw, i, j, k)
    !                   jj = jj + 1
    !                end if

    !                if (j < jl) then
    !                   j_U_fact(:, nw*kk+1:nw*(kk+1), i, k) = PCMat(:, 4*nw+1:5*nw, i, j, k)
    !                   kk =kk + 1
    !                end if
    !             end do

    !             ! ! Perform factorization
    !             ! call dgeblttrf(ny, nw, j_D_fact(:, :, i, k), j_L_fact(:, :, i, k), &
    !             !      j_U_fact(:, :, i, k), j_U2_fact(:, :, i, k), j_ipiv(:, :, i, k), info)

    !          end do
    !       end do

    !       ! ========= k-Lines ============
    !       do j=2, jl
    !          do i=2, il

    !             ! Copy data from PCMat
    !             ii = 0
    !             jj = 0
    !             kk = 0
    !             do k=2, kl
    !                k_D_fact(:, nw*ii+1:nw*(ii+1), i, j) = PCMat(:,      1:1*nw, i, j, k)
    !                ii = ii + 1

    !                if (k > 2) then
    !                   k_L_fact(:, nw*jj+1:nw*(jj+1), i, j) = PCMat(:, 5*nw+1:6*nw, i, j, k)
    !                   jj = jj + 1
    !                end if

    !                if (k < kl) then
    !                   k_U_fact(:, nw*kk+1:nw*(kk+1), i, j) = PCMat(:, 6*nw+1:7*nw, i, j, k)
    !                   kk =kk + 1
    !                end if
    !             end do

    !             ! ! Perform factorization
    !             ! call dgeblttrf(nz, nw, k_D_fact(:, :, i, j), k_L_fact(:, :, i,j), &
    !             !      k_U_fact(:, :, i, j), k_U2_fact(:, :, i, j), k_ipiv(:, :, i, j), info)

    !          end do
    !       end do
    !    end do
    ! end do

  end subroutine factorPCMatrix

  subroutine PCMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! PCMat with a vector.

    use constants
    use communication
    use blockPointers
    use iteration
    use flowVarRefState
    use inputAdjoint
    use ADjointVars
    use inputTimeSpectral
    use utils, only : EChk, setPointers
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr, i, j, k, l, nn, sps, ii
    real(kind=realType) :: sum1, sum2, sum3, sum4, sum5
    real(kind=realType) :: x1, x2, x3, x4, x5
    real(kind=realType), pointer :: yPtr(:)
    real(kind=realType) ::  RHS(nw*7)
    ! We first have to distribute xPtr into the halo cells for the local
    ! products.
    print *,'calling pcmatmult'
    call setPCVec(vecX)

    call VecGetArrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Now we can compute the acutal matrix vector product.
    ii = 1
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, 1, sps)

          do k=2, kl
             do j=2, jl
                do i=2, il

                   ! Fill up the RHS

                   rhs(0*nw+1:1*nw) = PCVec1(:, i  , j, k)
                   rhs(1*nw+1:2*nw) = PCVec1(:, i-1, j, k)
                   rhs(2*nw+1:3*nw) = PCVec1(:, i+1, j, k)
                   rhs(3*nw+1:4*nw) = PCVec1(:, i, j-1, k)
                   rhs(4*nw+1:5*nw) = PCVec1(:, i, j+1, k)
                   rhs(5*nw+1:6*nw) = PCVec1(:, i, j, k-1)
                   rhs(6*nw+1:7*nw) = PCVec1(:, i, j, k+1)

                   ! Call blass mat-vec. We can dump the result directly
                   ! into yPtr. There doesn't appear to be much
                   ! difference between the fortran matmul and blas for
                   ! these sized operations.

                   !yPtr(ii:ii+nw) = matmul(PCMat(:, :, i, j, k), rhs)
                   call DGEMV('n', nw, nw*7, one, PCMat(:, :, i, j, k), nw, rhs, 1, zero, yPtr(ii), 1)

                   ii = ii + nw

                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ierr = 0

  end subroutine PCMatMult

  subroutine myShellPCApply(pc, vecX, vecY, ierr)

    use precision
    use blockPointers
    use inputTimeSpectral
    use flowvarrefstate
    use communication
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Arguments
    PC   pc
    Vec   vecX, vecY
    integer(kind=intType) :: ierr, info, i, j, k, ii, nn, sps, ipiv(nw)
    real(kind=realType) :: blk(nw, nw), rhs(nw)
    real(kind=realType), pointer :: yPtr(:)


    ! First copy X to Y. This way we will continually transform vecY
    ! into the preconditioned vector.
    call vecCopy(vecX, vecY, ierr)

    call VecGetArrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Do something useful here....

    call VecRestoreArrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine myShellPCApply

  subroutine setPCVec(vecX)

    use constants
    use communication
    use blockPointers
    use inputTimeSpectral
    use flowVarRefState
    use utils, only : EChk, setPointers
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Arguments
    Vec   vecX
    ! integer(kind=intType) ::ierr, i, j, k, l, nn, mm, sps, nVar, size, index
    ! integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2, ii, jj, procID
    ! real(kind=realType), pointer :: xPtr(:)
    ! type(commType) :: commPattern
    ! type(internalCommType) :: internal
    ! integer, dimension(mpi_status_size) :: mpiStatus

    ! call VecGetArrayReadF90(vecX, xPtr, ierr)
    ! call EChk(ierr,__FILE__,__LINE__)

    ! ! First set all the owned cells...this is basically just a straight
    ! ! copy in order
    ! ii = 0
    ! do nn=1, nDom
    !    do sps=1, nTimeIntervalsSpectral
    !       call setPointers(nn, 1, sps)
    !       do k=2, kl
    !          do j=2, jl
    !             do i=2, il
    !                do l=1, nw
    !                   ii = ii + 1
    !                   flowDoms(nn, 1, sps)%PCVec1(l, i, j, k) = xPtr(ii)
    !                end do
    !             end do
    !          end do
    !       end do
    !    end do
    ! end do

    ! ! Done with the vecX.
    ! call VecRestorearrayReadF90(vecX, xPtr, ierr)
    ! call EChk(ierr,__FILE__,__LINE__)

    ! ! Now we do a custom halo exchange. We can use the same pattern as
    ! ! whlo, but the ordering of the unknowns is different so we do our
    ! ! own.

    ! nVar = nw
    ! internal = internalCell_1st(1)
    ! commPattern = commPatternCell_1st(1)

    ! spectralModes: do mm=1,nTimeIntervalsSpectral

    !    ! Send the variables. The data is first copied into
    !    ! the send buffer after which the buffer is sent asap.

    !    ii = 1
    !    sends: do i=1,commPattern%nProcSend

    !       ! Store the processor id and the size of the message
    !       ! a bit easier.

    !       procID = commPattern%sendProc(i)
    !       size    = nVar*commPattern%nsend(i)

    !       ! Copy the data in the correct part of the send buffer.

    !       jj = ii
    !       do j=1,commPattern%nsend(i)

    !          ! Store the block id and the indices of the donor
    !          ! a bit easier.

    !          d1 = commPattern%sendList(i)%block(j)
    !          i1 = commPattern%sendList(i)%indices(j,1)
    !          j1 = commPattern%sendList(i)%indices(j,2)
    !          k1 = commPattern%sendList(i)%indices(j,3)

    !          ! Copy the given range of the working variables for
    !          ! this cell in the buffer. Update the counter jj.

    !          do k=1, nw
    !             sendBuffer(jj) = flowDoms(d1,1,mm)%PCVec1(k, i1, j1, k1)
    !             jj = jj + 1
    !          enddo

    !       enddo

    !       ! Send the data.

    !       call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
    !            procID, ADflow_comm_world, sendRequests(i), &
    !            ierr)

    !       ! Set ii to jj for the next processor.

    !       ii = jj

    !    enddo sends

    !    ! Post the nonblocking receives.

    !    ii = 1
    !    receives: do i=1,commPattern%nProcRecv

    !       ! Store the processor id and the size of the message
    !       ! a bit easier.

    !       procID = commPattern%recvProc(i)
    !       size    = nVar*commPattern%nrecv(i)

    !       ! Post the receive.

    !       call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
    !            myID, ADflow_comm_world, recvRequests(i), ierr)

    !       ! And update ii.

    !       ii = ii + size

    !    enddo receives

    !    ! Copy the local data.

    !    localCopy: do i=1,internal%ncopy

    !       ! Store the block and the indices of the donor a bit easier.

    !       d1 = internal%donorBlock(i)
    !       i1 = internal%donorIndices(i,1)
    !       j1 = internal%donorIndices(i,2)
    !       k1 = internal%donorIndices(i,3)

    !       ! Idem for the halo's.

    !       d2 = internal%haloBlock(i)
    !       i2 = internal%haloIndices(i,1)
    !       j2 = internal%haloIndices(i,2)
    !       k2 = internal%haloIndices(i,3)

    !       ! Copy the given range of working variables.

    !       do k=1, nw
    !          flowDoms(d2,1,mm)%PCVec1(k, i2, j2, k2) = &
    !               flowDoms(d1,1,mm)%PCVec1(k, i1, j1, k1)
    !       enddo

    !    enddo localCopy

    !    ! Complete the nonblocking receives in an arbitrary sequence and
    !    ! copy the variables from the buffer into the halo's.

    !    size = commPattern%nProcRecv
    !    completeRecvs: do i=1,commPattern%nProcRecv

    !       ! Complete any of the requests.

    !       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

    !       ! Copy the data just arrived in the halo's.

    !       ii = index
    !       jj = nVar*commPattern%nrecvCum(ii-1)
    !       do j=1,commPattern%nrecv(ii)

    !          ! Store the block and the indices of the halo a bit easier.

    !          d2 = commPattern%recvList(ii)%block(j)
    !          i2 = commPattern%recvList(ii)%indices(j,1)
    !          j2 = commPattern%recvList(ii)%indices(j,2)
    !          k2 = commPattern%recvList(ii)%indices(j,3)

    !          do k=1, nw
    !             jj = jj + 1
    !             flowDoms(d2,1,mm)%PCVec1(k, i2, j2, k2) = recvBuffer(jj)
    !          enddo

    !       enddo

    !    enddo completeRecvs

    !    ! Complete the nonblocking sends.

    !    size = commPattern%nProcSend
    !    do i=1,commPattern%nProcSend
    !       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    !    enddo

    ! enddo spectralModes

  end subroutine setPCVec

end module fortranPC
