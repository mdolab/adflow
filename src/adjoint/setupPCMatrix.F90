subroutine setupPCMatrix(useAD,  useTranspose, frozenTurb, level)
#ifndef USE_NO_PETSC
  
  ! This routine generates a fortran form of the PCmatrix. It is
  ! currently not used anywhere, but it become useful in the future.

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
  use utils, only : setPointers, EChk
  use haloExchange, only : whalo2
  implicit none

  ! Input Variables
  logical, intent(in) :: useAD,  useTranspose, frozenTurb
  integer(kind=intType), intent(in) :: level

  ! Local variables.
  integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
  integer(kind=intType) :: nColor, iColor, jColor, irow, icol, fmDim, frow
  integer(kind=intType) :: nTransfer, nState, tmp, icount
  integer(kind=intType) :: n_stencil, i_stencil, ind1
  integer(kind=intType), dimension(:, :), pointer :: stencil
  real(kind=realType) :: delta_x, one_over_dx
  real(kind=realType) :: delta_x_turb, one_over_dx_turb

#ifdef USE_COMPLEX
  complex(kind=realType) :: alpha, beta, alphad, betad
  complex(kind=realType), dimension(:,:), allocatable :: blk
#else
  real(kind=realType) :: alpha, beta, alphad, betad
  real(kind=realType), dimension(:,:), allocatable :: blk
#endif
  integer(kind=intType) :: liftIndex
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, mm, colInd
  logical :: resetToRANS, secondOrdSave,  splitMat
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
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)


  ! Set a pointer to the correct set of stencil depending on if we are
  ! using the first order stencil or the full jacobian

  stencil => euler_pc_stencil
  n_stencil = N_euler_pc

  ! Very important to use only Second-Order dissipation for PC 
  lumpedDiss=.True.
  secondOrdSave = secondOrd
  secondOrd = .False.

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
        shockSensor => flowDoms(nn,1,sps)%shockSensor
        
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
                 call block_res_d(nn, sps, .False., &
                      alpha, alphad, beta, betad, liftIndex, frozenTurb)
#else
                 print *, 'Forward AD routines are not complexified'
                 stop
#endif
              else
                 call block_res(nn, sps, .False., alpha, beta, &
                      liftIndex, frozenTurb)
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

        ! Deallocate the shock sensor refernce if usePC 
        deallocate(flowDoms(nn, 1, sps)%shockSensor)
     end do
  end do

  ! Return dissipation Parameters to normal -> VERY VERY IMPORTANT
  lumpedDiss = .False.
  secondOrd = secondOrdSave
  
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
