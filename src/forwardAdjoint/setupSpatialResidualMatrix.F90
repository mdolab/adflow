subroutine setupSpatialResidualMatrix(matrix, useAD, useObjective, frozenTurb)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the spatial derivative matrix using a forward mode calc*
  !     * There is one flag to determine how this routine is run:        *
  !     *                                                                *
  !     * useAD: if True, AD is used for derivative calculation, if      *
  !     *        False, FD is used.                                      *
  !     *                                                                *
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
  use communication
  use adjointVars
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical, intent(in) :: useAD, useObjective, frozenTurb

  ! Local variables.
  integer(kind=intType) :: ierr,nn,sps,sps2,i,j,k,l,ll,ii,jj,kk, mm
  integer(kind=intType) :: irow, icol, level, fdim
  integer(kind=intType) :: n_stencil, i_stencil, n_force_stencil, nState
  integer(kind=intType), dimension(:,:), pointer :: stencil, force_stencil
  integer(kind=intType) :: nColor, iColor, jColor, ind, fmInd
  real(kind=realType) :: delta_x,one_over_dx, val
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, fmDim
  real(kind=realType) :: alpha, beta, sepSensor, Cavitation
  real(kind=realType) :: alphad, betad, sepSensord, Cavitationd
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment, forced, momentd
  integer(kind=intType) :: liftIndex
  logical :: resetToRANS

  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here
  alphad = zero
  betad  = zero
  machd  = zero
  machCoefd = zero
  machGridd = zero
  lengthRefd = zero
  pointRefd  = zero
  surfaceRefd = zero
  reynoldslengthd = zero
  reynoldsd = zero
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

! Setup number of state variable based on turbulence assumption
  if ( frozenTurb ) then
     nState = nwf
  else
     nState = nw
  endif

  ! Hardcode levels since assembling on coarser levels is not working
  level = 1
  rkStage = 0
  currentLevel =1 
  groundLevel = 1

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Run the initialize_stencils routine just in case
  call initialize_stencils

  ! Set a pointer to the correct set of stencil depending on if we are
  ! using the first order stencil or the full jacobian

  if(viscous) then
     stencil => visc_drdx_stencil
     n_stencil = N_visc_drdx
     force_stencil => visc_force_x_stencil
     n_force_stencil = N_visc_force_x
  else 
     stencil => euler_drdx_stencil
     n_stencil = N_euler_drdx
     force_stencil => euler_force_x_stencil
     n_force_stencil = N_euler_force_x
  end if

  ! Call the residual to make sure its up to date withe current w
  call whalo2(1_intType, 1_intType, nw, .True.,.True.,.True.)
  call computeResidualNK

  ! Set delta_x
  delta_x = 1e-5
  one_over_dx = 1.0/delta_x
     
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

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurb .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Allocate the memory we need for this block to do the forward
  ! mode derivatives and copy reference values
  if (.not. derivVarsAllocated) then 
     call alloc_derivative_values(level)
  end if
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        call zeroADSeeds(nn,level, sps)
     end do
  end do

  ! Master Domain Loop
  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn, level, 1)

     
     
     ! Setup the coloring for this block depending on if its
     ! drdw or a PC

     ! Debugging Colorings Below:
     !call setup_3x3x3_coloring(nn, level, nColor)
     !call setup_5x5x5_coloring(nn, level, nColor)
     !call setup_BF_coloring(nn, level, nColor)

     if (viscous) then
        call setup_dRdx_visc_coloring(nn, level, nColor)! Viscous/RANS
        !call setup_5x5x5_coloring(nn, level, nColor)! Viscous/RANS
        !call setup_BF_coloring(nn, level, nColor)! Viscous/RANS
     else 
        call setup_dRdx_euler_coloring(nn, level, nColor)
        !call setup_4x4x4_coloring(nn, level, nColor)
        !call setup_5x5x5_coloring(nn, level, nColor)
        !call setup_BF_Node_coloring(nn, level, nColor)
     end if

     spectralLoop: do sps=1,nTimeIntervalsSpectral
        ! Set pointers and derivative pointers
        call setPointers_d(nn, 1, sps)

        ! Set unknown sizes in diffSizes for AD routine
        ISIZE1OFDrfbcdata = nBocos
        ISIZE1OFDrfviscsubface = nViscBocos
   
        ! Do Coloring and perturb states
        colorLoop: do iColor = 1,nColor
           do sps2 = 1,nTimeIntervalsSpectral
              flowDomsd(nn,1,sps2)%dw_deriv = zero
           end do

           ! Master Node Loop
           dofLoop: do l = 1,3

              ! Reset All Coordinates and possibe AD seeds
              do sps2 = 1,nTimeIntervalsSpectral
                 do ll=1,3
                    do k=0,ke
                       do j=0,je
                          do i=0,ie
                             flowDoms(nn,1,sps2)%x(i,j,k,ll) = flowDomsd(nn,1,sps2)%xtmp(i,j,k,ll)
                          end do
                       end do
                    end do
                 end do
                 do ll=1,nw
                    do k=0,kb
                       do j=0,jb
                          do i=0,ib
                             flowDoms(nn,1,sps2)%w(i,j,k,ll) = flowDomsd(nn,1,sps2)%wtmp(i,j,k,ll)
                          end do
                       end do
                    end do
                 end do

                 if (useAD) then
                    flowdomsd(nn,1,sps2)%x = zero ! This is actually
                    flowdomsd(nn,1,sps2)%w = zero ! This is actually
                    ! the x seed
                 end if
              end do

              ! Peturb x or set AD Seed
              do k=0,ke
                 do j=0,je
                    do i=0,ie
                       if (flowdomsd(nn,1,1)%color(i,j,k) == icolor .and.&
                            globalnode(i,j,k) >= 0) then
                          if (useAD) then
                             flowdomsd(nn,1,sps)%x(i,j,k,l) = one
                          else
                             x(i,j,k,l) = x(i,j,k,l) + delta_x
                          end if
                       end if
                    end do
                 end do
              end do

              ! Block-based residual
              if (useAD) then
#ifndef USE_COMPLEX
                 call block_res_d(nn, sps, .True., &
                      alpha, alphad, beta, betad, liftIndex, force, forced, &
                      moment, momentd, sepSensor, sepSensord, Cavitation, Cavitationd, frozenTurb)
#else
                 print *,'Forward AD routines are not complexified!'
                 stop
#endif
              else
                 call block_res(nn, sps, .True., &
                      alpha, beta, liftIndex, force, moment, sepSensor, Cavitation, frozenTurb)
              end if

              ! Set the computed residual in dw_deriv. If using FD,
              ! actually do the FD calculation if AD, just copy out dw
              ! in flowdomsd

              ! Take all Derivatives
              do sps2 = 1,nTimeIntervalsSpectral
                 do ll=1,nState
                    do k=2,kl 
                       do j=2,jl
                          do i=2,il
                             if (useAD) then
                                flowDomsd(nn,1,sps2)%dw_deriv(i,j,k,ll,l) = &
                                     flowdomsd(nn,1,sps2)%dw(i,j,k,ll)
                             else
                                if (sps2 == sps) then
                                   ! If the peturbation is on this
                                   ! instance, we've computed the spatial
                                   ! contribution so subtract dwtmp
                                   flowDomsd(nn,1,sps2)%dw_deriv(i,j,k,ll,l) = &
                                        one_over_dx*(flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
                                        flowDomsd(nn,1,sps2)%dwtmp(i,j,k,ll))
                                else
                                   ! If the peturbation is on an off
                                   ! instance, only subtract dwtmp2
                                   ! which is the reference result
                                   ! after initres
                                   flowDomsd(nn,1,sps2)%dw_deriv(i,j,k,ll,l) = &
                                        one_over_dx*(flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
                                        flowDomsd(nn,1,sps2)%dwtmp2(i,j,k,ll))
                                end if
                             end if
                          end do
                       end do
                    end do
                 end do
              end do
           end do dofLoop

           ! Set derivatives by block in "matrix" after we've peturbed
           ! all states in "color"

           kLoop: do k=0,ke
              jLoop: do j=0,je
                 iLoop: do i=0,ie
                    icol = flowDoms(nn,1,sps)%globalNode(i,j,k)

                    colorCheck: if (flowdomsd(nn,1,1)%color(i,j,k) == icolor &
                         .and. icol >= 0) then

                       ! i,j,k are now the "Center" node that we
                       ! actually petrubed. From knowledge of the
                       ! stencil, we can simply take this node and
                       ! using the stencil, set the values around it
                       ! in PETSc

                       stencilLoop: do i_stencil=1,n_stencil
                          ii = stencil(i_stencil,1)
                          jj = stencil(i_stencil,2)
                          kk = stencil(i_stencil,3)

                          ! Check to see if the cell in this
                          ! sentcil is on a physical cell, not a
                          ! halo
                          onBlock: if ( i+ii >= 2 .and. i+ii <= il .and. &
                               j+jj >= 2 .and. j+jj <= jl .and. &
                               k+kk >= 2 .and. k+kk <= kl) then 

                             ! Eight of the cells around the node will
                             ! have off-time instance contributions so
                             ! we will have to loop over the time
                             ! instances for those:

                             if ( (ii == 0 .or. ii == 1) .and. &
                                  (jj == 0 .or. jj == 1) .and. &
                                  (kk == 0 .or. kk == 1)) then 

                                do sps2=1,nTimeIntervalsSpectral
                                   irow = flowDoms(nn,1,sps2)%&
                                        globalCell(i+ii,j+jj,k+kk)
                                   call setBlock(&
                                        flowDomsd(nn,1,sps2)%&
                                        dw_deriv(i+ii,j+jj,k+kk,1:nstate,1:nstate))
                                end do
                             else
                                irow = flowDoms(nn,1,sps)%globalCell(&
                                     i+ii,j+jj,k+kk)

                                call setBlock(flowDomsd(nn,1,sps)%&
                                     dw_deriv(i+ii,j+jj,k+kk,1:nstate,1:nstate))
                             end if
                          end if onBlock
                       end do stencilLoop
                    end if colorCheck
                 end do iLoop
              end do jLoop
           end do kLoop
        end do colorLoop
     end do spectralLoop
  end do domainLoopAD

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd  (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

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
     restrictEddyVis = .False.
     if( eddyModel ) restrictEddyVis = .True.
  end if

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol. Note that blk is actually
    ! (nState,nState) but since this is drdx, we're only using
    ! (nState,3) chunk of it. This should be ok, nState will always be
    ! greater than 3

    implicit none
#ifdef USE_COMPLEX
    complex(kind=realType), dimension(nState, nState) :: blk
#else
    real(kind=realType), dimension(nState, nState) :: blk
#endif
    integer(kind=intType) :: iii,jjj, nrows, ncols

    do jjj=1,3
       do iii=1,nState
          ! NOTE: We are setting the values in the tranpose
          ! sense. That's why icol is followed by irow. 
          if ( .not. blk(iii,jjj) == zero) then
             call MatSetValues(matrix, 1, icol*3+jjj-1, 1, irow*nState+iii-1, &
                  blk(iii, jjj), ADD_VALUES, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if
       end do
    end do

  end subroutine setBlock

#endif
end subroutine setupSpatialResidualMatrix
