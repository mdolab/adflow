subroutine setupSpatialResidualMatrix(matrix, useAD, useObjective)
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
  use ADjointPetsc, only : FMx
  use blockPointers      
  use BCTypes
  use inputDiscretization 
  USE inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState    
  use stencils
  use diffSizes

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical, intent(in) :: useAD, useObjective

  ! Local variables.
  integer(kind=intType) :: ierr,nn,sps,sps2,i,j,k,l,ll,ii,jj,kk, mm
  integer(kind=intType) :: irow,icol, level
  integer(kind=intType) :: n_stencil,i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: nColor, iColor, jColor
  real(kind=realType) :: delta_x,one_over_dx, val
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, fmDim
  real(kind=realType) :: alpha, beta, Lift, Drag, CL, CD
  real(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
  real(kind=realType) :: alphad, betad, Liftd, Dragd, CLd, CDd
  real(kind=realType), dimension(3) :: Forced, Momentd, cForced, cMomentd
  integer(kind=intType) :: liftIndex
  integer(kind=intType), dimension(:,:), pointer ::  colorPtr

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
     !stencil => visc_drdx_stencil
     !n_stencil = N_visc_drdx
  else 
     stencil => euler_drdx_stencil
     n_stencil = N_euler_drdx
  end if

  ! Call the residual to make sure its up to date withe current w
  call whalo2(1_intType, 1_intType, nw, .True.,.True.,.True.)
  call computeResidualNK

  ! Set delta_x
  delta_x = 1e-5
  one_over_dx = 1.0/delta_x
  if (useObjective .and. useAD) then
     do fmDim=1,6
        call VecZeroEntries(FMx(fmDim), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  ! Master Domain Loop
  domainLoopAD: do nn=1,nDom

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

     ! Debugging Colorings Below:
     !call setup_3x3x3_coloring(nn, level, nColor)
     !call setup_5x5x5_coloring(nn, level, nColor)
     !call setup_BF_coloring(nn, level, nColor)

     if(viscous ) then
        !call setup_dRdx_visc_coloring(nn, level, nColor)! Viscous/RANS
        print *,'not done yet'
        stop
     else 
        call setup_dRdx_euler_coloring(nn, level, nColor)
        !call setup_4x4x4_coloring(nn, level, nColor)
        !call setup_5x5x5_coloring(nn, level, nColor)
        !call setup_BF_Node_coloring(nn, level, nColor)
     end if

     spectralLoop: do sps=1,nTimeIntervalsSpectral
        ! Set pointers and derivative pointers
        call setPointers_d(nn, 1, sps)

        ! Do Coloring and perturb states
        colorLoop: do iColor = 1,nColor
           do sps2 = 1,nTimeIntervalsSpectral
              flowDomsd(nn,1,sps2)%dw_deriv = zero
           end do

           ! Master Node Loop
           dofLoop: do l = 1,3

              ! Reset All Coordinates and possibe AD seeds
              do sps2 = 1,nTimeIntervalsSpectral
                 flowDoms(nn,1,sps2)%x = flowDomsd(nn,1,sps2)%xtmp

                 if (useAD) then
                    flowdomsd(nn,1,sps2)%x = zero ! This is actually
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
                 call block_res_d(nn, sps, .True., .False., &
                      alpha, alphad, beta, betad, liftIndex, Force, Forced, &
                      Moment, Momentd, lift, liftd, drag, dragd, cForce, &
                      cForced, cMoment, cMomentd, CL, CLd, CD, CDd)
#else
                 print *,'Forward AD routines are not complexified!'
                 stop
#endif
              else
                 call block_res(nn, sps, .True., .False., &
                      alpha, beta, liftIndex, Force, Moment, Lift, Drag, &
                      cForce, cMoment, CL, CD)
              end if

              ! If required, set values in the 6 vectors defined in
              ! FMx. 
              if (useObjective .and. useAD) then
                 ! We need to loop over the faces on this block and
                 ! set values in FMx
                 
                 bocos: do mm=1,nBocos

                    if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic &
                         .or. BCType(mm) == NSWallIsothermal) then

                       ! Set the globalNodePtr depending on what face
                       ! we are on:

                       select case (BCFaceID(mm))
                       case (iMin)
                          colorPtr => flowDomsd(nn, 1, 1)%color(1, :, :)
                       case (iMax)
                          colorPtr => flowDomsd(nn, 1, 1)%color(il, :, :)
                       case (jMin)
                          colorPtr => flowDomsd(nn, 1, 1)%color(:, 1, :)
                       case (jMax)
                          colorPtr => flowDomsd(nn, 1, 1)%color(:, jl, :)
                       case (kMin)
                          colorPtr => flowDomsd(nn, 1, 1)%color(:, :, 1)
                       case (kMax)
                          colorPtr => flowDomsd(nn, 1, 1)%color(:, :, kl)
                       end select

                       ! These are the indices for the NODES!
                       jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
                       iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd
                 
                       do j=jBeg, jEnd ! This is a node loop
                          do i=iBeg, iEnd ! This is a node loop

                             ! The +1 in the indices are due to the
                             ! offset from the globalNodePtr
                             jcolor = colorPtr(i+1, j+1)
                             if (jColor == iColor) then
                                do fmDim = 1,3
                                   call VecSetValues(FMx(fmDim), 1, &
                                        BCData(mm)%FMIndex(i,j)*3 + l -1, &
                                        Forced(fmDim), ADD_VALUES, ierr) 
                                   call EChk(ierr, __FILE__, __LINE__)

                                   call VecSetValues(FMx(fmDim+3), 1, &
                                        BCData(mm)%FMIndex(i,j)*3 + l -1, &
                                        Momentd(fmDim), ADD_VALUES, ierr) 
                                   call EChk(ierr, __FILE__, __LINE__)
                                end do
                             end if
                          end do
                       end do
                    end if
                 end do bocos
              end if
              ! Set the computed residual in dw_deriv. If using FD,
              ! actually do the FD calculation if AD, just copy out dw
              ! in flowdomsd

              ! Take all Derivatives
              do sps2 = 1,nTimeIntervalsSpectral
                 do ll=1,nw
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
                                        dw_deriv(i+ii,j+jj,k+kk,:,:))
                                end do
                             else
                                irow = flowDoms(nn,1,sps)%globalCell(&
                                     i+ii,j+jj,k+kk)

                                call setBlock(flowDomsd(nn,1,sps)%&
                                     dw_deriv(i+ii,j+jj,k+kk,:,:))
                             end if
                          end if onBlock
                       end do stencilLoop
                    end if colorCheck
                 end do iLoop
              end do jLoop
           end do kLoop
        end do colorLoop
     end do spectralLoop

     ! Deallocate and reset Values
     call dealloc_derivative_values(nn, level)
  end do domainLoopAD

  if (useObjective .and. useAD) then
     do fmDim=1,6
        call VecAssemblyBegin(FMx(fmDim), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAssemblyEnd(FMx(fmDim), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd  (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol. Note that blk is actually (nw,nw) but
    ! since this is drdx, we're only using (nw,3) chunk of it. This
    ! should be ok, nw will always be greater than 3

    implicit none
#ifdef USE_COMPLEX
    complex(kind=realType), dimension(nw,nw) :: blk
#else
    real(kind=realType), dimension(nw,nw) :: blk
#endif
    integer(kind=intType) :: iii,jjj, nrows, ncols

    do jjj=1,3
       do iii=1,nw
         
          ! NOTE: We are setting the values in the tranpose
          ! sense. That's why icol is followed by irow. 

          call MatSetValues(matrix,1,icol*3+jjj-1,1,irow*nw+iii-1,&
               blk(iii,jjj),ADD_VALUES,ierr)
          call EChk(ierr,__FILE__,__LINE__)

       end do
    end do

  end subroutine setBlock

#endif
end subroutine setupSpatialResidualMatrix
