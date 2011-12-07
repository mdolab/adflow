subroutine setupSpatialResidualMatrix(matrix,useAD)
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
  use ADjointVars
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use inputDiscretization ! spaceDiscr
  USE inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputTimeSpectral   ! spaceDiscr
  use inputDiscretization
  use inputPhysics 
  use stencils

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix
  Mat mat_copy
  ! Input Variables
  logical :: useAD

  !     Local variables.
  integer(kind=intType) :: ierr,nn,sps,sps2,i,j,k,l,ll,ii,jj,kk
  integer(kind=intType) :: irow,icol,ilow,ihigh
  real(kind=realType) :: delta_x,one_over_dx

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               ::setupTime,trace,nrm
  integer(kind=intType) :: n_stencil,i_stencil, assembled
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: nColor,iColor
  logical :: secondHalo
  logical , dimension(:,:,:,:),allocatable :: x_peturb

  rkStage = 0
  currentLevel =1 
  groundLevel = 1
  ! Start Timer
  time(1) = mpi_wtime()

  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Run the  initialize_stencils routine just in case
  call initialize_stencils

  ! Set a pointer to the correct set of stencil depending on if we are
  ! using the first order stencil or the full jacobian

  if( .not. viscous ) then
     stencil => euler_drdx_stencil
     n_stencil = N_euler_drdx
  else 
     !stencil => visc_drdx_stencil
     !n_stencil = N_visc_drdx
  end if

  ! Call the residual to make sure its up to date withe current w
  call whalo2(1_intType, 1_intType, nw, .True.,.True.,.True.)
  call computeResidualNK

  ! Set delta_x
  delta_x = 1e-7
  one_over_dx = 1.0/delta_x
  rkStage = 0
  secondHalo = .True. 
  ! Master Domain Loop
  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn,1,1)

     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn)
     allocate(x_peturb(0:ie,0:je,0:ke,3))

     ! Setup the coloring for this block depending on if its
     ! drdw or a PC

     ! Debugging Colorings Below:
     !call setup_3x3x3_coloring(nn,nColor)
     !call setup_5x5x5_coloring(nn,nColor)
     !call setup_BF_coloring(nn,nColor)

     if( .not. viscous ) then
        call setup_dRdx_euler_coloring(nn,nColor) ! Euler Colorings
        !call setup_4x4x4_coloring(nn,nColor)
        !call setup_5x5x5_coloring(nn,nColor)
        !call setup_BF_Node_coloring(nn,nColor)
     else 
        !call setup_dRdx_visc_coloring(nn,nColor)! Viscous/RANS
        print *,'not done yet'
        stop
     end if

     spectralLoop: do sps=1,nTimeIntervalsSpectral

        ! Do Coloring and perturb states
        do iColor = 1,nColor
           do sps2 = 1,nTimeIntervalsSpectral
              flowDomsd(sps2)%dw_deriv(:,:,:,:,:) = 0.0
           end do

           ! Master Node Loop
           do l = 1,3
              call setPointersAdj(nn,1,sps)

              ! Reset All Coordinates and possibe AD seeds
              do sps2 = 1,nTimeIntervalsSpectral
                 flowDoms(nn,1,sps2)%x(:,:,:,:) =  flowDomsd(sps2)%xtmp

                 if (useAD) then
                    flowdomsd(sps2)%x = 0.0 ! This is actually the x seed
                 end if
              end do
              x_peturb = .False.
              ! Peturb x or set AD Seed
              do k=0,ke
                 do j=0,je
                    do i=0,ie
                       if (flowdomsd(1)%color(i,j,k) == icolor .and. globalnode(i,j,k) >= 0) then
                          if (useAD) then
                             flowdomsd(sps)%x(i,j,k,l) = 1.0
                          else
                             ! Save the peturbation
                             !flowdomsd(sps)%x(i,j,k,l) = one_over_dx * 1.0
                             x_peturb(i,j,k,l) = .True.
                          end if
                       end if
                    end do
                 end do
              end do

              ! Block-based residual
              if (useAD) then
                 call block_res_spatial_spatial_d(nn,sps)
              else
                 call block_res_spatial(nn,sps)
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
                                flowDomsd(sps2)%dw_deriv(i,j,k,ll,l) = &
                                     flowdomsd(sps2)%dw(i,j,k,ll)

                             else
                                if (sps2 == sps) then
                                   ! If the peturbation is on this
                                   ! instance, we've computed the spatial
                                   ! contribution so subtrace dwtmp

                                   flowDomsd(sps2)%dw_deriv(i,j,k,ll,l) = &
                                        one_over_dx*(flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
                                        flowDomsd(sps2)%dwtmp(i,j,k,ll))

                                else

                                   ! If the peturbation is on an off
                                   ! instance, only subtract dwtmp2
                                   ! which is the reference result
                                   ! after initres

                                   flowDomsd(sps2)%dw_deriv(i,j,k,ll,l) = &
                                        one_over_dx*(flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
                                        flowDomsd(sps2)%dwtmp2(i,j,k,ll))
                                end if

                             end if
                          end do
                       end do
                    end do
                 end do
              end do
           end do ! State Loop

           ! Set derivatives by block in "matrix" after we've peturbed
           ! all states in "color"

           call setPointersAdj(nn,1,sps)

           do k=0,ke
              do j=0,je
                 do i=0,ie
                    icol = flowDoms(nn,1,sps)%globalNode(i,j,k)

                    if (flowdomsd(1)%color(i,j,k) == icolor .and. icol >= 0) then

                       ! i,j,k are now the "Center" node that we
                       ! actually petrubed. From knowledge of the
                       ! stencil, we can simply take this cell and
                       ! using the stencil, set the values around it
                       ! in PETSc

                       do i_stencil=1,n_stencil
                          ii = stencil(i_stencil,1)
                          jj = stencil(i_stencil,2)
                          kk = stencil(i_stencil,3)

                          ! Check to see if the cell in this
                          ! sentcil is on a physical cell, not a
                          ! halo
                          if ( i+ii >= 2 .and. i+ii <= il .and. &
                               j+jj >= 2 .and. j+jj <= jl .and. &
                               k+kk >= 2 .and. k+kk <= kl) then 

                             ! Ignore TS Stuff... Fix later

                             irow = flowDoms(nn,1,sps)%globalCell(i+ii,j+jj,k+kk)
                             call setBlock(flowDomsd(sps)%dw_deriv(i+ii,j+jj,k+kk,:,:))

                          end if ! On block Check
                       end do ! Stencil Loop
                    end if ! Color If check
                 end do ! i loop
              end do ! j loop
           end do ! k loop
        end do ! Color Loop
     end do spectralLoop

     ! Deallocate and reset Values
     call dealloc_derivative_values(nn)
     deallocate(x_peturb)
  end do domainLoopAD

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatSetOption(matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol. Note that blk is actually (nw,nw) but
    ! since this is drdx, we're only using (nw,3) chunk of it. This
    ! should be ok, nw will always be greater than 3

    implicit none
    real(kind=realType), dimension(nw,nw) :: blk
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

  subroutine writeOutMatrix()

    integer(kind=intType) :: nrows,ncols,icell,jcell,kcell,inode,jnode,knode
    real(kind=realType) :: val1(nw,3),val2(nw,3),err,avgval

    call MatGetOwnershipRange(matrix,nrows,ncols,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    do nn=1,1!nDom
       call setPointersAdj(nn,1,1)


       do knode=1,kl
          do jnode=1,jl
             do inode=1,il

                icol = globalNode(inode,jnode,knode)

                do i_stencil=1,n_stencil
                   ii = stencil(i_stencil,1)
                   jj = stencil(i_stencil,2)
                   kk = stencil(i_stencil,3)


                   if ( inode+ii >= 2 .and. icell+ii <= il .and. &
                        jnode+jj >= 2 .and. jnode+jj <= jl .and. &
                        knode+kk >= 2 .and. knode+kk <= kl) then 

                      irow = globalCell(inode+ii,jnode+jj,knode+kk)

                      do i=1,nw
                         do j=1,3

                            call MatGetValues(matrix  ,1,irow*nw+i-1,1,icol*3+j-1,val1(i,j),ierr)
                            call MatGetValues(mat_copy,1,irow*nw+i-1,1,icol*3+j-1,val2(i,j),ierr)

                            write(13,30),nn,inode,jnode,knode,inode+ii,jnode+jj,knode+kk,i,j,val1(i,j)
                            write(14,30),nn,inode,jnode,knode,inode+ii,jnode+jj,knode+kk,i,j,val2(i,j)

                         end do
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

30  format(1x,I4,' | ', I4,' ',I4,'  ',I4,' | ',I4,' ',I4,' ',I4,' | ',I4,'  ',I4,' ',f20.4)
  end subroutine writeOutMatrix
#endif
end subroutine setupSpatialResidualMatrix
