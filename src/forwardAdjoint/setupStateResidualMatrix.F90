subroutine setupStateResidualMatrix(matrix,useAD,usePC,useTranspose)
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
  !     ******************************************************************
  !
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use inputDiscretization ! spaceDiscr
  use inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputAdjoint        ! useDiagTSPC
  use stencils

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix

  ! Input Variables
  logical :: useAD,usePC,useTranspose

  ! Local variables.
  integer(kind=intType) :: ierr, nn, sps, sps2, i, j, k, l, ll, ii, jj, kk
  integer(kind=intType) :: nColor, iColor, irow, icol
  integer(kind=intType) :: n_stencil, i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  real(kind=realType) :: delta_x, one_over_dx

  useDiagTSPC = .False.

  rkStage = 0
  currentLevel =1 
  groundLevel = 1
   
  ! Zero out the matrix before we start
  call MatZeroEntries(matrix,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Run the  initialize_stencils routine just in case
  call initialize_stencils

  ! Set a pointer to the correct set of stencil depending on if we are
  ! using the first order stencil or the full jacobian

  if (usePC) then
     if (viscous) then
        stencil => visc_pc_stencil
        n_stencil = N_visc_pc
     else
        stencil => euler_pc_stencil
        n_stencil = N_euler_pc
     end if

     ! Very important to use only second Order dissipation for PC 
     lumpedDiss=.True.
  else
     stencil => euler_drdw_stencil
     n_stencil = N_euler_drdw
  end if

  ! Exchange data and call the residual to make sure its up to date
  ! withe current w
  call whalo2(1_intType, 1_intType, nw, .True.,.True.,.True.)
  call computeResidualNK ! This is the easiest way to do this

  ! Set delta_x
  delta_x = 1e-6_realType
  one_over_dx = one/delta_x
  rkStage = 0

  ! Master Domain Loop
  domainLoopAD: do nn=1,nDom

     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn,1,1)

     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values(nn)

     ! Setup the coloring for this block depending on if its
     ! drdw or a PC
     
     ! List of all Coloring Routines:
     !   Debugging Colorings Below:
     !       call setup_3x3x3_coloring(nn,nColor)
     !       call setup_5x5x5_coloring(nn,nColor)
     !       call setup_BF_coloring(nn,nColor)
     !   Regular:
     !       call setup_PC_coloring(nn,nColor)
     !       call setup_dRdw_euler_coloring(nn,nColor)
     !       call setup_dRdw_visc_coloring(nn,nColor)
     
     if (usePC) then
        ! Note: The lumped dissipation doesn't quite result in a
        !3-cell stencil in each direction but we will still use PC
        !coloring. Not really a big deal.
     
        if (not (viscous)) then
           call setup_PC_coloring(nn,nColor) ! Euler Colorings
        else
           call setup_3x3x3_coloring(nn,nColor) ! dense 3x3x3 coloring
        end if
     else
        if( .not. viscous ) then
           call setup_dRdw_euler_coloring(nn,nColor) ! Euler Colorings
        else 
           call setup_dRdw_visc_coloring(nn,nColor)! Viscous/RANS
        end if
     end if
     
     spectralLoop: do sps=1,nTimeIntervalsSpectral

        ! Do Coloring and perturb states
        do iColor = 1,nColor
           do sps2 = 1,nTimeIntervalsSpectral
              flowDomsd(sps2)%dw_deriv(:,:,:,:,:) = zero
           end do

           ! Master State Loop
           do l = 1,nw

              call setPointersAdj(nn,1,sps)
              ! Reset All States and possibe AD seeds
              do sps2 = 1,nTimeIntervalsSpectral
                 flowDoms(nn,1,sps2)%w(:,:,:,:) =  flowDomsd(sps2)%wtmp
                 if (useAD) then
                    flowdomsd(sps2)%w = zero ! This is actually w seed
                 end if
              end do

              ! Peturb w or set AD Seed
              do k=0,kb
                 do j=0,jb
                    do i=0,ib
                       if (flowdomsd(1)%color(i,j,k) == icolor) then
                          if (useAD) then
                             flowdomsd(sps)%w(i,j,k,l) = one
                          else
                             w(i,j,k,l) = w(i,j,k,l) + delta_x
                          end if
                       end if
                    end do
                 end do
              end do

              ! Block-based residual
              if (useAD) then
                 call block_res_d(nn,sps)
              else
                 call block_res(nn,sps)
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
                                        one_over_dx*&
                                        (flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
                                        flowDomsd(sps2)%dwtmp(i,j,k,ll))
                                else
                                   ! If the peturbation is on an off
                                   ! instance, only subtract dwtmp2
                                   ! which is the reference result
                                   ! after initres

                                   flowDomsd(sps2)%dw_deriv(i,j,k,ll,l) = &
                                        one_over_dx*(&
                                        flowDoms(nn,1,sps2)%dw(i,j,k,ll) - &
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

           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    icol = flowDoms(nn,1,sps)%globalCell(i,j,k)
                    if (flowdomsd(1)%color(i,j,k) == icolor .and.&
                         icol >= 0) then

                       ! i,j,k are now the "Center" cell that we
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

                             ! If we're doing the PC and we want to
                             ! use the diagonal version we ONLY set
                             ! values from the ON-time instance
                             if (usePC .and. useDiagTSPC)then
                                irow = flowDoms(nn,1,sps)%&
                                     globalCell(i+ii,j+jj,k+kk)
                                call setBlock(flowDomsd(sps)%&
                                     dw_deriv(i+ii,j+jj,k+kk,:,:))
                             else
                                ! Center Stencil Cell has off
                                ! time-instance dependancies
                                if (ii == 0 .and. jj == 0 .and. kk == 0) then
                                   
                                   do sps2=1,nTimeIntervalsSpectral
                                      irow = flowDoms(nn,1,sps2)%&
                                           globalCell(i+ii,j+jj,k+kk)
                                      call setBlock(&
                                           flowDomsd(sps2)%&
                                           dw_deriv(i+ii,j+jj,k+kk,:,:))
                                   end do
                                   
                                else 
                                   irow = flowDoms(nn,1,sps)%&
                                        globalCell(i+ii,j+jj,k+kk)
                                   call setBlock(&
                                        flowDomsd(sps)%&
                                        dw_deriv(i+ii,j+jj,k+kk,:,:))
                                end if ! Center Cell Check
                             end if
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

  end do domainLoopAD

  ! Redo the complete residual to make sure all the halos/pressures
  ! are up to date
  call whalo2(1_intType, 1_intType, nw, .True., .True., .True.)
  call computeResidualNK()

  !Return dissipation Parameters to normal -> VERY VERY IMPORTANT
  if (usePC) then
     lumpedDiss = .False.
  end if

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (matrix,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatSetOption(matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)

contains

  subroutine setBlock(blk)
    ! Sets a block at irow,icol, if useTranspose is False
    ! Sets a block at icol,irow with transpose of blk if useTranspose is True

    implicit none
    real(kind=realType), dimension(nw,nw) :: blk

    if (useTranspose) then
       call MatSetValuesBlocked(matrix,1,icol,1,irow,transpose(blk),&
            ADD_VALUES,ierr)
       call EChk(ierr,__FILE__,__LINE__)
    else
       call MatSetValuesBlocked(matrix,1,irow,1,icol,blk,&
            ADD_VALUES,ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

  end subroutine setBlock


  subroutine writeOutMatrix()

    integer(kind=intType) :: nrows,ncols,icell,jcell,kcell
    real(kind=realType) :: val1(nw,nw),val2(nw,nw),err,avgval

    call MatGetOwnershipRange(matrix,nrows,ncols,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    do nn=1,nDom
       call setPointersAdj(nn,1,1)
      
     
       do kcell=2,kl
          do jcell=2,jl
             do icell=2,il
                
                irow = globalCell(icell,jcell,kcell)
               
                do i_stencil=1,n_stencil
                   ii = stencil(i_stencil,1)
                   jj = stencil(i_stencil,2)
                   kk = stencil(i_stencil,3)

                   
                   if ( icell+ii >= 2 .and. icell+ii <= il .and. &
                        jcell+jj >= 2 .and. jcell+jj <= jl .and. &
                        kcell+kk >= 2 .and. kcell+kk <= kl) then 

                      icol = globalCell(icell+ii,jcell+jj,kcell+kk)

                      do i=1,nw
                         do j=1,nw

                            call MatGetValues(matrix  ,1,icol*nw+j-1,1,irow*nw+i-1,val1(i,j),ierr)
                            call EChk(ierr,__FILE__,__LINE__)
!                            call MatGetValues(mat_copy,1,irow*nw+i-1,1,icol*nw+j-1,val2(i,j),ierr)
!                            call EChk(ierr,__FILE__,__LINE__)

                            if (useAD) then
                               write(18,30),nn,icell,jcell,kcell,icell+ii,jcell+jj,kcell+kk,i,j,val1(i,j)
                            else 
!!$                               if (val1(i,j) > 0.00001) then
                               write(16,30),nn,icell,jcell,kcell,icell+ii,jcell+jj,kcell+kk,i,j,val1(i,j)
!!$                               end if
                            end if
!!$                            if (val2(i,j) > 0.00001) then
                            write(17,30),nn,icell,jcell,kcell,icell+ii,jcell+jj,kcell+kk,i,j,val2(i,j)
!!$                            end if

                         end do
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

30  format(1x,I4,' | ', I4,' ',I4,'  ',I4,' | ',I4,' ',I4,' ',I4,' | ',I4,'  ',I4,' ',f20.6)
  end subroutine writeOutMatrix


#endif
end subroutine setupStateResidualMatrix
