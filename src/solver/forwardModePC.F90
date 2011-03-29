subroutine setupNK_KSP_PC3(dRdwPre)

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the dRdWPre matrix for the NK Solver using Finite      *
  !     * Difference and a color based stencil                           *
  !     ******************************************************************
  !
  use ADjointVars
  use ADjointPETSc , only:localInfo
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
  Mat dRdwPre

  !     Local variables.
  integer(kind=intType) :: ierr,nn,sps,sps2,i,j,k,l,ll,ii,jj,kk
  integer(kind=intType) :: iCell,jCell,kCell
  integer(kind=intType) :: idxngb,idxmgb
  integer(kind=intType) :: istencil,jstencil,kstencil
  real(kind=realType) :: delta_x,one_over_dx
  real(kind=realType) :: gm1,v2
  logical :: secondHalo
  real(kind=realType), dimension(2) :: time
  real(kind=realType)               ::setupTime
  integer(kind=intType) :: n_stencil,i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  

  ! Set the grid level of the current MG cycle, the value of the
  ! discretization and the logical correctForK.
  currentLevel = 1
  rkStage = 0

  time(1) = mpi_wtime()

  ! Very important to use only second Order dissipation!
  lumpedDiss=.True.

  ! Zero out the matrix before we start
  call MatZeroEntries(dRdwPre,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Run the  initialize_stencils routine just in case
  call initialize_stencils

  ! Set a pointer to the correct set of stencil depending on what time
  ! of equations we're solving
  if (viscous) then
     stencil => euler_pc_stencil
     n_stencil = N_euler_pc
  else
     stencil => euler_pc_stencil
     n_stencil = N_euler_pc
  end if
  print *,'n_stencil:',n_stencil
  print *,'stencil:',stencil
  ! Call the residual to make sure its up to date withe current w
  call computeResidualNK ! This is the easiest way to do this

  ! Allocate Memory and copy out w and dw for reference
  allocatedomains: do nn = 1,ndom
     allocspectralLoop: do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1,sps)
        allocate(flowDoms(nn,1,sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        allocate(flowDoms(nn,1,sps)%dw_FD(0:ib,0:jb,0:kb,1:nw,1:nw),stat=ierr)
        allocate(flowDoms(nn,1,sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        flowDoms(nn,1,sps)%dw_FD(:,:,:,:,:) = 0.0
        do l=1,nw
           do k=0,kb 
              do j=0,jb
                 do i=0,ib
                    flowdoms(nn,1,sps)%wtmp(i,j,k,l) = w(i,j,k,l)
                    flowdoms(nn,1,sps)%dwtmp(i,j,k,l) = &
                         flowdoms(nn,1,sps)%dw(i,j,k,l)/vol(i,j,k)
                 end do
              end do
           end do
        end do
     end do allocspectralLoop
  end do allocatedomains

  ! Set delta_x
  delta_x = 1e-8
  one_over_dx = 1.0/delta_x
  rkStage = 0
  currentLevel =1 
  secondHalo = .True. 

  ! Master Domain Loop
  domainLoopAD1: do nn=1,nDom

     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives
     call alloc_derivative_values(nn)

     spectralLoop1: do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1,sps)

        ! Master Stencil Loops
        do kstencil = 0,2
           do jstencil = 0,2
              do istencil = 0,2

                 ! Master State Loop
                 do l = 1,nw

                    ! Reset All States and Seeds
                    w(:,:,:,:) =  flowDoms(nn,1,sps)%wtmp
                    flowdomsd(sps)%w = 0.0 
                    
                    ! Do Coloring and perturb states
                    do kCell=kstencil,kb,3
                       do jCell=jstencil,jb,3
                          do iCell=istencil,ib,3
                             flowdoms(nn,1,sps)%w(iCell,jCell,kCell,l) = &
                                  flowDoms(nn,1,sps)%wtmp(iCell,jCell,kCell,l) &
                                  + delta_x
!                             flowdomsd(sps)%w(icell,jcell,kcell,l) = 1.0

                          end do
                       end do
                    end do
                    
                    ! Block-based residual
                    innerSpectralLoop: do sps2=1,nTimeIntervalsSpectral
                       call setPointersAdj(nn,1,sps2)
                       !call setPointersd(sps2)
                       !call block_res_d(nn,sps2)
                       call block_res(nn,sps2)
                    end do innerSpectralLoop
                    
                    ! This is actually the finite difference calc
                    do ll=1,nw
                       do k=2,kl 
                          do j=2,jl
                             do i=2,il
                                flowDoms(nn,1,sps)%dw_FD(i,j,k,ll,l) = &
                                     one_over_dx*&
                                     (dw(i,j,k,ll)/vol(i,j,k) - &
                                     flowDoms(nn,1,sps)%dwtmp(i,j,k,ll))
!                                  flowDoms(nn,1,sps)%dw_FD(i,j,k,ll,l) = &
!                                       flowdomsd(sps)%dw(i,j,k,ll)/vol(i,j,k)

                             end do
                          end do
                       end do
                    end do
                 end do ! State Loop

                 ! Set derivatives by block in dRdwPre after we've peturbed
                 ! all states in one block
                    
                 do kCell=kstencil,kb,3
                    do jCell=jstencil,jb,3
                       do iCell=istencil,ib,3

                          ! iCell,jCell,kCell are now the "Center" cell that we
                          ! actually petrubed. From knowledge of the
                          ! stencil, we can simply take this cell and
                          ! the extra 6 around it to set in PETSc

                          idxngb = flowDoms(nn,1,sps)% &
                               globalCell(iCell,jCell,kCell)

                          do i_stencil=1,n_stencil
                             ii = stencil(i_stencil,1)
                             jj = stencil(i_stencil,2)
                             kk = stencil(i_stencil,3)

                             idxmgb = flowDoms(nn,1,sps)% &
                                  globalCell(iCell+ii,jcell+jj,kcell+kk)

                             if (iCell+ii >= 2 .and. iCell+ii <= il .and. &
                                  jCell+jj >= 2 .and. jCell+jj <= jl .and. &
                                  kCell+kk >= 2 .and. kCell+kk <= kl) then 
                                
!                                 call checkBlock( flowDoms(nn,1,sps)%&
!                                      dw_FD(iCell+ii,jCell+jj,Kcell+kk,:,:),nn,&
!                                      icell+ii,jcell+jj,kcell+kk)
                               

                                call MatSetValuesBlocked(&
                                     dRdWPre,1,idxmgb,1,idxngb,&
                                     flowDoms(nn,1,sps)%&
                                     dw_FD(iCell+ii,jCell+jj,Kcell+kk,:,:),&
                                     ADD_VALUES,ierr)
                                call EChk(ierr,__FILE__,__LINE__)
                             end if
                          end do
                                        
                       end do ! istencil loop
                    end do ! jstencil loop
                 end do ! kstencil loop

              end do ! outer istencil loop
           end do ! outer jstencil loop
        end do ! outer k stencil loop

     end do spectralLoop1

     call dealloc_derivative_values(nn)
     
  end do domainLoopAD1

  ! Reset w and dw -> Its like nothing happened...
  deallocatedomains: do nn = 1,ndom
     deallocatespectral: do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1,sps)
        ! Reset w 
        w = flowDoms(nn,1,sps)%wtmp

        ! Set dw
        dw =flowDoms(nn,1,sps)%dwtmp

        ! Deallocate memtory
        deallocate(flowDoms(nn,1,sps)%dwtmp,stat=ierr)
        deallocate(flowDoms(nn,1,sps)%dw_FD,stat=ierr)
        deallocate(flowDoms(nn,1,sps)%wtmp,stat=ierr)

     end do deallocatespectral
  end do deallocatedomains

  !Return dissipation Parameters to normal -> VERY VERY IMPORTANT
  lumpedDiss = .False.

  ! PETSc Matrix Assembly and Options Set
  call MatAssemblyBegin(dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

#ifdef USE_PETSC_3
  call MatSetOption(dRdWPre,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call MatSetOption(dRdWPre,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  if (myid == 0) then
     print *,'Assembly time:',setupTime
  end if

!   call MatGetInfo(dRdwpre,MAT_LOCAL,localInfo,ierr)
!   call EChk(ierr,__FILE__,__LINE__)
!   call printLocal("dRdwT")

end subroutine setupNK_KSP_PC3


subroutine checkBlock(blk,nn,ii,jj,kk)
  use flowvarrefstate
  use communication
  implicit none

  real(kind=realType) :: blk(nw,nw)
  integer(kind=intType) :: i,ii,jj,kk,nn
  ! Just check to see if any of the diags on blk are zero
  !if (myid == 0) then
!   do i=1,nw
!      if (abs(blk(i,i)) < 1e-14) then
!         !print *,'Zero diag found on proc:',myid
!         !print *,nn,ii,jj,kk,i
!         !print *,blk(:,i)
!         !stop
!         blk(i,i) = 1.0
!      end if
!   end do
!   !end if

  
  if (myid == 0) then
     do i=1,nw
        if (abs(blk(i,i)) < 1e-14) then
           print *,'Zero diag found on proc 0'
           print *,'block i,j,k,row:'
           print *,nn,ii,jj,kk,i
           print *,blk(:,i)
           stop
           blk(i,i) = 1.0
        end if
     end do
  end if


end subroutine checkBlock

