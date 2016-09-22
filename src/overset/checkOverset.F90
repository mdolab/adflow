subroutine checkOverset (level, sps, totalOrphans)

  !
  !       CheckOverset checks the integrity of the overset connectivity  
  !       and holes. For every comptue cell (iblank = 1) it checks that  
  !       every cell in its stencil are not blanked. If even 1 cell is   
  !      * found with an incomplete stencil it is a fatal error. 
  use constants
  use blockPointers, only : il, jl, kl, iblank, flowDoms, nDom, orphans, &
       iBegOR, jBegOr, kBegOr, nbkGlobal, nOrphans
  !use overset, only : 
  use stencils, only : visc_drdw_stencil, N_visc_drdw
  use communication, only : myid, adflow_comm_world
  use utils, only : setPointers, EChk
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  integer(kind=intType), intent(out) :: totalOrphans

  ! Working
  integer(kind=intType) :: i, j, k, nn, ii, jj, kk, n, ierr
  integer(kind=intType) :: magic, localOrphans, i_stencil

  magic = 33
  localOrphans = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)

     ! On the first pass count up the total number of orphans for this block
     n = 0
     do k=2, kl
        do j=2, jl
           do i=2, il
              if (iblank(i,j,k) == 0 .or. iblank(i,j,k)==-2 .or. iblank(i,j,k)==-3) then 

                 stencilLoop: do i_stencil=1, N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (ii >= 2 .and. jj >= 2 .and. kk>=2 .and. &
                         ii <= il .and. jj <= jl .and. kk <= kl) then 

                       if (iBlank(ii, jj, kk) == 1) then 
                          ! This cell is an orphan:
                          print *,'Error in connectivity: ', nbkglobal, i+iBegOR-1, j+jBegOr-1, k+kBegOr-1
                          n = n + 1
                       end if
                    end if
                 end do stencilLoop
              end if
           end do
        end do
     end do

     localOrphans = localOrphans + n

     ! Remove any existing orphans 
     if (associated(flowDoms(nn, level, sps)%orphans)) then 
        deallocate(flowDoms(nn, level, sps)%orphans)
     end if
     allocate(flowDoms(nn, level, sps)%orphans(3, n))

     ! Save the total number of orphans on this block
     flowDoms(nn, level, sps)%nOrphans = n

     ! Manual set information from blockPointers that would be set
     ! with setPointers()
     orphans => flowDoms(nn, level, sps)%orphans
     nOrphans = n

     ! On the first pass count up the total number of orphans for this block
     n = 0
     do k=2, kl
        do j=2, jl
           do i=2, il
              if (iblank(i,j,k) == 0 .or. iblank(i,j,k)==-2 .or. iblank(i,j,k)==-3) then 

                 stencilLoop2: do i_stencil=1, N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (ii >= 2 .and. jj >= 2 .and. kk>=2 .and. &
                         ii <= il .and. jj <= jl .and. kk <= kl) then 

                       if (iBlank(ii, jj, kk) == 1) then 
                          ! This cell is an orphan:
                          n = n + 1
                          orphans(:, n) = (/ii, jj, kk/)
                       end if
                    end if
                 end do stencilLoop2
              end if
           end do
        end do
     end do
  end do

  ! Determine the total number of overset orphans
  call mpi_allreduce(localOrphans, totalOrphans, 1, adflow_integer, MPI_SUM, &
       adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  if (myid == 0) then 
     print *, 'Total number of orphans:', totalOrphans
  end if
 

end subroutine checkOverset

