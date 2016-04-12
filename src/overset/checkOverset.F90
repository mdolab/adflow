subroutine checkOverset (level, sps)

  !
  !      ******************************************************************
  !      *                                                                *
  !      * CheckOverset checks the integrity of the overset connectivity  *
  !      * and holes. For every comptue cell (iblank = 1) it checks that  *
  !      * every cell in its stencil are not blanked. If even 1 cell is   *
  !      * found with an incomplete stencil it is a fatal error. 
  !      *                                                                *
  !      ******************************************************************

  use overset
  use blockPointers
  use stencils
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(kind=intType) :: i, j, k, nn, ii, jj, kk, n
  integer(kind=intType) :: magic
  integer(kind=intType) :: i_stencil

  magic = 33
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
end subroutine checkOverset

