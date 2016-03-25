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
  integer(kind=intType) :: i, j, k, nn, ii, jj, kk
  integer(kind=intType) :: magic, ibval
  integer(kind=intType) :: i_stencil

  magic = 33
  do nn=1, nDom
     call setPointers(nn, level, sps)

     do k=2, kl
        do j=2, jl
           do i=2, il
              if (iblank(i,j,k) == 1) then 

                 ibval = 0

                 stencilLoop: do i_stencil=1, N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (iBlank(ii, jj, kk) == 1 .or. &
                         iBlank(ii, jj, kk) == -1) then 

                       ibval = ibval + 1
                    end if

                 end do stencilLoop
                 if (ibval /= magic) then
                    print *,'Error in connectivity: ', nbkglobal, i+iBegOR-1, j+jBegOr-1, k+kBegOr-1, ibval
                    !stop
                 end if
              end if
           end do
        end do
     end do
  end do
end subroutine checkOverset
