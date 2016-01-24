subroutine fringeReduction(level, sps)

  use blockPointers
  use stencils
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(kind=intType) :: i, j, k, nn, ii, jj, kk, i_stencil
  logical :: computeCellFound

  do nn=1, nDom
     call setPointers(nn, level, sps)

     do k=2, kl
        do j=2, jl
           do i=2, il

              ! Check if this cell is a fringe:
              if (fringes(i, j, k)%donorProc /= -1) then 

                 computeCellFound = .False.

                 stencilLoop2: do i_stencil=1, N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (fringes(ii, jj, kk)%isCompute) then 
                       ! This is a compute cell
                       computeCellFound = .True.
                    end if
                 end do stencilLoop2

                 if (.not. computeCellFound) then 
                    ! This cell is a hole no compute cell
                    ! surrounding a fringe, we can hard iblank it.
                    call emptyFringe(fringes(i, j, k))
                    fringes(i, j, k)%isHole = .True.
                    fringes(i, j, k)%isCompute = .False.
                 end if
              end if
           end do
        end do
     end do
  end do

end subroutine fringeReduction
