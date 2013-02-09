module stencils
  
  !      ******************************************************************
  !      *                                                                *
  !      * stencils defines indices for several types of stencils. These  *
  !      * are useful for setting block in dRdw,dRdwPre or dRdx depending *
  !      * on the type of equations  being solved                         *
  !      *                                                                *
  !      ******************************************************************

  use precision
  implicit none
  ! First define the sizes of the stncils
  integer(kind=intType), parameter :: N_euler_pc   = 7
  integer(kind=intType), parameter :: N_euler_drdw = 13
  integer(kind=intType), parameter :: N_euler_drdx = 32
  integer(kind=intType), parameter :: N_visc_pc     = 21

  integer(kind=intType), dimension(7 ,3),target :: euler_pc_stencil
  integer(kind=intType), dimension(13,3),target :: euler_drdw_stencil
  integer(kind=intType), dimension(32,3),target :: euler_drdx_stencil
  integer(kind=intType), dimension(21,3),target :: visc_pc_stencil
  
end module stencils
       
subroutine initialize_stencils

  use stencils
  integer(kind=intType) :: i,j,k,counter

  ! Euler PC Stencil

  euler_pc_stencil(1,:) = (/ 0, 0, 0 /)
  euler_pc_stencil(2,:) = (/-1, 0, 0 /)
  euler_pc_stencil(3,:) = (/ 1, 0, 0 /)
  euler_pc_stencil(4,:) = (/ 0,-1, 0 /)
  euler_pc_stencil(5,:) = (/ 0, 1, 0 /)
  euler_pc_stencil(6,:) = (/ 0, 0,-1 /)
  euler_pc_stencil(7,:) = (/ 0, 0, 1 /)

 ! Euler drdw Stencil

  euler_drdw_stencil(1 ,:) = (/ 0, 0, 0 /)
  euler_drdw_stencil(2 ,:) = (/-2, 0, 0 /)
  euler_drdw_stencil(3 ,:) = (/-1, 0, 0 /)
  euler_drdw_stencil(4 ,:) = (/ 1, 0, 0 /)
  euler_drdw_stencil(5 ,:) = (/ 2, 0, 0 /)
  euler_drdw_stencil(6 ,:) = (/ 0,-2, 0 /)
  euler_drdw_stencil(7 ,:) = (/ 0,-1, 0 /)
  euler_drdw_stencil(8 ,:) = (/ 0, 1, 0 /)
  euler_drdw_stencil(9 ,:) = (/ 0, 2, 0 /)
  euler_drdw_stencil(10,:) = (/ 0, 0,-2 /)
  euler_drdw_stencil(11,:) = (/ 0, 0,-1 /)
  euler_drdw_stencil(12,:) = (/ 0, 0, 1 /)
  euler_drdw_stencil(13,:) = (/ 0, 0, 2 /)

  ! Euler dRdx Stencil
  euler_drdx_stencil(1 ,:) = (/   -1,  0,  0 /)
  euler_drdx_stencil(2 ,:) = (/   -1,  0,  1 /)
  euler_drdx_stencil(3 ,:) = (/   -1,  1,  0 /)
  euler_drdx_stencil(4 ,:) = (/   -1,  1,  1 /)
  euler_drdx_stencil(5 ,:) = (/    0, -1,  0 /)
  euler_drdx_stencil(6 ,:) = (/    0, -1,  1 /)
  euler_drdx_stencil(7 ,:) = (/    0,  0, -1 /)
  euler_drdx_stencil(8 ,:) = (/    0,  0,  0 /)
  euler_drdx_stencil(9 ,:) = (/    0,  0,  1 /)
  euler_drdx_stencil(10,:) = (/    0,  0,  2 /)
  euler_drdx_stencil(11,:) = (/    0,  1, -1 /)
  euler_drdx_stencil(12,:) = (/    0,  1,  0 /)
  euler_drdx_stencil(13,:) = (/    0,  1,  1 /)
  euler_drdx_stencil(14,:) = (/    0,  1,  2 /)
  euler_drdx_stencil(15,:) = (/    0,  2,  0 /)
  euler_drdx_stencil(16,:) = (/    0,  2,  1 /)
  euler_drdx_stencil(17,:) = (/    1, -1,  0 /)
  euler_drdx_stencil(18,:) = (/    1, -1,  1 /)
  euler_drdx_stencil(19,:) = (/    1,  0, -1 /)
  euler_drdx_stencil(20,:) = (/    1,  0,  0 /)
  euler_drdx_stencil(21,:) = (/    1,  0,  1 /)
  euler_drdx_stencil(22,:) = (/    1,  0,  2 /)
  euler_drdx_stencil(23,:) = (/    1,  1, -1 /)
  euler_drdx_stencil(24,:) = (/    1,  1,  0 /)
  euler_drdx_stencil(25,:) = (/    1,  1,  1 /)
  euler_drdx_stencil(26,:) = (/    1,  1,  2 /)
  euler_drdx_stencil(27,:) = (/    1,  2,  0 /)
  euler_drdx_stencil(28,:) = (/    1,  2,  1 /)
  euler_drdx_stencil(29,:) = (/    2,  0,  0 /)
  euler_drdx_stencil(30,:) = (/    2,  0,  1 /)
  euler_drdx_stencil(31,:) = (/    2,  1,  0 /)
  euler_drdx_stencil(32,:) = (/    2,  1,  1 /)

  ! Viscous 3x3x3 stencil:
  counter = 0
!   do i=-1,1
!      do j=-1,1
!         do k=-1,1
!            counter = counter + 1
!            visc_pc_stencil(counter,:) = (/i,j,k/)
!         end do
!      end do
!   end do

  ! K = 0 plane
  visc_pc_stencil(1,:) = (/-1,-1, 0/)
  visc_pc_stencil(2,:) = (/ 0,-1, 0/)
  visc_pc_stencil(3,:) = (/ 1,-1, 0/)
  visc_pc_stencil(4,:) = (/-1, 0, 0/)
  visc_pc_stencil(5,:) = (/ 0, 0, 0/)
  visc_pc_stencil(6,:) = (/ 1, 0, 0/)
  visc_pc_stencil(7,:) = (/-1, 1, 0/)
  visc_pc_stencil(8,:) = (/ 0, 1, 0/)
  visc_pc_stencil(9,:) = (/ 1, 1, 0/)

  ! K=1 top star
  visc_pc_stencil(10,:) = (/-1, 0, 1/)
  visc_pc_stencil(11,:) = (/ 0, 0, 1/)
  visc_pc_stencil(12,:) = (/ 1, 0, 1/)
  visc_pc_stencil(13,:) = (/ 0,-1, 1/)
  visc_pc_stencil(14,:) = (/ 0, 1, 1/)

  ! K=-1 bottom star
  visc_pc_stencil(15,:) = (/-1, 0,-1/)
  visc_pc_stencil(16,:) = (/ 0, 0,-1/)
  visc_pc_stencil(17,:) = (/ 1, 0,-1/)
  visc_pc_stencil(18,:) = (/ 0,-1,-1/)
  visc_pc_stencil(19,:) = (/ 0, 1,-1/)

end subroutine initialize_stencils

