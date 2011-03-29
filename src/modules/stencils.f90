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
 
  integer(kind=intType), parameter :: N_rans_pc   = 19
  integer(kind=intType), parameter :: N_rans_drdw = 25

  integer(kind=intType), dimension(7,3),target :: euler_pc_stencil
  integer(kind=intType), dimension(19,3),target :: rans_pc_stencil

  integer(kind=intType), dimension(13,3),target :: euler_drdw_stencil

end module stencils
       

subroutine initialize_stencils

  use stencils

  ! Euler PC Stencil

  euler_pc_stencil(1,:) = (/ 0, 0, 0 /)
  euler_pc_stencil(2,:) = (/-1, 0, 0 /)
  euler_pc_stencil(3,:) = (/ 1, 0, 0 /)
  euler_pc_stencil(4,:) = (/ 0,-1, 0 /)
  euler_pc_stencil(5,:) = (/ 0, 1, 0 /)
  euler_pc_stencil(6,:) = (/ 0, 0,-1 /)
  euler_pc_stencil(7,:) = (/ 0, 0, 1 /)

 ! Euler PC Stencil

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

  ! NS PC Stencil - All of the 3x3 ECEPT corners

  rans_pc_stencil(1,:) = (/ 0, 0, 0 /)

  rans_pc_stencil(2,:) = (/-1, 0, 0 /)
  rans_pc_stencil(3,:) = (/ 1, 0, 0 /)
  rans_pc_stencil(4,:) = (/ 0,-1, 0 /)
  rans_pc_stencil(5,:) = (/ 0, 1, 0 /)
  rans_pc_stencil(6,:) = (/ 0, 0,-1 /)
  rans_pc_stencil(7,:) = (/ 0, 0, 1 /)

  rans_pc_stencil(8 ,:) = (/-1,-1, 0 /)
  rans_pc_stencil(9 ,:) = (/-1, 1, 0 /)
  rans_pc_stencil(10,:) = (/ 1,-1, 0 /)
  rans_pc_stencil(11,:) = (/ 1,-1, 0 /)

  rans_pc_stencil(12,:) = (/ 0,-1,-1 /)
  rans_pc_stencil(13,:) = (/ 0,-1, 1 /)
  rans_pc_stencil(14,:) = (/ 0, 1,-1 /)
  rans_pc_stencil(15,:) = (/ 0, 1, 1 /)
  rans_pc_stencil(16,:) = (/-1, 0,-1 /)
  rans_pc_stencil(17,:) = (/-1, 0, 1 /)
  rans_pc_stencil(18,:) = (/ 1, 0,-1 /)
  rans_pc_stencil(19,:) = (/ 1, 0, 1 /)

end subroutine initialize_stencils

