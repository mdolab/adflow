       subroutine inviscidDissFluxCusp
!
!       inviscidDissFluxCusp computes the artificial dissipation for   
!       the fine grid, i.e. second order, cusp scheme for a given      
!       block. Therefore it is assumed that the pointers in            
!       blockPointers already point to the correct block.              
!
       use blockPointers
       use constants
       use inputDiscretization
       use iteration
       implicit none

       call returnFail("inviscidDissFluxCusp", "not implemented yet")

       end subroutine inviscidDissFluxCusp
