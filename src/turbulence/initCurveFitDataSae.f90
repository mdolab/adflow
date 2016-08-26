       subroutine initCurveFitDataSae
!
!       initCurveFitDataSae contains the curve fit constants for       
!       the wall function data for the Spalart-Allmaras (Edwards       
!       modification) turbulence model.                                
!
       use flowVarRefState
       use paramTurb
       use utils, only : terminate
       implicit none
!
!      Local variables.
!
   !   integer :: ierr

       call terminate("initCurveFitDataSae", &
                      "Not implemented yet")

       end subroutine initCurveFitDataSae
