subroutine returnFail(routineName, errorMessage)
  !
  !       returnFail writes an error message to standard output and       
  !       returnFails the execution of the program.                       
  !
  use precision
  use communication
  use constants
  implicit none
  !
  !      Subroutine arguments
  !
  character(len=*), intent(in) :: routineName
  character(len=*), intent(in) :: errorMessage

end subroutine returnFail
