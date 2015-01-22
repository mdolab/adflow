module myPushPopLib

  ! This modules contains a fixed-sized integer stack that can be used
  ! to replace tapenade pushcontrol1b routines *AS LONG AS THEY ARE
  ! CONTAINED ONLY WITHIN A LOOP BODY* ie. only a fixed number of
  ! pushecontrol1b's are called before corresponding popcontrol1b. By
  ! modifiing the output code to make the code in-line, this is
  ! modification makes the reverse mode code faster. 

  use precision
  implicit none
  save

  integer(kind=intType) :: myIntStack(16)
  integer(kind=intType) :: myIntStackPtr=0

end module myPushPopLib
