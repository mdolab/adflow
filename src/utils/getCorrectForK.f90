function getCorrectForK()

  use flowVarRefState
  use inputPhysics
  use iteration
  implicit none

  logical :: getCorrectForK

  if( kPresent .and. currentLevel <= groundLevel) then
     getCorrectForK = .true.
  else
     getCorrectForK = .false.
  end if
end function getCorrectForK
