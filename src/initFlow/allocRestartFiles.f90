subroutine allocRestartFiles(nFiles)
  !
  !          Allocate memory for the restartfles                         *                                   
  !          The array is populated from Python using setRestartFiles    
  !          If memeory has been allocated for the array there exist at  
  !          least one element in the array.                             
  !
  use constants
  use inputIO, only : restartFiles
  use utils, only : terminate
  implicit none
  !
  !      Subroutine argument.
  !
  integer(kind=intType) :: nFiles
  !
  !      Local variables.
  !
  integer :: ierr

  if (allocated(restartFiles)) then
     deallocate(restartFiles)
  end if

  allocate(restartFiles(nFiles), stat=ierr)
  if(ierr /= 0)                          &
    call terminate("allocRestartFiles", &
            "Memory allocation failure for restartFiles")
  
  ! Zero Array with empty strings
  restartFiles = ""

end subroutine allocRestartFiles

