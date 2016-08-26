subroutine checkPartitioning(np,load_inbalance,face_inbalance)

  ! This subroutine runs the load balancing and partitioning algorithm
  ! to determine what the load balancing will be for a given number of
  ! procs np. The output is load_inbalance and face_inbalance. 

  use constants
  use communication, only : nProc
  use partitionMod, only : ubvec
  implicit none
 
  integer(kind=intType),intent(in) ::np
  integer(kind=intType) :: nproc_save
  real(kind=realType),intent(out) :: load_inbalance,face_inbalance

  ! Note: This file follows mostly partitionAndReadGrid. See
  ! partitionAndReadGrid.f90 for more infromation
  
  ! Trick it into thinking we have np processors:
  nproc_save = nproc
  nproc = np
  
  call blockDistribution
  
  ! Restore the number of procs
  nproc = nproc_save
  
  ! Extract the inbalance info:
  load_inbalance = ubvec(1)
  face_inbalance = ubvec(2)

end subroutine checkPartitioning
