!       This routine determines if there are any oveset boundaries     
!      * present in the mesh. 

function oversetPresent()

  use blockPointers
  use communication
  implicit none

  ! Function
  logical :: oversetPresent, local

  ! Working
  integer(Kind=intType) :: nn, mm, ierr

  local = .False.
  do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)

     do mm=1, nBocos
        if (BCType(mm) == OversetOuterBound) then
           local = .True.
        end if
     end do
  end do

  call mpi_allreduce(local, oversetPresent, 1, MPI_LOGICAL, MPI_LOR, ADflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

end function oversetPresent
