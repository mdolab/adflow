!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Reduce.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-11-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Reduce(sendbuf, recvbuf, count, datatype, op, &
                             root, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Reduce simply copies the send buffer into the receive      *
!      * buffer, which means that not all operations defined in MPI are *
!      * available in sequential mode.                                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer                 :: count, datatype, op, root, comm, error
       character, dimension(*) :: sendbuf, recvbuf
!
!      Local variables
!
       integer :: i, nbytes_copy
!
!      Function definition
!
       integer :: determine_size_entity
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the reduce operation to be performed.

       if(op == mpi_max .or. op == mpi_min .or. op == mpi_sum .or. &
          op == mpi_prod) then

         ! For sequential mode all these operations reduce to a copy of
         ! the send buffer in the receive buffer. Determine the number of
         ! bytes to be copied and do so.

         nbytes_copy = count*determine_size_entity(datatype)
         do i=1,nbytes_copy
           recvbuf(i) = sendbuf(i)
         enddo

       else
         ! Reduce operation not supported yet.

         call su_mpi_terminate("MPI_Reduce", &
                               "Reduce operation not supported &
                               &yet in sequential mode")
       endif

       ! Set error to mpi_success to indicate that everything
       ! went okay.

       error = mpi_success

       end subroutine MPI_Reduce
