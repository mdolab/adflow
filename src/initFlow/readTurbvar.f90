!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbvar.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-21-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbvar(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbvar controls the reading of the turbulent variables    *
!      * for a restart. It calls the routine, which corresponds to the  *
!      * turbulence model used.                                         *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputPhysics
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType), intent(inout) :: nTypeMismatch
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if the rans equations must be solved. If not return.

       if(equations /= RANSEquations) return

       ! Determine the turbulence model to be used and call the
       ! appropriate subroutine.

       select case (turbModel)

         case (baldwinLomax)

         !===============================================================

         case (spalartAllmaras, spalartAllmarasEdwards)
           call readTurbSA(nTypeMismatch)

         !===============================================================

         case (komegaWilcox, komegaModified, menterSST, ktau)
           call readTurbKwType(nTypeMismatch)

         !===============================================================

         case (v2f)
           call readTurbV2f(nTypeMismatch)

         !===============================================================

         case default
           if(myID == 0) &
             call terminate("readTurbvar", "Restart not implemented &
                            &for this turbulence model.")
           call mpi_barrier(SUmb_comm_world, ierr)

       end select

       end subroutine readTurbvar
