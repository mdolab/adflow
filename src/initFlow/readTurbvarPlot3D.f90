!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbvarPlot3D.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 07-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbvarPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readTurbvarPlot3D controls the reading of the turbulent        *
!      * variables for a restart. It calls the routine, which           *
!      * corresponds to the  turbulence model used.                     *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputPhysics
       implicit none
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
           call readTurbSAPlot3D

         !===============================================================

         case (komegaWilcox, komegaModified, menterSST, ktau)
           call readTurbKwTypePlot3D

         !===============================================================

         case (v2f)
           call readTurbV2fPlot3D

         !===============================================================

         case default
           if(myID == 0)                         &
             call terminate("readTurbvarPlot3D", &
                            "Restart not implemented for this &
                            &turbulence model.")
           call mpi_barrier(SUmb_comm_world, ierr)

       end select

       end subroutine readTurbvarPlot3D
