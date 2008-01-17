!
!      ******************************************************************
!      *                                                                *
!      * File:          convergenceHeader.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-27-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine convergenceHeader
!
!      ******************************************************************
!      *                                                                *
!      * convergenceHeader writes the convergence header to stdout.     *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use inputPhysics
       use inputUnsteady
       use flowVarRefState
       use monitor
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, nCharWrite

       logical :: writeIterations
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not the iterations must be written.

       writeIterations = .true.
       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) writeIterations = .false.

       ! Determine the number of characters to write.
       ! First initialize this number with the variables which are
       ! always written. This depends on the equation mode. For unsteady
       ! and spectral computations a bit more info is written.

       nCharWrite = 9
       if( writeIterations ) nCharWrite = nCharWrite + 7
       if(equationMode == unsteady) then
         nCharWrite = nCharWrite + 7 + fieldWidth + 1
       else if(equationMode == timeSpectral) then
         nCharWrite = nCharWrite + 11
       endif

       ! Add the number of characters needed for the actual variables.

       nCharWrite = nCharWrite + nMon*(fieldWidth+1)
       if( showCPU ) nCharWrite = nCharWrite + fieldWidth + 1

       ! Write the line of - signs. This line starts with a #, such
       ! that it is ignored by some plotting software.

       write(*,"(a)",advance="no") "#"
       do i=2,nCharWrite
         write(*,"(a)",advance="no") "-"
       enddo
       print "(1x)"

       ! Write the first line of the header. First the variables that
       ! will always be written. Some extra variables must be written
       ! for unsteady and time spectral problems.

       write(*,"(a)",advance="no") "#  Grid |"

       if(equationMode == unsteady) then
         write(*,"(a)",advance="no") " Time |    Time    |"
       else if(equationMode == timeSpectral) then
         write(*,"(a)",advance="no") " Spectral |"
       endif

       if( writeIterations ) write(*,"(a)",advance="no") " Iter |"
       if( showCPU )         write(*,"(a)",advance="no") "    Wall    |"

       ! Write the header for the variables to be monitored.

       do i=1,nMon

         ! Determine the variable name and write the
         ! corresponding text.

         select case (monNames(i))

           case (cgnsL2resRho)
             write(*,"(a)",advance="no") "   Res rho  |"

           case (cgnsL2resMomx)
             write(*,"(a)",advance="no") "  Res rhou  |"

           case (cgnsL2resMomy)
             write(*,"(a)",advance="no") "  Res rhov  |"

           case (cgnsL2resMomz)
             write(*,"(a)",advance="no") "  Res rhow  |"

           case (cgnsL2resRhoe)
             write(*,"(a)",advance="no") "  Res rhoE  |"

           case (cgnsL2resNu)
             write(*,"(a)",advance="no") " Res nuturb |"

           case (cgnsL2resK)
             write(*,"(a)",advance="no") "  Res kturb |"

           case (cgnsL2resOmega)
             write(*,"(a)",advance="no") "  Res wturb |"

           case (cgnsL2resTau)
             write(*,"(a)",advance="no") " Res tauturb|"

           case (cgnsL2resEpsilon)
             write(*,"(a)",advance="no") " Res epsturb|"

           case (cgnsL2resV2)
             write(*,"(a)",advance="no") "  Res v2turb|"

           case (cgnsL2resF)
             write(*,"(a)",advance="no") "  Res fturb |"

           case (cgnsCl)
             write(*,"(a)",advance="no") "   C_lift   |"

           case (cgnsClp)
             write(*,"(a)",advance="no") "  C_lift_p  |"

           case (cgnsClv)
             write(*,"(a)",advance="no") "  C_lift_v  |"

           case (cgnsCd)
             write(*,"(a)",advance="no") "   C_drag   |"

           case (cgnsCdp)
             write(*,"(a)",advance="no") "  C_drag_p  |"

           case (cgnsCdv)
             write(*,"(a)",advance="no") "  C_drag_v  |"

           case (cgnsCfx)
             write(*,"(a)",advance="no") "    C_Fx    |"

           case (cgnsCfy)
             write(*,"(a)",advance="no") "    C_Fy    |"

           case (cgnsCfz)
             write(*,"(a)",advance="no") "    C_Fz    |"

           case (cgnsCmx)
             write(*,"(a)",advance="no") "    C_Mx    |"

           case (cgnsCmy)
             write(*,"(a)",advance="no") "    C_My    |"

           case (cgnsCmz)
             write(*,"(a)",advance="no") "    C_Mz    |"

           case (cgnsHdiffMax)
             write(*,"(a)",advance="no") "  |H-H_inf| |"

           case (cgnsMachMax)
             write(*,"(a)",advance="no") "  Mach_max  |"

           case (cgnsYplusMax)
             write(*,"(a)",advance="no") "   Y+_max   |"

           case (cgnsEddyMax)
             write(*,"(a)",advance="no") "  Eddyv_max |"

         end select

       enddo

       print "(1x)"

       ! Write the second line of the header. Most of them are empty,
       ! but some variables require a second line.

       write(*,"(a)",advance="no")   "# level |"

       if(equationMode == unsteady) then
         write(*,"(a)",advance="no") " Step |            |"
       else if(equationMode == timeSpectral) then
         write(*,"(a)",advance="no") " Solution |"
       endif

       if( writeIterations ) write(*,"(a)",advance="no") "      |"
       if( showCPU )         write(*,"(a)",advance="no") " Clock (s)  |"

       ! Loop over the variables to be monitored and write the
       ! second line.

       do i=1,nMon

         ! Determine the variable name and write the
         ! corresponding text.

         select case (monNames(i))

           case (cgnsHdiffMax)
             write(*,"(a)",advance="no") "     max    |"

           case default
             write(*,"(a)",advance="no") "            |"

         end select

       enddo

       print "(1x)"

       ! Write again a line of - signs (starting with a #).

       write(*,"(a)",advance="no") "#"
       do i=2,nCharWrite
         write(*,"(a)",advance="no") "-"
       enddo
       print "(1x)"

       end subroutine convergenceHeader
