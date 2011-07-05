!
!      ******************************************************************
!      *                                                                *
!      * file:          tdbc.f90                                        *
!      * author:        Eran Arad                                       *
!      * starting date: 10-25-2004                                      *
!      * last modified: 04-20-2009                                      *
!      *                                                                *
!      ******************************************************************
!
      subroutine tdbc(timeUnsteady,rho,velM)
!
!      ******************************************************************
!      *                                                                *
!      * time dependent bc for subsonic inlet in unsteady simulation.   *
!      *                                                                *
!      *   U_i= V_i* ( C0Tdbc + sin(2*pi*oscilFreq*time + cPhaseR) )    *
!      *  NOTE:  inflow rho < 0 overites cPhase = pi                    *
!      * output here is velxM which multiplies bcdata%velx              *
!      ******************************************************************
!
        use constants
        use inputTDBC
        implicit none
!
! subroutine arguments.
!
        real(kind=realType) :: timeUnsteady, rho, velM
!
! local variables:
!
        real(kind=realType) :: omegaT
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
        if(rho < 0 )then
           cPhase = pi
        end if

        omegaT = 2.0*pi*oscillFreq*timeUnsteady
        velM = C0Tdbc  + sin(omegaT + cPhaseR)

        return
        end
