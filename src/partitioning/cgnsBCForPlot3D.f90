!
!      ******************************************************************
!      *                                                                *
!      * File:          cgnsBCForPlot3D.f90                             *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 02-22-2005                                      *
!      * Last modified: 08-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine cgnsBCForPlot3D(BCNamePlot3D, BCTypeCGNS, BCType)
!
!      ******************************************************************
!      *                                                                *
!      * cgnsBCForPlot3D sets the equivalent internal and CGNS boundary *
!      * condition for the given plot3D BC name.                        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use su_cgns
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer,               intent(out) :: BCTypeCGNS
       integer(kind=intType), intent(out) :: BCType
       character(len=*),      intent(in)  :: BCNamePlot3D
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the case we are having here and set the return
       ! value accordingly.

       select case (BCNamePlot3D)

         case ("BCAxisymmetricWedge")
          BCTypeCGNS = BCAxisymmetricWedge
          BCType     = Extrap

         case ("BCDegenerateLine")
           BCTypeCGNS = BCDegenerateLine
           BCType     = Extrap

         case ("BCDegeneratePoint")
           BCTypeCGNS = BCDegeneratePoint
           BCType     = Extrap

         case ("BCDirichlet")
           BCTypeCGNS = BCDirichlet
           BCType     = BCNull

         case ("BCExtrapolate")
           BCTypeCGNS = BCExtrapolate
           BCType     = Extrap

         case ("BCFarfield")
           BCTypeCGNS = BCFarfield
           BCType     = Farfield

         case ("BCGeneral")
           BCTypeCGNS = BCGeneral
           BCType     = BCNull

         case ("BCInflow")
           BCTypeCGNS = BCInflow
           BCType     = Farfield

         case ("BCInflowSubsonic")
           BCTypeCGNS = BCInflowSubsonic
           BCType     = SubsonicInflow

         case ("BCInflowSupersonic")
           BCTypeCGNS = BCInflowSupersonic
           BCType     = SupersonicInflow

         case ("BCNeumann")
           BCTypeCGNS = BCNeumann
           BCType     = BCNull

         case ("BCOutflow")
           BCTypeCGNS = BCOutflow
           BCType     = Farfield

         case ("BCOutflowSubsonic")
           BCTypeCGNS = BCOutflowSubsonic
           BCType     = SubsonicOutflow

         case ("BCOutflowSupersonic")
           BCTypeCGNS = BCOutflowSupersonic
           BCType     = SupersonicOutflow

         case ("BCSymmetryPlane")
           BCTypeCGNS = BCSymmetryPlane
           BCType     = Symm

         case ("BCSymmetryPolar")
           BCTypeCGNS = BCSymmetryPolar
           BCType     = SymmPolar

         case ("BCTunnelInflow")
           BCTypeCGNS = BCTunnelInflow
           BCType     = SubsonicInflow

         case ("BCTunnelOutflow")
           BCTypeCGNS = BCTunnelOutflow
           BCType     = SubsonicOutflow

         case ("BCWall")
           BCTypeCGNS = BCWall
           BCType     = NSWallAdiabatic
           if(equations == EulerEquations) BCType = EulerWall

         case ("BCWallInviscid")
           BCTypeCGNS = BCWallInviscid
           BCType     = EulerWall

         case ("BCWallViscous")
           BCTypeCGNS = BCWallViscous
           BCType     = NSWallAdiabatic
           if(equations == EulerEquations) BCType = EulerWall

         case ("BCWallViscousHeatFlux")
           BCTypeCGNS = BCWallViscousHeatFlux
           BCType     = NSWallAdiabatic
           if(equations == EulerEquations) BCType = EulerWall

         case ("BCWallViscousIsothermal")
           BCTypeCGNS = BCWallViscousIsothermal
           BCType     = NSWallIsothermal
           if(equations == EulerEquations) BCType = EulerWall

         case ("BCMassBleedInflow")
           BCTypeCGNS = UserDefined
           BCType     = MassBleedInflow

         case ("BCMassBleedOutflow")
           BCTypeCGNS = UserDefined
           BCType     = MassBleedOutflow

         case ("BCSlidingMesh")
           BCTypeCGNS = UserDefined
           BCType     = SlidingInterface

         case ("BCOverset")
           BCTypeCGNS = UserDefined
           BCType     = OversetOuterBound

         case ("BCDomainInterfaceAll")
           BCTypeCGNS = UserDefined
           BCType     = DomainInterfaceAll

         case ("BCDomainInterfaceRhoUVW")
           BCTypeCGNS = UserDefined
           BCType     = DomainInterfaceRhoUVW

         case ("BCDomainInterfaceP")
           BCTypeCGNS = UserDefined
           BCType     = DomainInterfaceP

         case ("BCDomainInterfaceRho")
           BCTypeCGNS = UserDefined
           BCType     = DomainInterfaceRho

         case ("BCDomainInterfaceTotal")
           BCTypeCGNS = UserDefined
           BCType     = DomainInterfaceTotal

         case default
           BCTypeCGNS = Null
           BCType     = BCNull

       end select

       end subroutine cgnsBCForPlot3D
