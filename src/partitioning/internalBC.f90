!
!      ******************************************************************
!      *                                                                *
!      * File:          internalBC.f90                                  *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 01-30-2003                                      *
!      * Last modified: 08-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function internalBC(cgnsBocoType, userDefinedName)
!
!      ******************************************************************
!      *                                                                *
!      * internalBC determines the corresponding internally used        *
!      * boundary condition type for the given CGNS boundary condition. *
!      * The flow equations to be solved are taken into account, e.g.   *
!      * a viscous wall BC for the Euler equations is set to an         *
!      * inviscid wall.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use constants
       use inputPhysics
       use su_cgns
       implicit none
!
!      Function type.
!
       integer(kind=intType) :: internalBC
!
!      Function argument.
!
       integer, intent(in) :: cgnsBocoType  ! Note integer and not
                                            ! integer(intType).
                                            ! Because of cgns.

       character(len=maxCGNSNameLen), intent(in) :: userDefinedName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the CGNS boundary condition type and set
       ! internalBC accordingly.

       select case (cgnsBocoType)
         case (BCWallInviscid)
           internalBC = EulerWall

         case (BCWall, BCWallViscous, BCWallViscousHeatFlux)
           internalBC = NSWallAdiabatic
           if(equations == EulerEquations) internalBC = EulerWall

         case (BCWallViscousIsothermal)
           internalBC = NSWallIsothermal
           if(equations == EulerEquations) internalBC = EulerWall

         case (BCSymmetryPlane)
           internalBC = Symm

         case (BCSymmetryPolar)
           internalBC = SymmPolar

         case (BCExtrapolate, BCDegenerateLine, BCDegeneratePoint, &
               BCAxisymmetricWedge)
           internalBC = Extrap

         case (BCFarfield, BCInflow, BCOutflow)
           internalBC = FarField

         case (BCInflowSubsonic)
           internalBC = SubsonicInflow

         case (BCInflowSupersonic)
           internalBC = SupersonicInflow

         case (BCOutflowSubsonic)
           internalBC = SubsonicOutflow

         case (BCOutflowSupersonic)
           internalBC = SupersonicOutflow

         case (BCTunnelInflow)
           internalBC = SubsonicInflow

         case (BCTunnelOutflow)
           internalBC = SubsonicOutflow

         case (UserDefined)

           ! Select the internal type base on the user defined name.

           select case (trim(adjustl(userDefinedName)))

             case ("BCMassBleedInflow")
               internalBC = MassBleedInflow

             case ("BCMassBleedOutflow")
               internalBC = MassBleedOutflow

             case ("BCSlidingMesh")
               internalBC = SlidingInterface

             case ("BCOverset")
               internalBC = OversetOuterBound

             case ("BCDomainInterfaceAll")
               internalBC = DomainInterfaceAll

             case ("BCDomainInterfaceRhoUVW")
               internalBC = DomainInterfaceRhoUVW

             case ("BCDomainInterfaceP")
               internalBC = DomainInterfaceP

             case ("BCDomainInterfaceRho")
               internalBC = DomainInterfaceRho

             case ("BCDomainInterfaceTotal")
               internalBC = DomainInterfaceTotal

             case default
               internalBC = bcNull

           end select

         case default
           internalBC = bcNull

       end select

       ! Check whether the boundary conditions is allowed for the
       ! type of flow problem to be solved.

       select case (flowType)

         case (internalFlow)

           ! Internal flow. Not allowed to specify a farfield
           ! boundary condition.

           if(internalBC == FarField) internalBC = bcNotValid

       end select

       end function internalBC
