!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbTreatment.f90                             *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 06-13-2003                                      *
!      * Last modified: 08-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbTreatment
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbTreatment sets the arrays bmti1, bvti1, etc, such that   *
!      * the physical boundary conditions are treated correctly.        *
!      * It is assumed that the variables in blockPointers already      *
!      * point to the correct block.                                    *
!      *                                                                *
!      * The turbulent variable in the halo is computed as follows:     *
!      * wHalo = -bmt*wInternal + bvt for every block facer. As it is   *
!      * possible to have a coupling in the boundary conditions bmt     *
!      * actually are matrices. If there is no coupling between the     *
!      * boundary conditions of the turbulence equations bmt is a       *
!      * diagonal matrix.                                               *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Local variable.
!
       integer(kind=intType) :: nn, i, j, k, l, m
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the arrays for the boundary condition treatment
       ! to zero, such that internal block boundaries are solved
       ! correctly (i.e. explicitly).

       do k=1,ke
         do j=1,je
           do l=nt1,nt2
             do m=nt1,nt2
               bmti1(j,k,l,m) = zero
               bmti2(j,k,l,m) = zero
             enddo
             bvti1(j,k,l) = zero
             bvti2(j,k,l) = zero
           enddo
         enddo
       enddo

       do k=1,ke
         do i=1,ie
           do l=nt1,nt2
             do m=nt1,nt2
               bmtj1(i,k,l,m) = zero
               bmtj2(i,k,l,m) = zero
             enddo
             bvtj1(i,k,l) = zero
             bvtj2(i,k,l) = zero
           enddo
         enddo
       enddo

       do j=1,je
         do i=1,ie
           do l=nt1,nt2
             do m=nt1,nt2
               bmtk1(i,j,l,m) = zero
               bmtk2(i,j,l,m) = zero
             enddo
             bvtk1(i,j,l) = zero
             bvtk2(i,j,l) = zero
           enddo
         enddo
       enddo
 
       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Determine the kind of boundary condition for this subface.

         typeBC: select case (BCType(nn))

           case (NSWallAdiabatic, NSWallIsothermal)

             ! Viscous wall. There is no difference between an adiabatic
             ! and an isothermal wall for the turbulent equations.
             ! Set the implicit treatment of the wall boundary conditions.

             call bcTurbWall(nn)

           !=============================================================

           case (SubsonicInflow, SupersonicInflow, MassBleedInflow)

             ! Inflow. Subsonic, supersonic or mass bleed inflow is
             ! identical for the turbulent transport equations.

             call bcTurbInflow(nn)

           !=============================================================

           case (SubsonicOutflow,  SupersonicOutflow, &
                 MassBleedOutflow, Extrap)

             ! Outflow. Subsonic, supersonic or mass bleed outflow is
             ! identical for the turbulent transport equations. The
             ! extrapolation boundary is also an outflow.

             call bcTurbOutflow(nn)

           !=============================================================

           case (Symm, SymmPolar, EulerWall)

             ! Symmetry, polar symmetry or inviscid wall. Treatment of
             ! the turbulent equations is identical.

             call bcTurbSymm(nn)

           !=============================================================

           case (FarField)

             ! Farfield. The kind of boundary condition to be applied,
             ! inflow or outflow, depends on the local conditions.

             call bcTurbFarfield(nn)

           !=============================================================

           case (SlidingInterface,   OversetOuterBound,     &
                 DomainInterfaceAll, DomainInterfaceRhoUVW, &
                 DomainInterfaceP,   DomainInterfaceRho,    &
                 DomainInterfaceTotal)

             ! Sliding mesh interface, overset outer boudaries, and 
             ! domain interface with another code are not really boundary
             ! condition and therefore the values are kept.

             call bcTurbInterface(nn)

           !=============================================================

           case default

             call terminate("bcTurbTreatment", &
                            "Unknown boundary condition")

         end select typeBC

       enddo bocos

       end subroutine bcTurbTreatment
