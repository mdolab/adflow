!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataSubsonicInflow.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2004                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataSubsonicInflow(boco, allTurbPresent)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataSubsonicInflow tries to extract the prescribed data      *
!      * for the currently active boundary face, which is a subsonic    *
!      * inflow. Either total conditions and velocity direction or the  *
!      * velocity and density can be prescribed. In the latter case the *
!      * mass flow is prescribed, which is okay as long as the flow is  *
!      * not choked.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use flowVarRefState
       use inputPhysics
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: boco
       logical,               intent(inout) :: allTurbPresent
!
!      Local variables.
!
       integer :: ierr, nn

       logical :: ptPresent,   ttPresent,   rhotPresent
       logical :: axPresent,   ayPresent,   azPresent
       logical :: xdirPresent, ydirPresent, zdirPresent
       logical :: rdirPresent, tdirPresent
       logical :: velxPresent, velyPresent, velzPresent
       logical :: rhoPresent,  velrPresent, veltPresent
       logical :: allTurbSubface
       logical :: totPresent, velPresent, dirPresent

       character(len=maxStringLen) :: errorMessage
!
!      Interfaces
!
       interface
         logical function setBcVarTurb(offset, boco, turbInlet)
         use paramTurb
         use BCDataMod
         implicit none

         integer(kind=intType), intent(in) :: offset, boco
         real(kind=realType), dimension(:,:,:), pointer :: turbInlet

         end function setBcVarTurb
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the buffer bcVarArray, which is used
       ! for the interpolation and set the cgns names.

       nbcVar = 17
       if(equations == RANSEquations) nbcVar = nbcVar + nwt

       allocate(bcVarArray(iBeg:iEnd,jBeg:jEnd,nbcVar), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("BCDataSubsonicInflow", &
                        "Memory allocation failure for bcVarArray")

       bcVarNames(1)  = cgnsPtot
       bcVarNames(2)  = cgnsTtot
       bcVarNames(3)  = cgnsRhotot
       bcVarNames(4)  = cgnsVelAnglex
       bcVarNames(5)  = cgnsVelAngley
       bcVarNames(6)  = cgnsVelAnglez
       bcVarNames(7)  = cgnsVelVecx
       bcVarNames(8)  = cgnsVelVecy
       bcVarNames(9)  = cgnsVelVecz
       bcVarNames(10) = cgnsVelVecr
       bcVarNames(11) = cgnsVelVectheta
       bcVarNames(12) = cgnsDensity
       bcVarNames(13) = cgnsVelx
       bcVarNames(14) = cgnsVely
       bcVarNames(15) = cgnsVelz
       bcVarNames(16) = cgnsVelr
       bcVarNames(17) = cgnsVeltheta

       call setBcVarNamesTurb(17_intType)

       ! Try to determine these variables.

       call extractFromDataSet(BCFaceID(boco))

       ! Store the logicals, which indicate succes or failure
       ! a bit more readable.

       ptPresent   = bcVarPresent(1)
       ttPresent   = bcVarPresent(2)
       rhotPresent = bcVarPresent(3)
       axPresent   = bcVarPresent(4)
       ayPresent   = bcVarPresent(5)
       azPresent   = bcVarPresent(6)
       xdirPresent = bcVarPresent(7)
       ydirPresent = bcVarPresent(8)
       zdirPresent = bcVarPresent(9)
       rdirPresent = bcVarPresent(10)
       tdirPresent = bcVarPresent(11)
       rhoPresent  = bcVarPresent(12)
       velxPresent = bcVarPresent(13)
       velyPresent = bcVarPresent(14)
       velzPresent = bcVarPresent(15)
       velrPresent = bcVarPresent(16)
       veltPresent = bcVarPresent(17)

       ! Check if the total conditions are present.

       nn = 0
       if( ptPresent )   nn = nn + 1
       if( ttPresent )   nn = nn + 1
       if( rhotPresent ) nn = nn + 1

       totPresent = .false.
       if(nn >= 2) totPresent = .true.

       ! Check if a velocity direction is present.

       dirPresent = .false.
       if(xdirPresent .and. rdirPresent) dirPresent = .true.
       if((axPresent .or. xdirPresent) .and. &
          (ayPresent .or. ydirPresent) .and. &
          (azPresent .or. zdirPresent))  dirPresent = .true.

       ! Check if a velocity vector is present.

       velPresent = .false.
       if(velxPresent .and. velrPresent) velPresent = .true.
       if(velxPresent .and. velyPresent .and. velzPresent) &
         velPresent = .true.

       ! Determine the situation we have here.

       if(totPresent .and. dirPresent) then

         ! Total conditions and velocity direction are prescribed.
         ! Determine the values for the faces of the subface.

         call totalSubsonicInlet

       else if(rhoPresent .and. velPresent) then

         ! Density and velocity vector are prescribed, i.e. mass flow.
         ! Determine the values for the faces of the subface.

         call massflowSubsonicInlet

       else

         ! Not enough data is prescribed. Print an error message
         ! and exit.

         write(errorMessage,100)                   &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 100     format("Zone ",a,", boundary subface ",a, &
                ": Not enough data specified for subsonic inlet")

         call terminate("BCDataSubsonicInflow", errorMessage)

       endif

       ! Set the turbulence variables and check if all of them are
       ! prescribed. If not set allTurbPresent to .false.

       allTurbSubface = setBcVarTurb(17_intType, boco, &
                                     BCData(boco)%turbInlet)

       if(.not. allTurbSubface) allTurbPresent = .false.

       ! Release the memory of the bcVarArray.

       deallocate(bcVarArray, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("BCDataSubsonicInflow", &
                        "Deallocation failure for bcVarArray")

       !=================================================================

       contains

         !===============================================================

         subroutine totalSubsonicInlet
!
!        ****************************************************************
!        *                                                              *
!        * TotalSubsonicInlet converts the prescribed total           *
!        * conditions and velocity direction into a useable format.     *
!        *                                                              *
!        ****************************************************************
!
         use section
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: i, j, nn

         real(kind=realType) :: rhot, mult, trans
         real(kind=realType) :: ax, r1, r2, var, wax, wrad, wtheta

         real(kind=realType), dimension(3) :: xc, dir
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Set the subsonic inlet treatment to totalConditions.

         BCData(boco)%subsonicInletTreatment = totalConditions

         ! If the total pressure is present, convert it to SI-units and
         ! store it.

         if( ptPresent ) then
           call siPressure(mass(1), length(1), time(1), mult, trans)

           do j=jBeg,jEnd
             do i=iBeg,iEnd
               BCData(boco)%ptInlet(i,j) = mult*bcVarArray(i,j,1) &
                                         + trans
             enddo
           enddo
         endif

         ! If the total temperature is present, convert it to SI-units
         ! and store it.

         if( ttPresent ) then
           call siTemperature(temp(2), mult, trans)

           do j=jBeg,jEnd
             do i=iBeg,iEnd
               BCData(boco)%ttInlet(i,j) = mult*bcVarArray(i,j,2) &
                                         + trans
             enddo
           enddo
         endif

         ! Check if the total density is present. If so, it may be used
         ! to determine the total temperature or pressure if one of these
         ! variables was not specified.

         if( rhotPresent ) then
           call siDensity(mass(3), length(3), mult, trans)

           if(ptPresent .and. (.not. ttPresent)) then

             ! Total pressure is present but total temperature is not.
             ! Convert the total density to SI-units and use the perfect
             ! gas law to obtain the total temperature.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 rhot = mult*bcVarArray(i,j,3) + trans

                 BCData(boco)%ttInlet(i,j) = &
                            BCData(boco)%ptInlet(i,j)/(RGasDim*rhot)
               enddo
             enddo

           else if(ttPresent .and. (.not. ptPresent)) then

             ! Total temperature is present but total pressure is not.
             ! Convert the total density to SI-units and use the perfect
             ! gas law to obtain the total pressure.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 rhot = mult*bcVarArray(i,j,3) + trans

                 BCData(boco)%ptInlet(i,j) = RGasDim*rhot &
                                           * BCData(boco)%ttInlet(i,j)
               enddo
             enddo

           endif
         endif

         ! Determine the velocity direction. There are multiple
         ! possibilities to specify this direction.

         radialTest: if( rdirPresent ) then

           ! Radial direction specified, i.e. a cylindrical coordinate
           ! system is used for the velocity direction.

           ! Determine the unit vectors, which define the cylindrical
           ! coordinate system aligned with the rotation axis.

           call unitVectorsCylSystem(boco)

           ! Initialize wtheta to zero. This value will be used if no
           ! theta velocity component was specified.

           wtheta = zero

           ! Loop over the faces of the subface.

           do j=jBeg,jEnd
             do i=iBeg,iEnd

               ! Determine the coordinates of the face center relative to
               ! the rotation point of this section. Normally this is an
               ! average of i-1, i, j-1, j, but due to the usage of the
               ! pointer xf and the fact that x originally starts at 0,
               ! an offset of 1 is introduced and thus the average should
               ! be taken of i, i+1, j and j+1.

               xc(1) = fourth*(xf(i,j,  1) + xf(i+1,j,  1)  &
                     +         xf(i,j+1,1) + xf(i+1,j+1,1)) &
                     - sections(sectionId)%rotCenter(1)
               xc(2) = fourth*(xf(i,j,  2) + xf(i+1,j,  2)  &
                     +         xf(i,j+1,2) + xf(i+1,j+1,2)) &
                     - sections(sectionId)%rotCenter(2)
               xc(3) = fourth*(xf(i,j,  3) + xf(i+1,j,  3)  &
                     +         xf(i,j+1,3) + xf(i+1,j+1,3)) &
                     - sections(sectionId)%rotCenter(3)

               ! Determine the coordinates in the local cartesian frame,
               ! i.e. the frame determined by axis, radVec1 and radVec2.

               ax = xc(1)*axis(1)    + xc(2)*axis(2)    &
                  + xc(3)*axis(3)
               r1 = xc(1)*radVec1(1) + xc(2)*radVec1(2) &
                  + xc(3)*radVec1(3)
               r2 = xc(1)*radVec2(1) + xc(2)*radVec2(2) &
                  + xc(3)*radVec2(3)

               ! Determine the weights of the unit vectors in the local
               ! cylindrical system.

               wax  = bcVarArray(i,j,7)
               wrad = bcVarArray(i,j,10)
               if( tdirPresent ) wtheta = bcVarArray(i,j,11)

               ! Determine the direction in the local cartesian frame,
               ! determined by axis, radVec1 and radVec2.

               var    = one/sqrt(max(eps,(r1*r1 + r2*r2)))
               dir(1) = wax
               dir(2) = var*(wrad*r1 - wtheta*r2)
               dir(3) = var*(wrad*r2 + wtheta*r1)

               ! Transform this direction to the global cartesian frame.

               BCData(boco)%flowXdirInlet(i,j) = dir(1)*axis(1)     &
                                               + dir(2)*radVec1(1) &
                                               + dir(3)*radVec2(1)

               BCData(boco)%flowYdirInlet(i,j) = dir(1)*axis(2)     &
                                               + dir(2)*radVec1(2) &
                                               + dir(3)*radVec2(2)

               BCData(boco)%flowZdirInlet(i,j) = dir(1)*axis(3)     &
                                               + dir(2)*radVec1(3) &
                                               + dir(3)*radVec2(3)
             enddo
           enddo

         else radialTest

           ! Cartesian direction specified. Either the angle or the
           ! direction should be present.

           ! X-direction.

           if( axPresent ) then

             ! Angle specified. Convert it to SI-units and determine
             ! the corresponding direction.

             call siAngle(angle(4), mult, trans)

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowXdirInlet(i,j) = &
                                  cos(mult*bcVarArray(i,j,4) + trans)
               enddo
             enddo

           else

             ! Direction specified. Simply copy it.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowXdirInlet(i,j) = bcVarArray(i,j,7)
               enddo
             enddo

           endif

           ! Y-direction.

           if( ayPresent ) then

             ! Angle specified. Convert it to SI-units and determine
             ! the corresponding direction.

             call siAngle(angle(5), mult, trans)

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowYdirInlet(i,j) = &
                                  cos(mult*bcVarArray(i,j,5) + trans)
               enddo
             enddo

           else

             ! Direction specified. Simply copy it.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowYdirInlet(i,j) = bcVarArray(i,j,8)
               enddo
             enddo

           endif

           ! Z-direction.

           if( azPresent ) then

             ! Angle specified. Convert it to SI-units and determine
             ! the corresponding direction.

             call siAngle(angle(6), mult, trans)

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowZdirInlet(i,j) = &
                                  cos(mult*bcVarArray(i,j,6) + trans)
               enddo
             enddo

           else

             ! Direction specified. Simply copy it.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%flowZdirInlet(i,j) = bcVarArray(i,j,9)
               enddo
             enddo

           endif

         endif radialTest

         ! Loop over the faces of the subface to compute some
         ! additional info.

         do j=jBeg,jEnd
           do i=iBeg,iEnd

             ! Compute the total enthalpy from the given
             ! total temperature.

             call computeHtot(BCData(boco)%ttInlet(i,j), &
                              BCData(boco)%htInlet(i,j))

             ! Determine the unit vector of the flow direction.

             dir(1) = BCData(boco)%flowXdirInlet(i,j)
             dir(2) = BCData(boco)%flowYdirInlet(i,j)
             dir(3) = BCData(boco)%flowZdirInlet(i,j)

             var = one/max(eps,sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2))

             BCData(boco)%flowXdirInlet(i,j) = var*dir(1)
             BCData(boco)%flowYdirInlet(i,j) = var*dir(2)
             BCData(boco)%flowZdirInlet(i,j) = var*dir(3)

           enddo
         enddo

         ! Check if the prescribed direction is an inflow. No halo's
         ! should be included here and therefore the nodal range
         ! (with an offset) must be used.

         nn = 0
         do j=(BCData(boco)%jnbeg+1), BCData(boco)%jnend
           do i=(BCData(boco)%inbeg+1), BCData(boco)%inend

             var = BCData(boco)%flowXdirInlet(i,j) &
                 * BCData(boco)%norm(i,j,1)          &
                 + BCData(boco)%flowYdirInlet(i,j) &
                 * BCData(boco)%norm(i,j,2)          &
                 + BCData(boco)%flowZdirInlet(i,j) &
                 * BCData(boco)%norm(i,j,3)

             if(var > zero) nn = nn + 1

           enddo
         enddo

         if(nn > 0) then
           write(errorMessage,200)                   &
                 trim(cgnsDoms(nbkGlobal)%zonename), &
                 trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 200       format("Zone ",a,", subsonic inlet boundary subface ",a, &
                  ": Flow direction points out of the domain for &
                  &some faces.")

           call terminate("totalSubsonicInlet", errorMessage)
         endif

         end subroutine totalSubsonicInlet

         !===============================================================

         subroutine massflowSubsonicInlet
!
!        ****************************************************************
!        *                                                              *
!        * MassflowSubsonicInlet converts the prescribed mass flow    *
!        * conditions (density and velocity) into a useable format.     *
!        *                                                              *
!        ****************************************************************
!
         use section
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: i, j, nn

         real(kind=realType) :: mult, trans
         real(kind=realType) :: ax, r1, r2, var, vax, vrad, vtheta

         real(kind=realType), dimension(3) :: xc, vloc
         real(kind=realType), dimension(3) :: multVel, transVel
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Set massflowPrescribed to .true. to indicate that this
         ! type of boundary condition is present. Set the subsonic inlet
         ! treatment accordingly.

         massflowPrescribed = .true.
         BCData(boco)%subsonicInletTreatment = massFlow

         ! Set the density. Take the conversion factor to SI-units
         ! into account.

         call siDensity(mass(12), length(12), mult, trans)

         do j=jBeg,jEnd
           do i=iBeg,iEnd
             BCData(boco)%rho(i,j) = mult*bcVarArray(i,j,12) + trans
           enddo
         enddo

         ! Check the situation we are having here for the velocity.

         testRadial: if( velrPresent ) then

           ! Radial velocity component prescribed. This must be converted
           ! to cartesian components.

           ! Determine the unit vectors, which define the cylindrical
           ! coordinate system aligned with the rotation axis.

           call unitVectorsCylSystem(boco)

           ! Determine the conversion factor to SI-units for the three
           ! components. Note that a test must be made whether the theta
           ! component is present.

           call siVelocity(length(13), time(13), multVel(1), transVel(1))
           call siVelocity(length(16), time(16), multVel(2), transVel(2))

           if( veltPresent ) &
            call siVelocity(length(17), time(17), multVel(3), transVel(3))

           ! Initialize vtheta to zero. This value will be used
           ! if no theta velocity component was specified.

           vtheta = zero

           ! Loop over the faces of the subface.

           do j=jBeg,jEnd
             do i=iBeg,iEnd

               ! Determine the coordinates of the face center relative to
               ! the rotation point of this section. Normally this is an
               ! average of i-1, i, j-1, j, but due to the usage of the
               ! pointer xf and the fact that x originally starts at 0,
               ! an offset of 1 is introduced and thus the average should
               ! be taken of i, i+1, j and j+1.

               xc(1) = fourth*(xf(i,j,  1) + xf(i+1,j,  1)  &
                     +         xf(i,j+1,1) + xf(i+1,j+1,1)) &
                     - sections(sectionId)%rotCenter(1)
               xc(2) = fourth*(xf(i,j,  2) + xf(i+1,j,  2)  &
                     +         xf(i,j+1,2) + xf(i+1,j+1,2)) &
                     - sections(sectionId)%rotCenter(2)
               xc(3) = fourth*(xf(i,j,  3) + xf(i+1,j,  3)  &
                     +         xf(i,j+1,3) + xf(i+1,j+1,3)) &
                     - sections(sectionId)%rotCenter(3)

               ! Determine the coordinates in the local cartesian frame,
               ! i.e. the frame determined by axis, radVec1 and radVec2.

               ax = xc(1)*axis(1)    + xc(2)*axis(2)    &
                  + xc(3)*axis(3)
               r1 = xc(1)*radVec1(1) + xc(2)*radVec1(2) &
                  + xc(3)*radVec1(3)
               r2 = xc(1)*radVec2(1) + xc(2)*radVec2(2) &
                  + xc(3)*radVec2(3)

               ! Determine the velocity components in the local
               ! cylindrical system. Take the conversion to SI units
               ! into account.

               vax  = multVel(1)*bcVarArray(i,j,13) + transVel(1)
               vrad = multVel(2)*bcVarArray(i,j,16) + transVel(2)
               if( veltPresent ) &
                 vtheta = multVel(3)*bcVarArray(i,j,17) + transVel(3)

               ! Determine the velocities in the local cartesian
               ! frame determined by axis, radVec1 and radVec2.

               var     = one/sqrt(max(eps,(r1*r1 + r2*r2)))
               vloc(1) = vax
               vloc(2) = var*(vrad*r1 - vtheta*r2)
               vloc(3) = var*(vrad*r2 + vtheta*r1)

               ! Transform vloc to the global cartesian frame and
               ! store the values.

               BCData(boco)%velx(i,j) = vloc(1)*axis(1)    &
                                      + vloc(2)*radVec1(1) &
                                      + vloc(3)*radVec2(1)

               BCData(boco)%vely(i,j) = vloc(1)*axis(2)    &
                                      + vloc(2)*radVec1(2) &
                                      + vloc(3)*radVec2(2)

               BCData(boco)%velz(i,j) = vloc(1)*axis(3)    &
                                      + vloc(2)*radVec1(3) &
                                      + vloc(3)*radVec2(3)
             enddo
           enddo

         else testRadial

           ! Cartesian components prescribed.

           ! Determine the conversion factor to SI-units for the three
           ! components.

           call siVelocity(length(13), time(13), multVel(1), transVel(1))
           call siVelocity(length(14), time(14), multVel(2), transVel(2))
           call siVelocity(length(15), time(15), multVel(3), transVel(3))

           ! Set the velocities.

           do j=jBeg,jEnd
             do i=iBeg,iEnd
               BCData(boco)%velx(i,j) = multVel(1)*bcVarArray(i,j,13) &
                                      + transVel(1)
               BCData(boco)%vely(i,j) = multVel(2)*bcVarArray(i,j,14) &
                                      + transVel(2)
               BCData(boco)%velz(i,j) = multVel(3)*bcVarArray(i,j,15) &
                                      + transVel(3)
             enddo
           enddo

         endif testRadial

         ! Check if the prescribed velocity is an inflow. No halo's
         ! should be included here and therefore the nodal range
         ! (with an offset) must be used.

         nn = 0
         do j=(BCData(boco)%jnbeg+1), BCData(boco)%jnend
           do i=(BCData(boco)%inbeg+1), BCData(boco)%inend

             var = BCData(boco)%velx(i,j)*BCData(boco)%norm(i,j,1) &
                 + BCData(boco)%vely(i,j)*BCData(boco)%norm(i,j,2) &
                 + BCData(boco)%velz(i,j)*BCData(boco)%norm(i,j,3)

             if(var > zero) nn = nn + 1

           enddo
         enddo

         if(nn > 0) then
           write(errorMessage,300)                   &
                 trim(cgnsDoms(nbkGlobal)%zonename), &
                 trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 300       format("Zone ",a,", subsonic inlet boundary subface ",a, &
                  ": Velocity points out of the domain for some faces.")

           call terminate("massflowSubsonicInlet", errorMessage)
         endif

         end subroutine massflowSubsonicInlet

       end subroutine BCDataSubsonicInflow
