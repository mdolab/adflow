!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataSupersonicInflow.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataSupersonicInflow(boco, allFlowPresent, &
                                         allTurbPresent)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataSupersonicInflow tries to extract the primitive state    *
!      * vector for the currently active boundary face, which is a      *
!      * supersonic inflow.                                             *
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
       integer(kind=intType), intent(in) :: boco
       logical,               intent(inout) :: allFlowPresent
       logical,               intent(inout) :: allTurbPresent
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, nn

       real(kind=realType) :: var

       character(len=maxStringLen) :: errorMessage

       logical :: rhoPresent,  pPresent,    velPresent
       logical :: velxPresent, velyPresent, velzPresent
       logical :: velrPresent, veltPresent
       logical :: allTurbSubface
!
!      Interfaces
!
       interface
         logical function setBCVarTurb(offset, boco, turbInlet)
         use paramTurb
         use BCDataMod
         implicit none

         integer(kind=intType), intent(in) :: offset, boco
         real(kind=realType), dimension(:,:,:), pointer :: turbInlet

         end function setBCVarTurb
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

       nbcVar = 7
       if(equations == RANSEquations) nbcVar = nbcVar + nwt

       allocate(bcVarArray(iBeg:iEnd,jBeg:jEnd,nbcVar), stat=ierr)
       if(ierr /= 0)                                &
         call terminate("BCDataSupersonicInflow", &
                        "Memory allocation failure for bcVarArray")

       bcVarNames(1) = cgnsDensity
       bcVarNames(2) = cgnsPressure
       bcVarNames(3) = cgnsVelx
       bcVarNames(4) = cgnsVely
       bcVarNames(5) = cgnsVelz
       bcVarNames(6) = cgnsVelr
       bcVarNames(7) = cgnsVeltheta

       call setBCVarNamesTurb(7_intType)

       ! Try to determine these variables.

       call extractFromDataSet(BCFaceID(boco))

       ! Store the logicals, which indicate success or failure
       ! a bit more readable.

       rhoPresent  = bcVarPresent(1)
       pPresent    = bcVarPresent(2)
       velxPresent = bcVarPresent(3)
       velyPresent = bcVarPresent(4)
       velzPresent = bcVarPresent(5)
       velrPresent = bcVarPresent(6)
       veltPresent = bcVarPresent(7)

       ! Check if a velocity vector is present.

       velPresent = .false.
       if(velxPresent .and. velrPresent) velPresent = .true.
       if(velxPresent .and. velyPresent .and. velzPresent) &
         velPresent = .true.

       ! Check if rho, p and the velocity vector are present.

       testPresent: if(rhoPresent .and. pPresent .and. velPresent) then

         ! All the variables needed are prescribed. Set them.

         call prescribedSupersonicInlet

       else testPresent

         ! Not all variables are present. Check what type of flow
         ! is to be solved.

         select case(flowType)

           case (internalFlow)

             ! Internal flow. Data at the inlet must be specified;
             ! no free stream data can be taken.

             write(errorMessage,100)               &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 100         format("Zone ",a,", boundary subface ",a, &
                    ": Not enough data specified for supersonic inlet")

             call terminate("BCDataSupersonicInflow", errorMessage)

           !=============================================================

           case (externalFlow)

             ! External flow. Free stream data is used. However the
             ! correct nondimensional data is not yet known. Therefore
             ! this subface is saved in freestreamSubfaces and the
             ! data will be set later on.

             call storeFreestreamSubface(boco)

             ! For the turbulence something similar must be done.
             ! Set bcVarPresent(8) to .false. to indicate that not all
             ! turbulent variables are present and call setBCVarTurb
             ! to do the job.

             bcVarPresent(8) = .false.

             allTurbSubface = setBCVarTurb(7_intType, boco, &
                                           BCData(boco)%turbInlet)

             ! Set the velocity to the free stream direction such that
             ! the inflow test below can be applied. Furthermore give
             ! the density and pressure a value, such that the
             ! interpolation to the coarse grid values can be done.
             ! The actual values do not matter, because they will
             ! be overwritten later on in setSupersonicInletFreeStream.

             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 BCData(boco)%rho(i,j)  = zero
                 BCData(boco)%velx(i,j) = velDirFreestream(1)
                 BCData(boco)%vely(i,j) = velDirFreestream(2)
                 BCData(boco)%velz(i,j) = velDirFreestream(3)
                 BCData(boco)%ps(i,j)   = zero
               enddo
             enddo

             ! Set allFlowPresent to .false.

             allFlowPresent = .false.

         end select

       endif testPresent

       ! Release the memory of the bcVarArray.

       deallocate(bcVarArray, stat=ierr)
       if(ierr /= 0)                              &
         call terminate("BCDataSupersonicInflow", &
                        "Deallocation failure for bcVarArray")

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
         write(errorMessage,102)                   &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 102     format("Zone ",a,", supersonic inlet boundary subface ",a, &
                ": Velocity points out of the domain for some faces.")

         call terminate("BCDataSupersonicInflow", errorMessage)
       endif

       !=================================================================

       contains

         !===============================================================

         subroutine prescribedSupersonicInlet
!
!        ****************************************************************
!        *                                                              *
!        * prescribedSupersonicInlet sets the variables for this        *
!        * supersonic inlet to prescribed values.                       *
!        *                                                              *
!        ****************************************************************
!
         use section
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: i, j

         real(kind=realType) :: mult, trans
         real(kind=realType) :: ax, r1, r2, var, vax, vrad, vtheta

         real(kind=realType), dimension(3) :: xc, vloc
         real(kind=realType), dimension(3) :: multVel, transVel

         logical :: allTurbSubface
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Set the density. Take the conversion factor to SI-units
         ! into account.

         call siDensity(mass(1), length(1), mult, trans)

         do j=jBeg,jEnd
           do i=iBeg,iEnd
             BCData(boco)%rho(i,j) = mult*bcVarArray(i,j,1) + trans
           enddo
         enddo

         ! Set the pressure. Take the conversion factor to SI-units
         ! into account.

         call siPressure(mass(1), length(2), time(2), mult, trans)

         do j=jBeg,jEnd
           do i=iBeg,iEnd
             BCData(boco)%ps(i,j) = mult*bcVarArray(i,j,2) + trans
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

           call siVelocity(length(3), time(3), multVel(1), transVel(1))
           call siVelocity(length(6), time(6), multVel(2), transVel(2))

           if( veltPresent ) &
            call siVelocity(length(7), time(7), multVel(3), transVel(3))

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
                     - sections(sectionID)%rotCenter(1)
               xc(2) = fourth*(xf(i,j,  2) + xf(i+1,j,  2)  &
                     +         xf(i,j+1,2) + xf(i+1,j+1,2)) &
                     - sections(sectionID)%rotCenter(2)
               xc(3) = fourth*(xf(i,j,  3) + xf(i+1,j,  3)  &
                     +         xf(i,j+1,3) + xf(i+1,j+1,3)) &
                     - sections(sectionID)%rotCenter(3)

               ! Determine the coordinates in the local cartesian frame,
               ! i.e. the frame determined by axis, radVec1 and radVec2.

               ax = xc(1)*axis(1)    + xc(2)*axis(2)    &
                  + xc(3)*axis(3)
               r1 = xc(1)*radVec1(1) + xc(2)*radVec1(2) &
                  + xc(3)*radVec1(3)
               r2 = xc(1)*radVec2(1) + xc(2)*radVec2(2) &
                  + xc(3)*radVec2(3)

               ! Determine the velocity components in the local
               ! cylindrical system. Take the conversion to si units
               ! into account.

               vax  = multVel(1)*bcVarArray(i,j,3) + transVel(1)
               vrad = multVel(2)*bcVarArray(i,j,6) + transVel(2)
               if( veltPresent ) &
                 vtheta = multVel(3)*bcVarArray(i,j,7) + transVel(3)

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

           call siVelocity(length(3), time(3), multVel(1), transVel(1))
           call siVelocity(length(4), time(4), multVel(2), transVel(2))
           call siVelocity(length(5), time(5), multVel(3), transVel(3))

           ! Set the velocities.

           do j=jBeg,jEnd
             do i=iBeg,iEnd
               BCData(boco)%velx(i,j) = multVel(1)*bcVarArray(i,j,3) &
                                      + transVel(1)
               BCData(boco)%vely(i,j) = multVel(2)*bcVarArray(i,j,4) &
                                      + transVel(2)
               BCData(boco)%velz(i,j) = multVel(3)*bcVarArray(i,j,5) &
                                      + transVel(3)
             enddo
           enddo

         endif testRadial

         ! Set the turbulence variables and check if all of them are
         ! prescribed. If not set allTurbPresent to .false.

         allTurbPresent = setBCVarTurb(7_intType, boco, &
                                       BCData(boco)%turbInlet)

         if(.not. allTurbSubface) allTurbPresent = .false.

         end subroutine prescribedSupersonicInlet

       end subroutine BCDataSupersonicInflow
