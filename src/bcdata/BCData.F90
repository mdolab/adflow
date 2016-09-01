module BCData
  !
  !       This local module contains the variables and subroutine to     
  !       handle the prescribed boundary data.                           
  !
  use cgnsGrid
  implicit none
  save

  ! nbcVarMax: Parameter, which defines the maximum number of
  !              prescribed variables for a boundary.

  integer, parameter :: nbcVarMax = 21

  ! mass(nbcVarMax):   Unit of mass for the prescribed data.
  ! length(nbcVarMax): Unit of length for the prescribed data.
  ! time(nbcVarMax):   Unit of time for the prescribed data.
  ! temp(nbcVarMax):   Unit of temperature for the prescribed data.
  ! angle(nbcVarMax):  Unit of angle for the prescribed data.

  integer, dimension(nbcVarMax) :: mass, length, time, temp, angle

  ! nDataSet: Number of data sets present for the active face.
  ! cgnsBoco: The corresponding boundary face in the cgns block.
  ! nbcVar:   Theoretically possible number of variables
  !           prescribed for this face.
  ! iBeg:     Starting index i-direction.
  ! jBeg:     Starting index j-direction.
  ! iEnd:     Ending index i-direction.
  ! jEnd:     Ending index j-direction.

  integer(kind=intType) :: nDataSet, cgnsBoco, nbcVar
  integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd

  ! nFreestreamSubfaces:     Number of supersonic inflow subfaces
  !                          for which the flow field variables must
  !                          be set to the free stream values.
  ! freestreamSubfaces(:,3): the corresponding block, subface and
  !                          spectral solution.

  integer(kind=intType) :: nFreestreamSubfaces
  integer(kind=intType), dimension(:,:), allocatable :: &
       freestreamSubfaces

  ! nTurbFreestreamSubfaces:     Number of subfaces for which the
  !                              turbulence variables must be set
  !                              to the free stream values.
  ! turbFreestreamSubfaces(:,3): the corresponding block, subface
  !                              and spectral solution.

  integer(kind=intType) :: nTurbFreestreamSubfaces
  integer(kind=intType), dimension(:,:), allocatable :: &
       turbFreestreamSubfaces

  ! bcVarArray(:,:,nbcVar): array to store the interpolated
  !                         values for the active face.

  real(kind=realType), dimension(:,:,:), allocatable :: bcVarArray

  ! xf(:,:,:): pointer to the coordinates of the block face to
  !            which the active subface belongs.

  real(kind=realType), dimension(:,:,:), pointer :: xf

  ! axis(3):    Axial unit vector in the local cylindrical
  !             coordinate system.
  ! radVec1(3): First radial unit vector in the local cylindrical
  !             coordinate system.
  ! radVec2(3): Second radial unit vector.

  real(kind=realType), dimension(3) :: axis, radVec1, radVec2

  ! bcVarNames(nbcVarMax): The cgns names of the possible
  !                        prescribed variables.

  character(len=maxCGNSNameLen), dimension(nbcVarMax) :: bcVarNames

  ! axAssumed:               Whether or not the x-axis is assumed
  !                          to be the axial direction.
  ! massflowPrescribed:      Whether or not subsonic inflow boundaries
  !                          are present with prescribed massflow.
  ! bcVarPresent(nbcVarMax): Whether or not the possible
  !                          variables are actually prescribed.

  logical :: axAssumed, massflowPrescribed
  logical, dimension(nbcVarMax) :: bcVarPresent

  ! dataSet: Pointer for the data sets of the corresponding cgns
  !          boundary face.

  type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet

contains

  subroutine BCDataIsothermalWall(boco)
    !
    !       BCDataIsothermalWall tries to extract the wall temperature     
    !       for the currently active boundary face, which is an isothermal 
    !       viscous wall.                                                  
    !
    use constants
    use cgnsNames
    use blockPointers, only : BCFaceID, BCData, nBKGlobal
    use utils, only : terminate, siTemperature
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType) :: boco
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j

    real(kind=realType) :: mult, trans

    character(len=maxStringLen) :: errorMessage

    ! Allocate the memory for the buffer bcVarArray, which is used
    ! for the interpolation and set the cgns names.

    nbcVar = 1
    allocate(bcVarArray(iBeg:iEnd,jBeg:jEnd,nbcVar), stat=ierr)
    if(ierr /= 0)                            &
         call terminate("BCDataIsothermalWall", &
         "Memory allocation failure for bcVarArray")

    bcVarNames(1) = cgnsTemp

    ! Try to determine the temperature from the data set.

    call extractFromDataSet(BCFaceID(boco))

    ! Write an error message and terminate if it was not
    ! possible to determine the temperature.

    if(.not. bcVarPresent(1)) then

       write(errorMessage,100)                    &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
100    format("Zone ",a,", boundary subface ",a, &
            ": Wall temperature not specified for isothermal wall")

       call terminate("BCDataIsothermalWall", errorMessage)

    endif

    ! Convert to si-units and store the temperature in TNS_Wall.

    call siTemperature(temp(1), mult, trans)

    do j=jBeg,jEnd
       do i=iBeg,iEnd
          BCData(boco)%TNS_Wall(i,j) = mult*bcVarArray(i,j,1) + trans
       enddo
    enddo

    ! Release the memory of the bcVarArray.

    deallocate(bcVarArray, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("BCDataIsothermalWall", &
         "Deallocation failure for bcVarArray")

  end subroutine BCDataIsothermalWall

  subroutine BCDataSubsonicInflow(boco, allTurbPresent)
    !
    !       BCDataSubsonicInflow tries to extract the prescribed data      
    !       for the currently active boundary face, which is a subsonic    
    !       inflow. Either total conditions and velocity direction or the  
    !       velocity and density can be prescribed. In the latter case the 
    !       mass flow is prescribed, which is okay as long as the flow is  
    !       not choked.                                                    
    !
    use constants 
    use cgnsNames
    use blockPointers, only : nbkGlobal, sectionID, BCFaceID, BCData
    use flowVarRefState, only : nwt
    use inputPhysics, only : equations
    use utils, only : siDensity, siVelocity, siPressure, siAngle, &
         siTemperature, terminate
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
100    format("Zone ",a,", boundary subface ",a, &
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
      !         TotalSubsonicInlet converts the prescribed total           
      !         conditions and velocity direction into a useable format.     
      !
      use inputPhysics, only : RGasDim
      use section, only : sections
      implicit none
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j, nn

      real(kind=realType) :: rhot, mult, trans
      real(kind=realType) :: ax, r1, r2, var, wax, wrad, wtheta

      real(kind=realType), dimension(3) :: xc, dir

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
200      format("Zone ",a,", subsonic inlet boundary subface ",a, &
              ": Flow direction points out of the domain for &
              &some faces.")

         call terminate("totalSubsonicInlet", errorMessage)
      endif

    end subroutine totalSubsonicInlet

    subroutine massflowSubsonicInlet
      !
      !         MassflowSubsonicInlet converts the prescribed mass flow    
      !         conditions (density and velocity) into a useable format.     
      !
      use section, only: sections
      implicit none
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j, nn

      real(kind=realType) :: mult, trans
      real(kind=realType) :: ax, r1, r2, var, vax, vrad, vtheta

      real(kind=realType), dimension(3) :: xc, vloc
      real(kind=realType), dimension(3) :: multVel, transVel

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
300      format("Zone ",a,", subsonic inlet boundary subface ",a, &
              ": Velocity points out of the domain for some faces.")

         call terminate("massflowSubsonicInlet", errorMessage)
      endif

    end subroutine massflowSubsonicInlet

  end subroutine BCDataSubsonicInflow

  subroutine BCDataSubsonicOutflow(boco)
    !
    !       BCDataSubsonicOutflow tries to extract the static pressure     
    !       for the currently active boundary face, which is a subsonic    
    !       outflow boundary.                                              
    !
    use constants
    use cgnsNames
    use blockPointers, only : BCData, nbkGlobal, BCFaceID
    use utils, only : terminate, siPressure
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType) :: boco
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j

    real(kind=realType) :: mult, trans

    character(len=maxStringLen) :: errorMessage

    ! Allocate the memory for the buffer bcVarArray, which is used
    ! for the interpolation and set the cgns names.

    nbcVar = 1
    allocate(bcVarArray(iBeg:iEnd,jBeg:jEnd,nbcVar), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("BCDataSubsonicOutflow", &
         "Memory allocation failure for bcVarArray")

    bcVarNames(1) = cgnsPressure

    ! Try to determine the static pressure from the data set.

    call extractFromDataSet(BCFaceID(boco))

    ! Write an error message and terminate if it was not
    ! possible to determine the static pressure.

    if(.not. bcVarPresent(1)) then

       write(errorMessage,100)            &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
100    format("Zone ",a,", boundary subface ",a, &
            ": Static pressure not specified for subsonic outlet")

       call terminate("BCDataSubsonicOutflow", errorMessage)

    endif

    ! Convert to SI-units and store the pressure in ps.

    call siPressure(mass(1), length(1), time(1), mult, trans)

    do j=jBeg,jEnd
       do i=iBeg,iEnd
          BCData(boco)%ps(i,j) = mult*bcVarArray(i,j,1) + trans
       enddo
    enddo

    ! Release the memory of the bcVarArray.

    deallocate(bcVarArray, stat=ierr)
    if(ierr /= 0)                               &
         call terminate("BCDataSubsonicOutflow", &
         "Deallocation failure for bcVarArray")

  end subroutine BCDataSubsonicOutflow

  subroutine BCDataSupersonicInflow(boco, allFlowPresent, &
       allTurbPresent)
    !
    !       BCDataSupersonicInflow tries to extract the primitive state    
    !       vector for the currently active boundary face, which is a      
    !       supersonic inflow.                                             
    !
    use constants
    use cgnsNames
    use blockPointers, only : BCData, nbkGlobal, BCFaceID, sectionID
    use flowVarRefState, only : nwt
    use inputPhysics, onlY : equations, flowType, velDirFreeStream
    use utils, only : siDensity, siPressure, siVelocity, siTemperature, terminate
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
100       format("Zone ",a,", boundary subface ",a, &
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
102    format("Zone ",a,", supersonic inlet boundary subface ",a, &
            ": Velocity points out of the domain for some faces.")

       call terminate("BCDataSupersonicInflow", errorMessage)
    endif


  contains

    subroutine prescribedSupersonicInlet
      !
      !         prescribedSupersonicInlet sets the variables for this        
      !         supersonic inlet to prescribed values.                       
      !
      use section, only : sections
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

  subroutine allocMemBCData
    !
    !       allocMemBCData allocates the memory for the prescribed         
    !       boundary data for all multigrid levels and all spectral        
    !       solutions for all blocks.                                      
    !
    use constants
    use blockPointers, only : BCData, flowDoms, nBocos, nDom, BCType
    use flowVarRefState, only : nt1, nt2
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : nALESteps
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: mm, nn, sps, level, nLevels
    integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd
    integer(kind=intType) :: inodeBeg, jnodeBeg, inodeEnd, jnodeEnd

    ! Determine the number of multigrid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of multigrid level, spectral solutions
    ! and local blocks.

    levelLoop: do level=1,nLevels
       spectralLoop: do sps=1,nTimeIntervalsSpectral
          domainsLoop: do nn=1,nDom

             ! Have the pointers in blockPointers point to the
             ! current block to make everything more readable.

             call setPointers(nn, level, sps)

             ! Loop over the number of boundary subfaces for this block.

             bocoLoop: do mm=1,nBocos

                ! Store the cell range of the boundary subface
                ! a bit easier.

                iBeg = BCData(mm)%icbeg; iEnd = BCData(mm)%icend
                jBeg = BCData(mm)%jcbeg; jEnd = BCData(mm)%jcend

                inodeBeg = BCData(mm)%inbeg; inodeEnd = BCData(mm)%inend
                jnodeBeg = BCData(mm)%jnbeg; jnodeEnd = BCData(mm)%jnend

                ! Determine the boundary condition we are having here
                ! and allocate the memory accordingly.

                select case (BCType(mm))

                case (NSWallAdiabatic)
                   allocate(BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%uSlipALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%T(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &                         
                        BCData(mm)%Tp(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%Tv(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), & 
                        BCData(mm)%fIndex(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                        BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                                ! Note: iBlank/delta are cell based, but uses the
                                ! node number to guarantee a halo exists. 
                        BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%area(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &an adiabatic wall")

                   !=======================================================

                case (NSWallIsothermal)

                   allocate(BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3),  &
                        BCData(mm)%uSlipALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%T(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &                         
                        BCData(mm)%Tp(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%Tv(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), & 
                        BCData(mm)%fIndex(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                        BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                                ! Note: iBlank/delta are cell based, but uses the
                                ! node number to guarantee a halo exists. 
                        BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%area(iBeg:iEnd,jBeg:jEnd), &
                        
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &an isothermal wall")

                   !=======================================================

                case (EulerWall)

                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%T(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &                         
                        BCData(mm)%Tp(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                        BCData(mm)%Tv(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), & 
                        BCData(mm)%fIndex(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                        BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                        BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                                ! Note: iBlank/delta are cell based, but uses the
                                ! node number to guarantee a halo exists. 
                        BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                        BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                        BCData(mm)%area(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &an Euler wall")

                   !=======================================================

                case (farField)

                   ! Just allocate the memory for the normal mesh
                   ! velocity.

                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a farfield")
                   !=======================================================

                case (symm, symmPolar)

                   ! Allocate for symm as well. This is not necessary
                   ! but we need it for the reverse AD.

                   ! Modified by HDN
                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a symm")

                   !=======================================================

                case (SupersonicInflow, DomainInterfaceAll)

                   ! Supersonic inflow or a domain interface with
                   ! all the data prescribed. Allocate the memory for
                   ! the entire state vector to be prescribed.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),   &
                        BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd),  &
                        BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd),  &
                        BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd),  &
                        BCData(mm)%ps(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a supersonic inflow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                      allocate(&
                           BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                           stat=ierr)
                      if(ierr /= 0)                      &
                           call terminate("allocMemBCData", &
                           "Memory allocation failure for &
                           &turbInlet for a supersonic &
                           &inflow")
                   endif

                   !=======================================================

                   ! Added by HDN
                case (SupersonicOutflow)
                   ! No state is needed for this boco


                   !=======================================================

                case (SubsonicInflow)

                   ! Subsonic inflow. Allocate the memory for the
                   ! variables needed. Note the there are two ways to
                   ! specify boundary conditions for a subsonic inflow.

                   allocate(BCData(mm)%flowXdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%flowYdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%flowZdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%ptInlet(iBeg:iEnd,jBeg:jEnd),       &
                        BCData(mm)%ttInlet(iBeg:iEnd,jBeg:jEnd),       &
                        BCData(mm)%htInlet(iBeg:iEnd,jBeg:jEnd),       &
                        BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),           &
                        BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd),          &
                        BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd),          &
                        BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a subsonic inflow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                      allocate(&
                           BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                           stat=ierr)
                      if(ierr /= 0)                      &
                           call terminate("allocMemBCData", &
                           "Memory allocation failure for &
                           &turbInlet for a subsonic inflow")
                   endif

                   !=======================================================

                case (SubsonicOutflow, MassBleedOutflow, &
                     DomainInterfaceP)

                   ! Subsonic outflow, outflow mass bleed or domain
                   ! interface with prescribed pressure. Allocate the
                   ! memory for the static pressure.

                   allocate(BCData(mm)%ps(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a subsonic outflow, outflow mass &
                        &bleed or domain interface with &
                        &prescribed pressure.")

                   ! Initialize the pressure to avoid problems for
                   ! the bleed flows.

                   BCData(mm)%ps = zero

                   !=======================================================

                case (DomainInterfaceRhoUVW)

                   ! Domain interface with prescribed density and 
                   ! velocities, i.e. mass flow is prescribed. Allocate
                   ! the memory for the variables needed.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),  &
                        BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a domain interface with a &
                        &prescribed mass flow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                      allocate(&
                           BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                           stat=ierr)
                      if(ierr /= 0)                      &
                           call terminate("allocMemBCData", &
                           "Memory allocation failure for &
                           &turbInlet for a domain interface &
                           &with a prescribed mass flow")
                   endif

                   !=======================================================

                case (DomainInterfaceTotal)

                   ! Domain interface with prescribed total conditions.
                   ! Allocate the memory for the variables needed.

                   allocate(BCData(mm)%flowXdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%flowYdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%flowZdirInlet(iBeg:iEnd,jBeg:jEnd), &
                        BCData(mm)%ptInlet(iBeg:iEnd,jBeg:jEnd),       &
                        BCData(mm)%ttInlet(iBeg:iEnd,jBeg:jEnd),       &
                        BCData(mm)%htInlet(iBeg:iEnd,jBeg:jEnd),       &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a domain interface with total &
                        &conditions")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                      allocate(&
                           BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                           stat=ierr)
                      if(ierr /= 0)                      &
                           call terminate("allocMemBCData", &
                           "Memory allocation failure for &
                           &turbInlet for a domain interface &
                           &with a prescribed mass flow")
                   endif

                   !=======================================================

                case (domainInterfaceRho)

                   ! Domain interface with prescribed density. 
                   ! Allocate the memory for the density.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd), &
                        stat=ierr)
                   if(ierr /= 0)                      &
                        call terminate("allocMemBCData", &
                        "Memory allocation failure for &
                        &a domain interface")

                end select

             enddo bocoLoop

          enddo domainsLoop
       enddo spectralLoop
    enddo levelLoop

  end subroutine allocMemBCData

  subroutine cart1D_InterpolSubface(iBeg, iEnd, jBeg, jEnd, &
       nbcVar, cgnsBoco,       &
       blockFaceId, coorId,    &
       indCoor, ind,           &
       bcVarPresent,           &
       bcVarArray)
    !
    !       cart1D_InterpolSubface interpolates the prescribed variables   
    !       in the data set of the given cgns subface onto the subface     
    !       indicated by iBeg, iEnd, jBeg, jEnd and blockFaceId.           
    !       This routine performs a 1d cartesian interpolation for the     
    !       coordinate coorId and it is assumed that there is no           
    !       variation in the other directions.                             
    !       The variables in blockPointers are already set.                
    !
    use constants
    use blockPointers, only : x, il, jl, kl, nBKGlobal
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
    integer(kind=intType), intent(in) :: nbcVar, cgnsBoco
    integer(kind=intType), intent(in) :: blockFaceId, coorId

    integer(kind=intType), dimension(2,3), intent(in) :: indCoor
    integer(kind=intType), dimension(2,nbcVar), intent(in) :: ind

    real(kind=realType), dimension(iBeg:iEnd,jBeg:jEnd,nbcVar), &
         intent(out) :: bcVarArray

    logical, dimension(nbcVar), intent(in) :: bcVarPresent
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l, n1, n2, nn, ii
    integer(kind=intType) :: nDim, nPoints, start

    real(kind=realType) :: fact, ww1, ww2, xc

    real(kind=realType), dimension(:),   pointer :: xx
    real(kind=realType), dimension(:,:), pointer :: xf

    character(len=maxStringLen) :: errorMessage

    type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet


    ! Set the pointer for dataSet to make the code more readable.

    dataSet => cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%dataSet

    ! Check the number of dimensions of the specified data set.
    ! This should be 1, because a 1d interpolation in radial
    ! direction is performed.

    k    = indCoor(1,coorId)
    l    = indCoor(2,coorId)
    nDim = dataSet(k)%dirichletArrays(l)%nDimensions

    if(nDim > 1) then
       write(errorMessage,101) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
       call terminate("cart1D_InterpolSubface", errorMessage)
101    format("Zone",1X,A,", subface",1X,A,": Multidimensional &
            &radially varying data specified. Only 1d data possible")
    endif

    ! Set the pointer for the current coordinate and abbreviate the
    ! number of interpolation points a bit easier.

    xx     => dataSet(k)%dirichletArrays(l)%dataArr
    nPoints = dataSet(k)%dirichletArrays(l)%dataDim(1)

    ! Check if the data is specified for increasing values of
    ! the coordinate.

    do i=2,nPoints
       if(xx(i) < xx(i-1)) exit
    enddo

    if(i <= nPoints) then
       write(errorMessage,102) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
       call terminate("cart1D_InterpolSubface", errorMessage)
102    format("Zone",1X,A,", subface",1X,A,": Data should be &
            &specified for increasing coordinate values.")
    endif

    ! Set the pointer for the coordinate of the block face on which
    ! the boundary subface is located.

    select case (blockFaceId)
    case (iMin)
       xf => x(1, :,:,coorId)
    case (iMax)
       xf => x(il,:,:,coorId)
    case (jMin)
       xf => x(:,1, :,coorId)
    case (jMax)
       xf => x(:,jl,:,coorId)
    case (kMin)
       xf => x(:,:,1, coorId)
    case (kMax)
       xf => x(:,:,kl,coorId)
    end select

    ! Compute the factor needed to compute the coordinates in the
    ! original units. The fourth comes from the averaging of the 4
    ! nodal coordinates.

    fact = fourth/cgnsDoms(nbkGlobal)%LRef

    ! Loop over the range of the subface.

    jloop: do j=jBeg,jEnd
       iloop: do i=iBeg,iEnd

          ! Determine the coordinate of the face center. Normally this
          ! is an average of i-1, i, j-1, j, but due to the usage of
          ! the pointer xf and the fact that x originally starts at 0,
          ! an offset of 1 is introduced and thus the average should
          ! be taken of i, i+1, j and j+1.

          xc = fact*(xf(i,j) + xf(i+1,j) + xf(i,j+1) + xf(i+1,j+1))

          ! Determine the interpolation interval and set the
          ! interpolation weights. Take care of the exceptions.

          checkInterpol: if(xc <= xx(1)) then

             ! Coordinate is less than the minimum value specified.
             ! Use constant extrapolation.

             n1  = 1;   n2  = 1
             ww1 = one; ww2 = zero

          else if(xc >= xx(nPoints)) then checkInterpol

             ! Coordinate is larger than the maximum value specified.
             ! Use constant extrapolation.

             n1  = nPoints; n2  = nPoints
             ww1 = one;     ww2 = zero

          else checkInterpol

             ! Coordinate is in the range. Determine the correct
             ! interval using a binary search algorithm.

             ii    = nPoints - 1
             start = 1
             interval: do

                ! Next guess for the interval and determine the new
                ! situation.

                nn = start + ii/2
                if(xc > xx(nn+1)) then

                   ! Coordinate is larger than the upper boundary of the
                   ! current interval. Update the lower boundary.

                   start = nn + 1; ii = ii - 1

                else if(xc >= xx(nn)) then

                   ! This is the correct range. Exit the loop.

                   exit

                endif

                ! Modify ii for the next branch to search.

                ii = ii/2

             enddo interval

             ! xc is in the interval nn, nn+1. Store this and
             ! determine the interpolation weight.

             n1 = nn
             n2 = nn + 1
             ww1 = (xx(nn+1) - xc)/(xx(nn+1) - xx(nn))
             ww2 = one - ww1

          endif checkInterpol

          ! Interpolate the values the values present for this face.

          do nn=1,nbcVar
             if( bcVarPresent(nn) ) then

                ! Easier storage of the indices in the data set.

                k = ind(1,nn)
                l = ind(2,nn)

                ! Interpolate this variable.

                bcVarArray(i,j,nn) =                                 &
                     ww1*dataSet(k)%dirichletArrays(l)%dataArr(n1) &
                     + ww2*dataSet(k)%dirichletArrays(l)%dataArr(n2)
             endif
          enddo

       enddo iloop
    enddo jloop

  end subroutine cart1D_InterpolSubface

  subroutine computeHtot(tt, ht)
    !
    !       computeHtot computes the total enthalpy from the given total   
    !       temperature. The total enthalpy is the integral of cp, which   
    !       is a very simple expression for constant cp. For a variable cp 
    !       it is a bit more work.                                         
    !
    use constants
    use cpCurveFits
    use communication, only : myid
    use inputPhysics, only : cpModel, gammaConstant,rGasDim
    use flowVarRefState, only : PinfDim
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), intent(in)  :: tt
    real(kind=realType), intent(out) :: ht
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, nn, mm, start

    real(kind=realType) :: t2

    ! Determine the cp model used in the computation.

    select case (cpModel)

    case (cpConstant)

       ! Constant cp. The total enthalpy is simply cp*tt.

       ht = gammaConstant*RGasDim*tt/(gammaConstant - one)

       !        ================================================================

    case (cpTempCurveFits)

       ! Cp as function of the temperature is given via curve fits.
       ! The actual integral must be computed.

       ! Determine the case we are having here.

       if(tt < cpTrange(0)) then

          ! Temperature is less than the smallest value in the
          ! curve fits. Print a warning and use extrapolation using
          ! constant cp.

          if(myId == 0) then
             print "(a)", "#"
             print "(a)", "#                    Warning"
             print 100, tt, cpTrange(0)
             print "(a)", "# Extrapolation with constant cp is used."
             print "(a)", "#"
100          format("# Prescribed total temperature ",e12.5,          &
                  " is less than smallest curve fit value, ",e12.5, &
                  ".")
          endif

          ht = RGasDim*(cpEint(0) + tt + cv0*(tt - cpTrange(0)))

       else if(tt > cpTrange(cpNparts)) then

          ! Temperature is larger than the largest value in the
          ! curve fits. Print a warning and use extrapolation using
          ! constant cp.

          if(myId == 0) then
             print "(a)", "#"
             print "(a)", "#                    Warning"
             print 101, tt, cpTrange(cpNparts)
             print "(a)", "# Extrapolation with constant cp is used."
             print "(a)", "#"
101          format("# Prescribed total temperature ",e12.5,     &
                  " is larger than largest curve fit value, ", &
                  e12.5, ".")
          endif

          ht = RGasDim*(cpEint(cpNparts) + tt &
               +           cvn*(tt - cpTrange(cpNparts)))

       else

          ! Temperature is in the curve fit range.
          ! First find the correct range for this temperature.

          ii    = cpNparts
          start = 1
          interval: do

             ! Next guess for the interval.

             nn = start + ii/2

             ! Determine the situation we are having here.

             if(tt > cpTrange(nn)) then

                ! Temperature is larger than the upper boundary of
                ! the current interval. Update the lower boundary.

                start = nn + 1
                ii    = ii - 1

             else if(tt >= cpTrange(nn-1)) then

                ! This is the correct range. Exit the do-loop.

                exit

             endif

             ! Modify ii for the next branch to search.

             ii = ii/2

          enddo interval

          ! nn contains the correct curve fit interval.
          ! Integrate cp to get ht.

          ht = cpTempFit(nn)%eint0
          do ii=1,cpTempFit(nn)%nterm
             if(cpTempFit(nn)%exponents(ii) == -1_intType) then
                ht = ht + cpTempFit(nn)%constants(ii)*log(tt)
             else
                mm = cpTempFit(nn)%exponents(ii) + 1
                t2 = tt**mm
                ht = ht + cpTempFit(nn)%constants(ii)*t2/mm
             endif
          enddo

          ! Multiply ht by RGasDim to obtain the correct
          ! dimensional value.

          ht = RGasDim*ht

       endif

    end select

  end subroutine computeHtot

  subroutine extractFromDataSet(blockFaceID)
    !
    !       extractFromDataSet tries to extract and interpolate the        
    !       variables in bcVarNames from the cgns data set.                
    !       If successful the corresponding entry of bcVarPresent is       
    !       set to .true., otherwise it is set to .false.                  
    !
    use constants
    use cgnsNames
    use blockPointers, onlY : nbkGlobal
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: blockFaceID
    !
    !      Local variables.
    !
    integer(kind=intType) :: k, l, m, n
    integer(kind=intType) :: nInter, nDim, nVarPresent, nCoor

    integer(kind=intType), dimension(3) :: dataDim, coor
    integer(kind=intType), dimension(2,3) :: indCoor
    integer(kind=intType), dimension(2,nbcVar) :: ind

    character(len=maxStringLen) :: errorMessage

    logical :: xPresent, yPresent, zPresent, rPresent
    logical :: sameInterpol, firstVar

    ! Determine whether the variables are specified and if so,
    ! where they are located in the data set. As the number of
    ! variables specified is usually not so big, a linear search
    ! algorithm is perfectly okay. At the moment only the Dirichlet
    ! arrays are checked.

    nVarPresent = 0

    do m=1,nbcVar
       bcVarPresent(m) = .false.

       dataSetLoop: do k=1,nDataSet
          do l=1,dataSet(k)%nDirichletArrays
             if(dataSet(k)%dirichletArrays(l)%arrayName == &
                  bcVarNames(m)) then

                ! Variable is present. Store the indices, update
                ! nVarPresent and set bcVarPresent(m) to .True.

                ind(1,m) = k; ind(2,m) = l

                nVarPresent      = nVarPresent + 1
                bcVarPresent(m) = .true.

                ! Set the units for this variable.

                mass(m)   = dataSet(k)%dirichletArrays(l)%mass
                length(m) = dataSet(k)%dirichletArrays(l)%len
                time(m)   = dataSet(k)%dirichletArrays(l)%time
                temp(m)   = dataSet(k)%dirichletArrays(l)%temp
                angle(m)  = dataSet(k)%dirichletArrays(l)%angle

                ! Exit the search loop, as the variable was found.

                exit dataSetLoop

             endif
          enddo
       enddo dataSetLoop
    enddo

    ! Return if none of the variables are present.

    if(nVarPresent == 0) return

    ! Find out whether the given data points are equal for every
    ! variable or that every variable must be interpolated
    ! differently.

    sameInterpol = .true.
    firstVar     = .true.

    do m=1,nbcVar
       if( bcVarPresent(m) ) then
          k = ind(1,m)
          l = ind(2,m)

          if( firstVar ) then
             nDim = dataSet(k)%dirichletArrays(l)%nDimensions
             firstVar = .false.

             do n=1,nDim
                dataDim(n) = dataSet(k)%dirichletArrays(l)%dataDim(n)
             enddo
          else
             if(nDim == dataSet(k)%dirichletArrays(l)%nDimensions) then
                do n=1,nDim
                   if(dataSet(k)%dirichletArrays(l)%dataDim(n) /= &
                        dataDim(n)) sameInterpol = .false.
                enddo
             else
                sameInterpol = .false.
             endif
          endif

       endif
    enddo

    ! Determine the situation we are dealing with here.

    testSameInterpol: if( sameInterpol ) then

       ! The interpolation is the same for all variables.
       ! First determine the number of interpolation points.

       nInter = dataDim(1)
       do m=2,nDim
          nInter = nInter*dataDim(m)
       enddo

       ! If nInter == 1 then the prescribed data is constant
       ! everywhere and the variables can be determined easily.

       testConstant1: if(nInter == 1) then

          ! Data is constant for this subface. Set the data.

          do m=1,nbcVar
             if( bcVarPresent(m) ) then
                k = ind(1,m)
                l = ind(2,m)

                bcVarArray(:,:,m) = &
                     dataSet(k)%dirichletArrays(l)%dataArr(1)
             endif
          enddo

       else testConstant1

          ! Data varies over the interface and must be interpolated.
          ! Determine the indices of the coordinates in the dataset.

          rPresent = .false.
          xPresent = .false.
          yPresent = .false.
          zPresent = .false.

          do k=1,nDataSet
             do l=1,dataSet(k)%nDirichletArrays

                if(dataSet(k)%dirichletArrays(l)%arrayName == &
                     cgnsCoorr) then
                   indCoor(1,1) = k; indCoor(2,1) = l
                   rPresent = .true.
                   exit
                endif

                if(dataSet(k)%dirichletArrays(l)%arrayName == &
                     cgnsCoorx) then
                   indCoor(1,1) = k; indCoor(2,1) = l
                   xPresent = .true.
                endif

                if(dataSet(k)%dirichletArrays(l)%arrayName == &
                     cgnsCoory) then
                   indCoor(1,2) = k; indCoor(2,2) = l
                   yPresent = .true.
                endif

                if(dataSet(k)%dirichletArrays(l)%arrayName == &
                     cgnsCoorz) then
                   indCoor(1,3) = k; indCoor(2,3) = l
                   zPresent = .true.
                endif

             enddo
          enddo

          ! Check if a radial coordinate is present.

          if( rPresent ) then

             ! Radial coordinate is present. Use radial interpolation
             ! for the given variable.

             call radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                  nbcVar, cgnsBoco,       &
                  blockFaceID, indCoor,   &
                  ind, bcVarPresent,      &
                  bcVarArray, axAssumed)

          else if(xPresent .or. yPresent .or. zPresent) then

             ! Cartesian interpolation will be performed. Determine
             ! which coordinates are present.

             nCoor = 0
             if( xPresent ) then
                nCoor = nCoor + 1; coor(nCoor) = 1
             endif
             if( yPresent ) then
                nCoor = nCoor + 1; coor(nCoor) = 2
             endif
             if( zPresent ) then
                nCoor = nCoor + 1; coor(nCoor) = 3
             endif

             ! The number of dimensions cannot be larger than the
             ! number of coordinates. Check this.

             if(nDim > nCoor) then
                write(errorMessage,100) &
                     trim(cgnsDoms(nbkGlobal)%zonename), &
                     trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
                call terminate("extractFromDataSet", errorMessage)
             endif

             ! Check what kind of interpolation must be used.

             select case (nCoor)

             case (1_intType)

                ! 1D line interpolation.

                call cart1D_InterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                     nbcVar, cgnsBoco,       &
                     blockFaceID, coor(1),   &
                     indCoor, ind,           &
                     bcVarPresent,           &
                     bcVarArray)

                !=======================================================

             case (2_intType, 3_intType)

                call terminate("extractFromDataSet", &
                     "Multi-D Cartesian interpolation &
                     &not implemented yet")

             end select

          else

             ! Neither the radial nor the cartesian coordinates are
             ! present. So there is not enough information available
             ! for the interpolation. Print an error message and exit.

             write(errorMessage,101) &
                  trim(cgnsDoms(nbkGlobal)%zonename), &
                  trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
             call terminate("extractFromDataSet", errorMessage)

          endif

       endif testConstant1

    else testSameInterpol

       ! Different interpolation must be used for the different
       ! variables. Loop over the number of variables and test
       ! whether they are present.

       bcVarLoop: do m=1,nbcVar
          testBcVarPresent: if( bcVarPresent(m) ) then

             ! Store the indices of the corresponding data set and
             ! dirichlet array a bit easier. Abbreviate the number
             ! of dimensions a bit easier.

             k    = ind(1,m)
             l    = ind(2,m)
             nDim = dataSet(k)%dirichletArrays(l)%nDimensions

             ! Determine the number of interpolation points.

             nInter = dataSet(k)%dirichletArrays(l)%dataDim(1)
             do n=2,nDim
                nInter = nInter &
                     * dataSet(k)%dirichletArrays(l)%dataDim(n)
             enddo

             ! If nInter == 1 then the prescribed data is constant
             ! everywhere and the variable can be determined easily.

             testConstant2: if(nInter == 1) then

                ! Data is constant for this subface. Set it

                bcVarArray(:,:,m) = &
                     dataSet(k)%dirichletArrays(l)%dataArr(1)

             else testConstant2

                ! Data varies over the interface and must be
                ! interpolated. Determine the indices of the
                ! coordinates in the dirichlet arrays of the given
                ! data set. Note that the coordinates now have to be
                ! specified in the same dataset as the variable.

                rPresent = .false.
                xPresent = .false.
                yPresent = .false.
                zPresent = .false.

                do l=1,dataSet(k)%nDirichletArrays

                   if(dataSet(k)%dirichletArrays(l)%arrayName == &
                        cgnsCoorr) then
                      indCoor(1,1) = k; indCoor(2,1) = l
                      rPresent = .true.
                      exit
                   endif

                   if(dataSet(k)%dirichletArrays(l)%arrayName == &
                        cgnsCoorx) then
                      indCoor(1,1) = k; indCoor(2,1) = l
                      xPresent = .true.
                   endif

                   if(dataSet(k)%dirichletArrays(l)%arrayName == &
                        cgnsCoory) then
                      indCoor(1,2) = k; indCoor(2,2) = l
                      yPresent = .true.
                   endif

                   if(dataSet(k)%dirichletArrays(l)%arrayName == &
                        cgnsCoorz) then
                      indCoor(1,3) = k; indCoor(2,3) = l
                      zPresent = .true.
                   endif

                enddo

                ! Check if a radial coordinate is present.

                if( rPresent ) then

                   ! Radial coordinate is present. Use radial
                   ! interpolation for the given variable.

                   call radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                        1_intType, cgnsBoco,    &
                        blockFaceID, indCoor,   &
                        ind(1,m),               &
                        bcVarPresent(m),        &
                        bcVarArray(1,1,m),      &
                        axAssumed)

                else if(xPresent .or. yPresent .or. zPresent) then

                   ! Cartesian interpolation will be performed. Determine
                   ! which coordinates are present.

                   nCoor = 0
                   if( xPresent ) then
                      nCoor = nCoor + 1; coor(nCoor) = 1
                   endif
                   if( yPresent ) then
                      nCoor = nCoor + 1; coor(nCoor) = 2
                   endif
                   if( zPresent ) then
                      nCoor = ncoor + 1; coor(nCoor) = 3
                   endif

                   ! The number of dimensions cannot be larger than the
                   ! number of coordinates. Check this.

                   if(nDim > nCoor) then
                      write(errorMessage,200) &
                           trim(cgnsDoms(nbkGlobal)%zonename), &
                           trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName), &
                           trim(bcVarNames(m))
                      call terminate("extractFromDataSet", errorMessage)
                   endif

                   ! Check what kind of interpolation must be used.

                   select case (nCoor)

                   case (1_intType)

                      ! 1D line interpolation.

                      call cart1D_InterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                           1_intType, cgnsBoco,    &
                           blockFaceID, coor(1),   &
                           indCoor, ind,           &
                           bcVarPresent(m),        &
                           bcVarArray(1,1,m))

                      !===================================================

                   case (2_intType, 3_intType)

                      call terminate("extractFromDataSet", &
                           "Multi-D Cartesian interpolation &
                           &not implemented yet")

                   end select

                else

                   ! Neither the radial nor the cartesian coordinates
                   ! are present. So there is not enough information
                   ! available for the interpolation. Print an error
                   ! message and exit.

                   write(errorMessage,201) &
                        trim(cgnsDoms(nbkGlobal)%zonename), &
                        trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName), &
                        trim(bcVarNames(m))
                   call terminate("extractFromDataSet", &
                        errorMessage)
                endif

             endif testConstant2

          endif testBcVarPresent
       enddo bcVarLoop

    endif testSameInterpol

    ! Format statements.

100 format("Zone",1X,A,", subface",1X,A,": Number of dimensions &
         &is larger than number of coordinates.")
101 format("Zone",1X,A,", subface",1X,A,": No coordinates &
         &are present for the interpolation.")
200 format("Zone",1X,A,", subface",1X,A,", variable",1X,A, &
         ": Number of dimensions is larger than number of &
         &coordinates.")
201 format("Zone",1X,A,", subface",1X,A,": No coordinates &
         &are present for the interpolation of",1X,A,".")

  end subroutine extractFromDataSet

  subroutine initBCData
    !
    !       initBCData allocates and initializes the arrays BCData for     
    !       all boundary subfaces on all grid levels for all spectral      
    !       solutions.                                                     
    !
    use constants
    use blockPointers, only : flowDoms, BCData, nDom, nBocos, inBeg, inEnd, &
         jnBeg, jnEnd, knBeg, knEnd, icBeg, icEnd, jcBeg, jcBeg, jcEnd, kcBeg, &
         kcEnd, BCFaceID
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, sps
    integer(kind=intType) :: nLevels, level

    ! Determine the number of grid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of grid levels.

    levelLoop: do level=1,nLevels

       ! Loop over the number of spectral solutions and number of
       ! blocks stored on this processor.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
          domainsLoop: do i=1,nDom

             ! Allocate the memory for the array of the boundary
             ! condition data.

             j = flowDoms(i,level,sps)%nBocos
             allocate(flowDoms(i,level,sps)%BCData(j), stat=ierr)
             if(ierr /= 0)                   &
                  call terminate("initBCData", &
                  "Memory allocation failure for BCData")

             ! Set the pointers to make it more readable.

             call setPointers(i,level,sps)

             ! Copy the range of the subfaces in BCData and nullify its
             ! pointers.

             bocoLoop: do j=1,nBocos

                ! Determine the block face on which the subface is located
                ! and set the dimensions accordingly.

                select case (BCFaceID(j))

                case (iMin,iMax)
                   BCData(j)%inBeg = jnBeg(j)
                   BCData(j)%inEnd = jnEnd(j)
                   BCData(j)%jnBeg = knBeg(j)
                   BCData(j)%jnEnd = knEnd(j)

                   BCData(j)%icbeg = jcbeg(j)
                   BCData(j)%icend = jcend(j)
                   BCData(j)%jcbeg = kcbeg(j)
                   BCData(j)%jcend = kcend(j)

                case (jMin,jMax)
                   BCData(j)%inBeg = inBeg(j)
                   BCData(j)%inEnd = inEnd(j)
                   BCData(j)%jnBeg = knBeg(j)
                   BCData(j)%jnEnd = knEnd(j)

                   BCData(j)%icbeg = icbeg(j)
                   BCData(j)%icend = icend(j)
                   BCData(j)%jcbeg = kcbeg(j)
                   BCData(j)%jcend = kcend(j)

                case (kMin,kMax)
                   BCData(j)%inBeg = inBeg(j)
                   BCData(j)%inEnd = inEnd(j)
                   BCData(j)%jnBeg = jnBeg(j)
                   BCData(j)%jnEnd = jnEnd(j)

                   BCData(j)%icbeg = icbeg(j)
                   BCData(j)%icend = icend(j)
                   BCData(j)%jcbeg = jcbeg(j)
                   BCData(j)%jcend = jcend(j)

                end select

                ! Initialize the boundary condition treatment for
                ! subsonic inlet to noSubInlet.

                BCData(j)%subsonicInletTreatment = noSubInlet

                ! Nullify the pointers of BCData.
                ! Some compilers require this.

                nullify(BCData(j)%norm)
                nullify(BCData(j)%rface)
                nullify(BCData(j)%F)
                nullify(BCData(j)%Fv)
                nullify(BCData(j)%Fp)
                nullify(BCData(j)%T)
                nullify(BCData(j)%tv)
                nullify(BCData(j)%Tp)
                nullify(BCData(j)%area)
                nullify(BCData(j)%fIndex)
                nullify(BCData(j)%uSlip)
                nullify(BCData(j)%TNS_Wall)

                nullify(BCData(j)%normALE)
                nullify(BCData(j)%rfaceALE)
                nullify(BCData(j)%uSlipALE)
                nullify(BCData(j)%sHeatFlux)

                nullify(BCData(j)%ptInlet)
                nullify(BCData(j)%ttInlet)
                nullify(BCData(j)%htInlet)
                nullify(BCData(j)%flowXdirInlet)
                nullify(BCData(j)%flowYdirInlet)
                nullify(BCData(j)%flowZdirInlet)

                nullify(BCData(j)%turbInlet)

                nullify(BCData(j)%rho)
                nullify(BCData(j)%velx)
                nullify(BCData(j)%vely)
                nullify(BCData(j)%velz)
                nullify(BCData(j)%ps)
                bcData(j)%symNormSet = .False.
                bcData(j)%symNorm = zero
                nullify(BCData(j)%iblank)
                nullify(BCData(j)%delta)
                nullify(BCData(j)%deltaNode)

             enddo bocoLoop
          enddo domainsLoop
       enddo spectralLoop
    enddo levelLoop

  end subroutine initBCData

  subroutine nondimBoundData
    !
    !       nondimBoundData nondimensionalizes the boundary data           
    !       specified in the cgns file.                                    
    !
    use constants
    use block, only : flowDoms, BCDataType, nDom
    use flowVarRefState, only : pRef, rhoRef, tRef, muRef, rhoRef, nt1, nt2
    use inputPhysics, only : turbModel
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, sps, nn, mm, i
    real(kind=realType)   :: hRef, uRef

    type(BCDataType), dimension(:), pointer :: BCData

    ! Compute the number of multigrid levels and the values of hRef
    ! and uRef.

    nLevels = ubound(flowDoms,2)
    hRef    = pRef/rhoRef
    uRef    = sqrt(hRef)

    ! Loop over the number of multigrid levels, spectral solutions
    ! and local blocks.

    do mm=1,nLevels
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom

             ! Set the pointer for BCData to make the code more readable.

             BCData => flowDoms(nn,mm,sps)%BCData

             ! Loop over the number of boundary faces.

             do i=1,flowDoms(nn,mm,sps)%nBocos

                ! Nondimensionalize the data if the pointer is associated
                ! with data.

                if( associated(BCData(i)%TNS_Wall) ) &
                     BCData(i)%TNS_Wall = BCData(i)%TNS_Wall/TRef

                if( associated(BCData(i)%ptInlet) ) &
                     BCData(i)%ptInlet = BCData(i)%ptInlet/pRef

                if( associated(BCData(i)%ttInlet) ) &
                     BCData(i)%ttInlet = BCData(i)%ttInlet/TRef

                if( associated(BCData(i)%htInlet) ) &
                     BCData(i)%htInlet = BCData(i)%htInlet/HRef

                if( associated(BCData(i)%turbInlet) ) &
                     call nondimTurb(BCData(i)%turbInlet)

                if( associated(BCData(i)%rho) ) &
                     BCData(i)%rho = BCData(i)%rho/rhoRef

                if( associated(BCData(i)%velx) ) &
                     BCData(i)%velx = BCData(i)%velx/uRef

                if( associated(BCData(i)%vely) ) &
                     BCData(i)%vely = BCData(i)%vely/uRef

                if( associated(BCData(i)%velz) ) &
                     BCData(i)%velz = BCData(i)%velz/uRef

                if( associated(BCData(i)%ps) ) &
                     BCData(i)%ps = BCData(i)%ps/pRef

             enddo
          enddo
       enddo
    enddo

    !=================================================================

  contains

    !===============================================================

    subroutine nondimTurb(turbInlet)
      !
      !         NondimTurb nondimensionalizes the turbulent data for inlet  
      !         boundary conditions.                                         
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(:,:,:), pointer :: turbInlet
      !
      !        Local variables.
      !
      integer(kind=intType) :: nn
      real(kind=realType)   :: nuRef, tmp

      real(kind=realType), dimension(nt1:nt2) :: ref

      ! Set the reference values depending on the turbulence model.

      nuRef = muRef/rhoRef
      select case (turbModel)

      case (spalartAllmaras, spalartAllmarasEdwards)
         ref(itu1) = nuRef

      case (komegaWilcox, komegaModified, menterSST)
         ref(itu1) = pRef/rhoRef
         ref(itu2) = ref(itu1)/nuRef

      case (ktau)
         ref(itu1) = pRef/rhoRef
         ref(itu2) = nuRef/ref(itu1)

      case (v2f)
         ref(itu1) = pRef/rhoRef
         ref(itu4) = ref(itu1)/nuRef
         ref(itu2) = ref(itu1)*ref(itu4)
         ref(itu3) = ref(itu1)

      end select

      ! Loop over the number of turbulence variables and make
      ! them nondimensional.

      do nn=nt1,nt2
         tmp = one/ref(nn)
         turbInlet(:,:,nn) = turbInlet(:,:,nn)*tmp
      enddo

    end subroutine nondimTurb

  end subroutine nondimBoundData

  subroutine radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, nbcVar, &
       cgnsBoco, blockFaceId, indCoor, &
       ind, bcVarPresent, bcVarArray,  &
       axAssumed)
    !
    !       radialInterpolSubface interpolates the prescribed variables    
    !       in the data set of the given cgns subface onto the subface     
    !       indicated by iBeg, iEnd, jBeg, jEnd and blockFaceId.           
    !       This routine performs a 1d interpolation in radial direction   
    !       assuming that there is no variation in the other directions.   
    !       The variables in blockPointers are already set.                
    !
    use constants
    use blockPointers, only: BCData, x, il, jl, kl, nbkglobal, sectionID
    use section, only : sections
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
    integer(kind=intType), intent(in) :: nbcVar, cgnsBoco
    integer(kind=intType), intent(in) :: blockFaceId

    integer(kind=intType), dimension(2), intent(in) :: indCoor
    integer(kind=intType), dimension(2,nbcVar), intent(in) :: ind

    real(kind=realType), dimension(iBeg:iEnd,jBeg:jEnd,nbcVar), &
         intent(out) :: bcVarArray

    logical, intent(inout) :: axAssumed
    logical, dimension(nbcVar), intent(in) :: bcVarPresent
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l, n1, n2, nn, ii
    integer(kind=intType) :: nDim, nPoints, start

    real(kind=realType) :: t, length, fact, rad, ww1, ww2

    real(kind=realType), dimension(3) :: rotAxis, xc, xaxis

    real(kind=realType), dimension(:),     pointer :: rr
    real(kind=realType), dimension(:,:,:), pointer :: xf

    character(len=maxStringLen) :: errorMessage

    type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet

    ! Store the rotation axis of the section a bit easier.

    rotAxis = sections(sectionId)%rotAxis

    ! Check if a rotation axis could be constructed for the section
    ! to which this block belongs. If not it is assumed that the
    ! x-axis is the axial direction.

    length = rotAxis(1)**2 + rotAxis(2)**2 + rotAxis(3)**2

    if(length < half) then

       ! No axis could be constructed from the rotation info.
       ! Assume this is the x-axis and set axAssumed to .True.

       rotAxis(1) = one; rotAxis(2) = zero; rotAxis(3) = zero
       axAssumed = .true.

    endif

    ! Set the pointer for dataSet to make the code more readable.

    dataSet => cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%dataSet

    ! Check the number of dimensions of the specified data set.
    ! This should be 1, because a 1d interpolation in radial
    ! direction is performed.

    k    = indCoor(1)
    l    = indCoor(2)
    nDim = dataSet(k)%dirichletArrays(l)%nDimensions

    if(nDim > 1) then
       write(errorMessage,101) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
       call terminate("radialInterpolSubface", errorMessage)
101    format("Zone",1X,A,", subface",1X,A,": Multidimensional &
            &radially varying data specified. Only 1d data possible")
    endif

    ! Set the pointer for the radial coordinate and abbreviate the
    ! number of interpolation points a bit easier.

    rr     => dataSet(k)%dirichletArrays(l)%dataArr
    nPoints = dataSet(k)%dirichletArrays(l)%dataDim(1)

    ! Check if the data is specified for increasing radius.

    do i=2,nPoints
       if(rr(i) < rr(i-1)) exit
    enddo

    if(i <= nPoints) then
       write(errorMessage,102) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
       call terminate("radialInterpolSubface", errorMessage)
102    format("Zone",1X,A,", subface",1X,A,": Data should be &
            &specified for increasing radius.")
    endif

    ! Set the pointer for the coordinates of the block face on which
    ! the boundary subface is located.

    select case (blockFaceId)
    case (iMin)
       xf => x(1,:,:,:)
    case (iMax)
       xf => x(il,:,:,:)
    case (jMin)
       xf => x(:,1,:,:)
    case (jMax)
       xf => x(:,jl,:,:)
    case (kMin)
       xf => x(:,:,1,:)
    case (kMax)
       xf => x(:,:,kl,:)
    end select

    ! Compute the factor needed to compute the coordinates in the
    ! original units. The fourth comes from the averaging of the 4
    ! nodal coordinates.

    fact = fourth/cgnsDoms(nbkGlobal)%LRef

    ! Loop over the range of the subface.

    jloop: do j=jBeg,jEnd
       iloop: do i=iBeg,iEnd

          ! Determine the coordinates of the face center. Normally this
          ! is an average of i-1, i, j-1, j, but due to the usage of
          ! the pointer xf and the fact that x originally starts at 0,
          ! an offset of 1 is introduced and thus the average should
          ! be taken of i, i+1, j and j+1.

          xc(1) = fact*(xf(i,j,  1) + xf(i+1,j,  1) &
               +       xf(i,j+1,1) + xf(i+1,j+1,1))
          xc(2) = fact*(xf(i,j,  2) + xf(i+1,j,  2) &
               +       xf(i,j+1,2) + xf(i+1,j+1,2))
          xc(3) = fact*(xf(i,j,  3) + xf(i+1,j,  3) &
               +       xf(i,j+1,3) + xf(i+1,j+1,3))

          ! Determine the parameter, which defines the closest point
          ! on the rotation axis. Note that rotAxis is a unit-vector.

          t = rotAxis(1)*(xc(1) - sections(sectionId)%rotCenter(1)) &
               + rotAxis(2)*(xc(2) - sections(sectionId)%rotCenter(2)) &
               + rotAxis(3)*(xc(3) - sections(sectionId)%rotCenter(3))

          ! Determine the coordinates of this point.

          xaxis(1) = sections(sectionId)%rotCenter(1) + t*rotAxis(1)
          xaxis(2) = sections(sectionId)%rotCenter(2) + t*rotAxis(2)
          xaxis(3) = sections(sectionId)%rotCenter(3) + t*rotAxis(3)

          ! Determine the radius of this point.

          rad = sqrt((xc(1)-xaxis(1))**2 + (xc(2)-xaxis(2))**2 &
               +      (xc(3)-xaxis(3))**2)

          ! Determine the interpolation interval and set the
          ! interpolation weights. Take care of the exceptions.

          checkInterpol: if(rad <= rr(1)) then

             ! Radius is less than the minimum value specified.
             ! Use constant extrapolation.

             n1  = 1;   n2  = 1
             ww1 = one; ww2 = zero

          else if(rad >= rr(nPoints)) then checkInterpol

             ! Radius is larger than the maximum value specified.
             ! Use constant extrapolation.

             n1  = nPoints; n2  = nPoints
             ww1 = one;     ww2 = zero

          else checkInterpol

             ! Radius is in the range. Determine the correct interval
             ! using a binary search algorithm.

             ii    = nPoints - 1
             start = 1
             interval: do

                ! Next guess for the interval and determine the new
                ! situation.

                nn = start + ii/2
                if(rad > rr(nn+1)) then

                   ! Rad is larger than the upper boundary of the
                   ! current interval. Update the lower boundary.

                   start = nn + 1; ii = ii - 1

                else if(rad >= rr(nn)) then

                   ! This is the correct range. Exit the loop.

                   exit

                endif

                ! Modify ii for the next branch to search.

                ii = ii/2

             enddo interval

             ! Rad is in the interval nn, nn+1. Store this and
             ! determine the interpolation weight.

             n1 = nn
             n2 = nn + 1
             ww1 = (rr(nn+1) - rad)/(rr(nn+1) - rr(nn))
             ww2 = one - ww1

          endif checkInterpol

          ! Interpolate the values the values present for this face.

          do nn=1,nbcVar
             if( bcVarPresent(nn) ) then

                ! Easier storage of the indices in the data set.

                k = ind(1,nn)
                l = ind(2,nn)

                ! Interpolate this variable.

                bcVarArray(i,j,nn) =                                 &
                     ww1*dataSet(k)%dirichletArrays(l)%dataArr(n1) &
                     + ww2*dataSet(k)%dirichletArrays(l)%dataArr(n2)
             endif
          enddo

       enddo iloop
    enddo jloop

  end subroutine radialInterpolSubface

  subroutine setBCDataCoarseGrid
    !
    !       setBCDataCoarseGrid determines the boundary condition info     
    !       on the coarse grid from the known info on the fine grid. It    
    !       will be stored in the BCData arrays of flowDoms.               
    !
    use constants
    use blockPointers, only : BCFaceID, BCData, nDom, flowDoms, il, jl, kl, &
         mgIFine, mgJFine, mgKFine, nBocos, BCType
    use flowVarRefState, only : nt1, nt2
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : groundLevel
    use utils, only : setPointers
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l, sps
    integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd, iiMax, jjMax
    integer(kind=intType) :: nLevels, level, levm1

    integer(kind=intType), dimension(:,:), pointer :: iFine, jFine

    real(kind=realType) :: var

    real(kind=realType), dimension(3) :: dir

    ! Determine the number of grid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over the coarser grid levels. It is assumed that the
    ! bc data of groundLevel is set correctly.

    coarseLevelLoop: do level=(groundLevel+1),nLevels

       ! Store the fine grid level a bit easier.

       levm1 = level - 1

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
          domainsLoop: do i=1,nDom

             ! Set the pointers to the coarse block.

             call setPointers(i, level, sps)

             ! Loop over the boundary subfaces and interpolate the
             ! prescribed boundary data for this grid level.

             bocoLoop: do j=1,nBocos

                ! Determine the block face on which the subface is
                ! located and set some multigrid variables accordingly.

                select case (BCFaceID(j))

                case (iMin,iMax)
                   iiMax = jl; jjMax = kl
                   iFine => mgJFine; jFine => mgKFine

                case (jMin,jMax)
                   iiMax = il; jjMax = kl
                   iFine => mgIFine; jFine => mgKFine

                case (kMin,kMax)
                   iiMax = il; jjMax = jl
                   iFine => mgIFine; jFine => mgJFine

                end select

                ! Abbreviate the size of the subface a bit easier.

                iBeg = BCData(j)%icBeg; iEnd = BCData(j)%icEnd
                jBeg = BCData(j)%jcBeg; jEnd = BCData(j)%jcEnd

                ! Copy the subsonic boundary conditions treatment.

                BCData(j)%subsonicInletTreatment = &
                     flowDoms(i,levm1,sps)%BCData(j)%subsonicInletTreatment

                ! Interpolate the data for the possible prescribed boundary
                ! data.

                call interpolateBcData(BCData(j)%TNS_Wall, &
                     flowDoms(i,levm1,sps)%BCData(j)%TNS_Wall)

                call interpolateBcData(BCData(j)%ptInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%ptInlet)
                call interpolateBcData(BCData(j)%ttInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%ttInlet)
                call interpolateBcData(BCData(j)%flowXdirInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%flowXdirInlet)
                call interpolateBcData(BCData(j)%flowYdirInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%flowYdirInlet)
                call interpolateBcData(BCData(j)%flowZdirInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%flowZdirInlet)

                call interpolateBCVecData(BCData(j)%turbInlet, &
                     flowDoms(i,levm1,sps)%BCData(j)%turbInlet, &
                     nt1, nt2)

                call interpolateBcData(BCData(j)%rho,  &
                     flowDoms(i,levm1,sps)%BCData(j)%rho)
                call interpolateBcData(BCData(j)%velx, &
                     flowDoms(i,levm1,sps)%BCData(j)%velx)
                call interpolateBcData(BCData(j)%vely, &
                     flowDoms(i,levm1,sps)%BCData(j)%vely)
                call interpolateBcData(BCData(j)%velz, &
                     flowDoms(i,levm1,sps)%BCData(j)%velz)
                call interpolateBcData(BCData(j)%ps,   &
                     flowDoms(i,levm1,sps)%BCData(j)%ps)

                ! Some additional variables should be computed/corrected
                ! for some boundary conditions. Determine the type of
                ! boundary condition.

                if((BCType(j) == SubsonicInflow .and.                         &
                     BCData(j)%subsonicInletTreatment == totalConditions) .or. &
                     BCType(j) == DomainInterfaceTotal) then

                   ! Total conditions are specified for subsonic inflow
                   ! or domain interfaces.

                   ! Compute the total enthalpy and make
                   ! sure that the unit vector is a unit vector.

                   ! Loop over the faces of the subface.

                   do l=jBeg,jEnd
                      do k=iBeg,iEnd

                         ! Compute the total enthalpy.

                         call computeHtot(BCData(j)%ttInlet(k,l), &
                              BCData(j)%htInlet(k,l))

                         ! Flow direction.

                         dir(1) = BCData(j)%flowXdirInlet(k,l)
                         dir(2) = BCData(j)%flowYdirInlet(k,l)
                         dir(3) = BCData(j)%flowZdirInlet(k,l)

                         var = one/max(eps,sqrt(dir(1)**2 + dir(2)**2 &
                              +                  dir(3)**2))

                         BCData(j)%flowXdirInlet(k,l) = var*dir(1)
                         BCData(j)%flowYdirInlet(k,l) = var*dir(2)
                         BCData(j)%flowZdirInlet(k,l) = var*dir(3)

                      enddo
                   enddo

                endif

             enddo bocoLoop
          enddo domainsLoop
       enddo spectralLoop
    enddo coarseLevelLoop

  contains

    subroutine interpolateBcData(varCoarse, varFine)
      !
      !         InterpolateBcData interpolates the given data array from   
      !         the fine to the coarse grid. Of course only if the fine      
      !         array is associated with some data.                          
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(:,:), pointer :: varCoarse
      real(kind=realType), dimension(:,:), pointer :: varFine
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j, if1, if2, jf1, jf2

      ! Check if varFine is associated to data. If not return.

      if(.not. associated(varFine)) return

      ! Loop over the faces of the given subface.
      ! First the j-direction.

      do j=jBeg,jEnd

         ! Determine the two children in this direction. Take care of
         ! the halo's, as this info is only available for owned cells.

         if(j < 2) then
            jf1 = 1; jf2 = 1
         else if(j > jjMax) then
            jf1 = jFine(jjMax,2) +1; jf2 = jf1
         else
            jf1 = jFine(j,1); jf2 = jFine(j,2)
         endif

         ! Loop in the i-direction.

         do i=iBeg,iEnd

            ! Determine the two children in this direction.
            ! Same story as in j-direction.

            if(i < 2) then
               if1 = 1; if2 = 1
            else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
            else
               if1 = iFine(i,1); if2 = iFine(i,2)
            endif

            ! Compute the coarse grid data as the average of the
            ! 4 fine grid values.

            varCoarse(i,j) = fourth*(varFine(if1,jf1) &
                 +         varFine(if2,jf1) &
                 +         varFine(if1,jf2) &
                 +         varFine(if2,jf2))
         enddo
      enddo

    end subroutine interpolateBcData

    subroutine interpolateBCVecData(varCoarse, varFine, &
         nstart, nend)
      !
      !         interpolateBCVecData interpolates the given data array       
      !         from the fine to the coarse grid. Of course only if the fine 
      !         array is associated with some data.                          
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      integer(kind=intType), intent(in) :: nstart, nend

      real(kind=realType), dimension(:,:,:), pointer :: varCoarse
      real(kind=realType), dimension(:,:,:), pointer :: varFine
      !
      !        Local variables.
      !
      integer(kind=intType) :: nn, i, j, if1, if2, jf1, jf2
      ! Check if varFine is associated to data. If not return.

      if(.not. associated(varFine)) return

      ! Loop over the faces of the given subface.
      ! First the j-direction.

      do j=jBeg,jEnd

         ! Determine the two children in this direction. Take care of
         ! the halo's, as this info is only available for owned cells.

         if(j < 2) then
            jf1 = 1; jf2 = 1
         else if(j > jjMax) then
            jf1 = jFine(jjMax,2) +1; jf2 = jf1
         else
            jf1 = jFine(j,1); jf2 = jFine(j,2)
         endif

         ! Loop in the i-direction.

         do i=iBeg,iEnd

            ! Determine the two children in this direction.
            ! Same story as in j-direction.

            if(i < 2) then
               if1 = 1; if2 = 1
            else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
            else
               if1 = iFine(i,1); if2 = iFine(i,2)
            endif

            ! Compute the coarse grid data as the average of the
            ! 4 fine grid values.

            do nn=nstart,nend
               varCoarse(i,j,nn) = fourth*(varFine(if1,jf1,nn) &
                    +         varFine(if2,jf1,nn) &
                    +         varFine(if1,jf2,nn) &
                    +         varFine(if2,jf2,nn))
            enddo
         enddo
      enddo

    end subroutine interpolateBCVecData

  end subroutine setBCDataCoarseGrid

  subroutine setBCDataFineGrid(initializationPart)
    !
    !       setBCDataFineGrid extracts the boundary condition data from    
    !       the cgnsGrid and stores it in useable form in the BCData       
    !       arrays of the currently finest grid, i.e. groundLevel.         
    !
    use constants
    use blockPointers, only : BCData, BCType, nBKGlobal, nBocos, nDom, cgnsSubFace
    use communication, only : sumb_comm_world, myid
    use inputTimeSpectral, only :nTimeIntervalsSpectral
    use iteration, only : groundLevel
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: initializationPart
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, sps

    logical :: allTurbMassBleedInflow,  allTurbSubsonicInflow
    logical :: allFlowSupersonicInflow, allTurbSupersonicInflow

    ! Initialize axAssumed and massflowPrescribed to .false.,
    ! indicating that no assumption is made about the axial direction
    ! and no subsonic inflow boundaries with prescribed mass flow
    ! are present.

    axAssumed          = .false.
    massflowPrescribed = .false.

    ! Initialize all the prescribed turbulence as well as flow
    ! variables for supersonic inlet to .true.

    allTurbMassBleedInflow  = .true.
    allTurbSubsonicInflow   = .true.
    allTurbSupersonicInflow = .true.

    allFlowSupersonicInflow = .true.

    ! Loop over the number of spectral solutions and local blocks.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domainsLoop: do i=1,nDom

          ! Set the pointers to this block on groundLevel to make
          ! the code readable.

          call setPointers(i,groundLevel,sps)

          ! Loop over the number of boundary condition subfaces.

          bocoLoop: do j=1,nBocos

             ! Store the cgns boundary subface number, the number of
             ! boundary condition data sets and the data sets a bit easier.

             cgnsBoco = cgnsSubface(j)
             nDataSet = &
                  cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%nDataSet
             dataSet => &
                  cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%dataSet

             ! Store the range of the boundary subface a bit easier.

             iBeg = BCData(j)%icBeg; iEnd = BCData(j)%icEnd
             jBeg = BCData(j)%jcBeg; jEnd = BCData(j)%jcEnd

             ! Determine the boundary condition we are having here and
             ! call the appropriate routine.

             select case (BCType(j))


             case (NSWallIsothermal)
                call BCDataIsothermalWall(j)

             case (SupersonicInflow)
                call BCDataSupersonicInflow(j, allFlowSupersonicInflow, &
                     allTurbSupersonicInflow)

             case (SubsonicInflow)
                call BCDataSubsonicInflow(j, allTurbSubsonicInflow)

             case (SubsonicOutflow)
                call BCDataSubsonicOutflow(j)

             case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                  DomainInterfaceP,   DomainInterfaceRho,    &
                  DomainInterfaceTotal)
                call terminate('setBCDataFineGrid', &
                     'Domain interface BCs are not fully implemented')
             end select

          enddo bocoLoop
       enddo domainsLoop
    enddo spectralLoop

    ! If this is the initialization part perform some checks
    ! to see if certain assumptions were made.

    checkInit: if( initializationPart ) then

       ! Check whether or not an assumption was made on the axial
       ! direction. If so, processor 0 prints a warning.

       i = 0
       if( axAssumed ) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Radial boundary data given while no &
               &rotation axis is present."
          print "(a)", "# It is assumed that the X-axis is the axial &
               &direction."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

       ! Check whether or not subsonic inflow boundaries are present with
       ! a prescribed mass flow. If so print a warning that the flow
       ! problem should not be a choked one.

       i = 0
       if( massflowPrescribed ) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Subsonic inflow boundaries present with &
               &prescribed mass flow."
          print "(a)", "# This is only a well posed problem if the &
               &flow is not choked."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

       ! Check whether or not mass bleed inflow regions are present
       ! for which the free stream turbulence is used.

       i = 0
       if(.not. allTurbMassBleedInflow) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Inflow bleed regions present for which the &
               &turbulence"
          print "(a)", "# quantities are not or insufficiently &
               &prescribed."
          print "(a)", "# Using free stream values instead."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

       ! Check whether or not subsonic inflow regions are present
       ! for which the free stream turbulence is used.

       i = 0
       if(.not. allTurbSubsonicInflow) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Subsonic inflow regions present for which &
               &the turbulence"
          print "(a)", "# quantities are not or insufficiently &
               &prescribed."
          print "(a)", "# Using free stream values instead."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

       ! Check whether or not supersonic inflow regions are present
       ! for which the free stream variables is used.

       i = 0
       if(.not. allFlowSupersonicInflow) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Supersonic inflow regions present for which &
               &the flow variables"
          print "(a)", "# are not or insufficiently prescribed."
          print "(a)", "# Using free stream values instead."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

       ! Check whether or not supersonic inflow regions are present
       ! for which the free stream turbulence is used.

       i = 0
       if(.not. allTurbSupersonicInflow) i = 1
       call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
            SUmb_comm_world, ierr)

       if(myID == 0 .and. j == 1) then

          print "(a)", "#"
          print "(a)", "#*==================== !!! Warning !!! &
               &======================"
          print "(a)", "# Supersonic inflow regions present for which &
               &the turbulence"
          print "(a)", "# quantities are not or insufficiently &
               &prescribed."
          print "(a)", "# Using free stream values instead."
          print "(a)", "#*=====================================&
               &======================"
          print "(a)", "#"

       endif

    endif checkInit

  end subroutine setBCDataFineGrid

  subroutine setBCVarNamesTurb(offset)
    !
    !       setBCVarNamesTurb sets the names for the turbulence            
    !       variables to be determined. This depends on the turbulence     
    !       model. If not the RANS equations are solved an immediate       
    !       return is made.                                                
    !
    use constants
    use cgnsNames
    use inputPhysics, only : equations, turbModel
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: offset

    ! Return immediately if not the RANS equations are solved.

    if(equations /= RANSEquations) return

    ! Determine the turbulence model and set the names accordingly.

    select case (turbModel)
    case (spalartAllmaras, spalartAllmarasEdwards)
       bcVarNames(offset+1) = cgnsTurbSaNu

    case (komegaWilcox, komegaModified, menterSST)
       bcVarNames(offset+1) = cgnsTurbK
       bcVarNames(offset+2) = cgnsTurbOmega

    case (ktau)
       bcVarNames(offset+1) = cgnsTurbK
       bcVarNames(offset+2) = cgnsTurbTau

    case (v2f)
       bcVarNames(offset+1) = cgnsTurbK
       bcVarNames(offset+2) = cgnsTurbEpsilon
       bcVarNames(offset+3) = cgnsTurbV2
       bcVarNames(offset+4) = cgnsTurbF

    end select

  end subroutine setBCVarNamesTurb

  !=================================================================

  logical function setBCVarTurb(offset, boco, turbInlet)
    !
    !       SetBCVarTurb sets the array for the turbulent halo data        
    !       for inlet boundaries. This function returns .true. If all      
    !       turbulence variables could be interpolated and .false.         
    !       otherwise.                                                     
    !
    use constants
    use flowVarRefState, only : nt1, nt2
    use inputPhysics, only : equations
    use utils, only : terminate, siTurb

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: offset, boco
    real(kind=realType), dimension(:,:,:), pointer :: turbInlet
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, i, j
    real(kind=realType)   :: mult, trans

    ! Initialize setBCVarTurb to .true. And return immediately
    ! if not the rans equations are solved.

    setBCVarTurb = .true.
    if(equations /= RANSEquations) return

    ! Loop over the number of turbulent variables. mm is the counter
    ! in the arrays bcVarArray and bcVarPresent.

    mm = offset
    turbLoop: do nn=nt1,nt2
       mm = mm + 1

       ! Check if the variable is present. If so, use the
       ! interpolated data.

       if( bcVarPresent(mm) ) then

          ! Conversion to SI units if possible.

          call siTurb(mass(mm), length(mm), time(mm), temp(mm), &
               bcVarNames(mm), mult, trans)

          ! Set the turbulent variables.

          do j=jBeg,jEnd
             do i=iBeg,iEnd
                turbInlet(i,j,nn) = mult*bcVarArray(i,j,mm) + trans
             enddo
          enddo

       else

          ! Turbulent variable not present.
          ! Set setBCVarTurb to .false.

          setBCVarTurb = .false.

       endif
    enddo turbLoop

    ! If not all the turbulence variables could be interpolated
    ! reset them to zero and store this face in the array
    ! turbFreestreamSubfaces, which stores these subfaces.
    ! They will be set to free stream values later on.

    if(.not. setBCVarTurb) then

       call storeTurbFreestreamSubface

       do nn=nt1,nt2
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                turbInlet(i,j,nn) = zero
             enddo
          enddo
       enddo
    endif

  contains

    subroutine storeTurbFreestreamSubface
      !
      !         storeTurbFreestreamSubface stores the currently active       
      !         subface in the array turbFreestreamSubfaces, such that the   
      !         turbulence variables can be set to the free stream values    
      !         later on in setInletFreestreamTurb.                          
      !
      use blockPointers, only: nbkLocal, spectralSol
      implicit none
      !
      !        Local variables.
      !
      integer :: ierr

      integer(kind=intType) :: nn

      integer(kind=intType), dimension(:,:), allocatable :: tmp

      ! Determine the situation we are dealing with here.

      testAllocated: if( allocated(turbFreestreamSubfaces) ) then

         ! TurbFreestreamSubfaces has already been allocated and
         ! thus contains information. It must be reallocated and the
         ! current subface should be added.

         ! Allocate the memory for tmp and copy the data from
         ! turbFreestreamSubfaces.

         allocate(tmp(nTurbFreestreamSubfaces,3), stat=ierr)
         if(ierr /= 0)                                  &
              call terminate("storeTurbFreestreamSubface", &
              "Memory allocation failure for tmp")

         tmp = turbFreestreamSubfaces

         ! Release turbFreestreamSubfaces and allocate it again
         ! with an increased dimension.

         deallocate(turbFreestreamSubfaces, stat=ierr)
         if(ierr /= 0)                                  &
              call terminate("storeTurbFreestreamSubface", &
              "Deallocation failure for &
              &turbFreestreamSubfaces")

         nTurbFreestreamSubfaces = nTurbFreestreamSubfaces + 1

         allocate(turbFreestreamSubfaces(nTurbFreestreamSubfaces,3), &
              stat=ierr)
         if(ierr /= 0)                                  &
              call terminate("storeTurbFreestreamSubface", &
              "Memory allocation failure for &
              &turbFreestreamSubfaces")

         ! Copy the data back from tmp into turbFreestreamSubfaces
         ! and release tmp again.

         do nn=1,(nTurbFreestreamSubfaces-1)
            turbFreestreamSubfaces(nn,1) = tmp(nn,1)
            turbFreestreamSubfaces(nn,2) = tmp(nn,2)
            turbFreestreamSubfaces(nn,3) = tmp(nn,3)
         enddo

         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                                  &
              call terminate("storeTurbFreestreamSubface", &
              "Deallocation failure for tmp")

      else testAllocated

         ! TurbFreestreamSubfaces has not been allocated yet. This
         ! is the first subface to store in this array. Allocate the
         ! array and set nTurbFreestreamSubfaces to 1.

         nTurbFreestreamSubfaces = 1
         allocate(turbFreestreamSubfaces(nTurbFreestreamSubfaces,3), &
              stat=ierr)
         if(ierr /= 0)                                  &
              call terminate("storeTurbFreestreamSubface", &
              "Memory allocation failure for &
              &turbFreestreamSubfaces")

      endif testAllocated

      ! Store the current subface in turbFreestreamSubfaces.

      nn = nTurbFreestreamSubfaces
      turbFreestreamSubfaces(nn,1) = nbkLocal
      turbFreestreamSubfaces(nn,2) = boco
      turbFreestreamSubfaces(nn,3) = spectralSol

    end subroutine storeTurbFreestreamSubface

  end function setBCVarTurb

  subroutine setInletFreestreamTurb
    !
    !       setInletFreestreamTurb sets for all boundary subfaces          
    !       stored in turbFreestreamSubfaces the turbulence variables to   
    !       the free stream variables. This is done for all multigrid      
    !       levels starting from groundLevel.                              
    !
    use constants
    use block, only : BCDataType, flowDoms
    use flowVarRefState, only : nt1, nt2, wInf
    use iteration, only : groundLevel
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, i, j, l
    integer(kind=intType) :: mm, boco, sps
    integer(kind=intType) :: level, nLevels

    type(BCDataType), dimension(:), pointer :: BCData

    ! If turbFreestreamSubfaces has not been allocated there are
    ! no such subfaces present and an immediate return should be made.

    if(.not. allocated(turbFreestreamSubfaces) ) return

    ! Determine the maximum number of multigrid levels present
    ! in flowDoms.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of flagged subfaces.

    subfaceLoop: do nn=1,nTurbFreestreamSubfaces

       ! Store the local block ID, the boundary subface and the
       ! spectral solution a bit easier.

       mm   = turbFreestreamSubfaces(nn,1)
       boco = turbFreestreamSubfaces(nn,2)
       sps  = turbFreestreamSubfaces(nn,3)

       ! Loop over the number multigrid levels, starting at the
       ! current groundLevel.

       do level=groundLevel,nLevels

          ! Set the pointer for BCData to make the code more readable.

          BCData => flowDoms(mm,level,sps)%BCData

          ! Loop over the range of this subface.

          do j=BCData(boco)%jcBeg,BCData(boco)%jcEnd
             do i=BCData(boco)%icBeg,BCData(boco)%icEnd

                ! Set the turbulence variables.

                do l=nt1,nt2
                   BCData(boco)%turbInlet(i,j,l) = wInf(l)
                enddo
             enddo
          enddo
       enddo

    enddo subfaceLoop

    ! Release the memory of turbFreestreamSubfaces.

    deallocate(turbFreestreamSubfaces, stat=ierr)
    if(ierr /= 0) &
         call terminate("setInletFreestreamTurb", &
         "Deallocation failure for &
         &turbFreestreamSubfaces")

  end subroutine setInletFreestreamTurb

  subroutine setSupersonicInletFreeStream
    !
    !       setSupersonicInletFreeStream sets for all boundary subfaces    
    !       stored in freestreamSubfaces the primitive flow variables to   
    !       the free stream variables. This is done for all multigrid      
    !       levels starting from groundLevel.                              
    !
    use constants
    use block, only : BCDataType, flowDoms
    use flowVarRefState, only : nt1, nt2, wInf, pInfCorr
    use iteration, only :  groundLevel
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, i, j
    integer(kind=intType) :: mm, boco, sps
    integer(kind=intType) :: level, nLevels

    type(BCDataType), dimension(:), pointer :: BCData

    ! Return immediately if freestreamSubfaces has not been allocated,
    ! i.e. if no subfaces present for which the inflow data should
    ! be set to free stream values.

    if(.not. allocated(freestreamSubfaces)) return

    ! Determine the maximum number of multigrid levels present
    ! in flowDoms.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of flagged subfaces.

    subfaceLoop: do nn=1,nFreestreamSubfaces

       ! Store the local block ID, the boundary subface and the
       ! spectral solution a bit easier.

       mm   = freestreamSubfaces(nn,1)
       boco = freestreamSubfaces(nn,2)
       sps  = freestreamSubfaces(nn,3)

       ! Loop over the number multigrid levels, starting at the
       ! current groundLevel.

       do level=groundLevel,nLevels

          ! Set the pointer for BCData to make the code more readable.

          BCData => flowDoms(mm,level,sps)%BCData

          ! Loop over the range of this subface.

          do j=BCData(boco)%jcBeg,BCData(boco)%jcEnd
             do i=BCData(boco)%icBeg,BCData(boco)%icEnd

                ! Set the flow field variables.

                BCData(boco)%rho(i,j)  = wInf(irho)
                BCData(boco)%velx(i,j) = wInf(ivx)
                BCData(boco)%vely(i,j) = wInf(ivy)
                BCData(boco)%velz(i,j) = wInf(ivz)
                BCData(boco)%ps(i,j)   = pInfCorr

             enddo
          enddo
       enddo

    enddo subfaceLoop

    ! Release the memory of freestreamSubfaces.

    deallocate(freestreamSubfaces, stat=ierr)
    if(ierr /= 0)                                    &
         call terminate("setSupersonicInletFreeStream", &
         "Deallocation failure for freestreamSubfaces")

  end subroutine setSupersonicInletFreeStream

  subroutine storeFreestreamSubface(boco)
    !
    !       storeFreestreamSubface stores the currently active subface in  
    !       the array freestreamSubfaces, such that the primitive flow     
    !       field variables can be set to the free stream values later on  
    !       in setSupersonicInletFreeStream.                               
    !
    use constants
    use blockPointers, only : nbkLocal, spectralSol
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: boco
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn

    integer(kind=intType), dimension(:,:), allocatable :: tmp

    ! Determine the situation we are dealing with here.

    testAllocated: if( allocated(freestreamSubfaces) ) then

       ! freestreamSubfaces has already been allocated and thus
       ! contains information. It must be reallocated and the current
       ! subface should be added.

       ! Allocate the memory for tmp and copy the data from
       ! freestreamSubfaces.

       allocate(tmp(nFreestreamSubfaces,3), stat=ierr)
       if(ierr /= 0)                              &
            call terminate("storeFreestreamSubface", &
            "Memory allocation failure for tmp")

       tmp = freestreamSubfaces

       ! Release freestreamSubfaces and allocate it again
       ! with an increased dimension.

       deallocate(freestreamSubfaces, stat=ierr)
       if(ierr /= 0)                              &
            call terminate("storeFreestreamSubface", &
            "Deallocation failure for freestreamSubfaces")

       nFreestreamSubfaces = nFreestreamSubfaces + 1

       allocate(freestreamSubfaces(nFreestreamSubfaces,3), &
            stat=ierr)
       if(ierr /= 0)                              &
            call terminate("storeFreestreamSubface", &
            "Memory allocation failure for &
            &freestreamSubfaces")

       ! Copy the data back from tmp into freestreamSubfaces
       ! and release tmp again.

       do nn=1,(nFreestreamSubfaces-1)
          freestreamSubfaces(nn,1) = tmp(nn,1)
          freestreamSubfaces(nn,2) = tmp(nn,2)
          freestreamSubfaces(nn,3) = tmp(nn,3)
       enddo

       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                              &
            call terminate("storeFreestreamSubface", &
            "Deallocation failure for tmp")

    else testAllocated

       ! freestreamSubfaces has not been allocated yet. This is the
       ! first subface to store in this array. Allocate the array and
       ! set nFreestreamSubfaces to 1.

       nFreestreamSubfaces = 1
       allocate(freestreamSubfaces(nFreestreamSubfaces,3), &
            stat=ierr)
       if(ierr /= 0)                              &
            call terminate("storeFreestreamSubface", &
            "Memory allocation failure for &
            &freestreamSubfaces")

    endif testAllocated

    ! Store the current subface in freestreamSubfaces.

    nn = nFreestreamSubfaces
    freestreamSubfaces(nn,1) = nbkLocal
    freestreamSubfaces(nn,2) = boco
    freestreamSubfaces(nn,3) = spectralSol

  end subroutine storeFreestreamSubface

  subroutine unitVectorsCylSystem(boco)
    !
    !       unitVectorsCylSystem determines the unit vectors of the        
    !       local coordinate systen of the boundary face defined by the    
    !       data in BCDataMod. In that local system the axial direction    
    !       is rotation axis.                                              
    !
    use constants
    use blockPointers, only : BCFaceID, BCData, x, si, sj, sk, il, jl, kl, &
         sectionID
    use section, only : sections
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: boco
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j
    real(kind=realType)   :: factInlet, var

    real(kind=realType), dimension(3) :: dir

    real(kind=realType), dimension(:,:,:), pointer :: ss

    ! Set the pointers for coordinates and normals of the block
    ! face on which this subface is located. Set factInlet
    ! such that factInlet*normals points into the domain.

    select case (BCFaceID(boco))
    case (iMin)
       xf => x(1,:,:,:);  ss => si(1 ,:,:,:); factInlet =  one
    case (iMax)
       xf => x(il,:,:,:); ss => si(il,:,:,:); factInlet = -one
    case (jMin)
       xf => x(:,1,:,:);  ss => sj(:,1 ,:,:); factInlet =  one
    case (jMax)
       xf => x(:,jl,:,:); ss => sj(:,jl,:,:); factInlet = -one
    case (kMin)
       xf => x(:,:,1,:);  ss => sk(:,:,1 ,:); factInlet =  one
    case (kMax)
       xf => x(:,:,kl,:); ss => sk(:,:,kl,:); factInlet = -one
    end select

    ! Loop over the physical range of the subface to store the sum of
    ! the normals. Note that jBeg, jEnd, iBeg, iEnd cannot be used
    ! here, because they may include the halo faces. Instead the
    ! nodal range is used, which defines the original subface. The
    ! offset of +1 in the start index is there because you need
    ! the face id's.

    dir(1) = zero; dir(2) = zero; dir(3) = zero

    do j=(BCData(boco)%jnBeg+1), BCData(boco)%jnEnd
       do i=(BCData(boco)%inBeg+1), BCData(boco)%inEnd
          dir(1) = dir(1) + ss(i,j,1)
          dir(2) = dir(2) + ss(i,j,2)
          dir(3) = dir(3) + ss(i,j,3)
       enddo
    enddo

    ! Multiply by factInlet to make sure that the normal
    ! is inward pointing.

    dir(1) = dir(1)*factInlet
    dir(2) = dir(2)*factInlet
    dir(3) = dir(3)*factInlet

    ! Determine three unit vectors, which define the local cartesian
    ! coordinate system of the rotation axis. First the axial
    ! direction. If the axis cannot be determined from rotation info,
    ! it is assumed to be the x-axis.

    axis = sections(sectionId)%rotAxis
    var  = axis(1)**2 + axis(2)**2 + axis(3)**2
    if(var < half) then

       ! No rotation axis specified. Assume the x-axis
       ! and set the logical axAssumed to .True.

       axis(1) = one; axis(2) = zero; axis(3) = zero
       axAssumed = .true.
    endif

    ! The axial axis must be such that it points into the
    ! computational domain. If the dot product with dir is
    ! negative the direction of axis should be reversed.

    var = axis(1)*dir(1) + axis(2)*dir(2) + axis(3)*dir(3)
    if(var < zero) then
       axis(1) = -axis(1); axis(2) = -axis(2); axis(3) = -axis(3)
    endif

    ! Two unit vectors define the radial plane. These vectors are
    ! defined up to a constants. Just pick a direction for the second
    ! and create a unit vector normal to axis.

    if(abs(axis(2)) < 0.707107_realType) then
       radVec1(1) = zero; radVec1(2) = one;  radVec1(3) = zero
    else
       radVec1(1) = zero; radVec1(2) = zero; radVec1(3) = one
    endif

    var = radVec1(1)*axis(1) + radVec1(2)*axis(2) &
         + radVec1(3)*axis(3)
    radVec1(1) = radVec1(1) - var*axis(1)
    radVec1(2) = radVec1(2) - var*axis(2)
    radVec1(3) = radVec1(3) - var*axis(3)

    var = one/sqrt(radVec1(1)**2 + radVec1(2)**2 &
         +          radVec1(3)**2)
    radVec1(1) = radVec1(1)*var
    radVec1(2) = radVec1(2)*var
    radVec1(3) = radVec1(3)*var

    ! The second vector of the radial plane is obtained
    ! by taking the cross product of axis and radVec1.

    radVec2(1) = axis(2)*radVec1(3) - axis(3)*radVec1(2)
    radVec2(2) = axis(3)*radVec1(1) - axis(1)*radVec1(3)
    radVec2(3) = axis(1)*radVec1(2) - axis(2)*radVec1(1)

  end subroutine unitVectorsCylSystem


  subroutine initBCDataDomainInterfaces
    !
    !       initBCDataDomainInterfaces initializes the prescribed boundary 
    !       data for domain interfaces. It is just an initialization such  
    !       the initial halo computations behave normally. During the      
    !       actual computation this data should be renewed constantly by   
    !       the coupler.                                                   
    !
    use constants
    use blockPointers, only : BCData, nBocos, BCType, nDom
    use flowVarRefState, only : nt1, nt2
    use inputIteration, only : mgStartLevel
    use inputTimeSpectral, only :nTimeIntervalsSpectral
    use BCPointers, only : ww2, pp2
    use utils, only : setPointers, setBCPointers
    use flowUtils, only : computePtot, computeTtot
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, nn, mm, l, sps

    real(kind=realType) :: rho, vvx, vvy, vvz, pres, vtotInv

    ! Loop over the number of spectral solutions and blocks of
    ! the multigrid start level.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, mgStartlevel, sps)

          ! Loop over the boundary condition subfaces of this block.

          bocos: do mm=1,nBocos

             ! Check for interface boundary condition types.

             select case (BCType(mm))

             case (DomainInterfaceAll)

                ! All data must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the state vector to set the prescribed state.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                      BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                      BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                      BCData(mm)%velz(i,j) = ww2(i,j,ivz)
                      BCData(mm)%ps(i,j)   = pp2(i,j)

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceRhoUVW)

                ! Density, velocity and turbulent variables are
                ! prescribed. Nullify the pointers and set them to the
                ! correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the state vector to set the prescribed state.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                      BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                      BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                      BCData(mm)%velz(i,j) = ww2(i,j,ivz)

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceP)

                ! Pressure must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the pressure to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                      BCData(mm)%ps(i,j) = pp2(i,j)
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceRho)

                ! Density must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the density to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                      BCData(mm)%rho(i,j) = ww2(i,j,irho)
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceTotal)

                ! Total conditions must be prescribed. Nullify the
                ! pointers and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the total conditions and the turbulence variables
                ! to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      ! Store the variables a bit easier.

                      rho  = ww2(i,j,irho)
                      vvx   = ww2(i,j,ivx)
                      vvy   = ww2(i,j,ivy)
                      vvz   = ww2(i,j,ivz)
                      pres = pp2(i,j)

                      ! Compute the total pressure, total temperature
                      ! and total entahlpy.

                      call computePtot(rho, vvx, vvy, vvz, pres, &
                           BCData(mm)%ptInlet(i,j))
                      call computeTtot(rho, vvx, vvy, vvz, pres, &
                           BCData(mm)%ttInlet(i,j))
                      BCData(mm)%htInlet(i,j) = (ww2(i,j,irhoE) + pres) &
                           / rho

                      ! Determine the velocity direction.

                      vtotInv = one/max(eps,sqrt(vvx*vvx + vvy*vvy + vvz*vvz))
                      BCData(mm)%flowXdirInlet(i,j) = vvx*vtotInv
                      BCData(mm)%flowYdirInlet(i,j) = vvy*vtotInv
                      BCData(mm)%flowZdirInlet(i,j) = vvz*vtotInv

                      ! Simply extrapolate the turbulence variables.

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

             end select

          enddo bocos
       enddo domains
    enddo spectralLoop

  end subroutine initBCDataDomainInterfaces

end module BCData
