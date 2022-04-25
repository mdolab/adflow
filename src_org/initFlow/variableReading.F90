module variableReading

  use constants, only : intType, maxCGNSNameLen, cgnsRealType, realType, maxStringLen
  use su_cgns, only : cgsize_t
  ! halosRead:Determines if the halos where read or not.
  logical :: halosRead

  ! cgnsInd:  File index of the CGNS file.
  ! cgnsBase: Base of the CGNS file, always set to 1.
  ! cgnsZone: Zone ID in the CGNS file.
  ! cgnsSol:  Solution ID in the zone ID of the CGNS file.
  ! location: Location where the variables are stored in CGNS.
  !           Supported possibilities are Vertex and CellCentered.

  integer :: cgnsInd, cgnsBase, cgnsZone, cgnsSol, location

  ! zoneNumbers: Corresponding zoneNumbers of the sorted
  !              zoneNames.
  ! zoneNames:   Zone names, sorted in increasing order, of the
  !              zones in the CGNS restart file.
  ! varNames: Variable names, sorted in increasing order,
  !           of the variables.

  integer(kind=intType), allocatable, dimension(:) :: zoneNumbers
  character(len=maxCGNSNameLen), allocatable, dimension(:) :: zoneNames
  character(len=maxCGNSNameLen), allocatable, dimension(:) :: varNames

  ! rangeMin(3):    Lower index in i, j and k direction of the
  !                 range to be read.
  ! rangeMax(3):    Upper index in i, j and k direction of the
  !                 range to be read.
  integer(kind=cgsize_t), dimension(3) :: rangeMin, rangeMax

  ! nVar:     Number of variables stored in the solution file.
  ! solID:    Loop variables for the number of solutions to be read.

  integer :: nVar
  integer(kind=intType) :: solID

  ! interpolSpectral:     Whether or not to interpolate the
  !                       coordinates/solutions for the time
  !                       spectral mode.
  ! copySpectral:         Whether or not to copy the solutions
  !                       for the time spectral mode.
  logical ::               interpolSpectral, copySpectral



  ! rhoScale: Scale factor for the density.
  ! velScale: Scale factor for the velocity.
  ! pScale:   Scale factor for the pressure.
  ! muScale:  Scale factor for the molecular viscosity.

  real(kind=realType) :: rhoScale, velScale, pScale, muScale

  ! nSolsRead:            Number of solution files to read.
  ! solFiles(nSolsRead):  Names of the solution files to be read.

  integer(kind=intType) :: nSolsRead
  character(len=maxStringLen), dimension(:), allocatable :: solFiles

  ! varTypes(nVar): Variable types of the variables stored.
  integer, allocatable, dimension(:) :: varTypes


  ! buffer(2:il,2:jl,2:kl): Buffer to read and store the cell
  !                         centered values.
  ! bufferVertex(:):        Additional buffer needed to read
  !                         vertex data and transform them into
  !                         cell centered values.

  real(kind=cgnsRealType), dimension(:,:,:), allocatable :: buffer
  real(kind=cgnsRealType), dimension(:,:,:), allocatable :: bufferVertex


contains
  subroutine readDensity(nTypeMismatch)
    !
    !       readDensity reads the density from the given place in the
    !       cgns file. If the density itself is not stored (unlikely),
    !       then it is tried to construct the density from other
    !       variables. If this is not possible an error message is printed
    !       and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the density is present in the solution file.

    nn = bsearchStrings(cgnsDensity, varNames)
    if(nn > 0) then

       ! Density is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the density from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the
       ! scaling factor to obtain to correct nondimensional value and
       ! take the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,irho) = buffer(i,j,k)*rhoScale
             enddo
          enddo
       enddo

       ! Density is read, so a return can be made.

       return

    endif

    ! Not able to determine the density.
    ! Print an error message and exit.

    call terminate("readDensity", &
         "Not able to retrieve density from the &
         &variables in the restart file.")

  end subroutine readDensity

  subroutine readEnergy(nTypeMismatch)
    !
    !       readEnergy reads the energy variable from the given place in
    !       the cgns file. If the energy is not stored then it is tried to
    !       construct it from the pressure, density and velocities. If it
    !       is not possible to create the energy an error message is
    !       printed and the program will stop. It is assumed that the
    !       pointers in blockPointers already point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    use flowUtils, only : eTot
    use flowVarRefState, only : kPresent
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: vvx, vvy, vvz, dummyK, pres, rhoInv

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the total energy is present in the solution file.

    nn = bsearchStrings(cgnsEnergy, varNames)

    testRhoEPresent: if(nn > 0) then

       ! Total energy is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the energy from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scaling
       ! factor to obtain to correct non-dimensional value and take the
       ! possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,irhoE) = buffer(i,j,k)*pScale
             enddo
          enddo
       enddo

       ! Energy has been read, so a return can be made.

       return

    endif testRhoEPresent

    ! Total energy is not present. Check for the pressure.

    nn = bsearchStrings(cgnsPressure, varNames)

    testPressure: if(nn > 0) then

       ! Pressure is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the pressure from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Compute the total energy. This depends whether or not
       ! a turbulent kinetic energy is present. Take the possible
       ! pointer offset into account.
       ! As this routine is only called to construct the states in
       ! the past for a time accurate computation, the momentum is
       ! stored and not the velocity.

       if( kPresent ) then

          do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
                jp = j+po
                do i=iBeg,iEnd
                   ip = i+po
                   rhoInv = one/w(ip,jp,kp,irho)
                   vvx     = w(ip,jp,kp,imx)*rhoInv
                   vvy     = w(ip,jp,kp,imy)*rhoInv
                   vvz     = w(ip,jp,kp,imz)*rhoInv
                   pres   = buffer(i,j,k)*pScale
                   call etot(w(ip,jp,kp,irho), vvx, vvy, vvz, pres,  &
                        w(ip,jp,kp,itu1), w(ip,jp,kp,irhoE), &
                        kPresent)
                enddo
             enddo
          enddo

       else

          dummyK = zero

          do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
                jp = j+po
                do i=iBeg,iEnd
                   ip = i+po
                   rhoInv = one/w(ip,jp,kp,irho)
                   vvx     = w(ip,jp,kp,imx)*rhoInv
                   vvy     = w(ip,jp,kp,imy)*rhoInv
                   vvz     = w(ip,jp,kp,imz)*rhoInv
                   pres   = buffer(i,j,k)*pScale
                   call etot(w(ip,jp,kp,irho), vvx, vvy, vvz, pres, &
                        dummyK, w(ip,jp,kp,irhoE), kPresent)
                enddo
             enddo
          enddo

       endif

       ! Energy has been created. So a return can be made.

       return

    endif testPressure

    ! Energy could not be created. Terminate.

    call terminate("readEnergy", &
         "Energy could not be created")

  end subroutine readEnergy

  subroutine readPressure(nTypeMismatch)
    !
    !       readPressure reads the pressure variable from the given place
    !       in the cgns file. If the pressure itself is not present it is
    !       tried to construct if from other variables. In that case it is
    !       assumed that the density, velocity and turbulent variables are
    !       already stored in the pointer variable w.
    !       If it is not possible to create the pressure an error message
    !       is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use flowUtils, only : computePressure
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the pressure is present in the solution file.

    nn = bsearchStrings(cgnsPressure, varNames)
    if(nn > 0) then

       ! Pressure is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the pressure from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into the position of rhoE
       ! in w. Multiply by the pressure scale factor to obtain the
       ! correct nondimensional value and take the possible pointer
       ! offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,irhoE) = buffer(i,j,k)*pScale
             enddo
          enddo
       enddo

       ! Pressure is read, so a return can be made.

       return

    endif

    ! Pressure is not present. Check for the total energy.

    nn = bsearchStrings(cgnsEnergy, varNames)
    if(nn > 0) then

       ! Total energy is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the total energy from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the
       ! pressure scale factor to obtain the correct nondimensional
       ! value and take the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,irhoE) = buffer(i,j,k)*pScale
             enddo
          enddo
       enddo

       ! Compute the pressure from energy, density and velocities.
       ! This will still be stored in the irhoE position of w.

       call computePressure(iBeg,iEnd,jBeg,jEnd,kBeg,kEnd,po)

       ! Pressure is constructed, so a return can be made.

       return

    endif

    ! Not able to determine the pressure.
    ! Print an error message and exit.

    call terminate("readPressure", &
         "Not able to retrieve the pressure from &
         &the variables in the restart file.")

  end subroutine readPressure

  subroutine readTurbEddyVis(nTypeMismatch, eddyVisPresent)
    !
    !       readTurbEddyVis tries to read the eddy viscosity from the
    !       restart file.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal, il, jl, kl, rev
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    use flowVarRefState, only : muInf
    use inputPhysics, only : eddyVisInfRatio
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    logical, intent(out)                 :: eddyVisPresent
    !
    !      Local variables.
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn

    ! Set the cgns real type

    realTypeCGNS = setCGNSRealType()

    ! Check if the eddy viscosity is present. If so, read it.

    nn = bsearchStrings(cgnsEddy, varNames)

    if(nn > 0) then

       ! Eddy viscosity is present. Determine if a type mismatch
       ! occured, read the buffer from the file and set
       ! eddyVisPresent to .true.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       call readRestartVariable(varNames(nn))

       eddyVisPresent = .true.

       ! Scale the eddy viscosity such that it contains the
       ! correct nonDimensional value.

       do k=2,kl
          do j=2,jl
             do i=2,il
                rev(i,j,k) = muScale*buffer(i,j,k)
             enddo
          enddo
       enddo

       ! Eddy viscosity has been read, so make a return.

       return

    endif

    ! Eddy viscosity is not present. Check if the eddy viscosity
    ! ratio is present. If so read it and construct the eddy
    ! viscosity from it.

    nn = bsearchStrings(cgnsEddyRatio, varNames)

    if(nn > 0) then

       ! Eddy viscosity ratio is present. Determine if a type
       ! mismatch occured, read the buffer from the file and set
       ! eddyVisPresent to .true.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       call readRestartVariable(varNames(nn))

       eddyVisPresent = .true.

       ! Multiply the eddy viscosity by the laminar viscosity such
       ! that it contains the correct nonDimensional value.
       ! As the laminar viscosity is not yet know, use the free
       ! stream viscosity.

       do k=2,kl
          do j=2,jl
             do i=2,il
                rev(i,j,k) = muInf*buffer(i,j,k)
             enddo
          enddo
       enddo

       ! Eddy viscosity has been read, so make a return.

       return

    endif

    ! Eddy viscosity cannot be constructed. Set it to the
    ! free stream eddy viscosity.

    do k=2,kl
       do j=2,jl
          do i=2,il
             rev(i,j,k) = muInf*eddyVisInfRatio
          enddo
       enddo
    enddo

    ! Eddy viscosity is not present, so set it to .false.

    eddyVisPresent = .false.

  end subroutine readTurbEddyVis

  subroutine readTurbKwType(nTypeMismatch)
    !
    !       readTurbKwType reads or constructs the k and omega values
    !       for two-equations turbulence models of the k-omega type.
    !       If no information could be retrieved some engineering guess of
    !       the turbulent variables is made.
    !
    use constants
    use cgnsNames
    use communication, only : myid
    use blockPointers, only : w, nbklocal, il, jl, kl, rlv, rev
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    use flowVarRefState, only : muInf
    use inputPhysics, only : turbModel
    use turbUtils, only : initKOmega
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables.
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: nuScale, kScale, omegaScale, val

    logical :: turbKPresent, omegaPresent, eddyVisPresent

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Compute the scales for nu, k and omega.

    nuScale    = muScale/rhoScale
    kScale     = velScale**2
    omegaScale = kScale/nuScale

    ! First check if k is present.

    turbKPresent = .false.

    nn = bsearchStrings(cgnsTurbK, varNames)

    if(nn > 0) then

       ! K is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read k from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct non-dimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu1) = kScale*buffer(i,j,k)
             enddo
          enddo
       enddo

       ! Set turbKPresent to .true.

       turbKPresent = .true.

    endif

    ! Check if omega is present.

    omegaPresent = .false.

    nn = bsearchStrings(cgnsTurbOmega, varNames)

    if(nn > 0) then

       ! Omega is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read omega from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct non-dimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu2) = omegaScale*buffer(i,j,k)
             enddo
          enddo
       enddo

       ! Set omegaPresent to .true.

       omegaPresent = .true.

    endif

    ! If omega is not present, check if tau is present and
    ! initialize omega accordingly.

    if(.not. omegaPresent) then

       nn = bsearchStrings(cgnsTurbTau, varNames)

       if(nn > 0) then

          ! Tau is present. First determine whether or not a type
          ! mismatch occurs. If so, update nTypeMismatch.

          if(realTypeCGNS /= varTypes(nn)) &
               nTypeMismatch = nTypeMismatch + 1

          ! Read tau from the restart file and store it in buffer.

          call readRestartVariable(varNames(nn))

          ! Transform tau to omega and copy the variables from buffer
          ! into w. Multiply by the scale factor to obtain the correct
          ! non-dimensional value and take the possible pointer offset
          ! into account.

          do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
                jp = j+po
                do i=iBeg,iEnd
                   ip = i+po

                   val = buffer(i,j,k)
                   w(ip,jp,kp,itu2) = omegaScale/max(eps,val)
                enddo
             enddo
          enddo

          ! Set omegaPresent to .true.

          omegaPresent = .true.

       endif

    endif

    ! Check if both variables were present.
    ! If so go to the check to transform omega to tau.

    if(turbKPresent .and. omegaPresent) goto 1001

    ! K and omega are not both present. It is tried to construct
    ! their values with the information that is present.

    ! Try to read the eddy viscosity.

    call readTurbEddyVis(nTypeMismatch, eddyVisPresent)

    ! The eddy viscosity is either known or still initialized
    ! to the free stream value. In any case determine the
    ! situation we are dealing with and try to initialize k and
    ! omega accordingly.

    if( turbKPresent ) then

       ! K is present. Compute omega using the eddy viscosity.
       ! Assume that the standard k-omega formula is also valid
       ! for the SST-model.
       ! Take the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu2) = w(ip,jp,kp,irho)*w(ip,jp,kp,itu1) &
                     / rev(i,j,k)
             enddo
          enddo
       enddo

       ! Print a warning that omega was not present and has been
       ! constructed. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

          print "(a)", "#"
          print "(a)", "#                 Warning"
          print "(a)", "# Omega is not present in the restart file."
          if( eddyVisPresent ) then
             print "(a)", "# It is initialized using the turbulent &
                  &kinetic energy and eddy viscosity."
          else
             print "(a)", "# It is initialized using the turbulent &
                  &kinetic energy and free stream eddy &
                  &viscosity."
          endif
          print "(a)", "#"

       endif

       ! K and omega are initialized.
       ! Go to the check to transform omega to tau.

       goto 1001

    endif

    if( omegaPresent ) then

       ! Omega is present. Compute k using the eddy viscosity.
       ! Assume that the standard k-omega formula is also valid
       ! for the SST-model.
       ! Take the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu1) = rev(i,j,k)*w(ip,jp,kp,itu2) &
                     / w(ip,jp,kp,irho)
             enddo
          enddo
       enddo

       ! Print a warning that k was not present and has been
       ! constructed. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

          print "(a)", "#"
          print "(a)", "#                 Warning"
          print "(a)", "# Turbulent kinetic energy is not present &
               &in the restart file."
          if( eddyVisPresent ) then
             print "(a)", "# It is initialized using omega and &
                  &the eddy viscosity."
          else
             print "(a)", "# It is initialized using omega and &
                  &the free stream eddy viscosity."
          endif
          print "(a)", "#"

       endif

       ! K and omega are initialized.
       ! Go to the check to transform omega to tau.

       goto 1001

    endif

    ! Both k and omega are not present. Use a guess for omega
    ! and compute k using the known value of the eddy viscosity.
    ! As the laminar viscosity is not yet known, set it to the
    ! free-stream value.

    rlv = muInf
    call initKOmega(po)

    ! Print a warning that both k and omega are not present in
    ! the restart file. Only processor 0 does this for block 1.

    if((myID == 0) .and. (nbkLocal == 1)) then

       print "(a)", "#"
       print "(a)", "#                 Warning"
       print "(a)", "# The turbulent kinetic energy and omega are &
            &not present in the restart file."
       if( eddyVisPresent ) then
          print "(a)", "# They have been initialized using the &
               &eddy viscosity."
       else
          print "(a)", "# The default initialization has been used."
       endif

       print "(a)", "#"

    endif

    ! For the k-tau model omega must be transformed to tau.
    ! Take the possible pointer offset into account.

1001 select case (turbModel)

    case (ktau)

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu2) = one/w(ip,jp,kp,itu2)
             enddo
          enddo
       enddo

    end select

  end subroutine readTurbKwType

  subroutine readTurbSA(nTypeMismatch)
    !
    !       readTurbSA reads or constructs the nu tilde transport
    !       variable for the Spalart-Allmaras type turbulence models.
    !       If no information could be retrieved some engineering guess of
    !       the turbulent variables is made.
    !
    use constants
    use cgnsNames
    use communication, only : myid
    use blockPointers, only : w, nbklocal, rlv, rev
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    use flowVarRefState, only : muInf, wInf
    use turbUtils, only : saNuKnownEddyRatio
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables.
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: nuScale, ratio, nu

    logical :: eddyVisPresent

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.
    ! Also compute the kinematic viscosity scale.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    nuScale = muScale/rhoScale

    ! Check if the nu tilde variable is present.

    nn = bsearchStrings(cgnsTurbSANu, varNames)

    nuTildePresent: if(nn > 0) then

       ! Nu tilde is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read nu tilde from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w and take the possible
       ! pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu1) = nuScale*buffer(i,j,k)
             enddo
          enddo
       enddo

       ! Variable is read, so a return can be made.

       return

    endif nuTildePresent

    ! NuTilde is not present. Try to construct the eddy viscosity.

    call readTurbEddyVis(nTypeMismatch, eddyVisPresent)

    ! Check if the eddy viscosity has been constructed.

    eddyPresent: if( eddyVisPresent ) then

       ! Eddy viscosity is present. As the laminar viscosity is not
       ! yet known, set it to the free-stream value.

       rlv = muInf

       ! Compute nuTilde from the known ratio of eddy and laminar
       ! viscosity. Take the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po

                ! Compute the eddy viscosity ratio and the laminar
                ! kinematic viscosity and call the function to
                ! compute the nu tilde variable.

                ratio            = rev(i,j,k)/rlv(i,j,k)
                nu               = rlv(i,j,k)/w(ip,jp,kp,irho)
                w(ip,jp,kp,itu1) = saNuKnownEddyRatio(ratio, nu)

             enddo
          enddo
       enddo

       ! Print a warning that nu tilde has been constructed and
       ! not read. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

          print "(a)", "#"
          print "(a)", "#                 Warning"
          print "(a)", "# Nu tilde for Spalart-Allmaras model not &
               &present in the restart file."
          print "(a)", "# Variable has been reconstructed from &
               &the eddy viscosity ratio."
          print "(a)", "#"

       endif

       ! Variable is constructed, so a return can be made.

       return

    endif eddyPresent

    ! No turbulence info is present in the restart file.
    ! Initialize nu tilde to the free stream value.
    ! Take the possible pointer offset into account.

    do k=kBeg,kEnd
       kp = k+po
       do j=jBeg,jEnd
          jp = j+po
          do i=iBeg,iEnd
             ip = i+po
             w(ip,jp,kp,itu1) = wInf(itu1)
          enddo
       enddo
    enddo

    ! Print a warning that nu tilde has been set to the
    ! free stream values. Only processor 0 does this for block 1.

    if((myID == 0) .and. (nbkLocal == 1)) then

       print "(a)", "#"
       print "(a)", "#                 Warning"
       print "(a)", "# No turbulent info present in the restart file."
       print "(a)", "# Nu tilde for Spalart-Allmaras model has &
            &been set to the free stream value."
       print "(a)", "#"

    endif

  end subroutine readTurbSA
  subroutine readTurbV2f(nTypeMismatch)
    !
    !       readTurbV2f reads or constructs the four transport variables
    !       for the v2f model. If no information could be retrieved some
    !       engineering guess of the turbulent variables is made.
    !
    use constants
    use cgnsNames
    use communication, only : myid
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    use flowVarRefState, only : wInf, nt1, nt2
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, itu
    integer, dimension(4) :: indW

    integer(kind=intType) :: i, j, k, ii, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: nuScale, kScale, epsScale, fScale

    real(kind=realType), dimension(4) :: turbScale

    character(len=maxCGNSNameLen), dimension(4) :: namesVar

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Set the names and indices for the four variables.

    indW(1) = itu1; namesVar(1) = cgnsTurbK
    indW(2) = itu2; namesVar(2) = cgnsTurbEpsilon
    indW(3) = itu3; namesVar(3) = cgnsTurbV2
    indW(4) = itu4; namesVar(4) = cgnsTurbF

    ! Compute the scales for nu, k, epsilon and f; v2 has the same
    ! scaling as k.

    nuScale  = muScale/rhoScale
    kScale   = velScale**2
    fScale   = kScale/nuScale
    epsScale = kScale*fScale

    turbScale(1) = kScale
    turbScale(2) = epsScale
    turbScale(3) = kScale
    turbScale(4) = fScale

    ! Loop over the four variables of the v2f model.

    varLoop: do ii=1,4

       ! Find the index of the variable in the solution file and check
       ! if it is present. If not exit the loop.

       nn = bsearchStrings(namesVar(ii), varNames)

       if(nn == 0) exit

       ! Variable is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the variable from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w.
       ! Take the possible pointer offset into account.

       itu = indW(ii)

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,itu) = turbScale(ii)*buffer(i,j,k)
             enddo
          enddo
       enddo

    enddo varLoop

    ! Check if all variables were present. If not, set all turbulence
    ! variables to the free-stream values.

    testPresent: if(ii <= 4) then

       ! Not all variables are present. Set all 4 to the free-stream
       ! values. Take the possible pointer offset into account.

       do ii=nt1,nt2
          do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
                jp = j+po
                do i=iBeg,iEnd
                   ip = i+po
                   w(ip,jp,kp,ii) = wInf(ii)
                enddo
             enddo
          enddo
       enddo

       ! Print a warning that the turbulence has been initialized to
       ! the free-stream. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

          print "(a)", "#"
          print "(a)", "#                 Warning"
          print "(a)", "# Not all turbulence variables are present &
               &for the v2f model."
          print "(a)", "# They have been initialized to the free &
               &stream values."
          print "(a)", "#"

       endif

    endif testPresent

  end subroutine readTurbV2f

  subroutine readTurbvar(nTypeMismatch)
    !
    !       readTurbvar controls the reading of the turbulent variables
    !       for a restart. It calls the routine, which corresponds to the
    !       turbulence model used.
    !
    use constants
    use communication, only : myid, adflow_comm_world
    use inputPhysics, only : equations, turbModel
    use utils, only : terminate
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables.
    !
    integer :: ierr

    ! Check if the rans equations must be solved. If not return.

    if(equations /= RANSEquations) return

    ! Determine the turbulence model to be used and call the
    ! appropriate subroutine.

    select case (turbModel)

    case (spalartAllmaras, spalartAllmarasEdwards)
       call readTurbSA(nTypeMismatch)

       ! !===============================================================

       ! case (komegaWilcox, komegaModified, menterSST, ktau)
       !   call readTurbKwType(nTypeMismatch)

       ! !===============================================================

       ! case (v2f)
       !   call readTurbV2f(nTypeMismatch)

       !===============================================================

    case default
       if(myID == 0) &
            call terminate("readTurbvar", "Restart not implemented &
            &for this turbulence model.")
       call mpi_barrier(ADflow_comm_world, ierr)

    end select

  end subroutine readTurbvar

  subroutine readXmomentum(nTypeMismatch)
    !
    !       readXmomentum reads the x-momentum variable from the given
    !       place in the cgns file. If the x-momentum itself is not stored
    !       then it is tried to construct it from the x-velocity and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the x-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings

    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, mm, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: momScale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Compute the momentum scaling factor, set the cgns real type and
    ! abbreviate the solution variable and the pointer offset to
    ! improve readability.

    momScale     = rhoScale*velScale
    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the X-momentum is present in the solution file.

    nn = bsearchStrings(cgnsMomX, varNames)

    testMxPresent: if(nn > 0) then

       ! X-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the x-momentum from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct non-dimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imx) = buffer(i,j,k)*momScale
             enddo
          enddo
       enddo

       ! X-momentum is read, so a return can be made.

       return

    endif testMxPresent

    ! X-momentum is not present. Check for x-velocity.

    nn = bsearchStrings(cgnsVelX, varNames)

    testVxPresent: if(nn > 0) then

       ! X-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the x-velocity from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the
       ! density and velocity scaling factor to obtain to correct
       ! non-dimensional value. Take the possible pointer offset
       ! into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imx) = buffer(i,j,k)*w(ip,jp,kp,irho)*velScale
             enddo
          enddo
       enddo

       ! X-momentum is constructed, so a return can be made.

       return

    endif testVxPresent

    ! X-momentum could not be created. Terminate.

    call terminate("readXmomentum", &
         "X-Momentum could not be created")

  end subroutine readXmomentum

  subroutine readXvelocity(nTypeMismatch)
    !
    !       readXvelocity reads the x-velocity variable from the given
    !       place in the cgns file. If the x-velocity itself is not stored
    !       then it is tried to construct it from the x-momentum and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the x-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: scale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the x-velocity is present in the solution file.

    nn = bsearchStrings(cgnsVelX, varNames)
    if(nn > 0) then

       ! X-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the x-velocity from the restart file and store it
       ! in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct nondimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivx) = buffer(i,j,k)*velScale
             enddo
          enddo
       enddo

       ! X-velocity is read, so a return can be made.

       return

    endif

    ! X-velocity not present. Check for x-momentum.

    nn = bsearchStrings(cgnsMomX, varNames)
    if(nn > 0) then

       ! X-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the x-momentum from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Construct the x-velocity; it is assumed that the density is
       ! already stored in w. Multiply by the momentum scale factor
       ! to obtain the correct nondimensional value and take the
       ! possible pointer offset into account.

       scale = rhoScale*velScale

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivx) = buffer(i,j,k)*scale/w(ip,jp,kp,irho)
             enddo
          enddo
       enddo

       ! X-velocity is constructed, so a return can be made.

       return

    endif

    ! Not able to determine the x-velocity.
    ! Print an error message and exit.

    call terminate("readXvelocity", &
         "Not able to retrieve x-velocity from the &
         &variables in the restart file.")

  end subroutine readXvelocity

  subroutine readYmomentum(nTypeMismatch)
    !
    !       readYmomentum reads the y-momentum variable from the given
    !       place in the cgns file. If the y-momentum itself is not stored
    !       then it is tried to construct it from the y-velocity and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the y-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: momScale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Compute the momentum scaling factor, set the cgns real type and
    ! abbreviate the solution variable and the pointer offset to
    ! improve readability.

    momScale     = rhoScale*velScale
    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the Y-momentum is present in the solution file.

    nn = bsearchStrings(cgnsMomY, varNames)

    testMyPresent: if(nn > 0) then

       ! Y-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the y-momentum from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct non-dimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imy) = buffer(i,j,k)*momScale
             enddo
          enddo
       enddo

       ! Y-momentum is read, so a return can be made.

       return

    endif testMyPresent

    ! Y-momentum is not present. Check for y-velocity.

    nn = bsearchStrings(cgnsVelY, varNames)

    testVyPresent: if(nn > 0) then

       ! Y-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the y-velocity from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the
       ! density and velocity scaling factor to obtain to correct
       ! non-dimensional value. Take the possible pointer offset
       ! into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imy) = buffer(i,j,k)*w(ip,jp,kp,irho)*velScale
             enddo
          enddo
       enddo

       ! Y-momentum is constructed, so a return can be made.

       return

    endif testVyPresent

    ! Y-momentum could not be created. Terminate.

    call terminate("readYmomentum", &
         "Y-Momentum could not be created")

  end subroutine readYmomentum

  subroutine readYvelocity(nTypeMismatch)
    !
    !       readYvelocity reads the y-velocity variable from the given
    !       place in the cgns file. If the y-velocity itself is not stored
    !       then it is tried to construct it from the y-momentum and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the y-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings

    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: scale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the y-velocity is present in the solution file.

    nn = bsearchStrings(cgnsVelY, varNames)
    if(nn > 0) then

       ! Y-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the y-velocity from the restart file and store it
       ! in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct nondimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivy) = buffer(i,j,k)*velScale
             enddo
          enddo
       enddo

       ! Y-velocity is read, so a return can be made.

       return

    endif

    ! Y-velocity not present. Check for y-momentum.

    nn = bsearchStrings(cgnsMomY, varNames)
    if(nn > 0) then

       ! Y-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the y-momentum from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Construct the y-velocity; it is assumed that the density is
       ! already stored in w. Multiply by the momentum scale factor
       ! to obtain the correct nondimensional value and take the
       ! possible pointer offset into account.

       scale = rhoScale*velScale

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivy) = buffer(i,j,k)*scale/w(ip,jp,kp,irho)
             enddo
          enddo
       enddo

       ! Y-velocity is constructed, so a return can be made.

       return

    endif

    ! Not able to determine the y-velocity.
    ! Print an error message and exit.

    call terminate("readYvelocity", &
         "Not able to retrieve y-velocity from the &
         &variables in the restart file.")

  end subroutine readYvelocity
  subroutine readZmomentum(nTypeMismatch)
    !
    !       readZmomentum reads the z-momentum variable from the given
    !       place in the cgns file. If the z-momentum itself is not stored
    !       then it is tried to construct it from the z-velocity and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the z-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: momScale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Compute the momentum scaling factor, set the cgns real type and
    ! abbreviate the solution variable and the pointer offset to
    ! improve readability.

    momScale     = rhoScale*velScale
    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the Z-momentum is present in the solution file.

    nn = bsearchStrings(cgnsMomZ, varNames)

    testMzPresent: if(nn > 0) then

       ! Z-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the z-momentum from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct non-dimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imz) = buffer(i,j,k)*momScale
             enddo
          enddo
       enddo

       ! Z-momentum is read, so a return can be made.

       return

    endif testMzPresent

    ! Z-momentum is not present. Check for z-velocity.

    nn = bsearchStrings(cgnsVelZ, varNames)

    testVzPresent: if(nn > 0) then

       ! Z-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the z-velocity from the restart file and store it in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the
       ! density and velocity scaling factor to obtain to correct
       ! non-dimensional value. Take the possible pointer offset
       ! into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,imz) = buffer(i,j,k)*w(ip,jp,kp,irho)*velScale
             enddo
          enddo
       enddo

       ! Z-momentum is constructed, so a return can be made.

       return

    endif testVzPresent

    ! Z-momentum could not be created. Terminate.

    call terminate("readZmomentum", &
         "Z-Momentum could not be created")

  end subroutine readZmomentum
  subroutine readZvelocity(nTypeMismatch)
    !
    !       readZvelocity reads the z-velocity variable from the given
    !       place in the cgns file. If the z-velocity itself is not stored
    !       then it is tried to construct it from the z-momentum and
    !       density; it is assumed that the latter is already stored in
    !       the pointer variable w.
    !       If it is not possible to create the z-velocity an error
    !       message is printed and the program will stop.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use cgnsNames
    use blockPointers, only : w, nbklocal
    use IOModule, only : IOVar
    use utils, only : setCGNSRealType, terminate
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: nTypeMismatch
    !
    !      Local variables
    !
    integer :: realTypeCGNS

    integer(kind=intType) :: i, j, k, nn, po, ip, jp, kp
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    real(kind=realType) :: scale

    ! Set the cell range to be copied from the buffer.

    iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
    jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
    kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

    ! Set the cgns real type and abbreviate the solution variable and
    ! the pointer offset to improve readability.

    realTypeCGNS = setCGNSRealType()

    po = IOVar(nbkLocal,solID)%pointerOffset
    w => IOVar(nbkLocal,solID)%w

    ! Find out if the z-velocity is present in the solution file.

    nn = bsearchStrings(cgnsVelZ, varNames)
    if(nn > 0) then

       ! Z-velocity is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the z-velocity from the restart file and store it
       ! in buffer.

       call readRestartVariable(varNames(nn))

       ! Copy the variables from buffer into w. Multiply by the scale
       ! factor to obtain the correct nondimensional value and take
       ! the possible pointer offset into account.

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivz) = buffer(i,j,k)*velScale
             enddo
          enddo
       enddo

       ! Z-velocity is read, so a return can be made.

       return

    endif

    ! Z-velocity not present. Check for z-momentum.

    nn = bsearchStrings(cgnsMomZ, varNames)
    if(nn > 0) then

       ! Z-momentum is present. First determine whether or not a type
       ! mismatch occurs. If so, update nTypeMismatch.

       if(realTypeCGNS /= varTypes(nn)) &
            nTypeMismatch = nTypeMismatch + 1

       ! Read the z-momentum from the restart file and store
       ! it in buffer.

       call readRestartVariable(varNames(nn))

       ! Construct the z-velocity; it is assumed that the density is
       ! already stored in w. Multiply by the momentum scale factor
       ! to obtain the correct nondimensional value and take the
       ! possible pointer offset into account.

       scale = rhoScale*velScale

       do k=kBeg,kEnd
          kp = k+po
          do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
                ip = i+po
                w(ip,jp,kp,ivz) = buffer(i,j,k)*scale/w(ip,jp,kp,irho)
             enddo
          enddo
       enddo

       ! Z-velocity is constructed, so a return can be made.

       return

    endif

    ! Not able to determine the z-velocity.
    ! Print an error message and exit.

    call terminate("readZvelocity", &
         "Not able to retrieve z-velocity from the &
         &variables in the restart file.")

  end subroutine readZvelocity

  subroutine readTimeHistory(fileIDs)
    !
    !       readTimeHistory attempts to read the time history of an
    !       unsteady computation from the given cgns restart file.
    !       If present it will be stored in the arrays timeArray and
    !       timeDataArray, for which memory is allocated.
    !
    use constants
    use cgnsNames
    use su_cgns
    use inputUnsteady, only : nTimeStepsFine
    use monitor, only : timeDataArray, nMon, nTimeStepsRestart, &
         timeUnsteadyRestart, monNames, timeArray
    use sorting, only : qsortStrings, bsearchStrings
    use utils, only : setCGNSRealType, terminate, allocTimeArrays
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, dimension(nSolsRead), intent(in) :: fileIDs

    !
    !      Local variables.
    !
    integer :: ierr, realTypeCGNS, dummyInt
    integer :: i, nConv, nDim
    integer(kind=cgsize_t) :: nSize(1)

    integer(kind=intType) :: j, ii, nn

    integer(kind=intType), dimension(:), allocatable :: ind

    character(len=maxCGNSNameLen) :: cgnsName
    character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
         convNames, tmpNames

    logical :: allConvInfo

    ! Store the file ID and the base a bit easier. Note that the time
    ! history only needs to be present in the first solution file.

    cgnsInd  = fileIDs(1)
    cgnsBase = 1

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Check if the time history is present by trying to read it.

    call cg_biter_read_f(cgnsInd, cgnsBase, cgnsName, &
         dummyInt, ierr)
    if(ierr /= all_ok) then

       ! No time history present. Set nTimeStepsRestart and
       ! timeUnsteadyRestart to zero, allocate the memory for the
       ! time history of the monitoring variables, print a warning
       ! and return.

       nTimeStepsRestart   = 0
       timeUnsteadyRestart = zero

       call allocTimeArrays(nTimeStepsFine)

       print "(a)", "#"
       print "(a)", "#                 Warning"
       print "(a)", "# No time history found in restart file."
       print "(a)", "# Starting at timestep 1 on the finest level."
       print "(a)", "#"

       return

    endif

    ! Store the number of old time levels.

    nTimeStepsRestart = dummyInt

    ! Go to the place in the cgns file where the time history
    ! should be stored.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "BaseIterativeData_t", 1, "end")
    if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
         "Something wrong when calling cg_goto_f")

    ! Find out how many convergence variables are stored.

    call cg_narrays_f(nConv, ierr)
    if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
         "Something wrong when calling cg_narrays_f")

    ! Allocate the memory for convNames, tmpNames and ind.

    allocate(convNames(nConv), tmpNames(nConv), ind(nConv), &
         stat=ierr)
    if(ierr /= 0)                       &
         call terminate("readTimeHistory", &
         "Memory allocation failure for convNames, etc.")

    ! Read the names of the convergence variables. Store them in
    ! tmpNames as well. Furthermore check the dimension of the
    ! data stored.

    do i=1,nConv
       call cg_array_info_f(i, convNames(i), dummyInt, nDim, &
            nSize, ierr)
       if(ierr /= all_ok)                  &
            call terminate("readConvHistory", &
            "Something wrong when calling cg_array_info_f")

       if(nDim /= 1) then
          print "(a)", "#"
          print "(a)", "#                 Warning"
          print 100, trim(convNames(i))
          print "(a)", "# Information is ignored."
          print "(a)", "#"
100       format("# Dimension of time history for",1X,A,1X, &
               "is not 1.")

          ! Screw up the string such that it does not correspond to
          ! a legal name. It is appended, because it is important that
          ! all strings differ.

          convNames(i) = convNames(i)//"#$@&^!#$%!"
       endif

       if(nSize(1) /= nTimeStepsRestart) then
          print "(a)", "#"
          print "(a)", "#                 Warning"
          print 110, trim(convNames(i))
          print "(a)", "# Displayed information might be incorrect."
          print "(a)", "#"
110       format("# Inconsistent time history for",1X,A,".")
       endif

       ! Copy the name in tmpNames for the sorting.

       tmpNames(i) = convNames(i)
    enddo

    ! Sort convNames in increasing order.

    nn = nConv
    call qsortStrings(convNames, nn)

    ! Find the numbers for the just sorted convergence names.

    do i=1,nConv
       ii      = bsearchStrings(tmpNames(i), convNames)
       ind(ii) = i
    enddo

    ! Find out whether the old time values are present.
    ! If not the time history stored will be ignored.

    ii = bsearchStrings(cgnsTimeValue, convNames)
    if(ii == 0) then
       print "(a)", "#"
       print "(a)", "#                 Warning"
       print "(a)", "# No time values found in the time history &
            &in the restart file."
       print "(a)", "# The rest of the time history is ignored."
       print "(a)", "# Starting at timestep 1 on the finest level."
       print "(a)", "#"

       ! Set nTimeStepsRestart and timeUnsteadyRestart to 0.

       nTimeStepsRestart   = 0
       timeUnsteadyRestart = zero
    endif

    ! Determine the total number of time levels and allocate the
    ! memory for the time history arrays.

    j = nTimeStepsRestart + nTimeStepsFine
    call allocTimeArrays(j)

    ! If the time values are not present, jump to the place where
    ! the memory of the variables used in this routine is released.

    if(ii == 0) goto 99

    ! Read the time values.

    i = ind(ii)
    call cg_array_read_as_f(i, realTypeCGNS, timeArray, ierr)
    if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
         "Something wrong when calling &
         &cg_array_read_as_f")

    ! Set the value of timeUnsteadyRestart to the last value in
    ! timeArray.

    timeUnsteadyRestart = timeArray(nTimeStepsRestart)

    ! Initialize allConvInfo to .true. and perform the loop over
    ! the number of monitoring variables.

    allConvInfo = .true.

    do j=1,nMon

       ! Search for the monitoring name in the sorted
       ! convergence names present in the restart file.

       ii = bsearchStrings(monNames(j), convNames)

       ! Check if the name was found.

       if(ii == 0) then

          ! Name not present in the restart file. Set allConvInfo
          ! to .false. and the corresponding entries in timeDataArray
          ! to zero

          allConvInfo = .false.
          do i=1,nTimeStepsRestart
             timeDataArray(i,j) = zero
          enddo

       else

          ! Name is present in the restart file. Read the corresponding
          ! time history.

          i = ind(ii)
          call cg_array_read_as_f(i, realTypeCGNS, &
               timeDataArray(1,j), ierr)
          if(ierr /= all_ok)                  &
               call terminate("readTimeHistory", &
               "Something wrong when calling &
               &cg_array_read_as_f")
       endif

    enddo

    ! Print a warning in case not all the time history could
    ! be retrieved from the restart file.

    if(.not. allConvInfo) then

       print "(a)", "#"
       print "(a)", "#                 Warning"
       print "(a)", "# Not all the time history could be &
            &retrieved from the restart file."
       print "(a)", "# Missing information is initialized to zero."
       print "(a)", "#"

    endif

99  continue

    ! Release the memory of the variables allocated in this routine.

    deallocate(convNames, tmpNames, ind, stat=ierr)
    if(ierr /= 0)                       &
         call terminate("readTimeHistory", &
         "Deallocation error for convNames, etc.")

  end subroutine readTimeHistory

  subroutine scaleFactors(fileIDs)
    !
    !       scaleFactors determines the scale factors for the density,
    !       pressure and velocity from either the reference state in the
    !       given cgns file or they are simply set to 1.0; the latter
    !       occurs if the input parameter checkRestartSol is .false.
    !       If no reference state is present checkRestartSol is .true.
    !       An error message will be printed and the program terminates.
    !
    use constants
    use cgnsNames
    use su_cgns
    use communication, only : myID, adflow_comm_world
    use flowVarRefState, only : pRef, muRef, rhoRef
    use inputIO, only : checkRestartSol
    use sorting, only : bsearchStrings, qsortStrings
    use utils, only : setCGNSRealType, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, dimension(nSolsRead), intent(in) :: fileIDs

    !      Local variables.
    !
    integer :: ierr, realTypeCGNS, typeCGNS
    integer :: i, nDim, nRef
    integer :: nsize
    integer(kind=cgsize_t) :: nsize2(1)
    integer(kind=intType) :: ii, nn

    integer(kind=intType), dimension(:), allocatable :: ind

    real(kind=cgnsRealType) :: tmpScale

    character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
         refNames, tmpNames

    logical :: muScalePresent

    ! Store the file ID and the base a bit easier. Note that the
    ! reference state only needs to be present in the first file.

    cgnsInd  = fileIDs(1)
    cgnsBase = 1

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Initialize the scale factors to 1.0, i.e. assume that the
    ! correct non-dimensional solution is stored in the restart file.

    rhoScale = one
    velScale = one
    pScale   = one
    muScale  = one

    ! Go to the place in the cgns file where the reference state
    ! should be stored.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "end")
    if(ierr /= all_ok)                &
         call terminate("scaleFactors", &
         "Something wrong when calling cg_goto_f")

    ! Try going to the reference state node. If we get an error code,
    ! it doesn't exist.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "ReferenceState_t", 1, "end")
    if(ierr /= all_ok) then

       ! Reference state does not exist. Check if the restart solution
       ! must be checked. If not, return; otherwise print an error
       ! message and terminate the execution. This error message is
       ! only printed by processor 0 to avoid a messy output.

       if(.not. checkRestartSol) return

       if(myId == 0)                    &
            call terminate("scaleFactors", &
            "Reference state not presented in restart &
            &file. Scaling factors cannot be determined.")

       ! The other processors will wait until they are killed.

       call mpi_barrier(ADflow_comm_world, ierr)

    endif

    ! Go to the reference state node.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "ReferenceState_t", 1, "end")
    if(ierr /= all_ok)                &
         call terminate("scaleFactors", &
         "Something wrong when calling cg_goto_f")

    ! Found out how many reference variables are stored.

    call cg_narrays_f(nRef, ierr)
    if(ierr /= all_ok)               &
         call terminate("scaleFactors", &
         "Something wrong when calling cg_narrays_f")

    ! Allocate the memory for refNames, tmpNames and ind.

    allocate(refNames(nRef), tmpNames(nRef), ind(nRef), &
         stat=ierr)
    if(ierr /= 0)                     &
         call terminate("scaleFactors", &
         "Memory allocation failure for refNames, etc.")

    ! Read the names of the reference variables. Store them in
    ! tmpNames as well.

    do i=1,nRef
       call cg_array_info_f(i, refNames(i), typeCGNS, nDim, &
            nsize2, ierr)
       if(ierr /= all_ok)                &
            call terminate("scaleFactors", &
            "Something wrong when calling cg_array_info_f")

       ! Check the dimension and the size of the array.
       ! Both should be 1. If not, screw up the name such that it
       ! will never be found in the search later on.

       if(nDim /= 1 .or. nsize2(1) /= 1) &
            refNames(i) = refNames(i)//"#$@&^!#$%!"

       ! And copy it in tmpNames.

       tmpNames(i) = refNames(i)
    enddo

    ! Sort refNames in increasing order.

    nn = nRef
    call qsortStrings(refNames, nn)

    ! Find the numbers for the just sorted reference names.

    do i=1,nRef
       ii      = bsearchStrings(tmpNames(i), refNames)
       ind(ii) = i
    enddo

    ! Determine the scale factors if these must be determined.

    if( checkRestartSol ) then

       ! Read the reference density from the restart file.

       ii = bsearchStrings(cgnsDensity, refNames)
       if(ii == 0) then
          if(myId == 0)                    &
               call terminate("scaleFactors", &
               "No reference density found in restart file")

          ! The other processors will wait until they are killed.

          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       i = ind(ii)
       call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
       if(ierr /= all_ok)               &
            call terminate("scaleFactors", &
            "Something wrong when calling &
            &cg_array_read_as_f")
       rhoScale = tmpScale

       ! Read the reference pressure from the restart file.

       ii = bsearchStrings(cgnsPressure, refNames)
       if(ii == 0) then
          if(myId == 0)                    &
               call terminate("scaleFactors", &
               "No reference pressure found in &
               &restart file")

          ! The other processors will wait until they are killed.

          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       i = ind(ii)
       call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
       if(ierr /= all_ok)                &
            call terminate("scaleFactors", &
            "Something wrong when calling &
            &cg_array_read_as_f")
       pScale = tmpScale

       ! Read the reference velocity from the restart file.

       ii = bsearchStrings(cgnsVelocity, refNames)
       if(ii == 0) then
          if(myId == 0)                    &
               call terminate("scaleFactors", &
               "No reference velocity found in &
               &restart file")

          ! The other processors will wait until they are killed.

          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       i = ind(ii)
       call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
       if(ierr /= all_ok)                &
            call terminate("scaleFactors", &
            "Something wrong when calling &
            &cg_array_read_as_f")
       velScale = tmpScale

       ! Set muScalePresent to .true. to indicate that it is present
       ! and read or construct the reference molecular viscosity.

       muScalePresent = .true.

       ii = bsearchStrings(cgnsViscMol, refNames)
       if(ii > 0) then

          ! Scale is present; read the value.

          i = ind(ii)
          call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
          if(ierr /= all_ok)                &
               call terminate("scaleFactors", &
               "Something wrong when calling &
               &cg_array_read_as_f")
          muScale = tmpScale

       else

          ! Try to read the kinematic viscosity.

          ii = bsearchStrings(cgnsViscKin, refNames)
          if(ii > 0) then

             ! Scale is present; read the value and multiply it by the
             ! density.

             i = ind(ii)
             call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
             if(ierr /= all_ok)                &
                  call terminate("scaleFactors", &
                  "Something wrong when calling &
                  &cg_array_read_as_f")

             muScale = tmpScale*rhoScale

          else

             ! Final possibility. Try to read the length scale.

             ii = bsearchStrings(cgnsLength, refNames)
             if(ii > 0) then

                ! Scale is present; read the value and create the
                ! reference viscosity.

                i = ind(ii)
                call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
                if(ierr /= all_ok)               &
                     call terminate("scaleFactors", &
                     "Something wrong when calling &
                     &cg_array_read_as_f")

                muScale = tmpScale*sqrt(pScale*rhoScale)

             else

                ! Set the logical muScalePresent to .false.

                muScalePresent = .false.

             endif
          endif
       endif

       ! Create the correct scaling factors for density, pressure,
       ! velocity and possibly viscosity.

       rhoScale = rhoScale/rhoRef
       pScale   = pScale/pRef
       velScale = velScale/sqrt(pRef/rhoRef)

       if( muScalePresent ) muScale = muScale/muRef

    endif

    ! Release the memory of refNames, tmpNames and ind.

    deallocate(refNames, tmpNames, ind, stat=ierr)
    if(ierr /= 0)                    &
         call terminate("scaleFactors", &
         "Deallocation error for convNames, etc.")

  end subroutine scaleFactors



  subroutine readRestartVariable(cgnsVarName)

    !      readRestartVariable reads the given variable name from the
    !      cgns restart file.
    !
    use constants
    use blockPointers, only : il, jl, kl
    use utils, only : terminate, setCGNSRealType
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(in) :: cgnsVarName
    !
    !      Local variables.
    !
    integer :: ierr, realTypeCGNS

    integer(kind=intType) :: i, j, k

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Check where the solution variables are stored.

    locationTest: if(location == CellCenter) then

       ! Cell centered values. Read the values directly in the buffer.

       call cg_field_read_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
            cgnsVarName, realTypeCGNS, rangeMin,  &
            rangeMax, buffer, ierr)
       if(ierr /= all_ok)                        &
            call terminate("readRestartVariable", &
            "Something wrong when calling cg_field_read_f")
    else locationTest

       ! Vertex centered values. First read the solution in the
       ! array bufferVertex.

       call cg_field_read_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
            cgnsVarName, realTypeCGNS, rangeMin, &
            rangeMax, bufferVertex, ierr)
       if(ierr /= all_ok)                        &
            call terminate("readRestartVariable", &
            "Something wrong when calling cg_field_read_f")

       ! Create the cell centered values by averaging the vertex values.

       do k=2,kl
          do j=2,jl
             do i=2,il
                buffer(i,j,k) = eighth*(bufferVertex(i-1,j-1,k-1) &
                     +         bufferVertex(i,  j-1,k-1) &
                     +         bufferVertex(i-1,j,  k-1) &
                     +         bufferVertex(i,  j,  k-1) &
                     +         bufferVertex(i-1,j-1,k)   &
                     +         bufferVertex(i,  j-1,k)   &
                     +         bufferVertex(i-1,j,  k)   &
                     +         bufferVertex(i,  j,  k))
             enddo
          enddo
       enddo

    endif locationTest
  end subroutine readRestartVariable

end module variableReading
