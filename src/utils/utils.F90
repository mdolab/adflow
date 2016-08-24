module utils

contains

  function TSbeta(degreePolBeta,   coefPolBeta,       &
       degreeFourBeta,  omegaFourBeta,     &
       cosCoefFourBeta, sinCoefFourBeta, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSbeta computes the angle of attack for a given Time interval *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics, only : equationMode
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSbeta
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolBeta
    integer(kind=intType), intent(in) :: degreeFourBeta

    real(kind=realType), intent(in) :: omegaFourBeta, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolBeta
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourBeta
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourBeta
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: beta, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSBeta = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    beta = coefPolBeta(0)
    do nn=1,degreePolBeta
       beta = beta + coefPolBeta(nn)*(t**nn)
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    beta = beta + cosCoefFourBeta(0)
    do nn=1,degreeFourBeta
       val = nn*omegaFourBeta*t
       beta = beta + cosCoefFourbeta(nn)*cos(val) &
            + sinCoefFourbeta(nn)*sin(val)
    enddo

    ! Set TSBeta to phi.

    TSBeta = beta

  end function TSbeta

  function TSbetadot(degreePolBeta,   coefPolBeta,       &
       degreeFourBeta,  omegaFourBeta,     &
       cosCoefFourBeta, sinCoefFourBeta, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSbeta computes the angle of attack for a given Time interval  *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSbetadot
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolBeta
    integer(kind=intType), intent(in) :: degreeFourBeta

    real(kind=realType), intent(in) :: omegaFourBeta, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolBeta
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourBeta
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourBeta
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: betadot, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSBetadot = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    betadot = zero
    do nn=1,degreePolBeta
       betadot = betadot + nn*coefPolBeta(nn)*(t**(nn-1))
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    do nn=1,degreeFourBeta
       val = nn*omegaFourBeta
       betadot = betadot -val* cosCoefFourbeta(nn)*sin(val*t) &
            +val* sinCoefFourbeta(nn)*cos(val*t)
    enddo

    ! Set TSBeta to phi.

    TSBetadot = betadot

  end function TSbetadot

  function TSMach(degreePolMach,   coefPolMach,       &
       degreeFourMach,  omegaFourMach,     &
       cosCoefFourMach, sinCoefFourMach, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSMach computes the Mach Number for a given time interval      *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSmach
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolMach
    integer(kind=intType), intent(in) :: degreeFourMach

    real(kind=realType), intent(in) :: omegaFourMach, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolMach
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourMach
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourMach
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: intervalMach, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSMach = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    intervalMach = coefPolMach(0)
    do nn=1,degreePolMach
       intervalMach = intervalMach + coefPolMach(nn)*(t**nn)
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    intervalMach = intervalMach + cosCoefFourMach(0)
    do nn=1,degreeFourMach
       val = nn*omegaFourMach*t
       intervalMach = intervalMach + cosCoefFourmach(nn)*cos(val) &
            + sinCoefFourmach(nn)*sin(val)
    enddo
    print *,'inTSMach',intervalMach,nn,val,t
    ! Set TSMach to phi.

    TSMach = intervalMach

  end function TSmach

  function TSMachdot(degreePolMach,   coefPolMach,       &
       degreeFourMach,  omegaFourMach,     &
       cosCoefFourMach, sinCoefFourMach, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSmach computes the angle of attack for a given Time interval *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSmachdot
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolMach
    integer(kind=intType), intent(in) :: degreeFourMach

    real(kind=realType), intent(in) :: omegaFourMach, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolMach
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourMach
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourMach
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: machdot, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSMachdot = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    machdot = zero
    do nn=1,degreePolMach
       machdot = machdot + nn*coefPolMach(nn)*(t**(nn-1))
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    do nn=1,degreeFourMach
       val = nn*omegaFourMach
       machdot = machdot -val* cosCoefFourmach(nn)*sin(val*t) &
            +val* sinCoefFourmach(nn)*cos(val*t)
    enddo

    ! Set TSMach to phi.

    TSMachdot = machdot

  end function TSmachdot

  function TSalpha(degreePolAlpha,   coefPolAlpha,       &
       degreeFourAlpha,  omegaFourAlpha,     &
       cosCoefFourAlpha, sinCoefFourAlpha, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSalpha computes the angle of attack for a given Time interval *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics, only : equationMode
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSalpha
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolAlpha
    integer(kind=intType), intent(in) :: degreeFourAlpha

    real(kind=realType), intent(in) :: omegaFourAlpha, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolAlpha
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourAlpha
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourAlpha
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: alpha, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSAlpha = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.
    alpha = coefPolAlpha(0)
    do nn=1,degreePolAlpha
       alpha = alpha + coefPolAlpha(nn)*(t**nn)
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    alpha = alpha + cosCoefFourAlpha(0)
    do nn=1,degreeFourAlpha
       val = nn*omegaFourAlpha*t
       alpha = alpha + cosCoefFouralpha(nn)*cos(val) &
            + sinCoefFouralpha(nn)*sin(val)
    enddo
    !print *,'inTSalpha',alpha,nn,val,t
    ! Set TSAlpha to phi.

    TSAlpha = alpha

  end function TSalpha

  function TSalphadot(degreePolAlpha,   coefPolAlpha,       &
       degreeFourAlpha,  omegaFourAlpha,     &
       cosCoefFourAlpha, sinCoefFourAlpha, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * TSalpha computes the angle of attack for a given Time interval *
    !      * in a time spectral solution.                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: TSalphadot
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolAlpha
    integer(kind=intType), intent(in) :: degreeFourAlpha

    real(kind=realType), intent(in) :: omegaFourAlpha, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolAlpha
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourAlpha
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourAlpha
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: alphadot, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       TSAlphadot = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    alphadot = zero
    do nn=1,degreePolAlpha
       alphadot = alphadot + nn*coefPolAlpha(nn)*(t**(nn-1))
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    do nn=1,degreeFourAlpha
       val = nn*omegaFourAlpha
       alphadot = alphadot -val* cosCoefFouralpha(nn)*sin(val*t) &
            +val* sinCoefFouralpha(nn)*cos(val*t)
    enddo

    ! Set TSAlpha to phi.

    TSAlphadot = alphadot

  end function TSalphadot


  function derivativeRigidRotAngle(degreePolRot,   &
       coefPolRot,     &
       degreeFourRot,  &
       omegaFourRot,   &
       cosCoefFourRot, &
       sinCoefFourRot, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * derivativeRigidRotAngle computes the time derivative of the    *
    !      * rigid body rotation angle at the given time for the given      *
    !      * arguments. The angle is described by a combination of a        *
    !      * polynomial and fourier series.                                 *
    !      *                                                                *
    !      ******************************************************************
    !
    use flowVarRefState
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: derivativeRigidRotAngle
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolRot
    integer(kind=intType), intent(in) :: degreeFourRot

    real(kind=realType), intent(in) :: omegaFourRot, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourRot
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: dPhi, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       derivativeRigidRotAngle = zero
       return
    endif

    ! Compute the polynomial contribution.

    dPhi = zero
    do nn=1,degreePolRot
       dPhi = dPhi + nn*coefPolRot(nn)*(t**(nn-1))
    enddo

    ! Compute the fourier contribution.

    do nn=1,degreeFourRot
       val = nn*omegaFourRot
       dPhi = dPhi - val*cosCoefFourRot(nn)*sin(val*t)
       dPhi = dPhi + val*sinCoefFourRot(nn)*cos(val*t)
    enddo

    ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
    ! to obtain the correct non-dimensional value.

    derivativeRigidRotAngle = timeRef*dPhi

  end function derivativeRigidRotAngle

  function myDim (x,y)

    use precision

    real(kind=realType) x,y
    real(kind=realType) :: myDim

    myDim = x - y 
    if (myDim < 0.0) then
       myDim = 0.0
    end if

  end function myDim

  function getCorrectForK()

    use flowVarRefState
    use inputPhysics
    use iteration
    implicit none

    logical :: getCorrectForK

    if( kPresent .and. currentLevel <= groundLevel) then
       getCorrectForK = .true.
    else
       getCorrectForK = .false.
    end if
  end function getCorrectForK
  subroutine terminate(routineName, errorMessage)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * terminate writes an error message to standard output and       *
    !      * terminates the execution of the program.                       *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    use communication
    use constants
#ifndef USE_TAPENADE
    use killSignals
#endif
    implicit none
    !
    !      Subroutine arguments
    !
    character(len=*), intent(in) :: routineName
    character(len=*), intent(in) :: errorMessage
#ifndef USE_TAPENADE

    !
    !      Local parameter
    !
    integer, parameter :: maxCharLine = 55
    !
    !      Local variables
    !
    integer :: ierr, len, i2
    logical :: firstTime

    character(len=len_trim(errorMessage)) :: message
    character(len=8) :: integerString
 
   !
    ! Copy the errorMessage into message. It is not possible to work
    ! with errorMessage directly, because it is modified in this
    ! routine. Sometimes a constant string is passed to this routine
    ! and some compilers simply fail then.

    message = errorMessage

    ! Print a nice error message. In case of a parallel executable
    ! also the processor id is printed.

    print "(a)", "#"
    print "(a)", "#--------------------------- !!! Error !!! &
         &----------------------------"

    write(integerString,"(i8)") myID
    integerString = adjustl(integerString)

    print "(2a)", "#* Terminate called by processor ", &
         trim(integerString)

    ! Write the header of the error message.

    print "(2a)", "#* Run-time error in procedure ", &
         trim(routineName)

    ! Loop to write the error message. If the message is too long it
    ! is split over several lines.

    firstTime = .true.
    do
       ! Determine the remaining error message to be written.
       ! If longer than the maximum number of characters allowed
       ! on a line, it is attempted to split the message.

       message = adjustl(message)
       len = len_trim(message)
       i2  = min(maxCharLine,len)

       if(i2 < len) i2 = index(message(:i2), " ", .true.) - 1
       if(i2 < 0)   i2 = index(message, " ") - 1
       if(i2 < 0)   i2 = len

       ! Write this part of the error message. If it is the first
       ! line of the message some additional stuff is printed.

       if( firstTime ) then
          print "(2a)", "#* Error message: ", &
               trim(message(:i2))
          firstTime = .false.
       else
          print "(2a)", "#*                ", &
               trim(message(:i2))
       endif

       ! Exit the loop if the entire message has been written.

       if(i2 == len) exit

       ! Adapt the string for the next part to be written.

       message = message(i2+1:)

    enddo

    ! Write the trailing message.

    print "(a)", "#*"
    print "(a)", "#* Now exiting"
    print "(a)", "#------------------------------------------&
         &----------------------------"
    print "(a)", "#"

    ! Call abort and stop the program. This stop should be done in
    ! abort, but just to be sure.

    call mpi_abort(SUmb_comm_world, 1, ierr)
    stop

#endif

  end subroutine terminate

  subroutine rotMatrixRigidBody(tNew, tOld, rotationMatrix, &
       rotationPoint)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * rotMatrixRigidBody determines the rotation matrix and the      *
    !      * rotation point to determine the coordinates of the new time    *
    !      * level starting from the coordinates of the old time level.     *
    !      *                                                                *
    !      ******************************************************************
    !
    use flowVarRefState
    use inputMotion
    use monitor
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), intent(in) :: tNew, tOld

    real(kind=realType), dimension(3),   intent(out) :: rotationPoint
    real(kind=realType), dimension(3,3), intent(out) :: rotationMatrix
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j

    real(kind=realType) :: phi
    real(kind=realType) :: cosX, cosY, cosZ, sinX, sinY, sinZ

    real(kind=realType), dimension(3,3) :: mNew, mOld

    ! Determine the rotation angle around the x-axis for the new
    ! time level and the corresponding values of the sine and cosine.

    phi = rigidRotAngle(degreePolXRot,   coefPolXRot,     &
         degreeFourXRot,  omegaFourXRot,   &
         cosCoefFourXRot, sinCoefFourXRot, tNew)
    sinX = sin(phi)
    cosX = cos(phi)

    ! Idem for the y-axis.

    phi = rigidRotAngle(degreePolYRot,   coefPolYRot,     &
         degreeFourYRot,  omegaFourYRot,   &
         cosCoefFourYRot, sinCoefFourYRot, tNew)
    sinY = sin(phi)
    cosY = cos(phi)

    ! Idem for the z-axis.

    phi = rigidRotAngle(degreePolZRot,   coefPolZRot,     &
         degreeFourZRot,  omegaFourZRot,   &
         cosCoefFourZRot, sinCoefFourZRot, tNew)
    sinZ = sin(phi)
    cosZ = cos(phi)

    ! Construct the transformation matrix at the new time level.
    ! It is assumed that the sequence of rotation is first around the
    ! x-axis then around the y-axis and finally around the z-axis.

    mNew(1,1) =  cosY*cosZ
    mNew(2,1) =  cosY*sinZ
    mNew(3,1) = -sinY

    mNew(1,2) = sinX*sinY*cosZ - cosX*sinZ
    mNew(2,2) = sinX*sinY*sinZ + cosX*cosZ
    mNew(3,2) = sinX*cosY

    mNew(1,3) = cosX*sinY*cosZ + sinX*sinZ
    mNew(2,3) = cosX*sinY*sinZ - sinX*cosZ
    mNew(3,3) = cosX*cosY

    ! Determine the rotation angle around the x-axis for the old
    ! time level and the corresponding values of the sine and cosine.

    phi = rigidRotAngle(degreePolXRot,   coefPolXRot,     &
         degreeFourXRot,  omegaFourXRot,   &
         cosCoefFourXRot, sinCoefFourXRot, tOld)
    sinX = sin(phi)
    cosX = cos(phi)

    ! Idem for the y-axis.

    phi = rigidRotAngle(degreePolYRot,   coefPolYRot,     &
         degreeFourYRot,  omegaFourYRot,   &
         cosCoefFourYRot, sinCoefFourYRot, tOld)
    sinY = sin(phi)
    cosY = cos(phi)

    ! Idem for the z-axis.

    phi = rigidRotAngle(degreePolZRot,   coefPolZRot,     &
         degreeFourZRot,  omegaFourZRot,   &
         cosCoefFourZRot, sinCoefFourZRot, tOld)
    sinZ = sin(phi)
    cosZ = cos(phi)

    ! Construct the transformation matrix at the old time level.

    mOld(1,1) =  cosY*cosZ
    mOld(2,1) =  cosY*sinZ
    mOld(3,1) = -sinY

    mOld(1,2) = sinX*sinY*cosZ - cosX*sinZ
    mOld(2,2) = sinX*sinY*sinZ + cosX*cosZ
    mOld(3,2) = sinX*cosY

    mOld(1,3) = cosX*sinY*cosZ + sinX*sinZ
    mOld(2,3) = cosX*sinY*sinZ - sinX*cosZ
    mOld(3,3) = cosX*cosY

    ! Construct the transformation matrix between the new and the
    ! old time level. This is mNew*inverse(mOld). However the
    ! inverse of mOld is the transpose.

    do j=1,3
       do i=1,3
          rotationMatrix(i,j) = mNew(i,1)*mOld(j,1) &
               + mNew(i,2)*mOld(j,2) &
               + mNew(i,3)*mOld(j,3)
       enddo
    enddo

    ! Determine the rotation point at the old time level; it is
    ! possible that this value changes due to translation of the grid.

    !  aInf = sqrt(gammaInf*pInf/rhoInf)

    !  rotationPoint(1) = LRef*rotPoint(1) &
    !                   + MachGrid(1)*aInf*tOld/timeRef
    !  rotationPoint(2) = LRef*rotPoint(2) &
    !                   + MachGrid(2)*aInf*tOld/timeRef
    !  rotationPoint(3) = LRef*rotPoint(3) &
    !                   + MachGrid(3)*aInf*tOld/timeRef

    rotationPoint(1) = LRef*rotPoint(1)
    rotationPoint(2) = LRef*rotPoint(2)
    rotationPoint(3) = LRef*rotPoint(3)

  end subroutine rotMatrixRigidBody

  function secondDerivativeRigidRotAngle(degreePolRot,   &
       coefPolRot,     &
       degreeFourRot,  &
       omegaFourRot,   &
       cosCoefFourRot, &
       sinCoefFourRot, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * 2ndderivativeRigidRotAngle computes the 2nd time derivative of *
    !      * the rigid body rotation angle at the given time for the given  *
    !      * arguments. The angle is described by a combination of a        *
    !      * polynomial and fourier series.                                 *
    !      *                                                                *
    !      ******************************************************************
    !
    use flowVarRefState
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: secondDerivativeRigidRotAngle
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolRot
    integer(kind=intType), intent(in) :: degreeFourRot

    real(kind=realType), intent(in) :: omegaFourRot, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourRot
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: dPhi, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       secondDerivativeRigidRotAngle = zero
       return
    endif

    ! Compute the polynomial contribution.

    dPhi = zero
    do nn=2,degreePolRot
       dPhi = dPhi + (nn-1)*nn*coefPolRot(nn)*(t**(nn-2))
    enddo

    ! Compute the fourier contribution.

    do nn=1,degreeFourRot
       val = nn*omegaFourRot
       dPhi = dPhi - val**2*sinCoefFourRot(nn)*sin(val*t)
       dPhi = dPhi - val**2*cosCoefFourRot(nn)*cos(val*t)
    enddo

    ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
    ! to obtain the correct non-dimensional value.

    secondDerivativeRigidRotAngle = timeRef**2*dPhi

  end function secondDerivativeRigidRotAngle

  function rigidRotAngle(degreePolRot,   coefPolRot,       &
       degreeFourRot,  omegaFourRot,     &
       cosCoefFourRot, sinCoefFourRot, t)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * rigidRotAngle computes the rigid body rotation angle at the    *
    !      * given time for the given arguments. The angle is described by  *
    !      * a combination of a polynomial and fourier series.              *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputPhysics
    implicit none
    !
    !      Function type
    !
    real(kind=realType) :: rigidRotAngle
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in) :: degreePolRot
    integer(kind=intType), intent(in) :: degreeFourRot

    real(kind=realType), intent(in) :: omegaFourRot, t

    real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
    real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
    real(kind=realType), dimension(*),   intent(in) :: sinCoefFourRot
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    real(kind=realType) :: phi, val

    ! Return immediately if this is a steady computation.

    if(equationMode == steady) then
       rigidRotAngle = zero
       return
    endif

    ! Compute the polynomial contribution. If no polynomial was
    ! specified, the value of index 0 is set to zero automatically.

    phi = coefPolRot(0)
    do nn=1,degreePolRot
       phi = phi + coefPolRot(nn)*(t**nn)
    enddo

    ! Compute the fourier contribution. Again the cosine coefficient
    ! of index 0 is defaulted to zero if not specified.

    phi = phi + cosCoefFourRot(0)
    do nn=1,degreeFourRot
       val = nn*omegaFourRot*t
       phi = phi + cosCoefFourRot(nn)*cos(val) &
            + sinCoefFourRot(nn)*sin(val)
    enddo

    ! Set rigidRotAngle to phi.

    rigidRotAngle = phi

  end function rigidRotAngle

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef  USE_TAPENADE

  function delta(val1,val2)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * delta is a function used to determine the contents of the full *
    !      * transformation matrix from the shorthand form. It returns 1    *
    !      * if the absolute value of the two arguments are identical.      *
    !      * Otherwise it returns 0.                                        *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Function type.
    !
    integer(kind=intType) :: delta
    !
    !      Function arguments.
    !
    integer(kind=intType) :: val1, val2

    if(abs(val1) == abs(val2)) then
       delta = 1_intType
    else
       delta = 0_intType
    endif

  end function delta


  logical function myIsNAN(val)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * myIsNAN determines whether or not the given value is a NAN and *
    !      * returns the according logical.                                 *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Function arguments.
    !
    real(kind=realType), intent(in) :: val
    !
    !      Local variable.
    !
    integer(kind=intType) :: res

    call myIsNaNC(val, res)
    if(res == 1) then
       myIsNAN = .true.
    else
       myIsNAN = .false.
    endif

  end function myIsNAN
  !
  subroutine nullifyCGNSDomPointers(nn)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * nullifyCGNSDomPointers nullifies all the pointers of the       *
    !      * given CGNS block.                                              *
    !      *                                                                *
    !      ******************************************************************
    !
    use cgnsGrid
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    nullify(cgnsDoms(nn)%procStored)
    nullify(cgnsDoms(nn)%conn1to1)
    nullify(cgnsDoms(nn)%connNonMatchAbutting)
    nullify(cgnsDoms(nn)%bocoInfo)

  end subroutine nullifyCGNSDomPointers

  subroutine nullifyFlowDomPointers(nn,level,sps)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * nullifyFlowDomPointers nullifies all the pointers of the       *
    !      * given block.                                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn, level, sps

    nullify(flowDoms(nn,level,sps)%BCType)
    nullify(flowDoms(nn,level,sps)%BCFaceID)
    nullify(flowDoms(nn,level,sps)%cgnsSubface)

    nullify(flowDoms(nn,level,sps)%inBeg)
    nullify(flowDoms(nn,level,sps)%jnBeg)
    nullify(flowDoms(nn,level,sps)%knBeg)
    nullify(flowDoms(nn,level,sps)%inEnd)
    nullify(flowDoms(nn,level,sps)%jnEnd)
    nullify(flowDoms(nn,level,sps)%knEnd)

    nullify(flowDoms(nn,level,sps)%dinBeg)
    nullify(flowDoms(nn,level,sps)%djnBeg)
    nullify(flowDoms(nn,level,sps)%dknBeg)
    nullify(flowDoms(nn,level,sps)%dinEnd)
    nullify(flowDoms(nn,level,sps)%djnEnd)
    nullify(flowDoms(nn,level,sps)%dknEnd)

    nullify(flowDoms(nn,level,sps)%icBeg)
    nullify(flowDoms(nn,level,sps)%jcBeg)
    nullify(flowDoms(nn,level,sps)%kcBeg)
    nullify(flowDoms(nn,level,sps)%icEnd)
    nullify(flowDoms(nn,level,sps)%jcEnd)
    nullify(flowDoms(nn,level,sps)%kcEnd)

    nullify(flowDoms(nn,level,sps)%neighBlock)
    nullify(flowDoms(nn,level,sps)%neighProc)
    nullify(flowDoms(nn,level,sps)%l1)
    nullify(flowDoms(nn,level,sps)%l2)
    nullify(flowDoms(nn,level,sps)%l3)
    nullify(flowDoms(nn,level,sps)%groupNum)

    nullify(flowDoms(nn,level,sps)%iblank)
    nullify(flowDoms(nn,level,sps)%fringes)
    nullify(flowDoms(nn,level,sps)%orphans)

    nullify(flowDoms(nn,level,sps)%BCData)
    nullify(flowDoms(nn,level,sps)%viscSubface)

    nullify(flowDoms(nn,level,sps)%viscIminPointer)
    nullify(flowDoms(nn,level,sps)%viscImaxPointer)
    nullify(flowDoms(nn,level,sps)%viscJminPointer)
    nullify(flowDoms(nn,level,sps)%viscJmaxPointer)
    nullify(flowDoms(nn,level,sps)%viscKminPointer)
    nullify(flowDoms(nn,level,sps)%viscKmaxPointer)

    nullify(flowDoms(nn,level,sps)%x)
    nullify(flowDoms(nn,level,sps)%xOld)
    nullify(flowDoms(nn,level,sps)%si)
    nullify(flowDoms(nn,level,sps)%sj)
    nullify(flowDoms(nn,level,sps)%sk)
    nullify(flowDoms(nn,level,sps)%vol)
    nullify(flowDoms(nn,level,sps)%volRef)
    nullify(flowDoms(nn,level,sps)%volOld)

    nullify(flowDoms(nn,level,sps)%pori)
    nullify(flowDoms(nn,level,sps)%porj)
    nullify(flowDoms(nn,level,sps)%pork)

    nullify(flowDoms(nn,level,sps)%indFamilyI)
    nullify(flowDoms(nn,level,sps)%indFamilyJ)
    nullify(flowDoms(nn,level,sps)%indFamilyK)

    nullify(flowDoms(nn,level,sps)%factFamilyI)
    nullify(flowDoms(nn,level,sps)%factFamilyJ)
    nullify(flowDoms(nn,level,sps)%factFamilyK)

    nullify(flowDoms(nn,level,sps)%rotMatrixI)
    nullify(flowDoms(nn,level,sps)%rotMatrixJ)
    nullify(flowDoms(nn,level,sps)%rotMatrixK)

    nullify(flowDoms(nn,level,sps)%sFaceI)
    nullify(flowDoms(nn,level,sps)%sFaceJ)
    nullify(flowDoms(nn,level,sps)%sFaceK)

    nullify(flowDoms(nn,level,sps)%w)
    nullify(flowDoms(nn,level,sps)%wOld)
    nullify(flowDoms(nn,level,sps)%p)
    nullify(flowDoms(nn,level,sps)%aa)
    nullify(flowDoms(nn,level,sps)%gamma)
    nullify(flowDoms(nn,level,sps)%rlv)
    nullify(flowDoms(nn,level,sps)%rev)
    nullify(flowDoms(nn,level,sps)%s)

    nullify(flowDoms(nn,level,sps)%ux)
    nullify(flowDoms(nn,level,sps)%uy)
    nullify(flowDoms(nn,level,sps)%uz)

    nullify(flowDoms(nn,level,sps)%vx)
    nullify(flowDoms(nn,level,sps)%vy)
    nullify(flowDoms(nn,level,sps)%vz)

    nullify(flowDoms(nn,level,sps)%wx)
    nullify(flowDoms(nn,level,sps)%wy)
    nullify(flowDoms(nn,level,sps)%wz)

    nullify(flowDoms(nn,level,sps)%qx)
    nullify(flowDoms(nn,level,sps)%qy)
    nullify(flowDoms(nn,level,sps)%qz)

    nullify(flowDoms(nn,level,sps)%dw)
    nullify(flowDoms(nn,level,sps)%fw)
    nullify(flowDoms(nn,level,sps)%scratch)

    nullify(flowDoms(nn,level,sps)%dwOldRK)

    nullify(flowDoms(nn,level,sps)%p1)
    nullify(flowDoms(nn,level,sps)%w1)
    nullify(flowDoms(nn,level,sps)%wr)

    nullify(flowDoms(nn,level,sps)%mgIFine)
    nullify(flowDoms(nn,level,sps)%mgJFine)
    nullify(flowDoms(nn,level,sps)%mgKFine)

    nullify(flowDoms(nn,level,sps)%mgIWeight)
    nullify(flowDoms(nn,level,sps)%mgJWeight)
    nullify(flowDoms(nn,level,sps)%mgKWeight)

    nullify(flowDoms(nn,level,sps)%mgICoarse)
    nullify(flowDoms(nn,level,sps)%mgJCoarse)
    nullify(flowDoms(nn,level,sps)%mgKCoarse)

    nullify(flowDoms(nn,level,sps)%ico)
    nullify(flowDoms(nn,level,sps)%jco)
    nullify(flowDoms(nn,level,sps)%kco)

    nullify(flowDoms(nn,level,sps)%wn)
    nullify(flowDoms(nn,level,sps)%pn)
    nullify(flowDoms(nn,level,sps)%dtl)
    nullify(flowDoms(nn,level,sps)%radI)
    nullify(flowDoms(nn,level,sps)%radJ)
    nullify(flowDoms(nn,level,sps)%radK)

    nullify(flowDoms(nn,level,sps)%d2Wall)

    nullify(flowDoms(nn,level,sps)%bmti1)
    nullify(flowDoms(nn,level,sps)%bmti2)
    nullify(flowDoms(nn,level,sps)%bmtj1)
    nullify(flowDoms(nn,level,sps)%bmtj2)
    nullify(flowDoms(nn,level,sps)%bmtk1)
    nullify(flowDoms(nn,level,sps)%bmtk2)

    nullify(flowDoms(nn,level,sps)%bvti1)
    nullify(flowDoms(nn,level,sps)%bvti2)
    nullify(flowDoms(nn,level,sps)%bvtj1)
    nullify(flowDoms(nn,level,sps)%bvtj2)
    nullify(flowDoms(nn,level,sps)%bvtk1)
    nullify(flowDoms(nn,level,sps)%bvtk2)

    nullify(flowDoms(nn,level,sps)%globalCell)
    nullify(flowDoms(nn,level,sps)%globalNode)
    nullify(flowDOms(nn,level,sps)%surfNodeIndices)
    nullify(flowDOms(nn,level,sps)%uv)
    nullify(flowDoms(nn,level,sps)%wallInd)
    nullify(flowDoms(nn,level,sps)%xSeed)

    ! *******************************
    ! Added by HDN
    ! *******************************
    nullify(flowDoms(nn,level,sps)%xALE)
    nullify(flowDoms(nn,level,sps)%sIALE)
    nullify(flowDoms(nn,level,sps)%sJALE)
    nullify(flowDoms(nn,level,sps)%sKALE)
    nullify(flowDoms(nn,level,sps)%sFaceIALE)
    nullify(flowDoms(nn,level,sps)%sFaceJALE)
    nullify(flowDoms(nn,level,sps)%sFaceKALE)
    nullify(flowDoms(nn,level,sps)%dwALE)
    nullify(flowDoms(nn,level,sps)%fwALE)
#ifndef USE_TAPENADE
    nullify(flowDoms(nn,level,sps)%PCMat)
    nullify(flowDoms(nn,level,sps)%i_D_Fact)
    nullify(flowDoms(nn,level,sps)%i_L_Fact)
    nullify(flowDoms(nn,level,sps)%i_U_Fact)
    nullify(flowDoms(nn,level,sps)%i_U2_Fact)

    nullify(flowDoms(nn,level,sps)%j_D_Fact)
    nullify(flowDoms(nn,level,sps)%j_L_Fact)
    nullify(flowDoms(nn,level,sps)%j_U_Fact)
    nullify(flowDoms(nn,level,sps)%j_U2_Fact)

    nullify(flowDoms(nn,level,sps)%k_D_Fact)
    nullify(flowDoms(nn,level,sps)%k_L_Fact)
    nullify(flowDoms(nn,level,sps)%k_U_Fact)
    nullify(flowDoms(nn,level,sps)%k_U2_Fact)
#endif

  end subroutine nullifyFlowDomPointers

  subroutine reallocateInteger(intArray, newSize, oldSize, &
       alwaysFreeMem)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * reallocateInteger reallocates the given integer array to the   *
    !      * given new size. The old values of the array are copied. Note   *
    !      * that newSize can be both smaller and larger than oldSize.      *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), dimension(:), pointer :: intArray
    integer(kind=intType), intent(in) :: newSize, oldSize
    logical, intent(in) :: alwaysFreeMem
    !
    !      Local variables.
    !
    integer(kind=intType), dimension(:), pointer :: tmp

    integer(kind=intType) :: i, nn, ll

    integer :: ierr

    ! Determine the minimum of newSize and oldSize.

    nn = min(newSize, oldSize)

    ! Set the pointer for tmp to intArray.

    tmp => intArray

    ! Allocate the memory for intArray in case newSize is larger
    ! than 0 or if alwaysFreeMem is .true. And copy the old data
    ! into it. Preserve the lower bound.

    if(newSize > 0 .or. alwaysFreeMem) then

       ll = 1
       if (associated(intArray)) ll = lbound(intArray,1)

       allocate(intArray(ll:newSize+ll-1), stat=ierr)
       if(ierr /= 0)                         &
            call returnFail("reallocateInteger", &
            "Memory allocation failure for intArray")
       do i=ll,ll+nn-1
          intArray(i) = tmp(i)
       enddo
    endif

    ! Release the memory for tmp in case oldSize is larger than 0 or
    ! if alwaysFreeMem is .true.

    if(oldSize > 0 .or. alwaysFreeMem) then
       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                         &
            call returnFail("reallocateInteger", &
            "Deallocation error for tmp")
    endif

  end subroutine reallocateInteger

  !================================================================

  subroutine reallocateMpiOffsetKindInteger(intArray, newSize, &
       oldSize, alwaysFreeMem)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * reallocateMpiOffsetKindInteger reallocates the given           *
    !      * mpi_offset_kind integer array to the given new size. The old   *
    !      * values of the array are copied. Note that newSize can be both  *
    !      * smaller and larger than oldSize.                               *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=mpi_offset_kind), dimension(:), pointer :: intArray
    integer(kind=intType), intent(in) :: newSize, oldSize
    logical, intent(in) :: alwaysFreeMem
    !
    !      Local variables.
    !
    integer(kind=mpi_offset_kind), dimension(:), pointer :: tmp

    integer(kind=intType) :: i, nn, ll

    integer :: ierr

    ! Determine the minimum of newSize and oldSize.

    nn = min(newSize, oldSize)

    ! Set the pointer for tmp to intArray.

    tmp => intArray

    ! Allocate the memory for intArray in case newSize is larger
    ! than 0 or if alwaysFreeMem is .true. And copy the old data
    ! into it. Preserve the lower bound.

    if(newSize > 0 .or. alwaysFreeMem) then

       ll = 1
       if (associated(intArray)) ll = lbound(intArray,1)

       allocate(intArray(ll:newSize+ll-1), stat=ierr)
       if(ierr /= 0)                                      &
            call returnFail("reallocateMpiOffsetKindInteger", &
            "Memory allocation failure for intArray")
       do i=ll,ll+nn-1
          intArray(i) = tmp(i)
       enddo
    endif

    ! Release the memory for tmp in case oldSize is larger than 0 or
    ! if alwaysFreeMem is .true.

    if(oldSize > 0 .or. alwaysFreeMem) then
       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                                      &
            call returnFail("reallocateMpiOffsetKindInteger", &
            "Deallocation error for tmp")
    endif

  end subroutine reallocateMpiOffsetKindInteger

  !================================================================

  subroutine reallocateInteger2(intArray, newSize1, newSize2, &
       oldSize1, oldSize2,           &
       alwaysFreeMem)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * reallocateInteger2 reallocates the given 2D integer array to   *
    !      * the given new sizes. The old values of the array are copied.   *
    !      * Note that the newSizes can be both smaller and larger than     *
    !      * the oldSizes.                                                  *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), dimension(:,:), pointer :: intArray
    integer(kind=intType), intent(in) :: newSize1, newSize2, &
         oldSize1, oldSize2
    logical, intent(in) :: alwaysFreeMem
    !
    !      Local variables.
    !
    integer(kind=intType), dimension(:,:), pointer :: tmp

    integer(kind=intType) :: newSize, oldSize
    integer(kind=intType) :: nn1, nn2, nn

    integer(kind=intType) :: i, j

    integer :: ierr

    ! Determine the total new and old size.

    newSize = newSize1*newSize2
    oldSize = oldSize1*oldSize2

    ! Determine for each of the 2 components the minimum of the new
    ! and the old size. Multiply these values to obtain the total
    ! amount of data that must be copied.

    nn1 = min(newSize1, oldSize1)
    nn2 = min(newSize2, oldSize2)

    nn = nn1*nn2

    ! Set the pointer for tmp.

    tmp => intArray

    ! Allocate the memory for intArray in case newSize is larger
    ! than 0 or if alwaysFreeMem is .true. and copy the old data
    ! into it.

    if(newSize > 0 .or. alwaysFreeMem) then
       allocate(intArray(newSize1,newSize2), stat=ierr)
       if(ierr /= 0)                          &
            call returnFail("reallocateInteger2", &
            "Memory allocation failure for intArray")
       do j=1,nn2
          do i=1,nn1
             intArray(i,j) = tmp(i,j)
          enddo
       enddo
    endif

    ! Release the memory of tmp in case oldSize is larger than 0
    ! or if alwaysFreeMem is .true..

    if(oldSize > 0 .or. alwaysFreeMem) then
       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                          &
            call returnFail("reallocateInteger2", &
            "Deallocation error for tmp")
    endif

  end subroutine reallocateInteger2

  subroutine reallocateReal(realArray, newSize, oldSize, &
       alwaysFreeMem)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * ReallocateReal reallocates the given real array to the given   *
    !      * new size. The old values of the array are copied. Note that    *
    !      * newSize can be both smaller and larger than oldSize.           *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), dimension(:), pointer :: realArray
    integer(kind=intType), intent(in) :: newSize, oldSize
    logical, intent(in) :: alwaysFreeMem
    !
    !      Local variables.
    !
    real(kind=realType), dimension(:), pointer :: tmp

    integer(kind=intType) :: i, nn

    integer :: ierr

    ! Determine the minimum of newSize and oldSize.

    nn = min(newSize, oldSize)

    ! Set the pointer for tmp to realArray.

    tmp => realArray

    ! Allocate the memory for realArray in case newSize is larger
    ! than 0 or if alwaysFreeMem is .True. And copy the old data
    ! into it.

    if(newSize > 0 .or. alwaysFreeMem) then
       allocate(realArray(newSize), stat=ierr)
       if(ierr /= 0)                       &
            call returnFail("reallocateReal", &
            "Memory allocation failure for realArray")
       do i=1,nn
          realArray(i) = tmp(i)
       enddo
    endif

    ! Release the memory for tmp in case oldSize is larger than 0 or
    ! if alwaysFreeMem is .True.

    if(oldSize > 0 .or. alwaysFreeMem) then
       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                       &
            call returnFail("reallocateReal", &
            "Deallocation error for tmp")
    endif

  end subroutine reallocateReal

  !================================================================

  subroutine reallocateReal2(realArray, newSize1, newSize2, &
       oldSize1, oldSize2,            &
       alwaysFreeMem)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * ReallocateReal2 reallocates the given 2d integer array to      *
    !      * the given new sizes. The old values of the array are copied.   *
    !      * Note that the newSizes can be both smaller and larger than     *
    !      * the oldSizes.                                                  *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), dimension(:,:), pointer :: realArray
    integer(kind=intType), intent(in) :: newSize1, newSize2, &
         oldSize1, oldSize2
    logical, intent(in) :: alwaysFreeMem
    !
    !      Local variables.
    !
    real(kind=realType), dimension(:,:), pointer :: tmp

    integer(kind=intType) :: newSize, oldSize
    integer(kind=intType) :: nn1, nn2, nn

    integer(kind=intType) :: i, j

    integer :: ierr

    ! Determine the total new and old size.

    newSize = newSize1*newSize2
    oldSize = oldSize1*oldSize2

    ! Determine for each of the 2 components the minimum of the new
    ! and the old size. Multiply these values to obtain the total
    ! amount of data that must be copied.

    nn1 = min(newSize1, oldSize1)
    nn2 = min(newSize2, oldSize2)

    nn = nn1*nn2

    ! Set the pointer for tmp.

    tmp => realArray

    ! Allocate the memory for realArray in case newSize is larger
    ! than 0 or if alwaysFreeMem is .True. And copy the old data
    ! into it.

    if(newSize > 0 .or. alwaysFreeMem) then
       allocate(realArray(newSize1,newSize2), stat=ierr)
       if(ierr /= 0)                           &
            call returnFail("reallocateReal2", &
            "Memory allocation failure for realArray")
       do j=1,nn2
          do i=1,nn1
             realArray(i,j) = tmp(i,j)
          enddo
       enddo
    endif

    ! Release the memory of tmp in case oldSize is larger than 0
    ! or if alwaysFreeMem is .True..

    if(oldSize > 0 .or. alwaysFreeMem) then
       deallocate(tmp, stat=ierr)
       if(ierr /= 0)                           &
            call returnFail("reallocateReal2", &
            "Deallocation error for tmp")
    endif

  end subroutine reallocateReal2

  subroutine setBufferSizes(level, sps, determine1to1Buf,    &
       determineSlidingBuf, & 
       determineOversetBuf)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setBufferSizes determines the size of the send and receive     *
    !      * buffers for this grid level. After that the maximum value of   *
    !      * these sizes and the currently stored value is taken, such that *
    !      * for all mg levels the same buffer can be used. Normally the    *
    !      * size on the finest grid should be enough, but it is just as    *
    !      * safe to check on all mg levels. A distinction is made between  *
    !      * 1 to 1 and sliding mesh communication, because these happen    *
    !      * consecutively and not simultaneously. Consequently the actual  *
    !      * buffer size is the maximum of the two and not the sum.         *
    !      * For steady state computations a mixing plane boundary          *
    !      * condition is used instead of a sliding mesh. However, both     *
    !      * communication patterns are allocated and initialized.          *
    !      * Therefore the maximum of the two can be taken without checking *
    !      * the situation we are dealing with.                             *
    !      *                                                                *
    !      ******************************************************************
    !
    use commMixing
    use commSliding
    use communication
    use flowVarRefState
    use inputPhysics
    use interfaceGroups
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    logical, intent(in) :: determine1to1Buf, determineSlidingBuf
    logical, intent(in) :: determineOversetBuf
    !
    !      Local variables.
    !
    integer(kind=intType) :: i
    integer(kind=intType) :: sendSize, recvSize, nVarComm

    ! Determine the maximum number of variables to be communicated.

    nVarComm = nw + 1
    if(cpModel == cpTempCurveFits) nVarComm = nVarComm + 1
    if( viscous )   nVarComm = nVarComm + 1
    if( eddyModel ) nVarComm = nVarComm + 1

    ! Check if the 1 to 1 communication must be considered.

    if( determine1to1Buf ) then

       ! Store the send and receive buffer sizes needed for the nodal
       ! exchange. Determine the maximum for the number of send and
       ! receive processors.

       i = commPatternNode_1st(level)%nProcSend
       sendSize = commPatternNode_1st(level)%nsendCum(i)

       i = commPatternNode_1st(level)%nProcRecv
       recvSize = commPatternNode_1st(level)%nrecvCum(i)

       ! Determine the buffer sizes for the 2nd level cell exchange and
       ! set the size for this processor to the maximum needed. Note
       ! that it is not needed to test the 1st level cell halo, because
       ! it is entirely incorporated in the 2nd level.
       ! Determine the maximum for the number of send and receive
       ! processors as well.

       i = commPatternCell_2nd(level)%nProcSend
       sendSize = max(sendSize, &
            commPatternCell_2nd(level)%nsendCum(i))

       i = commPatternCell_2nd(level)%nProcRecv
       recvSize = max(recvSize, &
            commPatternCell_2nd(level)%nrecvCum(i))

       ! Multiply sendSize and recvSize with the number of variables to
       ! be communicated.

       sendSize = sendSize*nVarComm
       recvSize = recvSize*nVarComm

       ! Store the maximum of the current values and the old values
       ! in sendBufferSize1to1 and recvBufferSize1to1.

       sendBufferSize_1to1 = max(sendBufferSize_1to1, sendSize)
       recvBufferSize_1to1 = max(recvBufferSize_1to1, recvSize)

    endif

    ! Check if the sliding mesh communication must be considered.

    if( determineSlidingBuf ) then

       ! Only the second level cell halo communication pattern needs
       ! to be considered. Note that there is no nodal communication
       ! pattern for sliding meshes.

       i = commSlidingCell_2nd(level,sps)%nProcSend
       sendSize = commSlidingCell_2nd(level,sps)%nsendCum(i)

       i = commSlidingCell_2nd(level,sps)%nProcRecv
       recvSize = commSlidingCell_2nd(level,sps)%nrecvCum(i)

       ! Multiply sendSize and recvSize with the number of variables to
       ! be communicated.

       sendSize = sendSize*nVarComm
       recvSize = recvSize*nVarComm

       ! Store the maximum of the current values and the old values
       ! in sendBufferSizeSlide and recvBufferSizeSlide.

       sendBufferSizeSlide = max(sendBufferSizeSlide, sendSize)
       recvBufferSizeSlide = max(recvBufferSizeSlide, recvSize)

       ! Take possible mixing plane boundaries into account.

       sendSize = 0
       do i=1,nInterfaceGroups
          sendSize = max(sendSize,                            &
               commPatternMixing(level,i,1)%nInter, &
               commPatternMixing(level,i,2)%nInter)
       enddo

       sendSize = sendSize*nVarComm
       recvSize = sendSize

       ! Store the maximum value in sendBufferSizeSlide and
       ! recvBufferSizeSlide.

       sendBufferSizeSlide = max(sendBufferSizeSlide, sendSize)
       recvBufferSizeSlide = max(recvBufferSizeSlide, recvSize)

    endif

    ! Check if the overset communication must be considered.

    if( determineOversetBuf ) then

       ! Same deal for the overset communication.

       i = commPatternOverset(level,sps)%nProcSend
       sendSize = commPatternOverset(level,sps)%nsendCum(i)

       i = commPatternOverset(level,sps)%nProcRecv
       recvSize = commPatternOverset(level,sps)%nrecvCum(i)

       ! Multiply sendSize and recvSize with the number of variables to
       ! be communicated.

       sendSize = sendSize*nVarComm
       recvSize = recvSize*nVarComm

       ! Store the maximum of the current values and the old values.

       sendBufferSizeOver = max(sendBufferSizeOver, sendSize)
       recvBufferSizeOver = max(recvBufferSizeOver, recvSize)

    endif

    ! Take the maximum for of all the buffers to
    ! obtain the actual size to be allocated.

    sendBufferSize = max(sendBufferSize_1to1, &
         sendBufferSizeOver,  &
         sendBufferSizeSlide)
    recvBufferSize = max(recvBufferSize_1to1, &
         recvBufferSizeOver,  &
         recvBufferSizeSlide)

  end subroutine setBufferSizes

  subroutine setPointers(nn,mm,ll)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setPointers makes the variables in blockPointers point to      *
    !      * block nn for grid level mm and spectral solution ll.           *
    !      *                                                                *
    !      ******************************************************************
    !
    use blockPointers
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: nn, mm, ll

    ! Store the info of the current block, such that inside the
    ! module blockPointers it is known to which block the data
    ! belongs.

    sectionID   = flowDoms(nn,mm,ll)%sectionID
    nbkLocal    = nn
    nbkGlobal   = flowDoms(nn,mm,ll)%cgnsBlockID
    mgLevel     = mm
    spectralSol = ll

    ! Block dimensions.

    nx = flowDoms(nn,mm,ll)%nx
    ny = flowDoms(nn,mm,ll)%ny
    nz = flowDoms(nn,mm,ll)%nz

    il = flowDoms(nn,mm,ll)%il
    jl = flowDoms(nn,mm,ll)%jl
    kl = flowDoms(nn,mm,ll)%kl

    ie = flowDoms(nn,mm,ll)%ie
    je = flowDoms(nn,mm,ll)%je
    ke = flowDoms(nn,mm,ll)%ke

    ib = flowDoms(nn,mm,ll)%ib
    jb = flowDoms(nn,mm,ll)%jb
    kb = flowDoms(nn,mm,ll)%kb

    imaxDim = max(ie,je)
    jmaxDim = max(je,ke)

    rightHanded = flowDoms(nn,mm,ll)%righthanded

    ! Point range in the corresponding cgns block

    iBegor = flowDoms(nn,mm,ll)%iBegor
    iEndor = flowDoms(nn,mm,ll)%iEndor
    jBegor = flowDoms(nn,mm,ll)%jBegor
    jEndor = flowDoms(nn,mm,ll)%jEndor
    kBegor = flowDoms(nn,mm,ll)%kBegor
    kEndor = flowDoms(nn,mm,ll)%kEndor

    ! Subface info. Note that the pointers point to the 1st spectral
    ! mode, because this is the only one allocated. The info is the
    ! same for all modes.

    nSubface   = flowDoms(nn,mm,ll)%nSubface
    n1to1      = flowDoms(nn,mm,ll)%n1to1
    nBocos     = flowDoms(nn,mm,ll)%nBocos
    nViscBocos = flowDoms(nn,mm,ll)%nViscBocos

    BCType      => flowDoms(nn,mm,1)%BCType
    BCFaceID    => flowDoms(nn,mm,1)%BCFaceID
    cgnsSubface => flowDoms(nn,mm,1)%cgnsSubface

    inBeg => flowDoms(nn,mm,1)%inBeg
    jnBeg => flowDoms(nn,mm,1)%jnBeg
    knBeg => flowDoms(nn,mm,1)%knBeg
    inEnd => flowDoms(nn,mm,1)%inEnd
    jnEnd => flowDoms(nn,mm,1)%jnEnd
    knEnd => flowDoms(nn,mm,1)%knEnd

    dinBeg => flowDoms(nn,mm,1)%dinBeg
    djnBeg => flowDoms(nn,mm,1)%djnBeg
    dknBeg => flowDoms(nn,mm,1)%dknBeg
    dinEnd => flowDoms(nn,mm,1)%dinEnd
    djnEnd => flowDoms(nn,mm,1)%djnEnd
    dknEnd => flowDoms(nn,mm,1)%dknEnd

    icBeg => flowDoms(nn,mm,1)%icBeg
    jcBeg => flowDoms(nn,mm,1)%jcBeg
    kcBeg => flowDoms(nn,mm,1)%kcBeg
    icEnd => flowDoms(nn,mm,1)%icEnd
    jcEnd => flowDoms(nn,mm,1)%jcEnd
    kcEnd => flowDoms(nn,mm,1)%kcEnd

    neighBlock => flowDoms(nn,mm,1)%neighBlock
    neighProc  => flowDoms(nn,mm,1)%neighProc
    l1         => flowDoms(nn,mm,1)%l1
    l2         => flowDoms(nn,mm,1)%l2
    l3         => flowDoms(nn,mm,1)%l3
    groupNum   => flowDoms(nn,mm,1)%groupNum

    ! Overset boundary and hole info.
    iblank => flowDoms(nn,mm,ll)%iblank
    fringes => flowDoms(nn,mm,ll)%fringes
    orphans => flowDoms(nn,mm,ll)%orphans
    nOrphans = flowDoms(nn,mm,ll)%nOrphans

    ! The data for boundary subfaces.

    BCData => flowDoms(nn,mm,ll)%BCData

    ! The stress tensor and heat flux vector at viscous wall faces
    ! as well as the face pointers to these viscous wall faces.
    ! The latter point to the 1st spectral mode, because they are
    ! the only ones allocated. The info is the same for all modes.

    viscSubface => flowDoms(nn,mm,ll)%viscSubface

    viscIminPointer => flowDoms(nn,mm,1)%viscIminPointer
    viscImaxPointer => flowDoms(nn,mm,1)%viscImaxPointer
    viscJminPointer => flowDoms(nn,mm,1)%viscJminPointer
    viscJmaxPointer => flowDoms(nn,mm,1)%viscJmaxPointer
    viscKminPointer => flowDoms(nn,mm,1)%viscKminPointer
    viscKmaxPointer => flowDoms(nn,mm,1)%viscKmaxPointer

    ! Mesh related variables. The porosities point to the 1st
    ! spectral mode, because they are the only ones allocated.
    ! The info is the same for all modes.
    ! Note that xOld and volOld always point to the finest
    ! grid level.

    x    => flowDoms(nn,mm,ll)%x
    xOld => flowDoms(nn,1,ll)%xOld

    si     => flowDoms(nn,mm,ll)%si
    sj     => flowDoms(nn,mm,ll)%sj
    sk     => flowDoms(nn,mm,ll)%sk

    vol    => flowDoms(nn,mm,ll)%vol
    volRef => flowDoms(nn,mm,ll)%volRef
    volOld => flowDoms(nn,1,ll)%volOld

    porI => flowDoms(nn,mm,1)%porI
    porJ => flowDoms(nn,mm,1)%porJ
    porK => flowDoms(nn,mm,1)%porK

    indFamilyI => flowDoms(nn,mm,1)%indFamilyI
    indFamilyJ => flowDoms(nn,mm,1)%indFamilyJ
    indFamilyK => flowDoms(nn,mm,1)%indFamilyK

    factFamilyI => flowDoms(nn,mm,1)%factFamilyI
    factFamilyJ => flowDoms(nn,mm,1)%factFamilyJ
    factFamilyK => flowDoms(nn,mm,1)%factFamilyK

    rotMatrixI => flowDoms(nn,mm,ll)%rotMatrixI
    rotMatrixJ => flowDoms(nn,mm,ll)%rotMatrixJ
    rotMatrixK => flowDoms(nn,mm,ll)%rotMatrixK

    blockIsMoving     = flowDoms(nn,mm,ll)%blockIsMoving
    addGridVelocities = flowDoms(nn,mm,ll)%addGridVelocities

    sFaceI => flowDoms(nn,mm,ll)%sFaceI
    sFaceJ => flowDoms(nn,mm,ll)%sFaceJ
    sFaceK => flowDoms(nn,mm,ll)%sFaceK

    ! Flow variables. Note that wOld, gamma and the laminar viscosity
    ! point to the entries on the finest mesh. The reason is that
    ! they are computed from the other variables. For the eddy
    ! viscosity this is not the case because in a segregated solver
    ! its values are obtained from the fine grid level.

    w     => flowDoms(nn,mm,ll)%w
    wOld  => flowDoms(nn,1, ll)%wOld
    p     => flowDoms(nn,mm,ll)%p
    aa    => flowDoms(nn,mm,ll)%aa

    gamma => flowDoms(nn,1, ll)%gamma
    rlv   => flowDoms(nn,1, ll)%rlv
    rev   => flowDoms(nn,mm,ll)%rev
    s     => flowDoms(nn,mm,ll)%s

    ux => flowDoms(nn,mm,ll)%ux
    uy => flowDoms(nn,mm,ll)%uy
    uz => flowDoms(nn,mm,ll)%uz

    vx => flowDoms(nn,mm,ll)%vx
    vy => flowDoms(nn,mm,ll)%vy
    vz => flowDoms(nn,mm,ll)%vz

    wx => flowDoms(nn,mm,ll)%wx
    wy => flowDoms(nn,mm,ll)%wy
    wz => flowDoms(nn,mm,ll)%wz

    qx => flowDoms(nn,mm,ll)%qx
    qy => flowDoms(nn,mm,ll)%qy
    qz => flowDoms(nn,mm,ll)%qz


    ! Residual and multigrid variables. The residual point to the
    ! finest grid entry, the multigrid variables to their own level.

    dw => flowDoms(nn,1,ll)%dw
    fw => flowDoms(nn,1,ll)%fw
    dwOldRK => flowDoms(nn,1,ll)%dwOldRK
    scratch => flowDoms(nn,1,ll)%scratch

    p1 => flowDoms(nn,mm,ll)%p1
    w1 => flowDoms(nn,mm,ll)%w1
    wr => flowDoms(nn,mm,ll)%wr

    ! Variables, which allow a more flexible multigrid treatment.
    ! They are the same for all spectral modes and therefore they
    ! point to the 1st mode.

    mgIFine => flowDoms(nn,mm,1)%mgIFine
    mgJFine => flowDoms(nn,mm,1)%mgJFine
    mgKFine => flowDoms(nn,mm,1)%mgKFine

    mgIWeight => flowDoms(nn,mm,1)%mgIWeight
    mgJWeight => flowDoms(nn,mm,1)%mgJWeight
    mgKWeight => flowDoms(nn,mm,1)%mgKWeight

    mgICoarse => flowDoms(nn,mm,1)%mgICoarse
    mgJCoarse => flowDoms(nn,mm,1)%mgJCoarse
    mgKCoarse => flowDoms(nn,mm,1)%mgKCoarse

    ! Time-stepping variables and spectral radIi.
    ! They all point to the fine mesh entry.

    wn  => flowDoms(nn,1,ll)%wn
    pn  => flowDoms(nn,1,ll)%pn
    dtl => flowDoms(nn,1,ll)%dtl

    radI => flowDoms(nn,1,ll)%radI
    radJ => flowDoms(nn,1,ll)%radJ
    radK => flowDoms(nn,1,ll)%radK

    ! Wall distance for the turbulence models.

    d2Wall => flowDoms(nn,mm,ll)%d2Wall
    filterDES   => flowDoms(nn,mm,ll)%filterDES  ! eran-des

    ! Arrays used for the implicit treatment of the turbulent wall
    ! boundary conditions. As these variables are only allocated for
    ! the 1st spectral solution of the fine mesh, the pointers point
    ! to those arrays.

    bmti1 => flowDoms(nn,1,1)%bmti1
    bmti2 => flowDoms(nn,1,1)%bmti2
    bmtj1 => flowDoms(nn,1,1)%bmtj1
    bmtj2 => flowDoms(nn,1,1)%bmtj2
    bmtk1 => flowDoms(nn,1,1)%bmtk1
    bmtk2 => flowDoms(nn,1,1)%bmtk2

    bvti1 => flowDoms(nn,1,1)%bvti1
    bvti2 => flowDoms(nn,1,1)%bvti2
    bvtj1 => flowDoms(nn,1,1)%bvtj1
    bvtj2 => flowDoms(nn,1,1)%bvtj2
    bvtk1 => flowDoms(nn,1,1)%bvtk1
    bvtk2 => flowDoms(nn,1,1)%bvtk2

    ! Pointers for globalCell/Node
    globalCell =>flowDoms(nn,mm,ll)%globalCell
    globalNode =>flowDoms(nn,mm,ll)%globalNode
    ifaceptb  =>flowDoms(nn,1,1)%ifaceptb
    iedgeptb  =>flowDoms(nn,1,1)%iedgeptb

    xSeed => flowDoms(nn,mm,ll)%xSeed
    wallInd => flowDoms(nn,mm,ll)%wallInd

    ! *******************************
    ! Added by HDN
    ! Kept the same dim as their counterparts
    ! *******************************
    xALE      => flowDoms(nn,mm,ll)%xALE
    sVeloIALE => flowDoms(nn,mm,ll)%sVeloIALE
    sVeloJALE => flowDoms(nn,mm,ll)%sVeloJALE
    sVeloKALE => flowDoms(nn,mm,ll)%sVeloKALE
    sIALE     => flowDoms(nn,mm,ll)%sIALE
    sJALE     => flowDoms(nn,mm,ll)%sJALE
    sKALE     => flowDoms(nn,mm,ll)%sKALE
    sFaceIALE => flowDoms(nn,mm,ll)%sFaceIALE
    sFaceJALE => flowDoms(nn,mm,ll)%sFaceJALE
    sFaceKALE => flowDoms(nn,mm,ll)%sFaceKALE
    dwALE     => flowDoms(nn,1,ll)%dwALE
    fwALE     => flowDoms(nn,1,ll)%fwALE

#ifndef USE_TAPENADE
    ! Pointers for PC
    PCMat => flowDoms(nn,mm,ll)%pcMat

    i_D_Fact => flowDoms(nn,mm,ll)%i_D_fact
    i_L_Fact => flowDoms(nn,mm,ll)%i_L_fact
    i_U_Fact => flowDoms(nn,mm,ll)%i_U_fact
    i_U2_Fact => flowDoms(nn,mm,ll)%i_U2_fact

    j_D_Fact => flowDoms(nn,mm,ll)%j_D_fact
    j_L_Fact => flowDoms(nn,mm,ll)%j_L_fact
    j_U_Fact => flowDoms(nn,mm,ll)%j_U_fact
    j_U2_Fact => flowDoms(nn,mm,ll)%j_U2_fact

    k_D_Fact => flowDoms(nn,mm,ll)%k_D_fact
    k_L_Fact => flowDoms(nn,mm,ll)%k_L_fact
    k_U_Fact => flowDoms(nn,mm,ll)%k_U_fact
    k_U2_Fact => flowDoms(nn,mm,ll)%k_U2_fact

    PCVec1 => flowDoms(nn,mm,ll)%PCVec1
    PCVec2 => flowDoms(nn,mm,ll)%PCVec2

    i_ipiv => flowDoms(nn,mm,ll)%i_ipiv
    j_ipiv => flowDoms(nn,mm,ll)%j_ipiv
    k_ipiv => flowDoms(nn,mm,ll)%k_ipiv
#endif
  end subroutine setPointers

  subroutine siAngle(angle, mult, trans)

    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: angle
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    if(angle == Radian) then

       ! Angle is already given in radIans. No need for a conversion.

       mult  = one
       trans = zero

    else if(angle == Degree) then

       ! Angle is given in degrees. A multiplication must be performed.

       mult  = pi/180.0_realType
       trans = zero

    else

       call returnFail("siAngle", &
            "No idea how to convert this to SI units")

    endif

  end subroutine siAngle


  subroutine siDensity(mass, len, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siDensity computes the conversion from the given density       *
    !      * unit, which can be constructed from mass and length, to the    *
    !      * SI-unit kg/m^3. The conversion will look like:                 *
    !      * density in kg/m^3 = mult*(density in NCU) + trans.             *
    !      * NCU means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: mass, len
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    if(mass == Kilogram .and. len == Meter) then

       ! Density is given in kg/m^3, i.e. no need for a conversion.

       mult  = one
       trans = zero

    else

       call returnFail("siDensity", &
            "No idea how to convert this to SI units")

    endif

  end subroutine siDensity

  subroutine siLen(len, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siLen computes the conversion from the given length unit to    *
    !      * the SI-unit meter. The conversion will look like:              *
    !      * length in meter = mult*(length in NCU) + trans.                *
    !      * NCU means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: len
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    select case (len)

    case (Meter)
       mult = one; trans = zero

    case (CenTimeter)
       mult = 0.01_realType; trans = zero

    case (Millimeter)
       mult = 0.001_realType; trans = zero

    case (Foot)
       mult = 0.3048_realType; trans = zero

    case (Inch)
       mult = 0.0254_realType; trans = zero

    case default
       call returnFail("siLen", &
            "No idea how to convert this to SI units")

    end select

  end subroutine siLen

  subroutine siPressure(mass, len, time, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siPressure computes the conversion from the given pressure     *
    !      * unit, which can be constructed from mass, length and time, to  *
    !      * the SI-unit Pa. The conversion will look like:                 *
    !      * pressure in Pa = mult*(pressure in NCU) + trans.               *
    !      * NCU means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: mass, len, time
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    if(mass == Kilogram .and. len == Meter .and. time == Second) then

       ! Pressure is given in Pa, i.e. no need for a conversion.

       mult  = one
       trans = zero

    else

       call returnFail("siPressure", &
            "No idea how to convert this to SI units")

    endif

  end subroutine siPressure

  subroutine siTemperature(temp, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siTemperature computes the conversion from the given           *
    !      * temperature unit to the SI-unit kelvin. The conversion will    *
    !      * look like:                                                     *
    !      * temperature in K = mult*(temperature in NCU) + trans.          *
    !      * NCU means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: temp
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    select case (temp)

    case (Kelvin)

       ! Temperature is already given in Kelvin. No need to convert.

       mult  = one
       trans = zero

    case (Celsius)      ! is it Celcius or Celsius?

       ! Temperature is in Celsius. Only an offset must be applied.

       mult  = one
       trans = 273.16_realType

    case (Rankine)

       ! Temperature is in Rankine. Only a multiplication needs to
       ! be performed.

       mult  = 5.0_realType/9.0_realType
       trans = zero

    case (Fahrenheit)

       ! Temperature is in Fahrenheit. Both a multiplication and an
       ! offset must be applied.

       mult  = 5.0_realType/9.0_realType
       trans = 255.382

    case default

       ! Unknown temperature unit.

       call returnFail("siTemperature", &
            "No idea how to convert this to SI units")

    end select

  end subroutine siTemperature
  !
  !      ******************************************************************
  !      *                                                                *
  !      * File:          siTurb.f90                                      *
  !      * Author:        Edwin van der Weide                             *
  !      * Starting date: 09-14-2004                                      *
  !      * Last modified: 06-12-2005                                      *
  !      *                                                                *
  !      ******************************************************************
  !
  subroutine siTurb(mass, len, time, temp, turbName, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siTurb computes the conversion from the given turbulence       *
    !      * unit, which can be constructed from mass, len, time and temp,  *
    !      * to the SI-unit for the given variable. The conversion will     *
    !      * look like: var in SI = mult*(var in NCU) + trans.              *
    !      * NCU means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: mass, len, time, temp
    character(len=*), intent(in)     :: turbName
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.

    if(mass == Kilogram .and. len  == Meter .and. &
         time == Second   .and. temp == Kelvin) then

       ! Everthing is already in SI units. No conversion needed.

       mult  = one
       trans = zero

    else

       call returnFail("siTurb", &
            "No idea how to convert this to SI units")

    endif

  end subroutine siTurb

  subroutine siVelocity(length, time, mult, trans)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * siVelocity computes the conversion from the given velocity     *
    !      * unit, which can be constructed from length and time, to the    *
    !      * SI-unit m/s. The conversion will look like:                    *
    !      * velocity in m/s = mult*(velocity in ncu) + trans.              *
    !      * Ncu means non-christian units, i.e. everything that is not SI. *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use su_cgns
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)              :: length, time
    real(kind=realType), intent(out) :: mult, trans

    ! Determine the situation we are having here.
    ! First the length.

    select case (length)

    case (Meter)
       mult = one; trans = zero

    case (CenTimeter)
       mult = 0.01_realType; trans = zero

    case (Millimeter)
       mult = 0.001_realType; trans = zero

    case (Foot)
       mult = 0.3048_realType; trans = zero

    case (Inch)
       mult = 0.0254_realType; trans = zero

    case default
       call returnFail("siVelocity", &
            "No idea how to convert this length to &
            &SI units")

    end select

    ! And the time.

    select case (time)

    case (Second)
       mult = mult

    case default
       call returnFail("siVelocity", &
            "No idea how to convert this time to &
            &SI units")

    end select

  end subroutine siVelocity

  subroutine spectralInterpolCoef(nsps, t, alpScal, alpMat)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * spectralInterpolCoef determines the scalar and matrix          *
    !      * spectral interpolation coefficients for the given number of    *
    !      * spectral solutions for the given t, where t is the ratio of    *
    !      * the time and the periodic interval time. Note that the index   *
    !      * of the spectral solutions of both alpScal and alpMat start     *
    !      * at 0. In this way these coefficients are easier to determine.  *
    !      *                                                                *
    !      ******************************************************************
    !
    use constants
    use inputTimeSpectral
    use section
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nsps
    real(kind=realType),   intent(in) :: t

    real(kind=realType), dimension(0:nsps-1), intent(out) :: alpScal
    real(kind=realType), dimension(nSections,0:nsps-1,3,3), &
         intent(out) :: alpMat
    !
    !      Local variables.
    !
    integer(kind=intType) :: jj, nn, j, p, r, nhalfM1, m, mhalfM1

    real(kind=realType) :: nspsInv, mInv, tm, alp

    real(kind=realType), dimension(3,3) :: rp, tmp

    !      ******************************************************************
    !      *                                                                *
    !      * Scalar coefficients.                                           *
    !      *                                                                *
    !      ******************************************************************
    !
    ! Loop over the number of spectral solutions to compute the
    ! coefficients. Note that the loop starts at 0.

    if (mod(nsps,2).eq.0) then
       nhalfM1 = nsps/2 - 1
    else
       nhalfM1 = (nsps-1)/2
    endif

    nspsInv = one/real(nsps,realType)

    do j=0,(nsps-1)
       if (mod(nsps,2).eq.0) then
          alpScal(j) = one + cos(j*pi)*cos(nsps*pi*t)
       else 
          alpScal(j) = one + cos(j*pi*(nsps+1)/nsps)*cos((nsps+1)*pi*t)
       endif

       do r=1,nhalfM1
          alpScal(j) = alpScal(j)                                  &
               + two*cos(r*j*two*pi*nspsInv)*cos(r*two*pi*t) &
               + two*sin(r*j*two*pi*nspsInv)*sin(r*two*pi*t)
       enddo

       alpScal(j) = alpScal(j)*nspsInv

    enddo
    !
    !      ******************************************************************
    !      *                                                                *
    !      * Matrix coefficients. These are (can be) different for every    *
    !      * section and they must therefore be determined for every        *
    !      * section.                                                       *
    !      *                                                                *
    !      ******************************************************************
    !
    ! Loop over the number of sections in the grid.

    sectionLoop: do nn=1,nSections

       ! Compute the numbers for the entire wheel for this section.
       ! Note that also t must be adapted, because t is a ratio between
       ! the actual time and the periodic time.

       m        = nsps*sections(nn)%nSlices
       if (mod(m,2).eq.0) then
          mhalfM1 = m/2 - 1
       else 
          mhalfM1 = (m-1)/2 
       endif
       mInv    = one/real(m,realType)
       tm       = t/real(sections(nn)%nSlices,realType)

       ! Loop over the number of spectral solutions.

       spectralLoop: do jj=0,(nsps-1)

          ! Initialize the matrix coefficients to zero and the matrix
          ! rp to the identity matrix. Rp is the rotation matrix of this
          ! section to the power p, which starts at 0, i.e. rp = i.

          alpMat(nn,jj,1,1) = zero
          alpMat(nn,jj,1,2) = zero
          alpMat(nn,jj,1,3) = zero

          alpMat(nn,jj,2,1) = zero
          alpMat(nn,jj,2,2) = zero
          alpMat(nn,jj,2,3) = zero

          alpMat(nn,jj,3,1) = zero
          alpMat(nn,jj,3,2) = zero
          alpMat(nn,jj,3,3) = zero

          rp(1,1) = one
          rp(1,2) = zero
          rp(1,3) = zero

          rp(2,1) = zero
          rp(2,2) = one
          rp(2,3) = zero

          rp(3,1) = zero
          rp(3,2) = zero
          rp(3,3) = one

          ! Loop over the number of slices of this section. Note that
          ! this loop starts at zero, which simplifies the formulas.

          slicesLoop: do p=0,(sections(nn)%nSlices-1)

             ! Determine the index j, the index of alp in the entire
             ! wheel.

             j = jj + p*nsps

             ! Compute the scalar coefficient alp of the index j in
             ! the entire wheel.

             if (mod(m,2).eq.0) then
                alp = one + cos(j*pi)*cos(m*pi*tm)
             else
                alp = one + cos(j*pi*(m+1)/m)*cos((m+1)*pi*tm)
             endif
             do r=1,mhalfM1
                alp = alp + two*cos(r*j*two*pi*mInv)*cos(r*two*pi*tm) &
                     +       two*sin(r*j*two*pi*mInv)*sin(r*two*pi*tm)
             enddo

             alp = alp*mInv

             ! Update the matrix coefficient.

             do r=1,3
                do j=1,3
                   alpMat(nn,jj,r,j) = alpMat(nn,jj,r,j) + alp*rp(r,j)
                enddo
             enddo

             ! Multiply rp by the rotation matrix to obtain the correct
             ! matrix for the next slice. Use tmp as temporary storage.

             do r=1,3
                do j=1,3
                   tmp(r,j) = rp(r,1)*rotMatrixSpectral(nn,1,j) &
                        + rp(r,2)*rotMatrixSpectral(nn,2,j) &
                        + rp(r,3)*rotMatrixSpectral(nn,3,j)
                enddo
             enddo

             rp = tmp

          enddo slicesLoop
       enddo spectralLoop
    enddo sectionLoop

  end subroutine spectralInterpolCoef

  subroutine deallocateTempMemory(resNeeded)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * deallocateTempMemory deallocates memory used in the solver,    *
    !      * but which is not needed to store the actual solution. In this  *
    !      * way the memory can be used differently, e.g. when writing the  *
    !      * solution or computing the wall distances.                      *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    use communication
    use inputIteration
    use inputTimeSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: resNeeded
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm

    ! Deallocate the communication buffers

    deallocate(sendBuffer, recvBuffer, stat=ierr)
    if(ierr /= 0)                              &
         call returnFail("deallocateTempMemory", &
         "Deallocation error for communication buffers")

    ! Loop over the spectral modes and domains. Note that only memory
    ! on the finest grid is released, because a) most of these
    ! variables are only allocated on the fine grid and b) the coarser
    ! grids do not contribute that much in the memory usage anyway.

    spectralModes: do mm=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Check if the residual, time step, etc. Is needed.

          if(.not. resNeeded) then

             ! Residual, etc. Not needed.
             ! Deallocate residual, the time step and the spectral radii
             ! of the fine level.

             deallocate(flowDoms(nn,1,mm)%dw,   flowDoms(nn,1,mm)%fw,   &
                  flowDoms(nn,1,mm)%dtl,  flowDoms(nn,1,mm)%radI, &
                  flowDoms(nn,1,mm)%radJ, flowDoms(nn,1,mm)%radK, &
                  stat=ierr)
             if(ierr /= 0)                            &
                  call returnFail("deallocateTempMemory", &
                  "Deallocation error for dw, fw, dtl and &
                  &spectral radii.")
          endif

          ! The memory for the zeroth Runge Kutta stage
          ! if a Runge Kutta scheme is used.

          if(smoother == RungeKutta) then

             deallocate(flowDoms(nn,1,mm)%wn, flowDoms(nn,1,mm)%pn, &
                  stat=ierr)
             if(ierr /= 0)                            &
                  call returnFail("deallocateTempMemory", &
                  "Deallocation error for wn and pn")
          endif

       enddo domains
    enddo spectralModes

  end subroutine deallocateTempMemory

  subroutine allocateTempMemory(resNeeded)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * AllocateTempMemory allocates the memory again that was         *
    !      * temporarily deallocted by deallocateTempMemory.                *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    use communication
    use constants
    use flowVarRefState
    use inputIteration
    use inputTimeSpectral
    use iteration
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: resNeeded
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn,mm
    integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb

    ! The memory for the receive buffers.

    allocate(sendBuffer(sendBufferSize), &
         recvBuffer(recvBufferSize), stat=ierr)
    if(ierr /= 0)                          &
         call returnFail("allocateTempMemory", &
         "Memory allocation failure for comm buffers")

    ! Loop over the spectral modes and domains. Note that only memory
    ! on the finest mesh level needs to be reallocated, because the
    ! memory on the coarser levels has not been released or is not
    ! needed .

    spectralModes: do mm=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Store some dimensions a bit easier.

          il = flowDoms(nn,1,mm)%il
          jl = flowDoms(nn,1,mm)%jl
          kl = flowDoms(nn,1,mm)%kl

          ie = flowDoms(nn,1,mm)%ie
          je = flowDoms(nn,1,mm)%je
          ke = flowDoms(nn,1,mm)%ke

          ib = flowDoms(nn,1,mm)%ib
          jb = flowDoms(nn,1,mm)%jb
          kb = flowDoms(nn,1,mm)%kb

          ! Check if the residual, time step, etc. was deallocated.

          if(.not. resNeeded) then

             ! Allocate the residual, the time step and
             ! the spectral radii.

             allocate(flowDoms(nn,1,mm)%dw(0:ib,0:jb,0:kb,1:nw),  &
                  flowDoms(nn,1,mm)%fw(0:ib,0:jb,0:kb,1:nwf), &
                  flowDoms(nn,1,mm)%dtl(1:ie,1:je,1:ke),      &
                  flowDoms(nn,1,mm)%radI(1:ie,1:je,1:ke),     &
                  flowDoms(nn,1,mm)%radJ(1:ie,1:je,1:ke),     &
                  flowDoms(nn,1,mm)%radK(1:ie,1:je,1:ke), stat=ierr)
             if(ierr /= 0)                            &
                  call returnFail("allocateTempMemory", &
                  "Memory allocation failure for dw, fw, &
                  &dtl and the spectral radii.")

             ! Initialize dw and fw to zero to avoid possible overflows
             ! of the halo's.

             flowDoms(nn,1,mm)%dw = zero
             flowDoms(nn,1,mm)%fw = zero

          endif

          ! The memory for the zeroth runge kutta stage
          ! if a runge kutta scheme is used.

          if(smoother == RungeKutta) then

             allocate(flowDoms(nn,1,mm)%wn(2:il,2:jl,2:kl,1:nwf), &
                  flowDoms(nn,1,mm)%pn(2:il,2:jl,2:kl), stat=ierr)
             if(ierr /= 0)                            &
                  call returnFail("allocateTempMemory", &
                  "Memory allocation failure for wn and pn")
          endif

       enddo domains
    enddo spectralModes

  end subroutine allocateTempMemory

  subroutine getLiftDirFromSymmetry(liftDir)

    ! The purpose of this function is to determine what coordinate
    ! direction the mirror plane is in. It does NOT handle multiple mirror
    ! planes. It is used just to determine what the lift direction is. 

    use constants
    use bcTypes
    use blockPointers
    use communication
    implicit none

    ! Output
    integer(kind=intType), intent(out) :: liftDir
    integer(kind=intType), dimension(3) :: sym_local, sym

    ! Working
    integer(kind=intType) :: nn, i_index(1), mm, ierr
    real(kind=realType), dimension(:, :, :), pointer :: xx
    real(kind=realType) :: cp(3), v1(3), v2(3)
    ! Loop over each block and each subFace

    sym_local = 0_intType
    sym       = 0_intType
    liftDir   = 0_intType
    do nn=1, nDom
       call setPointers(nn, 1, 1)
       do mm=1,nBocos
          if (bcType(mm) == symm) then

             select case (BCFaceID(mm))
             case (iMin)
                xx => x(1, :, :, :)
             case (iMax)
                xx => x(il, :, :, :)
             case (jMin)
                xx => x(:, 1, :, :)
             case (jMax)
                xx => x(:, jl, :, :)
             case (kMin)
                xx => x(:, :, 1, :)
             case (kMax)
                xx => x(:, :, kl, :)
             end select

             ! Take the cross product
             v1(:) = xx(bcData(mm)%inEnd, bcData(mm)%jnEnd, :) - &
                  xx(bcData(mm)%inBeg, bcData(mm)%jnBeg, :)
             v2(:) = xx(bcData(mm)%inBeg, bcData(mm)%jnEnd, :) - &
                  xx(bcData(mm)%inEnd, bcData(mm)%jnBeg, :)

             ! Cross Product
             cp(1) = (v1(2)*v2(3) - v1(3)*v2(2))
             cp(2) = (v1(3)*v2(1) - v1(1)*v2(3))
             cp(3) = (v1(1)*v2(2) - v1(2)*v2(1))  

             ! Only interesed in abs values
             cp = abs(cp)

             ! Location, ie coordiante direction of dominate direction
             i_index = maxloc(real(cp))

             sym_local(i_index(1)) = 1_intType
          end if
       end do
    end do

    ! Now we have a bunch of sym_locals, mpi_allreduce them and SUM

    call MPI_Allreduce (sym_local, sym, 3, sumb_integer, &
         MPI_SUM, sumb_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Now we should make sure that only ONE of the values is
    ! non-zero. If more than one value is zero, it means we have
    ! multiple symmetry planes which we can't support.
    if (sym(1) == 0 .and. sym(2) == 0 .and. sym(3) == 0) then
       ! Pass - no sym, can't determine lift dir:
    else if(sym(1) .ne. 0 .and. sym(2) == 0 .and. sym(3) == 0) then
       ! Pass - x dir can't be symmetry
    else if(sym(1) == 0 .and. sym(2) .ne. 0 .and. sym(3) == 0) then
       liftDir = 3
    else if(sym(1) == 0 .and. sym(2) == 0 .and. sym(3) .ne. 0) then
       liftDir = 2
    else
       ! Multiple orientations...can't do anything
    end if

  end subroutine getLiftDirFromSymmetry


  subroutine writeIntroMessage
    !
    !      ******************************************************************
    !      *                                                                *
    !      * writeIntroMessage writes a message to stdout with              *
    !      * information how the executable was built, e.g. whether single  *
    !      * or double precision is used for the integers and reals, etc.   *
    !      * To avoid a messy output only processor 0 prints this info.     *
    !      *                                                                *
    !      ******************************************************************
    !
    use communication
    use constants
    implicit none
    !
    !      Local variables
    !
    character(len=7) :: integerString

    ! Return if this is not processor 0.

    if(myID > 0) return

    ! I'm processor 0. Write the info to stdout.

    print "(a)", "#"
    print "(a)", "# SUmb, multiblock structured flow solver"
    print "(a)", "#"
    print "(a)", "# This code solves the 3D RANS, laminar NS or &
         &Euler equations"
    print "(a)", "# on multiblock structured hexahedral grids."


    write(integerString,"(i7)") nProc
    integerString = adjustl(integerString)
    print "(3a)", "# This is a parallel executable running on ", &
         trim(integerString), " processors."
    print "(a)", "# It has been compiled with the &
         &following options:"

    if( debug ) then
       print "(a)", "# - Debug mode."
    else
       print "(a)", "# - Optimized mode."
    endif

#ifdef USE_LONG_INT
    print "(a)", "# - Size of standard integers: 8 bytes."
#else
    print "(a)", "# - Size of standard integers: 4 bytes."
#endif

#ifdef USE_SINGLE_PRECISION
    print "(a)", "# - Size of standard floating point types: &
         &4 bytes."

#elif  USE_QUADRUPLE_PRECISION
    print "(a)", "# - Size of standard floating point types: &
         &16 bytes."
#else
    print "(a)", "# - Size of standard floating point types: &
         &8 bytes."
#endif

#ifdef USE_NO_CGNS
    print "(a)", "# - Without cgns support"
#else
    print "(a)", "# - With cgns support"
#endif

#ifdef USE_NO_SIGNALS
    print "(a)", "# - Without support for signals."
#else
    print "(a)", "# - With support for signals."
#endif

    print "(a)", "#"

  end subroutine writeIntroMessage

  subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

    ! Given a list of N points (pts) in three space, with possible
    ! duplicates, (to within tol) return a list of the nUnique
    ! uniquePoints of points and a link array of length N, that points
    ! into the unique list
    use precision
    use kdtree2_module
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: N
    real(kind=realType), intent(in), dimension(3, N) :: pts
    real(kind=realType), intent(in) :: tol

    ! Output Parametres
    real(kind=realType), intent(out), dimension(3, N) :: uniquePts
    integer(kind=intType), intent(out), dimension(N) :: link
    integer(kind=intType), intent(out) :: nUnique

    ! Working paramters
    type(kdtree2), pointer :: mytree
    real(kind=realType) :: tol2, timeb, timea
    integer(kind=intType) :: nFound, i, j, nAlloc
    type(kdtree2_result), allocatable, dimension(:) :: results

    if (N==0) then 
       nUnique = 0
       return 
    end if

    ! We will use the KD_tree to do most of the heavy lifting here:

    mytree => kdtree2_create(pts, sort=.True.)

    ! KD tree works with the square of the tolerance
    tol2 = tol**2

    ! Unlikely we'll have more than 20 points same, but there is a
    ! safetly check anwyay.
    nalloc = 20
    allocate(results(nalloc))

    link = 0
    nUnique = 0

    ! Loop over all nodes
    do i=1, N
       if (link(i) == 0) then 
          call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)

          ! Expand if necesary and re-run
          if (nfound > nalloc) then 
             deallocate(results)
             nalloc = nfound
             allocate(results(nalloc))
             call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)
          end if

          if (nFound == 1) then 
             ! This one is easy, it is already a unique node
             nUnique = nUnique + 1
             link(i) = nUnique
             uniquePts(:, nUnique) = pts(:, i)
          else
             if (link(i) == 0) then 
                ! This node hasn't been assigned yet:
                nUnique = nUnique + 1
                uniquePts(:, nUnique) = pts(:, i)

                do j=1, nFound
                   link(results(j)%idx) = nUnique
                end do
             end if
          end if
       end if
    end do

    ! Done with the tree and the result vector
    call kdtree2destroy(mytree)
    deallocate(results)

  end subroutine pointReduce
  !
  !      ******************************************************************
  !      *                                                                *
  !      * File:          releaseMemory.f90                               *
  !      * Author:        Edwin van der Weide, Steve Repsher              *
  !      * Starting date: 08-16-2004                                      *
  !      * Last modified: 11-27-2007                                      *
  !      *                                                                *
  !      ******************************************************************
  !
  subroutine releaseMemoryPart1
    !
    !      ******************************************************************
    !      *                                                                *
    !      * releaseMemoryPart1 releases all the memory on the coarser      *
    !      * grids of flowDoms and the fine grid memory which is not needed *
    !      * for the possible interpolation of the spectral solution.       *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    use inputIteration
    use inputTimeSpectral
    use inputPhysics
    use inputUnsteady
    use monitor
    use BCTypes
    use cgnsGrid
    use communication
    use iteration
    use cgnsGrid
    use bleedFlows
    use section
    use interfaceGroups
    use commSliding
    use commMixing
    use wallDistanceData
    use adjointVars
    use ADJointPETSc
    use surfaceFamilies
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    integer(kind=intType) :: sps, nLevels, level, nn, l, i, j

    ! Determine the number of grid levels present in flowDoms.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the coarser grid levels and local blocks and
       ! deallocate all the memory.

       do level=2,nLevels
          do nn=1,nDom
             call deallocateBlock(nn, level, sps)
          enddo
       enddo

       ! Release some memory of the fine grid, which is not needed
       ! anymore.

       do nn=1,nDom
          ! *******************************
          ! Modified by HDN
          ! Added dwALE, fwALE
          ! *******************************
          deallocate( &
               flowDoms(nn,1,sps)%dw,    flowDoms(nn,1,sps)%fw,    &
               flowDoms(nn,1,sps)%dwALE, flowDoms(nn,1,sps)%fwALE, &
               flowDoms(nn,1,sps)%dtl,   flowDoms(nn,1,sps)%radI,  &
               flowDoms(nn,1,sps)%radJ,  flowDoms(nn,1,sps)%radK,  &
               stat=ierr)
          if(ierr /= 0)                          &
               call returnFail("releaseMemoryPart1", &
               "Deallocation error for dw, fw, dwALE, fwALE, dtl and &
               &spectral radii.")

          ! Nullify the pointers, such that no attempt is made to
          ! release the memory again.

          nullify(flowDoms(nn,1,sps)%dw)
          nullify(flowDoms(nn,1,sps)%fw)
          nullify(flowDoms(nn,1,sps)%dwALE) ! Added by HDN
          nullify(flowDoms(nn,1,sps)%fwALE) ! Added by HDN
          nullify(flowDoms(nn,1,sps)%dtl)
          nullify(flowDoms(nn,1,sps)%radI)
          nullify(flowDoms(nn,1,sps)%radJ)
          nullify(flowDoms(nn,1,sps)%radK)
          nullify(flowDoms(nn,1,sps)%scratch)

          ! Check if the zeroth stage runge kutta memory has been
          ! allocated. If so deallocate it and nullify the pointers.

          if(smoother == RungeKutta) then

             deallocate(flowDoms(nn,1,sps)%wn, flowDoms(nn,1,sps)%pn, &
                  stat=ierr)
             if(ierr /= 0)                          &
                  call returnFail("releaseMemoryPart1", &
                  "Deallocation error for wn and pn")

             nullify(flowDoms(nn,1,sps)%wn)
             nullify(flowDoms(nn,1,sps)%pn)

          endif

          ! Release the memory of the old residuals for the time
          ! accurate Runge-Kutta schemes.

          if(equationMode          == unsteady .and. &
               timeIntegrationScheme == explicitRK) then

             deallocate(flowDoms(nn,1,sps)%dwOldRK, stat=ierr)
             if(ierr /= 0)                          &
                  call returnFail("releaseMemoryPart1", &
                  "Deallocation error for dwOldRK,")

             nullify(flowDoms(nn,1,sps)%dwOldRK)
          endif

       enddo

    enddo spectralLoop

    ! derivative values
    if (derivVarsAllocated) then 
       call deallocDerivativeValues(1)
    end if

    ! Bunch of extra sutff that hasn't been deallocated
    if (allocated(cycleStrategy)) then
       deallocate(cycleStrategy)
    end if

    if (allocated(monNames)) then 
       deallocate(monNames)
    end if

    if (allocated(monLoc)) then 
       deallocate(monLoc)
    end if

    if (allocated(monGlob)) then 
       deallocate(monGlob)
    end if

    if (allocated(monRef)) then 
       deallocate(monRef)
    end if

    if (allocated(cgnsFamilies)) then 
       deallocate(cgnsFamilies)
    end if

    deallocate(cgnsDomsd)
    deallocate(famIDsDomainInterfaces, &
         bcIDsDomainInterfaces,  &
         famIDsSliding)
    deallocate(inflowBleeds, outflowBleeds)
    deallocate(sections)
    deallocate(myinterfaces)

    ! Destroy wall distance stuff if necessary
    do l=1,nLevels
       call destroyWallDistanceData(l)
    end do
    deallocate(xSurfVec, xVolumeVec, wallScatter)

    ! Destroy the traction force stuff
    do j=1, size(familyExchanges, 2)
       do i=1, size(familyExchanges, 1)
          call destroyFamilyExchange(familyExchanges(i, j))
       end do
    end do
    deallocate(familyExchanges)

    do i=1, size(wallExchange)
       call destroyFamilyExchange(wallExchange(i))
    end do
    deallocate(wallExchange)

    ! From Communication Stuff
    do l=1,nLevels

       call deallocateCommType(commPatternCell_1st(l))
       call deallocateCommType(commPatternCell_2nd(l))
       call deallocateCommType(commPatternNode_1st(l))

       call deallocateInternalCommType(internalCell_1st(l))
       call deallocateInternalCommType(internalCell_2nd(l))
       call deallocateInternalCommType(internalNode_1st(l))

       do sps=1,nTimeIntervalsSpectral
          call deallocateSlidingCommType(commslidingCell_1st(l,sps))
          call deallocateSlidingCommType(commslidingCell_2nd(l,sps))
       end do

    end do
    deallocate(nCellGlobal)

    ! Now deallocate the containers
    deallocate(&
         commPatternCell_1st, commPatternCell_2nd, commPatternNode_1st, &
         internalCell_1st, internalCell_2nd, internalNode_1st)

    ! The remainder of the comms are just deallocated...these still need
    ! to be treated properly
    deallocate(commSlidingCell_1st, commSlidingCell_2nd, &
         intSlidingCell_1st, intSlidingCell_2nd, commPatternMixing)

    ! Send/recv buffer
    if (allocated(sendBuffer)) then
       deallocate(sendBuffer)
    end if

    if (allocated(recvBuffer)) then 
       deallocate(recvBuffer)
    end if

    ! massFlow stuff from setFamilyInfoFaces.f90
    deallocate(massFLowFamilyInv, massFlowFamilyDiss)

  end subroutine releaseMemoryPart1

  subroutine deallocateCommType(comm)
    use communication
    implicit none
    integer(kind=intType) :: ierr, i

    type(commType) :: comm
    ! Deallocate memory in comm

    ! Deallocate the sendLists
    do i=1, comm%nProcSend
       deallocate(comm%sendList(i)%block, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(comm%sendList(i)%indices, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(comm%sendList(i)%interp, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end do

    ! Deallocate the recvLists
    do i=1, comm%nProcRecv
       deallocate(comm%recvList(i)%block, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(comm%recvList(i)%indices, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end do

    deallocate(comm%sendProc, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%nsend, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%nsendcum, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%sendlist, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%recvProc, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%nrecv, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%nrecvcum, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%recvlist, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%indexsendproc, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%indexrecvproc, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (comm%nPeriodic > 0) then 
       do i=1,comm%nPeriodic
          deallocate(comm%periodicData(i)%block, stat=ierr)
          call EChk(ierr, __FILE__, __LINE__)

          deallocate(comm%periodicData(i)%indices)
          call EChk(ierr, __FILE__, __LINE__)
       end do

       deallocate(comm%periodicData, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

  end subroutine deallocateCommType

  subroutine deallocateInternalCommType(comm)
    use communication
    implicit none
    integer(kind=intType) :: ierr, i

    type(internalCommType) :: comm
    ! Deallocate memory in comm
    deallocate(comm%donorBlock, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%donorIndices, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%donorInterp, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%haloBlock, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%haloIndices, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (comm%nPeriodic > 0) then 
       do i=1,comm%nPeriodic
          deallocate(comm%periodicData(i)%block, stat=ierr)
          call EChk(ierr, __FILE__, __LINE__)

          deallocate(comm%periodicData(i)%indices)
          call EChk(ierr, __FILE__, __LINE__)
       end do
       deallocate(comm%periodicData, stat=ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

  end subroutine deallocateInternalCommType

  subroutine deallocateslidingCommType(comm)
    use communication
    use commSliding
    implicit none
    type(slidingCommType) :: comm
    integer(kind=intType) :: ierr
    deallocate(comm%nSendCum, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate(comm%nRecvCum, stat=ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine deallocateslidingCommType

  !      ==================================================================

  subroutine releaseMemoryPart2
    !
    !      ******************************************************************
    !      *                                                                *
    !      * releaseMemoryPart2 releases all the memory of flowDoms on the  *
    !      * finest grid as well as the memory allocated in the other       *
    !      * modules.                                                       *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    use inputTimeSpectral
    use ADjointPETSc
    use cgnsGrid
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    integer(kind=intType) :: nn, sps

    ! Release the memory of flowDoms of the finest grid and of the
    ! array flowDoms afterwards.

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call deallocateBlock(nn, 1_intType, sps)
       enddo
    enddo
    deallocate(flowDoms, stat=ierr)
    if(ierr /= 0)                          &
         call returnFail("releaseMemoryPart2", &
         "Deallocation failure for flowDoms")

    ! Some more memory should be deallocated if this code is to
    ! be used in combination with adaptation.

    ! Destroy variables allocated in preprocessingAdjoint

    call vecDestroy(w_like1,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    call vecDestroy(w_like2,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    call vecDestroy(psi_like1,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    call vecDestroy(psi_like2,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    call vecDestroy(psi_like3,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    call vecDestroy(x_like,PETScIerr)
    call EChk(PETScIerr, __FILE__, __LINE__)

    ! Finally delete cgnsDoms...but there is still more
    ! pointers that need to be deallocated...
    do nn=1,cgnsNDom
       if (associated(cgnsDoms(nn)%procStored)) &
            deallocate(cgnsDoms(nn)%procStored)

       if (associated(cgnsDoms(nn)%conn1to1)) &
            deallocate(cgnsDoms(nn)%conn1to1)

       if (associated(cgnsDoms(nn)%connNonMatchAbutting)) &
            deallocate(cgnsDoms(nn)%connNonMatchAbutting)

       if (associated(cgnsDoms(nn)%bocoInfo)) &
            deallocate(cgnsDoms(nn)%bocoInfo)

       deallocate(&
            cgnsDoms(nn)%iBegOr, cgnsDoms(nn)%iEndOr, &
            cgnsDoms(nn)%jBegOr, cgnsDoms(nn)%jEndOr, &
            cgnsDoms(nn)%kBegOr, cgnsDoms(nn)%kEndOr, &
            cgnsDoms(nn)%localBlockID)
    end do

  end subroutine releaseMemoryPart2

  subroutine deallocateBlock(nn, level, sps)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * deallocateBlock deallocates all the allocated memory of the    *
    !      * given block.                                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use block
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn, level, sps
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i

    type(viscSubfaceType), dimension(:), pointer :: viscSubface
    type(BCDataType),      dimension(:), pointer :: BCData

    logical :: deallocationFailure

    ! Initialize deallocationFailure to .false.

    deallocationFailure = .false.

    ! Set the pointer for viscSubface and deallocate the memory
    ! stored in there. Initialize ierr to 0, such that the returnFail
    ! routine is only called at the end if a memory deallocation
    ! failure occurs.
    ierr = 0
    viscSubface => flowDoms(nn,level,sps)%viscSubface
    do i=1,flowDoms(nn,level,sps)%nViscBocos
       deallocate(viscSubface(i)%tau,  viscSubface(i)%q, &
            viscSubface(i)%utau, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       nullify(viscSubface(i)%tau)
       nullify(viscSubface(i)%q)
       nullify(viscSubface(i)%utau)
    enddo

    ! Set the pointer for BCData and deallocate the memory
    ! stored in there.
    BCData => flowDoms(nn,level,sps)%BCData
    do i=1,flowDoms(nn,level,sps)%nBocos

       if( associated(BCData(i)%norm) ) &
            deallocate(BCData(i)%norm, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%area) ) &
            deallocate(BCData(i)%area, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%fIndex) ) &
            deallocate(BCData(i)%fIndex, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%F) ) &
            deallocate(BCData(i)%F, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%Fv) ) &
            deallocate(BCData(i)%Fv, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%Fp) ) &
            deallocate(BCData(i)%Fp, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%T) ) &
            deallocate(BCData(i)%T, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%Tv) ) &
            deallocate(BCData(i)%Tv, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%Tp) ) &
            deallocate(BCData(i)%Tp, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%rface) ) &
            deallocate(BCData(i)%rface, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%uSlip) ) &
            deallocate(BCData(i)%uSlip, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%TNS_Wall) ) &
            deallocate(BCData(i)%TNS_Wall, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%ptInlet) ) &
            deallocate(BCData(i)%ptInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%ttInlet) ) &
            deallocate(BCData(i)%ttInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%htInlet) ) &
            deallocate(BCData(i)%htInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%flowXdirInlet) ) &
            deallocate(BCData(i)%flowXdirInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%flowYdirInlet) ) &
            deallocate(BCData(i)%flowYdirInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%flowZdirInlet) ) &
            deallocate(BCData(i)%flowZdirInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%rho) ) &
            deallocate(BCData(i)%rho, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%velx) ) &
            deallocate(BCData(i)%velx, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%vely) ) &
            deallocate(BCData(i)%vely, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%velz) ) &
            deallocate(BCData(i)%velz, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%ps) ) &
            deallocate(BCData(i)%ps, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%turbInlet) ) &
            deallocate(BCData(i)%turbInlet, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%normALE) ) &
            deallocate(BCData(i)%normALE, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.
       if( associated(BCData(i)%rFaceALE) ) &
            deallocate(BCData(i)%rFaceALE, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.
       if( associated(BCData(i)%uSlipALE) ) &
            deallocate(BCData(i)%uSlipALE, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.
       if( associated(BCData(i)%sHeatFlux) ) &
            deallocate(BCData(i)%sHeatFlux, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(BCData(i)%iBlank) ) &
            deallocate(BCData(i)%iBlank, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       nullify(BCData(i)%norm)
       nullify(BCData(i)%rface)
       nullify(BCData(i)%F)
       nullify(BCData(i)%Fv)
       nullify(BCData(i)%Fp)
       nullify(BCData(i)%T)
       nullify(BCData(i)%Tv)
       nullify(BCData(i)%Tp)

       nullify(BCData(i)%uSlip)
       nullify(BCData(i)%TNS_Wall)

       nullify(BCData(i)%normALE)
       nullify(BCData(i)%rfaceALE)
       nullify(BCData(i)%uSlipALE)
       nullify(BCData(i)%sHeatFlux)

       nullify(BCData(i)%ptInlet)
       nullify(BCData(i)%ttInlet)
       nullify(BCData(i)%htInlet)
       nullify(BCData(i)%flowXdirInlet)
       nullify(BCData(i)%flowYdirInlet)
       nullify(BCData(i)%flowZdirInlet)

       nullify(BCData(i)%turbInlet)

       nullify(BCData(i)%rho)
       nullify(BCData(i)%velx)
       nullify(BCData(i)%vely)
       nullify(BCData(i)%velz)
       nullify(BCData(i)%ps)
       nullify(BCData(i)%iblank)

    enddo

    if( associated(flowDoms(nn,level,sps)%BCType) ) &
         deallocate(flowDoms(nn,level,sps)%BCType, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%BCFaceID) ) &
         deallocate(flowDoms(nn,level,sps)%BCFaceID, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%cgnsSubface) ) &
         deallocate(flowDoms(nn,level,sps)%cgnsSubface, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%inBeg) ) &
         deallocate(flowDoms(nn,level,sps)%inBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%inEnd) ) &
         deallocate(flowDoms(nn,level,sps)%inEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%jnBeg) ) &
         deallocate(flowDoms(nn,level,sps)%jnBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%jnEnd) ) &
         deallocate(flowDoms(nn,level,sps)%jnEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%knBeg) ) &
         deallocate(flowDoms(nn,level,sps)%knBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%knEnd) ) &
         deallocate(flowDoms(nn,level,sps)%knEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dinBeg) ) &
         deallocate(flowDoms(nn,level,sps)%dinBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dinEnd) ) &
         deallocate(flowDoms(nn,level,sps)%dinEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%djnBeg) ) &
         deallocate(flowDoms(nn,level,sps)%djnBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%djnEnd) ) &
         deallocate(flowDoms(nn,level,sps)%djnEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dknBeg) ) &
         deallocate(flowDoms(nn,level,sps)%dknBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dknEnd) ) &
         deallocate(flowDoms(nn,level,sps)%dknEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%icBeg) ) &
         deallocate(flowDoms(nn,level,sps)%icBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%icEnd) ) &
         deallocate(flowDoms(nn,level,sps)%icEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%jcBeg) ) &
         deallocate(flowDoms(nn,level,sps)%jcBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%jcEnd) ) &
         deallocate(flowDoms(nn,level,sps)%jcEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%kcBeg) ) &
         deallocate(flowDoms(nn,level,sps)%kcBeg, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%kcEnd) ) &
         deallocate(flowDoms(nn,level,sps)%kcEnd, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%neighBlock) ) &
         deallocate(flowDoms(nn,level,sps)%neighBlock, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%neighProc) ) &
         deallocate(flowDoms(nn,level,sps)%neighProc, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%l1) ) &
         deallocate(flowDoms(nn,level,sps)%l1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%l2) ) &
         deallocate(flowDoms(nn,level,sps)%l2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%l3) ) &
         deallocate(flowDoms(nn,level,sps)%l3, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%groupNum) ) &
         deallocate(flowDoms(nn,level,sps)%groupNum, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%iblank) ) &
         deallocate(flowDoms(nn,level,sps)%iblank, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%BCData) ) &
         deallocate(flowDoms(nn,level,sps)%BCData, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscSubface) ) &
         deallocate(flowDoms(nn,level,sps)%viscSubface, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.


    if( associated(flowDoms(nn,level,sps)%viscIminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscIminPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscImaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscImaxPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscJminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscJminPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscJmaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscJmaxPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscKminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscKminPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%viscKmaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscKmaxPointer, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%x) ) &
         deallocate(flowDoms(nn,level,sps)%x, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%xOld) ) &
         deallocate(flowDoms(nn,level,sps)%xOld, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%si) ) &
         deallocate(flowDoms(nn,level,sps)%si, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sj) ) &
         deallocate(flowDoms(nn,level,sps)%sj, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sk) ) &
         deallocate(flowDoms(nn,level,sps)%sk, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%vol) ) &
         deallocate(flowDoms(nn,level,sps)%vol, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%volRef) ) &
         deallocate(flowDoms(nn,level,sps)%volRef, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%volOld) ) &
         deallocate(flowDoms(nn,level,sps)%volOld, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%pori) ) &
         deallocate(flowDoms(nn,level,sps)%pori, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%porj) ) &
         deallocate(flowDoms(nn,level,sps)%porj, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%pork) ) &
         deallocate(flowDoms(nn,level,sps)%pork, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%indFamilyI) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyI, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%indFamilyJ) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyJ, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%indFamilyK) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%factFamilyI) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyI, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%factFamilyJ) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyJ, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%factFamilyK) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%rotMatrixI) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixI, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%rotMatrixJ) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixJ, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%rotMatrixK) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceI) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceI, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceJ) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceJ, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceK) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%w) ) &
         deallocate(flowDoms(nn,level,sps)%w, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wOld) ) &
         deallocate(flowDoms(nn,level,sps)%wOld, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%p) ) &
         deallocate(flowDoms(nn,level,sps)%p, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%aa) ) &
         deallocate(flowDoms(nn,level,sps)%aa, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%gamma) ) &
         deallocate(flowDoms(nn,level,sps)%gamma, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%ux) ) &
         deallocate(flowDoms(nn,level,sps)%ux, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%uy) ) &
         deallocate(flowDoms(nn,level,sps)%uy, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%uz) ) &
         deallocate(flowDoms(nn,level,sps)%uz, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%vx) ) &
         deallocate(flowDoms(nn,level,sps)%vx, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%vy) ) &
         deallocate(flowDoms(nn,level,sps)%vy, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%vz) ) &
         deallocate(flowDoms(nn,level,sps)%vz, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wx) ) &
         deallocate(flowDoms(nn,level,sps)%wx, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wy) ) &
         deallocate(flowDoms(nn,level,sps)%wy, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wz) ) &
         deallocate(flowDoms(nn,level,sps)%wz, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%qx) ) &
         deallocate(flowDoms(nn,level,sps)%qx, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%qy) ) &
         deallocate(flowDoms(nn,level,sps)%qy, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%qz) ) &
         deallocate(flowDoms(nn,level,sps)%qz, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%rlv) ) &
         deallocate(flowDoms(nn,level,sps)%rlv, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%rev) ) &
         deallocate(flowDoms(nn,level,sps)%rev, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%s) ) &
         deallocate(flowDoms(nn,level,sps)%s, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%p1) ) &
         deallocate(flowDoms(nn,level,sps)%p1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dw) ) &
         deallocate(flowDoms(nn,level,sps)%dw, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%fw) ) &
         deallocate(flowDoms(nn,level,sps)%fw, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dwOldRK) ) &
         deallocate(flowDoms(nn,level,sps)%dwOldRK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%w1) ) &
         deallocate(flowDoms(nn,level,sps)%w1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wr) ) &
         deallocate(flowDoms(nn,level,sps)%wr, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgIFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgIFine, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgJFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgJFine, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgKFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgKFine, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgIWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgIWeight, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgJWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgJWeight, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgKWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgKWeight, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgICoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgICoarse, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgJCoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgJCoarse, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%mgKCoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgKCoarse, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%iCo) ) &
         deallocate(flowDoms(nn,level,sps)%iCo, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%jCo) ) &
         deallocate(flowDoms(nn,level,sps)%jCo, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%kCo) ) &
         deallocate(flowDoms(nn,level,sps)%kCo, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%wn) ) &
         deallocate(flowDoms(nn,level,sps)%wn, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%pn) ) &
         deallocate(flowDoms(nn,level,sps)%pn, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dtl) ) &
         deallocate(flowDoms(nn,level,sps)%dtl, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%radI) ) &
         deallocate(flowDoms(nn,level,sps)%radI, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%radJ) ) &
         deallocate(flowDoms(nn,level,sps)%radJ, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%radK) ) &
         deallocate(flowDoms(nn,level,sps)%radK, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.


    if( associated(flowDoms(nn,level,sps)%d2Wall) ) &
         deallocate(flowDoms(nn,level,sps)%d2Wall, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.


    if( associated(flowDoms(nn,level,sps)%bmti1) ) &
         deallocate(flowDoms(nn,level,sps)%bmti1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bmti2) ) &
         deallocate(flowDoms(nn,level,sps)%bmti2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bmtj1) ) &
         deallocate(flowDoms(nn,level,sps)%bmtj1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bmtj2) ) &
         deallocate(flowDoms(nn,level,sps)%bmtj2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bmtk1) ) &
         deallocate(flowDoms(nn,level,sps)%bmtk1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bmtk2) ) &
         deallocate(flowDoms(nn,level,sps)%bmtk2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.


    if( associated(flowDoms(nn,level,sps)%bvti1) ) &
         deallocate(flowDoms(nn,level,sps)%bvti1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bvti2) ) &
         deallocate(flowDoms(nn,level,sps)%bvti2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bvtj1) ) &
         deallocate(flowDoms(nn,level,sps)%bvtj1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bvtj2) ) &
         deallocate(flowDoms(nn,level,sps)%bvtj2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bvtk1) ) &
         deallocate(flowDoms(nn,level,sps)%bvtk1, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%bvtk2) ) &
         deallocate(flowDoms(nn,level,sps)%bvtk2, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%globalCell) ) &
         deallocate(flowDoms(nn,level,sps)%globalCell, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%globalNode) ) &
         deallocate(flowDoms(nn,level,sps)%globalNode, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.


    ! *******************************
    ! Added by HDN
    ! *******************************
    if( associated(flowDoms(nn,level,sps)%xALE) ) &
         deallocate(flowDoms(nn,level,sps)%xALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sIALE) ) &
         deallocate(flowDoms(nn,level,sps)%sIALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sJALE) ) &
         deallocate(flowDoms(nn,level,sps)%sJALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sKALE) ) &
         deallocate(flowDoms(nn,level,sps)%sKALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sVeloIALE) ) &
         deallocate(flowDoms(nn,level,sps)%sVeloIALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sVeloJALE) ) &
         deallocate(flowDoms(nn,level,sps)%sVeloJALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sVeloKALE) ) &
         deallocate(flowDoms(nn,level,sps)%sVeloKALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceIALE) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceIALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceJALE) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceJALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%sFaceKALE) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceKALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%dwALE) ) &
         deallocate(flowDoms(nn,level,sps)%dwALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.

    if( associated(flowDoms(nn,level,sps)%fwALE) ) &
         deallocate(flowDoms(nn,level,sps)%fwALE, stat=ierr)
    if(ierr /= 0) deallocationFailure = .true.




    ! Check for errors in the deallocation.

    if( deallocationFailure ) &
         call returnFail("deallocateBlock", &
         "Something went wrong when deallocating memory")

    ! Nullify the pointers of this block.
    call nullifyFlowDomPointers(nn,level,sps)

  end subroutine deallocateBlock

  subroutine setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
       rev1, rev2, offset)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setBCPointers sets the pointers needed for the boundary        *
    !      * condition treatment on a general face, such that the boundary  *
    !      * routines are only implemented once instead of 6 times.         *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn, offset

    real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
    real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
    real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
    real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
    !
    !      Local variables
    !
    integer(kind=intType) :: id, ih

    ! Determine the face id on which the subface is located and set
    ! the pointers accordinly.

    select case (BCFaceID(nn))

    case (iMin)

       id = 2 + offset;   ih = 1 - offset
       ww1 => w(ih,1:,1:,:); ww2 => w(id,1:,1:,:)
       pp1 => p(ih,1:,1:);   pp2 => p(id,1:,1:)

       if( viscous ) then
          rlv1 => rlv(ih,1:,1:); rlv2 => rlv(id,1:,1:)
       endif

       if( eddyModel ) then
          rev1 => rev(ih,1:,1:); rev2 => rev(id,1:,1:)
       endif

       !===============================================================

    case (iMax)

       id = il - offset;  ih = ie + offset
       ww1 => w(ih,1:,1:,:); ww2 => w(id,1:,1:,:)
       pp1 => p(ih,1:,1:);   pp2 => p(id,1:,1:)

       if( viscous ) then
          rlv1 => rlv(ih,1:,1:); rlv2 => rlv(id,1:,1:)
       endif

       if( eddyModel ) then
          rev1 => rev(ih,1:,1:); rev2 => rev(id,1:,1:)
       endif

       !===============================================================

    case (jMin)

       id = 2 + offset;   ih = 1 - offset
       ww1 => w(1:,ih,1:,:); ww2 => w(1:,id,1:,:)
       pp1 => p(1:,ih,1:);   pp2 => p(1:,id,1:)

       if( viscous ) then
          rlv1 => rlv(1:,ih,1:); rlv2 => rlv(1:,id,1:)
       endif

       if( eddyModel ) then
          rev1 => rev(1:,ih,1:); rev2 => rev(1:,id,1:)
       endif

       !===============================================================

    case (jMax)

       id = jl - offset;  ih = je + offset
       ww1 => w(1:,ih,1:,:); ww2 => w(1:,id,1:,:)
       pp1 => p(1:,ih,1:);   pp2 => p(1:,id,1:)

       if( viscous ) then
          rlv1 => rlv(1:,ih,1:); rlv2 => rlv(1:,id,1:)
       endif

       if( eddyModel ) then
          rev1 => rev(1:,ih,1:); rev2 => rev(1:,id,1:)
       endif

       !===============================================================

    case (kMin)

       id = 2 + offset;   ih = 1 - offset
       ww1 => w(1:,1:,ih,:); ww2 => w(1:,1:,id,:)
       pp1 => p(1:,1:,ih);   pp2 => p(1:,1:,id)

       if( viscous ) then
          rlv1 => rlv(1:,1:,ih); rlv2 => rlv(1:,1:,id)
       endif

       if( eddyModel ) then
          rev1 => rev(1:,1:,ih); rev2 => rev(1:,1:,id)
       endif

       !===============================================================

    case (kMax)

       id = kl - offset;  ih = ke + offset
       ww1 => w(1:,1:,ih,:); ww2 => w(1:,1:,id,:)
       pp1 => p(1:,1:,ih);   pp2 => p(1:,1:,id)

       if( viscous ) then
          rlv1 => rlv(1:,1:,ih); rlv2 => rlv(1:,1:,id)
       endif

       if( eddyModel ) then
          rev1 => rev(1:,1:,ih); rev2 => rev(1:,1:,id)
       endif

    end select

  end subroutine setBCPointers

  subroutine setBCPointersBwd(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
       rev1, rev2, offset)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setBCPointers sets the pointers needed for the boundary        *
    !      * condition treatment on a general face, such that the boundary  *
    !      * routines are only implemented once instead of 6 times.         *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn, offset

    real(kind=realType), dimension(imaxDim,jmaxDim,nw) :: ww1, ww2
    real(kind=realType), dimension(imaxDim,jmaxDim) :: pp1, pp2
    real(kind=realType), dimension(imaxDim,jmaxDim) :: rlv1, rlv2
    real(kind=realType), dimension(imaxDim,jmaxDim) :: rev1, rev2
    !
    !      Local variables
    !
    integer(kind=intType) :: id, ih, ierr, i, j, k

    ww1 = zero
    ww2 = zero
    pp1 = zero
    pp2 = zero
    rlv1 = zero
    rlv2 = zero
    rev1 = zero
    rev2 = zero
    ! Determine the face id on which the subface is located and set
    ! the pointers accordinly.

    select case (BCFaceID(nn))

    case (iMin)

       id = 2 + offset;   ih = 1 - offset
       do k = 1,ke
          do j = 1,je
             ww1(j,k,:) = w(ih,j,k,:)
             ww2(j,k,:) = w(id,j,k,:)
             pp1(j,k)   = p(ih,j,k)
             pp2(j,k)   = p(id,j,k)

             if( viscous ) then
                rlv1(j,k) = rlv(ih,j,k)
                rlv2(j,k) = rlv(id,j,k)
             endif

             if( eddyModel ) then
                rev1(j,k) = rev(ih,j,k)
                rev2(j,k) = rev(id,j,k)
             endif
          end do
       end do
       !===============================================================

    case (iMax)

       id = il - offset;  ih = ie + offset
       do k = 1,ke
          do j = 1,je
             ww1(j,k,:) = w(ih,j,k,:)
             ww2(j,k,:) = w(id,j,k,:)
             pp1(j,k)   = p(ih,j,k)
             pp2(j,k)   = p(id,j,k)

             if( viscous ) then
                rlv1(j,k) = rlv(ih,j,k)
                rlv2(j,k) = rlv(id,j,k)
             endif

             if( eddyModel ) then
                rev1(j,k) = rev(ih,j,k)
                rev2(j,k) = rev(id,j,k)
             endif
          end do
       end do
       !===============================================================

    case (jMin)

       id = 2 + offset;   ih = 1 - offset
       do k = 1,ke
          do i = 1,ie
             ww1(i,k,:) = w(i,ih,k,:)
             ww2(i,k,:) = w(i,id,k,:)
             pp1(i,k)   = p(i,ih,k)
             pp2(i,k)   = p(i,id,k)

             if( viscous ) then
                rlv1(i,k) = rlv(i,ih,k)
                rlv2(i,k) = rlv(i,id,k)
             endif

             if( eddyModel ) then
                rev1(i,k) = rev(i,ih,k)
                rev2(i,k) = rev(i,id,k)
             endif
          end do
       end do
       !===============================================================

    case (jMax)

       id = jl - offset;  ih = je + offset
       do k = 1,ke
          do i = 1,ie
             ww1(i,k,:) = w(i,ih,k,:)
             ww2(i,k,:) = w(i,id,k,:)
             pp1(i,k)   = p(i,ih,k)
             pp2(i,k)   = p(i,id,k)

             if( viscous ) then
                rlv1(i,k) = rlv(i,ih,k)
                rlv2(i,k) = rlv(i,id,k)
             endif

             if( eddyModel ) then
                rev1(i,k) = rev(i,ih,k)
                rev2(i,k) = rev(i,id,k)
             endif
          end do
       end do
       !===============================================================

    case (kMin)

       id = 2 + offset;   ih = 1 - offset
       do j = 1,je
          do i = 1,ie
             ww1(i,j,:) = w(i,j,ih,:)
             ww2(i,j,:) = w(i,j,id,:)
             pp1(i,j)   = p(i,j,ih)
             pp2(i,j)   = p(i,j,id)

             if( viscous ) then
                rlv1(i,j) = rlv(i,j,ih)
                rlv2(i,j) = rlv(i,j,id)
             endif

             if( eddyModel ) then
                rev1(i,j) = rev(i,j,ih)
                rev2(i,j) = rev(i,j,id)
             endif
          end do
       end do
       !===============================================================

    case (kMax)

       id = kl - offset;  ih = ke + offset
       do j = 1,je
          do i = 1,ie
             ww1(i,j,:) = w(i,j,ih,:)
             ww2(i,j,:) = w(i,j,id,:)
             pp1(i,j)   = p(i,j,ih)
             pp2(i,j)   = p(i,j,id)

             if( viscous ) then
                rlv1(i,j) = rlv(i,j,ih)
                rlv2(i,j) = rlv(i,j,id)
             endif

             if( eddyModel ) then
                rev1(i,j) = rev(i,j,ih)
                rev2(i,j) = rev(i,j,id)
             endif
          end do
       end do
    end select

  end subroutine setBCPointersBwd

  integer function setCGNSRealType()
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setCGNSRealType sets the cgns real type, depending on the      *
    !      * compiler options. Note that quadrupole precision is not        *
    !      * supported by CGNS; double precision is used instead for the    *
    !      * CGNS IO.                                                       *
    !      *                                                                *
    !      ******************************************************************
    !
    use su_cgns
    implicit none

#ifdef USE_NO_CGNS

    call returnFail("setCGNSRealType", &
         "Function should not be called if no cgns support &
         &is selected.")

#else

# ifdef USE_SINGLE_PRECISION
    setCGNSRealType = RealSingle
# else
    setCGNSRealType = RealDouble
# endif

#endif

  end function setCGNSRealType



  subroutine returnFail(routineName, errorMessage)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * returnFail writes an error message to standard output and      *
    !      * sets fail flags to be returned to python.                      *
    !      *                                                                *
    !      ******************************************************************
    !
    use precision
    use communication
    use constants
#ifndef USE_TAPENADE
    use killSignals
#endif
    implicit none
    !
    !      Subroutine arguments
    !
    character(len=*), intent(in) :: routineName
    character(len=*), intent(in) :: errorMessage
#ifndef USE_TAPENADE

    !
    !      Local parameter
    !
    integer, parameter :: maxCharLine = 55
    !
    !      Local variables
    !
    integer :: ierr, len, i2
    logical :: firstTime

    character(len=len_trim(errorMessage)) :: message
    character(len=8) :: integerString

    ! Copy the errorMessage into message. It is not possible to work
    ! with errorMessage directly, because it is modified in this
    ! routine. Sometimes a constant string is passed to this routine
    ! and some compilers simply fail then.

    message = errorMessage

    ! Print a nice error message. In case of a parallel executable
    ! also the processor id is printed.

    print "(a)", "#"
    print "(a)", "#--------------------------- !!! Error !!! &
         &----------------------------"

    write(integerString,"(i8)") myID
    integerString = adjustl(integerString)

    print "(2a)", "#* returnFail called by processor ", &
         trim(integerString)

    ! Write the header of the error message.

    print "(2a)", "#* Run-time error in procedure ", &
         trim(routineName)

    ! Loop to write the error message. If the message is too long it
    ! is split over several lines.

    firstTime = .true.
    do
       ! Determine the remaining error message to be written.
       ! If longer than the maximum number of characters allowed
       ! on a line, it is attempted to split the message.

       message = adjustl(message)
       len = len_trim(message)
       i2  = min(maxCharLine,len)

       if(i2 < len) i2 = index(message(:i2), " ", .true.) - 1
       if(i2 < 0)   i2 = index(message, " ") - 1
       if(i2 < 0)   i2 = len

       ! Write this part of the error message. If it is the first
       ! line of the message some additional stuff is printed.

       if( firstTime ) then
          print "(2a)", "#* Error message: ", &
               trim(message(:i2))
          firstTime = .false.
       else
          print "(2a)", "#*                ", &
               trim(message(:i2))
       endif

       ! Exit the loop if the entire message has been written.

       if(i2 == len) exit

       ! Adapt the string for the next part to be written.

       message = message(i2+1:)

    enddo

    ! Write the trailing message.

    print "(a)", "#*"
    print "(a)", "#------------------------------------------&
         &----------------------------"
    print "(a)", "#"

    ! Call abort and stop the program. This stop should be done in
    ! abort, but just to be sure.

    if (fromPython)then
       routineFailed=.True.
       fatalFail = .True.
    else
       call mpi_abort(SUmb_comm_world, 1, ierr)
       stop
    end if
#endif

  end subroutine returnFail


  subroutine EChk(ierr, file, line)

    ! Check if ierr that resulted from a petsc or MPI call is in fact an
    ! error. 
    use precision
    use communication
    implicit none

    integer(kind=intType),intent(in) :: ierr
    character*(*),intent(in) :: file
    integer(kind=intType),intent(in) :: line

    if (ierr == 0) then
       return ! No error, return immediately
    else
#ifndef USE_TAPENADE
#ifndef USE_COMPLEX
       print *,'================================================================='
       write(*,900) "PETSc or MPI Error. Error Code ",ierr,". Detected on Proc ",myid
       write(*,901) "Error at line: ",line," in file: ",file
       print *,'================================================================='
#else
       print *,'-----------------------------------------------------------------'
       write(*,900) "PETSc or MPI Error. Error Code ",ierr,". Detected on Proc ",myid
       write(*,901) "Error at line: ",line," in file: ",file
       print *,'-----------------------------------------------------------------'
#endif
       call MPI_Abort(sumb_comm_world,ierr)
       stop ! Just in case
#else
       stop
#endif
    end if

900 format(A,I2,A,I2)
901 format(A,I5,A,A)
  end subroutine EChk

  subroutine convertToLowerCase(string)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * convertToLowerCase converts the given string to lower case.    *
    !      *                                                                *
    !      ******************************************************************
    !
    implicit none
    !
    !      Subroutine arguments
    !
    character (len=*), intent(inout) :: string
    !
    !      Local variables
    !
    integer, parameter :: upperToLower = iachar("a") - iachar("A")

    integer :: i, lenString

    ! Determine the length of the given string and convert the upper
    ! case characters to lower case.

    lenString = len_trim(string)
    do i=1,lenString
       if("A" <= string(i:i) .and. string(i:i) <= "Z")    &
            string(i:i) = achar(iachar(string(i:i)) + upperToLower)
    enddo

  end subroutine convertToLowerCase
#endif
end module utils
