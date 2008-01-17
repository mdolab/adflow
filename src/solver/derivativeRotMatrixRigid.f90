!
!      ******************************************************************
!      *                                                                *
!      * File:          derivativeRotMatrixRigid.f90                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-01-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine derivativeRotMatrixRigid(rotationMatrix, &
                                           rotationPoint, t)
!
!      ******************************************************************
!      *                                                                *
!      * derivativeRotMatrixRigid determines the derivative of the      *
!      * rotation matrix at the given time for the rigid body rotation, *
!      * such that the grid velocities can be determined analytically.  *
!      * Also the rotation point of the current time level is           *
!      * determined. This value can change due to translation of the    *
!      * entire grid.                                                   *
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
       real(kind=realType), intent(in) :: t

       real(kind=realType), dimension(3),   intent(out) :: rotationPoint
       real(kind=realType), dimension(3,3), intent(out) :: rotationMatrix
!
!      Local variables.
!
       integer(kind=intType) :: i, j

       real(kind=realType) :: phi, dphiX, dphiY, dphiZ
       real(kind=realType) :: cosX, cosY, cosZ, sinX, sinY, sinZ

       real(kind=realType), dimension(3,3) :: dm, m
!
!      Function definitions.
!
       real(kind=realType) :: rigidRotAngle
       real(kind=realType) :: derivativeRigidRotAngle
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the rotation angle around the x-axis for the new
       ! time level and the corresponding values of the sine and cosine.

       phi = rigidRotAngle(degreePolXRot,   coefPolXRot,      &
                           degreeFourXRot,  omegaFourXRot,    &
                           cosCoefFourXRot, sinCoefFourXRot, t)
       sinX = sin(phi)
       cosX = cos(phi)

       ! Idem for the y-axis.

       phi = rigidRotAngle(degreePolYRot,   coefPolYRot,      &
                           degreeFourYRot,  omegaFourYRot,    &
                           cosCoefFourYRot, sinCoefFourYRot, t)
       sinY = sin(phi)
       cosY = cos(phi)

       ! Idem for the z-axis.

       phi = rigidRotAngle(degreePolZRot,   coefPolZRot,      &
                           degreeFourZRot,  omegaFourZRot,    &
                           cosCoefFourZRot, sinCoefFourZRot, t)
       sinZ = sin(phi)
       cosZ = cos(phi)

       ! Compute the time derivative of the rotation angles around the
       ! x-axis, y-axis and z-axis.

       dphiX = derivativeRigidRotAngle(degreePolXRot,   &
                                       coefPolXRot,     &
                                       degreeFourXRot,  &
                                       omegaFourXRot,   &
                                       cosCoefFourXRot, &
                                       sinCoefFourXRot, t)

       dphiY = derivativeRigidRotAngle(degreePolYRot,   &
                                       coefPolYRot,     &
                                       degreeFourYRot,  &
                                       omegaFourYRot,   &
                                       cosCoefFourYRot, &
                                       sinCoefFourYRot, t)

       dphiZ = derivativeRigidRotAngle(degreePolZRot,   &
                                       coefPolZRot,     &
                                       degreeFourZRot,  &
                                       omegaFourZRot,   &
                                       cosCoefFourZRot, &
                                       sinCoefFourZRot, t)

       ! Compute the time derivative of the rotation matrix applied to
       ! the coordinates at t == 0.

       ! Part 1. Derivative of the z-rotation matrix multiplied by the
       ! x and y rotation matrix, i.e. dmz * my * mx

       dm(1,1) = -cosY*sinZ*dphiZ
       dm(1,2) = (-cosX*cosZ - sinX*sinY*sinZ)*dphiZ
       dm(1,3) = ( sinX*cosZ - cosX*sinY*sinZ)*dphiZ

       dm(2,1) = cosY*cosZ*dphiZ
       dm(2,2) = (sinX*sinY*cosZ - cosX*sinZ)*dphiZ
       dm(2,3) = (cosX*sinY*cosZ + sinX*sinZ)*dphiZ

       dm(3,1) = zero
       dm(3,2) = zero
       dm(3,3) = zero

       ! Part 2: mz * dmy * mx.

       dm(1,1) = dm(1,1) - sinY*cosZ*dphiY
       dm(1,2) = dm(1,2) + sinX*cosY*cosZ*dphiY
       dm(1,3) = dm(1,3) + cosX*cosY*cosZ*dphiY

       dm(2,1) = dm(2,1) - sinY*sinZ*dphiY
       dm(2,2) = dm(2,2) + sinX*cosY*sinZ*dphiY
       dm(2,3) = dm(2,3) + cosX*cosY*sinZ*dphiY

       dm(3,1) = dm(3,1) - cosY*dphiY
       dm(3,2) = dm(3,2) - sinX*sinY*dphiY
       dm(3,3) = dm(3,3) - cosX*sinY*dphiY

       ! Part 3: mz * my * dmx

       dm(1,2) = dm(1,2) + (sinX*sinZ + cosX*sinY*cosZ)*dphiX
       dm(1,3) = dm(1,3) + (cosX*sinZ - sinX*sinY*cosZ)*dphiX

       dm(2,2) = dm(2,2) + (cosX*sinY*sinZ - sinX*cosZ)*dphiX
       dm(2,3) = dm(2,3) - (sinX*sinY*sinZ + cosX*cosZ)*dphiX

       dm(3,2) = dm(3,2) + cosX*cosY*dphiX
       dm(3,3) = dm(3,3) - sinX*cosY*dphiX

       ! Determine the rotation matrix at t == t.

       m(1,1) =  cosY*cosZ
       m(2,1) =  cosY*sinZ
       m(3,1) = -sinY

       m(1,2) = sinX*sinY*cosZ - cosX*sinZ
       m(2,2) = sinX*sinY*sinZ + cosX*cosZ
       m(3,2) = sinX*cosY

       m(1,3) = cosX*sinY*cosZ + sinX*sinZ
       m(2,3) = cosX*sinY*sinZ - sinX*cosZ
       m(3,3) = cosX*cosY

       ! Determine the matrix product dm * inverse(m), which will give
       ! the derivative of the rotation matrix when applied to the
       ! current coordinates. Note that inverse(m) == transpose(m).

       do j=1,3
         do i=1,3
           rotationMatrix(i,j) = dm(i,1)*m(j,1) + dm(i,2)*m(j,2) &
                               + dm(i,3)*m(j,3)
         enddo
       enddo

       ! Determine the rotation point at the new time level; it is
       ! possible that this value changes due to translation of the grid.

    !  aInf = sqrt(gammaInf*pInf/rhoInf)

    !  RotationPoint(1) = LRef*rotPoint(1) &
    !                    + MachGrid(1)*aInf*t/timeRef
    !  rotationPoint(2) = LRef*rotPoint(2) &
    !                    + MachGrid(2)*aInf*t/timeRef
    !  rotationPoint(3) = LRef*rotPoint(3) &
    !                    + MachGrid(3)*aInf*t/timeRef

       rotationPoint(1) = LRef*rotPoint(1)
       rotationPoint(2) = LRef*rotPoint(2)
       rotationPoint(3) = LRef*rotPoint(3)

       end subroutine derivativeRotMatrixRigid
