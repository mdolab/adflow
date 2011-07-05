!
!      ******************************************************************
!      *                                                                *
!      * File:          rotMatrixRigidBody.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-15-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
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
!
!      Function definition.
!
       real(kind=realType) :: rigidRotAngle
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
