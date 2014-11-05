   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of adjustinflowangle in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: veldirfreestream dragdirection
   !                liftdirection alpha
   !   with respect to varying inputs: alpha beta
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          adjustInflowAngle.f90                           *
   !      * Author:        C.A.(Sandy) Mader                               *
   !      * Starting date: 07-13-2011                                      *
   !      * Last modified: 07-13-2011                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE ADJUSTINFLOWANGLE_B(alpha, alphab, beta, betab, liftindex)
   USE CONSTANTS
   USE INPUTPHYSICS
   IMPLICIT NONE
   !Subroutine Vars
   REAL(kind=realtype), INTENT(IN) :: alpha, beta
   REAL(kind=realtype) :: alphab, betab
   INTEGER(kind=inttype), INTENT(IN) :: liftindex
   !Local Vars
   REAL(kind=realtype), DIMENSION(3) :: refdirection
   ! Velocity direction given by the rotation of a unit vector
   ! initially aligned along the positive x-direction (1,0,0)
   ! 1) rotate alpha radians cw about y or z-axis
   ! 2) rotate beta radians ccw about z or y-axis
   refdirection(:) = zero
   refdirection(1) = one
   ! Drag direction given by the rotation of a unit vector
   ! initially aligned along the positive x-direction (1,0,0)
   ! 1) rotate alpha radians cw about y or z-axis
   ! 2) rotate beta radians ccw about z or y-axis
   CALL PUSHREAL8ARRAY(refdirection, 3)
   refdirection(:) = zero
   refdirection(1) = one
   ! Lift direction given by the rotation of a unit vector
   ! initially aligned along the positive z-direction (0,0,1)
   ! 1) rotate alpha radians cw about y or z-axis
   ! 2) rotate beta radians ccw about z or y-axis
   CALL PUSHREAL8ARRAY(refdirection, 3)
   refdirection(:) = zero
   refdirection(liftindex) = one
   betab = 0.0_8
   CALL GETDIRVECTOR_B(refdirection, alpha, alphab, beta, betab, &
   &               liftdirection, liftdirectionb, liftindex)
   CALL POPREAL8ARRAY(refdirection, 3)
   CALL GETDIRVECTOR_B(refdirection, alpha, alphab, beta, betab, &
   &               dragdirection, dragdirectionb, liftindex)
   CALL POPREAL8ARRAY(refdirection, 3)
   CALL GETDIRVECTOR_B(refdirection, alpha, alphab, beta, betab, &
   &               veldirfreestream, veldirfreestreamb, liftindex)
   END SUBROUTINE ADJUSTINFLOWANGLE_B
