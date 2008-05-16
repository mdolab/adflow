!
!     ******************************************************************
!     *                                                                *
!     * File:          getDirVector.f90                                *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 10-25-2005                                      *
!     * Last modified: 10-26-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine getDirVector(xb,yb,zb,alpha,beta,xw,yw,zw)
!
!     ******************************************************************
!     *                                                                *
!     * Convert the angle of attack and side slip angle to wind axes.  *
!     * The components of the wind direction vector (xw,yw,zw) are     *
!     * computed given the direction angles in radians and the body    *
!     * direction by performing two rotations on the original          *
!     * direction vector:                                              *
!     *   1) Rotation about the yb-axis: alpha clockwise (CW)          *
!     *      (xb,yb,zb) -> (x1,y1,z1)                                  *
!     *                                                                *
!     *   2) Rotation about the z1-axis: beta counter-clockwise (CCW)  *
!     *      (x1,y1,z1) -> (xw,yw,zw)                                  *
!     *                                                                *
!     *    input arguments:                                            *
!     *       alpha    = angle of attack in radians                    *
!     *       beta     = side slip angle in radians                    *
!     *    output arguments:                                           *
!     *       xw,yw,zw = unit wind vector in body axes                 *
!     *                                                                *
!     ******************************************************************
!
      use constants
      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), intent(in)  :: xb, yb, zb
      real(kind=realType), intent(in)  :: alpha, beta
      real(kind=realType), intent(out) :: xw, yw, zw
!
!     Local variables.
!
      real(kind=realType) :: rnorm, xbn, ybn, zbn, x1, y1, z1
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Normalize the input vector.

      rnorm = sqrt( xb**2 + yb**2 + zb**2 )
      xbn = xb / rnorm
      ybn = yb / rnorm
      zbn = zb / rnorm

!!$      ! Compute the wind direction vector.
!!$
!!$      ! 1) rotate alpha radians cw about y-axis
!!$      !    ( <=> rotate y-axis alpha radians ccw)
!!$
!!$      call vectorRotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
!!$
!!$      ! 2) rotate beta radians ccw about z-axis
!!$      !    ( <=> rotate z-axis -beta radians ccw)
!!$
!!$      call vectorRotation(xw, yw, zw, 3, -beta, x1, y1, z1)

      ! Compute the wind direction vector.Aerosurf axes different!!

      ! 1) rotate alpha radians cw about z-axis
      !    ( <=> rotate z-axis alpha radians ccw)

      call vectorRotation(x1, y1, z1, 3, -alpha, xbn, ybn, zbn)

      ! 2) rotate beta radians ccw about y-axis
      !    ( <=> rotate z-axis -beta radians ccw)

      call vectorRotation(xw, yw, zw, 2, -beta, x1, y1, z1)

      end subroutine getDirVector
