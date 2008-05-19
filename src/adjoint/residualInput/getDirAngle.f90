!
!     ******************************************************************
!     *                                                                *
!     * File:          getDirAngle.f90                                 *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 10-25-2005                                      *
!     * Last modified: 10-26-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine getDirAngle(xw,yw,zw,alpha,beta)
!
!     ******************************************************************
!     *                                                                *
!     * Convert the wind axes to angle of attack and side slip angle.  *
!     * The direction angles alpha and beta are computed given the     *
!     * components of the wind direction vector (xw,yw,zw) and the     *
!     * body direction (xb,yb,zb), where the former was obtained from  *
!     * the latter by performing two rotations on the original         *
!     * direction vector:                                              *
!     *   1) Rotation about the yb-axis: alpha clockwise (CW)          *
!     *      (xb,yb,zb) -> (x1,y1,z1)                                  *
!     *                                                                *
!     *   2) Rotation about the z1-axis: beta counter-clockwise (CCW)  *
!     *      (x1,y1,z1) -> (xw,yw,zw)                                  *
!     *                                                                *
!     *    input arguments:                                            *
!     *       xw,yw,zw = wind vector in body axes                      *
!     *    output arguments:                                           *
!     *       alpha    = angle of attack in radians                    *
!     *       beta     = side slip angle in radians                    *
!     *                                                                *
!     ******************************************************************
!
      use constants
      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), intent(in)  :: xw, yw, zw
      real(kind=realType), intent(out) :: alpha, beta
!
!     Local variables.
!
      real(kind=realType) :: rnorm, xwn, ywn, zwn
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Normalize the input vector.

      rnorm = sqrt( xw**2 + yw**2 + zw**2 )
      xwn = xw / rnorm
      ywn = yw / rnorm
      zwn = zw / rnorm

!!$      ! Compute angle of attack alpha.
!!$
!!$      alpha = asin(zwn)
!!$
!!$      ! Compute side-slip angle beta.
!!$
!!$      beta  = atan2(ywn,xwn)
      
      !different coordinate system for aerosurf
      ! Compute angle of attack alpha.
      
      alpha = asin(ywn)

      ! Compute side-slip angle beta.

      beta  = atan2(zwn,xwn)

      end subroutine getDirAngle
