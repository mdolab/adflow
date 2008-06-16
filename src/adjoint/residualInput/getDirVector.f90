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
      subroutine getDirVector(refDirection, alpha, beta,&
           windDirection,liftIndex)
        !(xb,yb,zb,alpha,beta,xw,yw,zw)
!
!     ******************************************************************
!     *                                                                *
!     * Convert the angle of attack and side slip angle to wind axes.  *
!     * The components of the wind direction vector (xw,yw,zw) are     *
!     * computed given the direction angles in radians and the body    *
!     * direction by performing two rotations on the original          *
!     * direction vector:                                              *
!     *   1) Rotation about the zb or yb-axis: alpha clockwise (CW)    *
!     *      (xb,yb,zb) -> (x1,y1,z1)                                  *
!     *                                                                *
!     *   2) Rotation about the yl or z1-axis: beta counter-clockwise  *
!     *      (CCW)  (x1,y1,z1) -> (xw,yw,zw)                           *
!     *                                                                *
!     *    input arguments:                                            *
!     *       alpha    = angle of attack in radians                    *
!     *       beta     = side slip angle in radians                    *
!     *       refDirection = reference direction vector                *
!     *    output arguments:                                           *
!     *       windDirection = unit wind vector in body axes            *
!     *                                                                *
!     ******************************************************************
!
      use constants
      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType),dimension(3), intent(in)  :: refDirection
      real(kind=realType), intent(in)  :: alpha, beta
      real(kind=realType),dimension (3), intent(out) :: windDirection
      integer(kind=intType)::liftIndex
!
!     Local variables.
!
      real(kind=realType) :: rnorm,x1,y1,z1,xbn,ybn,zbn,xw,yw,zw
      real(kind=realType),dimension (3)::refDirectionNorm
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Normalize the input vector.

      rnorm = sqrt( refDirection(1)**2 + refDirection(2)**2 + refDirection(3)**2 )
      xbn = refDirection(1) / rnorm
      ybn = refDirection(2) / rnorm
      zbn = refDirection(3) / rnorm

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

      if (liftIndex==2)then
         ! Compute the wind direction vector.Aerosurf axes different!!
         
         ! 1) rotate alpha radians cw about z-axis
         !    ( <=> rotate z-axis alpha radians ccw)
         
         
         call vectorRotation(x1, y1, z1, 3, -alpha, xbn, ybn, zbn)
         
         ! 2) rotate beta radians ccw about y-axis
         !    ( <=> rotate z-axis -beta radians ccw)
         
         call vectorRotation(xw, yw, zw, 2, -beta, x1, y1, z1)

      elseif(liftIndex==3)then
         ! Compute the wind direction vector.Aerosurf axes different!!
         
         ! 1) rotate alpha radians cw about z-axis
         !    ( <=> rotate z-axis alpha radians ccw)
         
         call vectorRotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
         
         ! 2) rotate beta radians ccw about y-axis
         !    ( <=> rotate z-axis -beta radians ccw)
         
         call vectorRotation(xw, yw, zw, 3, beta, x1, y1, z1)
         
         
      else
         call terminate('getDirVector', 'Invalid Lift Direction')
         
      endif
      
      windDirection(1) = xw
      windDirection(2) = yw
      windDirection(3) = zw
      
      end subroutine getDirVector
