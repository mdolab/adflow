!
!     ******************************************************************
!     *                                                                *
!     * File:          vectorRotation.f90                              *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 06-23-2006                                      *
!     * Last modified: 07-28-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine vectorRotation(xp, yp, zp, iaxis, angle, x, y, z)
!
!     ****************************************************************
!     *                                                              *
!     * vectorRotation rotates a given vector with respect to a      *
!     * specified axis by a given angle.                             *
!     *                                                              *
!     *    input arguments:                                          *
!     *       iaxis      = rotation axis (1-x, 2-y, 3-z)             *
!     *       angle      = rotation angle (measured ccw in radians)  *
!     *       x, y, z    = coordinates in original system            *
!     *    output arguments:                                         *
!     *       xp, yp, zp = coordinates in rotated system             *
!     *                                                              *
!     ****************************************************************
!
      use precision
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: iaxis
      real(kind=realType), intent(in)   :: angle, x, y, z
      real(kind=realType), intent(out)  :: xp, yp, zp
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution                                                *
!     *                                                                *
!     ******************************************************************
!
      ! rotation about specified axis by specified angle

      select case(iaxis)

        ! rotation about the x-axis

        case(1)
          xp =        1.    * x +     0.     * y +     0.     * z
          yp =        0.    * x + cos(angle) * y + sin(angle) * z
          zp =        0.    * x - sin(angle) * y + cos(angle) * z

        ! rotation about the y-axis

        case(2)
          xp =   cos(angle) * x +     0.     * y - sin(angle) * z
          yp =        0.    * x +     1.     * y +     0.     * z
          zp =   sin(angle) * x +     0.     * y + cos(angle) * z

        ! rotation about the z-axis

        case(3)
          xp =   cos(angle) * x + sin(angle) * y +     0.     * z
          yp = - sin(angle) * x + cos(angle) * y +     0.     * z
          zp =       0.     * x +     0.     * y +     1.     * z

        case default
          write(*,*) "vectorRotation called with invalid arguments"
          stop

      end select

      end subroutine vectorRotation
