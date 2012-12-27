!
!     ******************************************************************
!     *                                                                *
!     * File:          gradD2hexa.f90                                  *
!     * Author:        Juan J. Alonso                                  *
!     * Starting date: 04-20-2006                                      *
!     * Last modified: 04-25-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine gradD2Hexa(xP,x1,x2,x3,x4,x5,x6,x7,x8,chi,x0,y0,z0,grad,iErr)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the gradient of the square of the distance between the *
!     * point xP, and the actual point in the hexahedron represented   *
!     * by chi(1:3).                                                   *
!     *                                                                *
!     ******************************************************************
!
      use precision

      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), intent(out)  :: x0, y0, z0
      real(kind=realType), dimension(3), intent(in)  :: xP, chi
      real(kind=realType), dimension(3), intent(in)  :: x1, x2, x3, x4
      real(kind=realType), dimension(3), intent(in)  :: x5, x6, x7, x8
      real(kind=realType), dimension(3), intent(out) :: grad

      integer(kind=intType), intent(out) :: iErr

!
!     Local variables.
!
      integer(kind=intType) :: i

      real(kind=realType) :: ksi, eta, zeta
      real(kind=realType), dimension(8,3) :: alpha
      real(kind=realType) :: dxdksi, dxdeta, dxdzeta
      real(kind=realType) :: dydksi, dydeta, dydzeta
      real(kind=realType) :: dzdksi, dzdeta, dzdzeta

!
!     ******************************************************************
!     *                                                                *
!     * Compute the gradient components.                               *
!     *                                                                *
!     ******************************************************************
!

      ! Initialize the parametric coordinates with a more recognizable
      ! name.

      ksi  = chi(1)
      eta  = chi(2)
      zeta = chi(3)

      ! Initialize the alpha array (obtained by regrouping terms in the
      ! parametric expansion of the (ksi,eta,zeta)->(x,y,z) mapping.
      ! The second index of the alpha array corresponds to the x, y, or
      ! z coordinate mapping.

      do i=1,3
         alpha(1,i) = x1(i)
         alpha(2,i) = x2(i) -x1(i)
         alpha(3,i) = x4(i) -x1(i)
         alpha(4,i) = x5(i) -x1(i)
         alpha(5,i) = x3(i) -x2(i) +x1(i) -x4(i)
         alpha(6,i) = x6(i) -x5(i) +x1(i) -x2(i)
         alpha(7,i) = x8(i) -x5(i) +x1(i) -x4(i)
         alpha(8,i) = x7(i) -x8(i) +x5(i) -x6(i) +x4(i) -x3(i) +x2(i) -x1(i)
      end do

      ! Compute the value of the partial derivatives of x, y, and z with
      ! respect to ksi, eta, and zeta

      dxdksi = alpha(2,1) +alpha(5,1)*eta +alpha(6,1)*zeta +alpha(8,1)*eta*zeta
      dydksi = alpha(2,2) +alpha(5,2)*eta +alpha(6,2)*zeta +alpha(8,2)*eta*zeta
      dzdksi = alpha(2,3) +alpha(5,3)*eta +alpha(6,3)*zeta +alpha(8,3)*eta*zeta

      dxdeta = alpha(3,1) +alpha(5,1)*ksi +alpha(7,1)*zeta +alpha(8,1)*ksi*zeta
      dydeta = alpha(3,2) +alpha(5,2)*ksi +alpha(7,2)*zeta +alpha(8,2)*ksi*zeta
      dzdeta = alpha(3,3) +alpha(5,3)*ksi +alpha(7,3)*zeta +alpha(8,3)*ksi*zeta

      dxdzeta = alpha(4,1) +alpha(6,1)*ksi +alpha(7,1)*eta +alpha(8,1)*ksi*eta
      dydzeta = alpha(4,2) +alpha(6,2)*ksi +alpha(7,2)*eta +alpha(8,2)*ksi*eta
      dzdzeta = alpha(4,3) +alpha(6,3)*ksi +alpha(7,3)*eta +alpha(8,3)*ksi*eta

      ! Compute the actual location of the point

      x0 = alpha(1,1) +alpha(2,1)*ksi +alpha(3,1)*eta +alpha(4,1)*zeta            &
          +alpha(5,1)*ksi*eta +alpha(6,1)*ksi*zeta +alpha(7,1)*eta*zeta           &
          +alpha(8,1)*ksi*eta*zeta
      y0 = alpha(1,2) +alpha(2,2)*ksi +alpha(3,2)*eta +alpha(4,2)*zeta            &
          +alpha(5,2)*ksi*eta +alpha(6,2)*ksi*zeta +alpha(7,2)*eta*zeta           &
          +alpha(8,2)*ksi*eta*zeta
      z0 = alpha(1,3) +alpha(2,3)*ksi +alpha(3,3)*eta +alpha(4,3)*zeta            &
          +alpha(5,3)*ksi*eta +alpha(6,3)*ksi*zeta +alpha(7,3)*eta*zeta           &
          +alpha(8,3)*ksi*eta*zeta

      ! Compute the gradient components

      grad(1) = -2.0_realType*(xP(1) -x0)*dxdksi                               &
                -2.0_realType*(xP(2) -y0)*dydksi                               &
                -2.0_realType*(xP(3) -z0)*dzdksi
      grad(2) = -2.0_realType*(xP(1) -x0)*dxdeta                               &
                -2.0_realType*(xP(2) -y0)*dydeta                               &
                -2.0_realType*(xP(3) -z0)*dzdeta
      grad(3) = -2.0_realType*(xP(1) -x0)*dxdzeta                              &
                -2.0_realType*(xP(2) -y0)*dydzeta                              &
                -2.0_realType*(xP(3) -z0)*dzdzeta

      ! Return iErr = 0 for the time being.

      iErr = 0

      return

      end subroutine gradD2Hexa
