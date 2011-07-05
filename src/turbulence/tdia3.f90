!
!      ******************************************************************
!      *                                                                *
!      * File:          tdia3.f90                                       *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-30-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine tdia3(nb, ne, l, c, u, r)
!
!      ******************************************************************
!      *                                                                *
!      * tdia3 solves the tridiagonal linear system (l+c+u) v = r,      *
!      * where l is the lower, c the central and u the upper diagonal.  *
!      * Every entry in the matrix is a 2x2 block matrix, i.e.          *
!      *                                                                *
!      *              x x  x 0  0 0  ........        = c(nb) u(nb)      *
!      *              x x  0 x  0 0  ........                           *
!      *                                                                *
!      *              x 0  x x  x 0  ........        = l(i) c(i) u(i)   *
!      *              0 x  x x  0 x  ........                           *
!      *                                                                *
!      *                   ........  x 0  x x        = l(ne) c(ne)      *
!      *                   ........  0 x  x x                           *
!      *                                                                *
!      *         With c = x x     u,l = x 0                             *
!      *                  x x           0 x                             *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nb, ne

       real(kind=realType), dimension(2,nb:ne), intent(inout) :: l, u, r
       real(kind=realType), dimension(2,2,nb:ne), intent(inout) :: c
!
!      Local variables.
!
       integer(kind=intType) :: n
       real(kind=realType)   :: deti, f11, f12, f21, f22, r1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Perform the backward sweep to eliMinate the upper diagonal uu.
       ! f     = u(n)*c^-1(n+1),
       ! c'(n) = c(n) - f*l(n+1)
       ! r'(n) = r(n) - f*r(n+1)

       do n=ne-1,nb,-1

         deti =  one/(c(1,1,n+1)*c(2,2,n+1) - c(1,2,n+1)*c(2,1,n+1))
         f11  =  u(1,n)*c(2,2,n+1)*deti
         f12  = -u(1,n)*c(1,2,n+1)*deti
         f21  = -u(2,n)*c(2,1,n+1)*deti
         f22  =  u(2,n)*c(1,1,n+1)*deti

         c(1,1,n) = c(1,1,n) - f11*l(1,n+1)
         c(1,2,n) = c(1,2,n) - f12*l(2,n+1)
         c(2,1,n) = c(2,1,n) - f21*l(1,n+1)
         c(2,2,n) = c(2,2,n) - f22*l(2,n+1)

         r(1,n)  = r(1,n) - f11*r(1,n+1) - f12*r(2,n+1)
         r(2,n)  = r(2,n) - f21*r(1,n+1) - f22*r(2,n+1)

       enddo

       ! The matrix is now in low block bi-diagonal form and can thus
       ! be solved be a forward sweep. The solution is stored in r.

       deti    =  one/(c(1,1,nb)*c(2,2,nb) - c(1,2,nb)*c(2,1,nb))
       r1      =  r(1,nb)
       r(1,nb) =  deti*(c(2,2,nb)*r1 - c(1,2,nb)*r(2,nb))
       r(2,nb) = -deti*(c(2,1,nb)*r1 - c(1,1,nb)*r(2,nb))

       do n=nb+1,ne

         r(1,n) = r(1,n) - l(1,n)*r(1,n-1)
         r(2,n) = r(2,n) - l(2,n)*r(2,n-1)

         deti   =  one/(c(1,1,n)*c(2,2,n) - c(1,2,n)*c(2,1,n))
         r1     =  r(1,n)
         r(1,n) =  deti*(c(2,2,n)*r1 - c(1,2,n)*r(2,n))
         r(2,n) = -deti*(c(2,1,n)*r1 - c(1,1,n)*r(2,n))

       enddo

       end subroutine tdia3
