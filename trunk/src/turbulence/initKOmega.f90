!
!      ******************************************************************
!      *                                                                *
!      * File:          initKOmega.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-01-2003                                      *
!      * Last modified: 07-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initKOmega(pOffset)
!
!      ******************************************************************
!      *                                                                *
!      * initKOmega initializes the values of k and omega a bit more    *
!      * intelligently than just free-stream values.                    *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block on the correct level.               *
!      * The argument pOffset is present such that in case of restart   *
!      * a possible pointer offset is taken into account. For more      *
!      * details see the corresponding routines in initFlow.            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use flowVarRefState
       use paramTurb
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: pOffset
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ip, jp, kp
       real(kind=realType)   :: rhoi, tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the owned cells of the block.

       do k=2,kl
         kp = k + pOffset
         do j=2,jl
           jp = j + pOffset
           do i=2,il
             ip = i + pOffset

             ! Compute a value of omega based on the wall distance.

             rhoi = one/w(ip,jp,kp,irho)
             tmp  = six*rhoi*rlv(i,j,k)/(rkwBeta1*(d2Wall(i,j,k)**2))

             ! Initialize omega using the value just computed; make sure
             ! that the free stream value is the lowest possible value.
             ! After that initialize k using the current value of omega
             ! and the eddy viscosity. Again the free-stream value is
             ! the lower bound.

             w(ip,jp,kp,itu2) = max(tmp,wInf(itu2))
             tmp              = rhoi*rev(i,j,k)*w(ip,jp,kp,itu2)
             w(ip,jp,kp,itu1) = max(tmp,wInf(itu1))

           enddo
         enddo
       enddo

       end subroutine initKOmega
