!
!     ******************************************************************
!     *                                                                *
!     * File:          mind2hexa.f90                                   *
!     * Author:        Juan J. Alonso                                  *
!     * Starting date: 04-20-2006                                      *
!     * Last modified: 04-25-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine minD2Hexa(xP,x1,x2,x3,x4,x5,x6,x7,x8,d2,chi,iErr)
!
!     ******************************************************************
!     *                                                                *
!     * Subroutine to compute a fail-safe minimum distance computation *
!     * between a point, P, and a hexahedral element.  This subroutine *
!     * can provide the minimum distance whether P is inside or        *
!     * outside of the hexahedral element.  If P is inside the element *
!     * the distance returned is zero, while if P is outside of the    *
!     * element, the distance returned is the minimum distance between *
!     * P and any of the bounding faces of the hexahedral element.     *
!     *                                                                *
!     * The basic idea of this subroutine is to perform a              *
!     * bound-constrained minimization of the distance between the     *
!     * point P and a point inside the element given by parametric     *
!     * coordinates chi(3) using a modified Newton step bound          *
!     * constrained optimizer.                                         *
!     *                                                                *
!     * For more details of the optimization algorithms, a good        *
!     * reference is:                                                  *
!     *                                                                *
!     * Nocedal, J. and Wright, S. J., Numerical Optimization,         *
!     * Springer, 1999.                                                *
!     *                                                                *
!     * Subroutine arguments:                                          *
!     *                                                                *
!     * xp(3)        - The Cartesian coordinates of point P.           *
!     * x1(3)-x8(3)  - The Cartesian coordinates of the 8 nodes that   *
!     *                 make up the hexahedron, in the order specified *
!     *                 in the CHIMPS standard.                        *
!     * d2           - The squared of the minimum distance between the *
!     *                 point P and the hexahedral element.            *
!     * chi(3)       - Parametric coordinates of the point that        *
!     *                 belongs to the hexahedron where the minimum    *
!     *                 distance has been found.  If P is inside the   *
!     *                 element, 0<(ksi,eta,zeta)<1, while if P is     *
!     *                 strictly outside of the element, then one or   *
!     *                 more of the values of (ksi,eta,zeta) will be   *
!     *                 exactly zero or one.                           *
!     * iErr         - Output status of this subroutine:               *
!     *                 iErr =  0, Proper minimum found.               *
!     *                 iErr = -1, Distance minimization failed.       *
!     *                                                                *
!     ******************************************************************
!
      use precision

      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), dimension(3), intent(in) :: xP
      real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
      real(kind=realType), dimension(3), intent(in) :: x5, x6, x7, x8

      integer(kind=intType), intent(out) :: iErr

      real(kind=realType), intent(out) :: d2
      real(kind=realType), dimension(3), intent(out) :: chi
!
!     Local variables.
!
      integer(kind=intType) , parameter :: maxIt = 30
      integer(kind=intType) :: i, itCount = 0
      integer(kind=intType), dimension(3) :: actSet, chiGradConv

      real(kind=realType) :: inactGradNorm, normDeltaChi, x0, y0, z0
      real(kind=realType), parameter :: gradTol = 1.0e-14_realType
      real(kind=realType), parameter :: deltaChiTol = 1.0e-14_realType
      real(kind=realType), dimension(3) :: deltaChi, actualDeltaChi
      real(kind=realType), dimension(3) :: lwrBnd, uppBnd
      real(kind=realType), dimension(3) :: grad, oldChi
      real(kind=realType), dimension(3,3) :: hess

      logical :: convDeltaChi, convGradD2
!
!     ******************************************************************
!     *                                                                *
!     * Initialization section.                                        *
!     *                                                                *
!     ******************************************************************
!

      ! Setup initial values for the parametric coordinates of the
      ! minimum distance point.  One may be able to do better than this
      ! in the future, but this is probably a good guess.

      chi(1) = 0.5_realType
      chi(2) = 0.5_realType
      chi(3) = 0.5_realType

      ! Initialize the active set array, actSet(i), i=1,2,3
      ! actSet(i) =  0 if bound constraint i is inactive
      ! actSet(i) = +1 if bound constraint i is active at upper bound
      ! actSet(i) = -1 if bound constraint i is active at lower bound

      actSet(:) = 0

      ! Initialize the upper and lower bounds for all chi variables

      lwrBnd(:) = 0.0_realType
      uppBnd(:) = 1.0_realType

      ! Initialize actualDeltaChi (step size) to a large value so that the
      ! convergence criteria for the step size is guaranteed not to pass
      ! on the first iteration.

      actualDeltaChi(:) = 1.0_realType

!
!     ******************************************************************
!     *                                                                *
!     * Main iteration loop.                                           *
!     *                                                                *
!     ******************************************************************
!
      itCount = 0

      IterLoop: do while (itCount <= maxIt)

         ! Increment iteration counter

         itCount = itCount +1

         !
         !     *********************************************************
         !     *                                                       *
         !     * Compute the gradient of d2                            *
         !     *                                                       *
         !     *********************************************************
         !

         call gradD2Hexa(xP,x1,x2,x3,x4,x5,x6,x7,x8,chi,x0,y0,z0,grad,iErr)

         !
         !     *********************************************************
         !     *                                                       *
         !     * Convergence test                                      *
         !     *                                                       *
         !     *********************************************************
         !

         convDeltaChi = .false.
         convGradD2   = .false. 

         ! Convergence test for step size

         normDeltaChi = sqrt(actualDeltaChi(1)**2 +actualDeltaChi(2)**2 +actualDeltaChi(3)**2)

         if (normDeltaChi < deltaChiTol) convDeltaChi = .true.

         ! Convergence test for the gradient.  Note that this gradient
         ! test is such that, in order to pass it, the following must
         ! be true:
         !
         ! 1) \frac{\partial d2}{\partial \chi_i} > 0 for i in the active
         !     set at the lower bound, actSet(i) = -1
         ! 2) \frac{\partial d2}{\partial \chi_i} < 0 for i in the active
         !     set at the upper bound, actSet(i) = +1
         ! 3) norm of components of the gradient in the inactive set must
         !     be smaller than gradTol
         !

         inactGradNorm  = 0.0_realType
         chiGradConv(:) = 0

         do i=1,3

            ! If this component of chi is in the inactive set, accumulate
            ! the total value of the gradient and deal with it later.

            if (actSet(i) == 0) inactGradNorm = inactGradNorm +grad(i)**2

            ! If this component of chi is active at a lower bound, check
            ! for the convergence criterion.  Also deactivate the constraint
            ! if the gradient is pointing in the proper direction.

            if (actSet(i) == -1) then
               if (grad(i) > 0.0_realType) then
                  chiGradConv(i) = 1
               else
                  actSet(i) = 0
               end if
            end if

            ! If this component of chi is active at an upper bound, check
            ! for the convergence criterion.

            if (actSet(i) ==  1) then
               if (grad(i) < 0.0_realType) then
                  chiGradConv(i) = 1
               else
                  actSet(i) = 0
               end if
            end if

         end do

         ! Check for convergence on the accumulated values of the inactive
         ! components of the gradient.

         if (inactGradNorm < gradTol) then
            do i=1,3
               if (actSet(i) == 0) chiGradConv(i) = 1
            end do
         end if

         if (sum(chiGradConv(1:3)) == 3) convGradD2 = .true.

         ! Test for convergence using both criteria

         if (convDeltaChi .and. convGradD2) exit IterLoop
         
         !
         !     *********************************************************
         !     *                                                       *
         !     * Compute the Hessian of d2                             *
         !     *                                                       *
         !     *********************************************************
         !

         call hessD2Hexa(xP,x1,x2,x3,x4,x5,x6,x7,x8,chi,hess,iErr)
         
         !
         !     *********************************************************
         !     *                                                       *
         !     * Compute the Newton step                               *
         !     *                                                       *
         !     *********************************************************
         !

         call newtonStep(hess,grad,deltaChi,iErr)

         !
         !     *********************************************************
         !     *                                                       *
         !     * Update the current guess (appropriately clipped at    *
         !     * the bounds) and update the active set array as needed *
         !     *                                                       *
         !     *********************************************************
         !

         ! Loop over the components of chi

         oldChi(:) = chi(:)

         do i=1,3

            ! Update only the components of chi that were not active

            if (actSet(i) == 0) then

               chi(i) = chi(i) +deltaChi(i)

               ! Check to see this degree of freedom has become
               ! actively constrained at either an upper or
               ! lower bound and update active set.

               if (chi(i) > uppBnd(i)) then
                  chi(i)    = uppBnd(i)
                  actSet(i) = 1
               else if (chi(i) < lwrBnd(i)) then
                  chi(i)    = lwrBnd(i)
                  actSet(i) = -1
               end if
            end if
         end do

         actualDeltaChi(:) = chi(:) -oldChi(:)
      end do IterLoop

!
!     ******************************************************************
!     *                                                                *
!     * Return the results to the calling function and print info      *
!     * for debugging purposes.                                        *
!     *                                                                *
!     ******************************************************************
!

      ! Compute the minimum distance (squared) that the algorithm found

      d2 = (xP(1) -x0)**2 +(xP(2) -y0)**2 +(xP(3) -z0)**2

      ! Print some stuff out to the screen for debugging purposes

     !write(*,*)
     !write(*,*) 'Results of minD2Hexa'
     !write(*,20)'Point P                = (',xP(1),xP(2),xP(3),' )'
     !write(*,20)'Found point            = (',x0,y0,z0,' )'
     !write(*,20)'Parametric coordinates = (',chi(1),chi(2),chi(3),' )'
     !write(*,30)'Minimum distance       =',sqrt(d2)
     !write(*,10)'Number of iterations   =',itCount
      
10    format(a,1x,i3)
20    format(a,3f10.6,a)
30    format(a,f20.17,a)
      return

      end subroutine minD2Hexa
