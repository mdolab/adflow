!      ==================================================================

         subroutine leftRightStateAdj(du1, du2, du3, left, right,nwInt,omk,opk,factminmod,firstOrderK)
!
!        ****************************************************************
!        *                                                              *
!        * leftRightState computes the differences in the left and      *
!        * right state compared to the first order interpolation. For a *
!        * monotonic second order discretization the interpolations     *
!        * need to be nonlinear. The linear second order scheme can be  *
!        * stable (depending on the value of kappa), but it will have   *
!        * oscillations near discontinuities.                           *
!        *                                                              *
!        ****************************************************************
!
         use precision
         use constants
         use inputDiscretization

         implicit none
!
!        Local parameter.
!
         real(kind=realType), parameter :: epsLim = 1.e-10_realType
!
!        Subroutine arguments.
!
         integer(kind=intType) :: nwInt
         real(kind=realType) :: omk, opk, factMinmod
         real(kind=realType), dimension(*), intent(in)  :: du1, du2, du3
         real(kind=realType), dimension(*) :: left, right
         logical :: firstOrderK
!
!        Local variables.
!
         integer(kind=intType) :: l

         real(kind=realType) :: rl1, rl2, rr1, rr2, tmp
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
             ! Linear interpolation; no limiter.
             ! Loop over the number of variables to be interpolated.

         select case (limiter)

            case (noLimiter)

               do l=1,nwInt
                  left(l)  =  omk*du1(l) + opk*du2(l)
                  right(l) = -omk*du3(l) - opk*du2(l)
               enddo
            
!          ==============================================================

            case (vanAlbeda)
             ! Nonlinear interpolation using the van albeda limiter.
             ! Loop over the number of variables to be interpolated.
               
             do l=1,nwInt

                ! Compute the limiter argument rl1, rl2, rr1 and rr2.
                ! Note the cut off to 0.0.

                tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
                rl1 = max(zero, &
                     du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
                rl2 = max(zero,du1(l)*tmp)
                
               rr1 = max(zero,du3(l)*tmp)
               rr2 = max(zero, &
                    du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))
               
               ! Compute the corresponding limiter values.
               
               rl1 = rl1*(rl1 + one)/(rl1*rl1 + one)
               rl2 = rl2*(rl2 + one)/(rl2*rl2 + one)
               rr1 = rr1*(rr1 + one)/(rr1*rr1 + one)
               rr2 = rr2*(rr2 + one)/(rr2*rr2 + one)
               
               ! Compute the nonlinear corrections to the first order
               ! scheme.
               
               left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
               right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)
               
            enddo
            
!          ==============================================================

           case (minmod)

             ! Nonlinear interpolation using the minmod limiter.
             ! Loop over the number of variables to be interpolated.

             do l=1,nwInt

               ! Compute the limiter argument rl1, rl2, rr1 and rr2.
               ! Note the cut off to 0.0.

               tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
               rl1 = max(zero, &
                         du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
               rl2 = max(zero,du1(l)*tmp)

               rr1 = max(zero,du3(l)*tmp)
               rr2 = max(zero, &
                         du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

               ! Compute the corresponding limiter values.

               rl1 = min(one, factMinmod*rl1)
               rl2 = min(one, factMinmod*rl2)
               rr1 = min(one, factMinmod*rr1)
               rr2 = min(one, factMinmod*rr2)

               ! Compute the nonlinear corrections to the first order
               ! scheme.

               left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
               right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

             enddo

         end select
            
         ! In case only a first order scheme must be used for the
         ! turbulent transport equations, set the correction for the
         ! turbulent kinetic energy to 0.

         if( firstOrderK ) then
           left(itu1)  = zero
           right(itu1) = zero
         endif

         end subroutine leftRightStateAdj
