
!      ******************************************************************
!      *                                                                *
!      * File:          computeTtot.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-17-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeTtot(rho, u, v, w, p, Ttot, kk)
!
!      ******************************************************************
!      *                                                                *
!      * computeTtot computes the total temperature for the given       *
!      * pressures, densities and velocities.                           *
!      *                                                                *
!      ******************************************************************
!
       use cpCurveFits
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: kk

       real(kind=realType), dimension(*), intent(in)  :: rho, p, u, v, w
       real(kind=realType), dimension(*), intent(out) :: Ttot
!
!      Local variables.
!
       integer(kind=intType) :: i

       real(kind=realType) :: govgm1, T, kin
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the cp model used.

       select case (cpModel)

         case (cpConstant)

           ! Constant cp and thus constant gamma. The well-known
           ! formula is valid.

           govgm1 = gammainf/(gammainf-one)

           do i=1,kk
             T       = p(i)/(rho(i)*RGas)
             kin     = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
             Ttot(i) = T*(one + rho(i)*kin/(govgm1*p(i)))
           enddo

         !===============================================================

         case (cpTempCurveFits)

           ! Cp is a function of the temperature. The formula used for
           ! constant cp is not valid anymore and a more complicated
           ! procedure must be followed.

           call terminate("computeTtot", &
                          "Variable cp formulation not implemented yet")

       end select

       end subroutine computeTtot
