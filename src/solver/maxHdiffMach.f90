!
!      ******************************************************************
!      *                                                                *
!      * File:          maxHdiffMach.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-01-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine maxHdiffMach(hdiffMax, MachMax)
!
!      ******************************************************************
!      *                                                                *
!      * maxHdiffMach determines the maximum value of the Mach number   *
!      * and total enthalpy (or better the relative total enthalpy      *
!      * difference with the freestream).                               *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use flowVarRefState
       use inputPhysics
       use monitor
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(out) :: hdiffMax, MachMax
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: hdiff, hInf, Mach2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the maximum values to zero.

       hdiffMax = zero
       MachMax  = zero

       ! In case none of the two variables needs to be monitored,
       ! a return is made.

       if(.not. monMachOrHMax) return

       ! Set the free stream value of the total enthalpy.

       hInf = (wInf(irhoE) + pInfCorr)/rhoInf

       ! Loop over the owned cells of this block.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the local total enthalpy and Mach number squared.

             hdiff = abs((w(i,j,k,irhoE) + p(i,j,k))/w(i,j,k,irho) - hInf)
             Mach2 = (w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                   +  w(i,j,k,ivz)**2)*w(i,j,k,irho)/(gamma(i,j,k)*p(i,j,k))

             ! Determine the maximum of these values and the
             ! currently stored maximum values.

             hdiffMax = max(hdiffMax, hdiff)
             MachMax  = max(MachMax,  Mach2)

           enddo
         enddo
       enddo

       ! Currently the maximum Mach number squared is stored in
       ! MachMax. Take the square root. Also create a relative
       ! total enthalpy difference.

       MachMax  = sqrt(MachMax)
       hdiffMax = hdiffMax/hInf

       end subroutine maxHdiffMach
