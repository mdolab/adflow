!
!      ******************************************************************
!      *                                                                *
!      * File:          readCpTempCurveFits.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-11-2003                                      *
!      * Last modified: 11-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readCpTempCurveFits
!
!      ******************************************************************
!      *                                                                *
!      * readCpTempCurveFits reads the curve fits for the cp as a       *
!      * function of the temperature from the file cpFile.              *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use cpCurveFits
       use inputIO
       implicit none

       ! Local variables.

       integer, parameter :: readUnit = 32

       integer :: ios, ierr

       integer(kind=intType) :: nn, mm, kk, ii
       real(kind=realType)   :: T1, T2, e0

       character(len=2*maxStringLen) :: errorMessage
       character(len=512)            :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Open the file for reading and check if it went okay. If the file
       ! is not found, processor 0 prints an error message.

       open(unit=readUnit, file=cpFile, status="old", &
            action="read", iostat=ios)

       if(ios /= 0) then

         write(errorMessage,*) "Cp curve fit file ", trim(cpFile), &
                               " not found."
         if(myID == 0) &
           call terminate("readCpTempCurveFits", errorMessage)

         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Skip the comment lines and read the number of parts.
       ! Check if a valid number is read.

       call findNextInfoLine(readUnit, string)
       read(string,*) cpNparts

       if(cpNparts <= 0) then
         if(myID == 0)                           &
           call terminate("readCpTempCurveFits", &
                          "Wrong number of temperature ranges in &
                          &Cp curve fit file.")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Allocate the memory for the variables to store the curve fit
       ! data.

       allocate(cpTrange(0:cpNparts), cpEint(0:cpNparts),  &
                cpHint(0:cpNparts),   cpTempFit(cpNparts), &
                stat=ierr)
       if(ierr /= 0)                           &
         call terminate("readCpTempCurveFits", &
                        "Memory allocation failure for cpTrange, &
                        &cpEint, cpHint and cpTempFit")

       ! Loop over the number of temperature ranges.

       nRanges: do nn=1,cpNparts

         ! Find the next line with information and read the temperature
         ! range.

         call findNextInfoLine(readUnit, string)
         read(string,*) T1, T2

         ! If this is the first range, set the temperature range;
         ! otherwise check if the lower boundary equals the upper
         ! boundary of the previous range.

         if(nn == 1) then
           cpTrange(0) = T1
           cpTrange(1) = T2
         else
           cpTrange(nn) = T2

           if(T1 /= cpTrange(nn-1)) then
             if(myID == 0)                           &
               call terminate("readCpTempCurveFits", &
                              "Curve fit boundary not continuous")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif
         endif

         ! Read the number of points in the fit.

         call findNextInfoLine(readUnit, string)
         read(string,*) cpTempFit(nn)%nterm

         ! Allocate the memory for the exponents and the constants.

         ii = cpTempFit(nn)%nterm
         allocate(cpTempFit(nn)%exponents(ii), &
                  cpTempFit(nn)%constants(ii), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("readCpTempCurveFits", &
                          "Memory allocation failure for exponents and &
                          &constants")

         ! Read the exponents from the file.

         call findNextInfoLine(readUnit, string)
         do ii=1,cpTempFit(nn)%nterm

           ! Read the exponent from the string.

           read(string,*) cpTempFit(nn)%exponents(ii)

           ! Remove this value from the string if this is not the
           ! last exponent to be read.

           if(ii < cpTempFit(nn)%nterm) then
             ios = index(string," ")
             if(ios > 0) then
               string = string(ios:)
               string = adjustl(string)
               string = trim(string)
             else
               if(myID == 0)                           &
                 call terminate("readCpTempCurveFits", &
                                "Not enough exponents on line; &
                                &Cp curve fit file not valid.")
               call mpi_barrier(SUmb_comm_world, ierr)
             endif
           endif
         enddo

         ! Read the constants from the file.

         call findNextInfoLine(readUnit, string)
         do ii=1,cpTempFit(nn)%nterm

           ! Read the constant from the string.

           read(string,*) cpTempFit(nn)%constants(ii)

           ! Remove this value from the string if this is not the
           ! last constant to be read.

           if(ii < cpTempFit(nn)%nterm) then
             ios = index(string," ")
             if(ios > 0) then
               string = string(ios:)
               string = adjustl(string)
               string = trim(string)
             else
               if(myID == 0)                           &
                 call terminate("readCpTempCurveFits", &
                                "Not enough constants on line; &
                                &Cp curve fit file not valid.")
               call mpi_barrier(SUmb_comm_world, ierr)
             endif
           endif
         enddo

       enddo nRanges

       ! Close the file

       close(unit=readUnit)
!
!      ******************************************************************
!      *                                                                *
!      * Compute the constants eint0, such that the internal energy is  *
!      * a continous function of the temperature.                       *
!      *                                                                *
!      ******************************************************************
!
       ! First for the first interval, such that at T = 0 Kelvin the
       ! energy is also zero.

       T1  =  cpTrange(0)
       cv0 = -one          ! cv/R = cp/R - 1.0
       e0  = -T1           ! e = integral of cv, not of cp.

       do ii=1,cpTempFit(1)%nterm

         ! Update cv0.

         T2  = T1**(cpTempFit(1)%exponents(ii))
         cv0 = cv0 + cpTempFit(1)%constants(ii)*T2

         ! Update e0, for which this contribution must be integrated.
         ! Take the exceptional case exponent is -1 into account.

         if(cpTempFit(1)%exponents(ii) == -1_intType) then
           e0 = e0 + cpTempFit(1)%constants(ii)*log(T1)
         else
           T2 = T1*T2
           e0 = e0 + cpTempFit(1)%constants(ii)*T2 &
              / (cpTempFit(1)%exponents(ii) + 1)
         endif

       enddo

       ! Set the value of the internal energy at the temperature T1.
       ! Cv is assumed to be constant in the temperature range 0 - T1.
       ! Idem for the internal enthalpy.

       cpEint(0) = cv0*T1
       cpHint(0) = cpEint(0) + T1

       ! Compute the integration constant for the energy.

       cpTempFit(1)%eint0 = cpEint(0) - e0

       ! Loop over the other temperature ranges to compute their
       ! integration constant and the energy at the curve fit boundary.

       nRanges2: do nn=2,cpNparts

         ! Store nn-1, the previous temperature range, in mm.

         mm = nn - 1

         ! Store the temperature at the interface a bit easier.

         T1 = cpTrange(mm)

         ! First compute the internal energy (scaled by r) from the
         ! previous range. Actually not the energy but the enthalpy is
         ! computed. This leads to the same integraton constant.
         ! Again check for exponent -1 when integrating.

         e0 = cpTempFit(mm)%eint0

         do ii=1,cpTempFit(mm)%nterm
           if(cpTempFit(mm)%exponents(ii) == -1_intType) then
             e0 = e0 + cpTempFit(mm)%constants(ii)*log(T1)
           else
             kk = cpTempFit(mm)%exponents(ii) + 1
             T2 = T1**kk
             e0 = e0 + cpTempFit(mm)%constants(ii)*T2/kk
           endif
         enddo

         ! Store the enthalpy and energy at the curve fit boundary.
         ! Remember that cp was integrated.

         cpHint(mm) = e0
         cpEint(mm) = e0 - T1

         ! Substract the part coming from the integration of cp/r of
         ! the range nn.

         do ii=1,cpTempFit(nn)%nterm
           if(cpTempFit(nn)%exponents(ii) == -1_intType) then
             e0 = e0 - cpTempFit(nn)%constants(ii)*log(T1)
           else
             kk = cpTempFit(nn)%exponents(ii) + 1
             T2 = T1**kk
             e0 = e0 - cpTempFit(nn)%constants(ii)*T2/kk
           endif
         enddo

         ! Store the integration constant for the range nn.

         cpTempFit(nn)%eint0 = e0

       enddo nRanges2

       ! Compute the values of cv and the internal energy at the upper
       ! boundary of the curve fit. This is needed for the extrapolation
       ! of the energy if states occur with a higher temperature than
       ! the validness of the curve fits.

       ! First initialize these values.

       nn  =  cpNparts
       T1  =  cpTrange(nn)
       cvn = -one                       ! cv/R = cp/R - 1.0
       e0  =  cpTempFit(nn)%eint0 - T1  ! e = integral of cv, not of cp.

       do ii=1,cpTempFit(nn)%nterm

         ! Update cvn.

         T2  = T1**(cpTempFit(nn)%exponents(ii))
         cvn = cvn + cpTempFit(nn)%constants(ii)*T2

         ! Update e0, for which this contribution must be integrated.
         ! Take the exceptional case exponent is -1 into account.

         if(cpTempFit(nn)%exponents(ii) == -1_intType) then
           e0 = e0 + cpTempFit(nn)%constants(ii)*log(T1)
         else
           e0 = e0 + cpTempFit(nn)%constants(ii)*T2*T1 &
              /     (cpTempFit(nn)%exponents(ii) + 1)
         endif

       enddo

       ! Store e0 correctly.

       cpEint(nn) = e0
       cpHint(nn) = e0 + T1

       ! Compute the values of the integrands of cp/(R*T) at the lower
       ! and upper curve fit boundary. This is needed to compute the
       ! total pressure. This cannot be done with a single integration
       ! constant, because of the singularity at T = 0.

       nRanges3: do nn=1,cpNparts

         ! Store the temperatures of the lower and upper boundary a
         ! bit easier.

         T1 = cpTrange(nn-1)
         T2 = cpTrange(nn)

         ! Initializes the integrands to zero.

         cpTempFit(nn)%intCpovrT_1 = zero
         cpTempFit(nn)%intCpovrT_2 = zero

         ! Loop over the number of terms of the curve fits and compute
         ! the integral cp/(r*t).

         do ii=1,cpTempFit(nn)%nterm

           ! Store the coefficient a bit easier in mm. As the integral
           ! of cp/(R*T) must be computed, this is also the exponent
           ! of the primitive function; except of course when the exponent
           ! is 0.

           mm = cpTempFit(nn)%exponents(ii)

           ! Update the integrands if the temperature is larger than
           ! 0 kelvin. In case the boundary is 0 kelvin the value is not
           ! needed anyway.

           if(T1 > zero) then
             if(mm == 0_intType) then
               cpTempFit(nn)%intCpovrT_1 = cpTempFit(nn)%intCpovrT_1 &
                      + cpTempFit(nn)%constants(ii)*log(T1)
             else
               cpTempFit(nn)%intCpovrT_1 = cpTempFit(nn)%intCpovrT_1 &
                      + (cpTempFit(nn)%constants(ii)*T1**mm)/mm
             endif
           endif

           if(T2 > zero) then
             if(mm == 0_intType) then
               cpTempFit(nn)%intCpovrT_2 = cpTempFit(nn)%intCpovrT_2 &
                      + cpTempFit(nn)%constants(ii)*log(T2)
             else
               cpTempFit(nn)%intCpovrT_2 = cpTempFit(nn)%intCpovrT_2 &
                      + (cpTempFit(nn)%constants(ii)*T2**mm)/mm
             endif
           endif

         enddo

       enddo nRanges3

       end subroutine readCpTempCurveFits

!      ==================================================================

       subroutine findNextInfoLine(readUnit, string)
!
!      ******************************************************************
!      *                                                                *
!      * findNextInfoLine skips the comment lines in the given unit     *
!      * and finds the first line containing information.               *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in)             :: readUnit
       character(len=512), intent(out) :: string
!
!      Local variables.
!
       integer :: ios, ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop to skip the comment lines.

       do
         read(unit=readUnit, fmt="(a512)", iostat=ios) string

         ! Test if everything went okay.

         if(ios /= 0) then
           if(myID == 0)                        &
             call terminate("findNextInfoLine", &
                            "Unexpected end of Cp curve fit file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Get rid of the leading and trailing spaces in string.

         string = adjustl(string)
         string = trim(string)

         ! Check if this is the correct line. If so, exit

         if((len_trim(string) > 0) .and. (string(:1) /= "#")) exit
       enddo

       end subroutine findNextInfoLine
