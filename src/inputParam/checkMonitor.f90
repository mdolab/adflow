!
!      ******************************************************************
!      *                                                                *
!      * File:          checkMonitor.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkMonitor
!
!      ******************************************************************
!      *                                                                *
!      * checkMonitor checks and possibly corrects the variables        *
!      * to be monitored during the convergence.  This depends on the   *
!      * governing equations to be solved. After the correction the     *
!      * sequence of the monitoring variable names is changed, such     *
!      * that the output is independent of the specified sequence.      *
!      * Furthermore memory is allocated for the arrays used to compute *
!      * the monitoring variables and it is checked whether or not the  *
!      * maximum Mach number of total enthalpy difference is to be      *
!      * monitored.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use monitor
       use cgnsNames
       use allInputParam
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, ii, nn
       integer(kind=intType), dimension(:), allocatable :: sortNumber
       integer(kind=intType), dimension(:), allocatable :: tmpNumber

       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                                tmpNames
       logical :: RKExplicit
!
!      Function definition
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Find out if an explicit RK scheme is used in unsteady mode and
       ! set the logical RKExplicit accordingly. For explicit RK schemes
       ! no residuals are monitored.

       RKExplicit = .false.
       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) RKExplicit = .true.

       ! If the turbulent residuals must be monitored add them to the
       ! list of monitoring names. Enough memory should have been
       ! allocated for this.

       if(monDturb .and. equations == RANSEquations .and. &
          (.not. RKExplicit)) then

         select case (turbModel)

           ! One equation models of the spalart-allmaras family.

           case (spalartAllmaras, spalartAllmarasEdwards)
             nMon = nMon + 1; nMonSum = nMonSum + 1
             monNames(nMon) = cgnsL2ResNu

           ! Two equation models of the k-w family.

           case (komegaWilcox, komegaModified, menterSST)
             nMon = nMon + 2; nMonSum = nMonSum + 2
             monNames(nMon-1) = cgnsL2ResK
             monNames(nMon)   = cgnsL2ResOmega

           ! Two equation k-tau model.

           case (ktau)
             nMon = nMon + 2; nMonSum = nMonSum + 2
             monNames(nMon-1) = cgnsL2ResK
             monNames(nMon)   = cgnsL2ResTau

           ! V2f model.

           case (v2f)
             nMon = nMon + 4; nMonSum = nMonSum + 4
             monNames(nMon-3) = cgnsL2ResK
             monNames(nMon-2) = cgnsL2ResEpsilon
             monNames(nMon-1) = cgnsL2ResV2
             monNames(nMon)   = cgnsL2ResF

         end select
       endif

       ! Allocate the memory for sortNumber tmpNumber and tmpNames.

       allocate(sortNumber(nMon), tmpNumber(nMon), tmpNames(nMon), &
                stat=ierr)
       if(ierr /= 0)                    &
         call terminate("checkMonitor", &
                        "Memory allocation failure for sortNumber, etc.")

       ! Loop over the monitoring variables, copy the name into tmpName
       ! and set a number to determine its place in the sequence. If the
       ! variable cannot be monitored for the governing equations, the
       ! priority is set to a high number, such that it will be at the
       ! end of the sorted numbers. At the same time the number of
       ! variables to be monitored, nMonSum, nMonMax and nMon, is
       ! corrected.

       nn = nMon
       do i=1,nn

         tmpNames(i) = monNames(i)

         ! Determine the place in the sequence for this string.

         select case (monNames(i))
           case (cgnsL2ResRho)
             sortNumber(i) = 1
             if( RKExplicit ) then
               sortNumber(i) = 10001
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResMomx)
             sortNumber(i) = 2
             if( RKExplicit ) then
               sortNumber(i) = 10002
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResMomy)
             sortNumber(i) = 3
             if( RKExplicit ) then
               sortNumber(i) = 10003
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResMomz)
             sortNumber(i) = 4
             if( RKExplicit ) then
               sortNumber(i) = 10004
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResRhoE)
             sortNumber(i) = 5
             if( RKExplicit ) then
               sortNumber(i) = 10005
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResNu)
             sortNumber(i) = 6
             if(equations /= RANSEquations) then
               sortNumber(i) = 10001
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResK)
             sortNumber(i) = 7
             if(equations /= RANSEquations) then
               sortNumber(i) = 10002
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResOmega)
             sortNumber(i) = 8
             if(equations /= RANSEquations) then
               sortNumber(i) = 10003
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResTau)
             sortNumber(i) = 9
             if(equations /= RANSEquations) then
               sortNumber(i) = 10004
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResEpsilon)
             sortNumber(i) = 10
             if(equations /= RANSEquations) then
               sortNumber(i) = 10005
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResV2)
             sortNumber(i) = 11
             if(equations /= RANSEquations) then
               sortNumber(i) = 10006
               nMonSum       = nMonSum - 1
             endif

           case (cgnsL2ResF)
             sortNumber(i) = 12
             if(equations /= RANSEquations) then
               sortNumber(i) = 10007
               nMonSum       = nMonSum - 1
             endif

           case (cgnsCl)
             sortNumber(i) = 101
             if(flowType == internalFlow) then
               sortNumber(i) = 11001
               nMonSum       = nMonSum - 1
             endif

           case (cgnsClp)
             sortNumber(i) = 102
             if(flowType == internalFlow) then
               sortNumber(i) = 11002
               nMonSum       = nMonSum - 1
             endif

           case (cgnsClv)
             sortNumber(i) = 103
             if(equations == EulerEquations .or. &
                flowType == internalFlow) then
               sortNumber(i) = 11003
               nMonSum       = nMonSum - 1
             endif

           case (cgnsCd)
             sortNumber(i) = 104
             if(flowType == internalFlow) then
               sortNumber(i) = 11004
               nMonSum       = nMonSum - 1
             endif

           case (cgnsCdp)
             sortNumber(i) = 105
             if(flowType == internalFlow) then
               sortNumber(i) = 11005
               nMonSum       = nMonSum - 1
             endif

           case (cgnsCdv)
             sortNumber(i) = 106
             if(equations == EulerEquations .or. &
                flowType == internalFlow) then
               sortNumber(i) = 11006
               nMonSum       = nMonSum - 1
             endif

           case (cgnsCfx)
             sortNumber(i) = 107

           case (cgnsCfy)
             sortNumber(i) = 108

           case (cgnsCfz)
             sortNumber(i) = 109

           case (cgnsCmx)
             sortNumber(i) = 110

           case (cgnsCmy)
             sortNumber(i) = 111

           case (cgnsCmz)
             sortNumber(i) = 112

           case (cgnsHdiffMax)
             sortNumber(i) = 201

           case (cgnsMachMax)
             sortNumber(i) = 202

           case (cgnsYplusMax)
             sortNumber(i) = 203
             if(equations /= RANSEquations) then
               sortNumber(i) = 12003
               nMonMax       = nMonMax - 1
             endif

           case (cgnsEddyMax)
             sortNumber(i) = 204
             if(equations /= RANSEquations) then
               sortNumber(i) = 12004
               nMonMax       = nMonMax - 1
             endif

           case default
             call terminate("checkMonitor", "This should not happen")
         end select

       enddo

       ! Set the new value of nMon, because this might have changed
       ! due to the corrections.

       nMon = nMonSum + nMonMax

       ! Copy sortNumber in tmpNumber and sort it in increasing order.
       ! Note that here nn must be used and not nMon.

       do i=1,nn
         tmpNumber(i) = sortNumber(i)
       enddo

       call qsortIntegers(sortNumber, nn)

       ! Loop over the the number of monitoring variables and store the
       ! new sequence in monNames.

       do i=1,nn
         ii = bsearchIntegers(tmpNumber(i), sortNumber, nn)
         monNames(ii) = tmpNames(i)
       enddo

       ! Release the memory of sortNumber, tmpNumber and tmpNames.

       deallocate(sortNumber, tmpNumber, tmpNames, stat=ierr)
       if(ierr /= 0)                    &
         call terminate("checkMonitor", &
                        "Deallocation error for sortNumber, etc.")

       ! Allocate the memory for the monitoring variables.

       allocate(monLoc(nMon), monGlob(nMon), monRef(nMon), stat=ierr)
       if(ierr /= 0)                    &
         call terminate("checkMonitor", &
                        "Memory allocation for monitoring variables")

       ! Check if the maximum Mach number or the maximum total enthalpy
       ! difference must be monitored.

       monMachOrHMax = .false.
       do i=(nMonSum+1),nMon
         if(monNames(i) == cgnsHdiffMax .or. &
            monNames(i) == cgnsMachMax) monMachOrHMax = .true.
       enddo

       end subroutine checkMonitor
