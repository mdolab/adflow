module inputParamRoutines

  use constants, only : maxStringLen, intType

  ! Set the parameter none, which is used as a check to see
  ! whether or not some key parameters were specified.

  integer(kind=intType), parameter :: none = 0

  ! monDturb:             Whether or not the turbulent residuals
  !                        must be monitored. This must be done via
  !                        this construction, because during the
  !                        reading of the monitoring variables the
  !                        turbulence model might not be known.
  ! monitorSpecified:     Whether or not the monitoring variables
  !                        were specified.
  ! surfaceOutSpecified: Whether or not the surface output
  !                        variables were specified.
  ! volumeOutSpecified:  Whether or not the volume output
  !                        variables were specified.
  ! isoOutSpecified:      Wheter or not the isosurface output
  !                        variables were specified
  logical :: monDturb
  logical :: monitorSpecified
  logical :: surfaceOutSpecified
  logical :: volumeOutSpecified
  logical :: isoOutSpecified

  ! liftDirSpecified: Whether or not the lift direction was
  !                     specified.

  logical :: liftDirSpecified

contains

  subroutine checkMonitor
    !
    !       checkMonitor checks and possibly corrects the variables
    !       to be monitored during the convergence.  This depends on the
    !       governing equations to be solved. After the correction the
    !       sequence of the monitoring variable names is changed, such
    !       that the output is independent of the specified sequence.
    !       Furthermore memory is allocated for the arrays used to compute
    !       the monitoring variables and it is checked whether or not the
    !       maximum Mach number of total enthalpy difference is to be
    !       monitored.
    !
    use constants
    use cgnsNames
    use monitor, only : monNames, monGlob, monLoc, nMonMax, nMonSum, &
         nMon, monMachOrHMax, monRef
    use inputPhysics, only : equations, flowType, turbModel, equationMode
    use inputUnsteady, only : timeIntegrationScheme
    use sorting, only :qsortIntegers, bsearchIntegers
    use utils, only : terminate
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

       case('totalR')
          sortNumber(i) = 113

       case(cgnsSepSensor)
          sortNumber(i) = 114

       case (cgnsCavitation)
          sortNumber(i) = 115

       case(cgnsAxisMoment)
          sortNumber(i) = 116

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
       ii = bsearchIntegers(tmpNumber(i), sortNumber)
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
  subroutine checkOutput
    !
    !       checkOutput checks and possibly corrects the and output
    !       variables. This depends on the set of governing equations to
    !       be solved.
    !
    use constants
    use extraOutput
    use inputPhysics, only : equations, equationMode, wallDistanceNeeded
    use inputUnsteady, only : timeIntegrationScheme
    use flowVarRefState, only : kPresent, eddyModel
    implicit none

    ! Determine the governing equations to be solved and set the
    ! variables accordingly.

    select case (equations)
    case (EulerEquations)
       surfWriteCf    = .false.
       surfWriteCh    = .false.
       surfWriteYplus = .false.
       surfWriteCfx   = .false.
       surfWriteCfy   = .false.
       surfWriteCfz   = .false.

       volWriteMachTurb     = .false.
       volWriteEddyVis      = .false.
       volWriteRatioEddyVis = .false.
       volWriteDist         = .false.
       volWriteResTurb      = .false.

    case (NSEquations)
       surfWriteYplus = .false.

       volWriteMachTurb     = .false.
       volWriteEddyVis      = .false.
       volWriteRatioEddyVis = .false.
       volWriteDist         = .false.
       volWriteResTurb      = .false.

    case (RANSEquations)

       ! Check if it is possible to write a turbulent Mach
       ! number and eddy viscosity; this depends on the turbulence
       ! model used.

       if(.not. kPresent) volWriteMachTurb = .false.
       if(.not. eddyModel) then
          volWriteEddyVis      = .false.
          volWriteRatioEddyVis = .false.
       endif

       ! If a wall distance free turbulence model is used,
       ! set volWriteDist to .false.

       if(.not. wallDistanceNeeded) volWriteDist = .false.

    end select

    if(equationMode          == unsteady .and. &
         timeIntegrationScheme == explicitRK) then
       volWriteResRho  = .false.
       volWriteResMom  = .false.
       volWriteResRhoE = .false.
       volWriteResTurb = .false.
    endif

  end subroutine checkOutput
  subroutine defaultIsoOut
    !
    !       defaultIsoOut sets the default set of additional
    !       variables to be written to the solution file; the primitive
    !       variables are always written. This additional set depends on
    !       the governing equations to be solved.
    !
    use constants
    use extraOutput
    use inputPhysics, only : equations
    implicit none

    ! First set the variables, which are independent from the
    ! governing equations to be solved.
    isoWriteRho      = .false.
    isoWriteVx       = .false.
    isoWriteVy       = .false.
    isoWriteVz       = .false.
    isoWriteP        = .false.
    isoWriteMx       = .false.
    isoWriteMy       = .false.
    isoWriteMz       = .false.
    isoWriteRhoe     = .false.
    isoWriteTemp     = .false.
    isoWriteCp       = .false.
    isoWriteMach     = .false.
    isoWriteMachTurb = .false.
    isoWriteDist     = .false.
    isoWriteVort     = .false.
    isoWriteVortx    = .false.
    isoWriteVorty    = .false.
    isoWriteVortz    = .false.
    isoWritePtotloss = .false.
    isoWriteResRho   = .false.
    isoWriteResMom   = .false.
    isoWriteResRhoe  = .false.
    isoWriteShock    = .false.
    isoWriteFilteredShock = .false.
    ! Set the values which depend on the equations to be solved.

    select case (equations)
    case (EulerEquations)
       isoWriteEddyVis      = .false.
       isoWriteRatioEddyVis = .false.
       isoWriteResTurb      = .false.

    case (NSEquations)
       isoWriteEddyVis      = .false.
       isoWriteRatioEddyVis = .false.
       isoWriteResTurb      = .false.

    case (RANSEquations)
       isoWriteEddyVis      = .false.
       isoWriteRatioEddyVis = .false.
       isoWriteResTurb      = .false.
    end select

  end subroutine defaultIsoOut
  subroutine defaultMonitor
    !
    !       defaultMonitor sets the default set of variables to be
    !       monitored during the convergence. This set depends on the
    !       governing equations to be solved.
    !
    use constants
    use cgnsNames
    use inputPhysics, only : equations, flowType
    use monitor, only : nMOnSum, nMonMax, nMon, monNames, showCPU
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! CPU time is written to stdout.

    showCPU = .true.

    ! Determine the governing equations to be solved.

    select case (equations)
    case (EulerEquations)

       ! Set the number of summation and maximum monitor variables
       ! and allocate the memory for the monitoring names.
       ! A distinction is made between internal and external flows,
       ! because cl and cd do not make a lot of sense for the former.

       if(flowType == internalFlow) then

          ! Internal flow; only the density residual is monitored.

          nMonSum = 1; nMonMax = 0; nMon = 1
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho

       else

          ! External; also lift and drag is monitored.

          nMonSum = 3; nMonMax = 0; nMon = 3
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho
          monNames(2) = cgnsCl
          monNames(3) = cgnsCd

       endif

    case (NSEquations)

       ! Set the number of summation and maximum monitor variables
       ! and allocate the memory for the monitoring names.
       ! A distinction is made between internal and external flows,
       ! because cl and cd do not make a lot of sense for the former.

       if(flowType == internalFlow) then

          ! Internal flow; only the density residual is monitored.

          nMonSum = 1; nMonMax = 0; nMon = 1
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho

       else

          ! External; also lift and drag (total and viscous)
          ! is monitored.

          nMonSum = 4; nMonMax = 0; nMon = 4
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho
          monNames(2) = cgnsCl
          monNames(3) = cgnsCd
          monNames(4) = cgnsCdv

       endif

    case (RANSEquations)

       ! Set the number of summation and maximum monitor variables
       ! and allocate the memory for the monitoring names.
       ! A distinction is made between internal and external flows,
       ! because cl and cd do not make a lot of sense for the former.

       if(flowType == internalFlow) then

          ! Internal flow; the density residual as well as the
          ! maximum values of yplus and the eddy viscosity ration
          ! are monitored.

          nMonSum = 1; nMonMax = 2; nMon = 3
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho
          monNames(2) = cgnsYplusMax
          monNames(3) = cgnsEddyMax

       else

          ! External; also lift and drag (total and viscous)
          ! is monitored.


          nMonSum = 4; nMonMax = 2; nMon = 6
          allocate(monNames(nMon), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
               "Memory allocation failure for monNames")

          ! Set the names for the variables to be monitored.

          monNames(1) = cgnsL2resRho
          monNames(2) = cgnsCl
          monNames(3) = cgnsCd
          monNames(4) = cgnsCdv
          monNames(5) = cgnsYplusMax
          monNames(6) = cgnsEddyMax

       endif

    end select

  end subroutine defaultMonitor
  subroutine defaultSurfaceOut
    !
    !       defaultSurfaceOut sets the default set of surface variables
    !       to be written to the solution file. This set depends on the
    !       governing equations to be solved.
    !
    use constants
    use extraOutput
    use inputPhysics, only : equations
    implicit none

    ! First set the variables, which are independent from the
    ! governing equations to be solved.

    surfWriteRho  = .true.
    surfWriteP    = .false.
    surfWriteTemp = .false.
    surfWriteVx   = .true.
    surfWriteVy   = .true.
    surfWriteVz   = .true.
    surfWriteCp   = .true.
    surfWriteMach = .true.

    ! Set the values which depend on the equations to be solved.

    select case (equations)
    case (EulerEquations)
       surfWritePtotloss = .true.
       surfWriteCf       = .false.
       surfWriteCh       = .false.
       surfWriteYplus    = .false.
       surfWriteCfx      = .false.
       surfWriteCfy      = .false.
       surfWriteCfz      = .false.

    case (NSEquations)
       surfWritePtotloss = .false.
       surfWriteCf       = .true.
       surfWriteCh       = .false.
       surfWriteYplus    = .false.
       surfWriteCfx      = .true.
       surfWriteCfy      = .true.
       surfWriteCfz      = .true.

    case (RANSEquations)
       surfWritePtotloss = .false.
       surfWriteCf       = .true.
       surfWriteCh       = .false.
       surfWriteYplus    = .true.
       surfWriteCfx      = .true.
       surfWriteCfy      = .true.
       surfWriteCfz      = .true.
    end select

  end subroutine defaultSurfaceOut
  subroutine defaultVolumeOut
    !
    !       defaultVolumeOut sets the default set of additional
    !       variables to be written to the solution file; the primitive
    !       variables are always written. This additional set depends on
    !       the governing equations to be solved.
    !
    use constants
    use extraOutput
    use inputPhysics, only : equations
    implicit none

    ! First set the variables, which are independent from the
    ! governing equations to be solved.

    volWriteMx       = .false.
    volWriteMy       = .false.
    volWriteMz       = .false.
    volWriteRhoe     = .false.
    volWriteTemp     = .false.
    volWriteCp       = .false.
    volWriteMach     = .false.
    volWriteMachTurb = .false.
    volWriteDist     = .false.
    volWriteVort     = .false.
    volWriteVortx    = .false.
    volWriteVorty    = .false.
    volWriteVortz    = .false.
    volWritePtotloss = .true.
    volWriteResRho   = .true.
    volWriteResMom   = .false.
    volWriteResRhoe  = .false.

    ! Set the values which depend on the equations to be solved.

    select case (equations)
    case (EulerEquations)
       volWriteEddyVis      = .false.
       volWriteRatioEddyVis = .false.
       volWriteResTurb      = .false.

    case (NSEquations)
       volWriteEddyVis      = .false.
       volWriteRatioEddyVis = .false.
       volWriteResTurb      = .false.

    case (RANSEquations)
       volWriteEddyVis      = .true.
       volWriteRatioEddyVis = .true.
       volWriteResTurb      = .true.
    end select

  end subroutine defaultVolumeOut

  subroutine dummyreadParamFile
    !
    !       This subroutine is the same as readParamFile EXCEPT it does not
    !       read the actual file. Values are set diectly from python for
    !       all the options and then this file is run.
    !
    use constants
    use inputPhysics, only : cpModel

    implicit none
    !
    !      Local variables
    !
    integer, parameter :: readUnit = 32


    ! Check if all the desired input parameters were specified and
    ! print warnings if some irrelevant ones are specified.

    call checkInputParam

    ! Read the cp curve fits from file if a variable cp model
    ! must be used.

    if(cpModel == cpTempCurveFits) call readCpTempCurveFits

    ! Determine the number of governing equations and set the
    ! corresponding parameters accordingly.

    call setEquationParameters

    ! Extract the multigrid info from the string.

    call extractMgInfo

    ! If no monitoring variables were specified, set the default set.
    ! Idem for the surface and volume output variables.

    if(.not. monitorSpecified)    call defaultMonitor
    if(.not. surfaceOutSpecified) call defaultSurfaceOut
    if(.not. volumeOutSpecified)  call defaultVolumeOut
    if(.not. isoOutSpecified)     call defaultIsoOut

    ! Check the monitoring and output variables.

    call checkMonitor
    call checkOutput


  end subroutine dummyreadParamFile
  subroutine extractMGInfo
    !
    !       extractMgInfo creates the integer array cycleStrategy from
    !       the string describing the multigrid strategy. This string
    !       either contains a predefined strategy, like sg, 2v, 4w, etc.,
    !       or a combination of -1's, 0's and 1's, which defines a user
    !       defined strategy. The integers -1, 0 and 1 have the following
    !       meaning:  0 -> perform an iteration step on the current grid.
    !                 1 -> go to next coarser grid.
    !                -1 -> go to next finer grid.
    !       For a valid cycling strategy the sum of the elements of the
    !       array should be 0.
    !
    use constants
    use inputIteration, only :cycleStrategy, nMGLevels, nMGSteps, mgStartLevel, mgDescription
    use inputPhysics, only : equationMode
    use inputUnsteady, only : timeIntegrationScheme
    use communication, only : myID
    use utils, only : convertToLowerCase, terminate
    implicit none
    !
    !      Local variables
    !
    integer                      :: stringLen, error
    integer(kind=intType)        :: i, ii, nMinus, nn
    character (len=maxStringLen) :: errorMessage

    ! For an unsteady computation using explicit Runge-Kutta schemes
    ! overrule mgDescription to sg.

    if(equationMode          == unsteady .and. &
         timeIntegrationScheme == explicitRK) mgDescription = "sg"

    ! Create a lower case version of mgDescription and determine the
    ! length of the string. Note that if the string contains a user
    ! defined cycling strategy the contents (0's, 1's, -1's) is not
    ! changed by the call to convertToLowerCase

    mgDescription = trim(mgDescription)
    call convertToLowerCase(mgDescription)
    stringLen = len_trim(mgDescription)

    ! Check for predefined cycling strategies.

    if(mgDescription == "sg") then

       ! Single grid computation. Set the values for the parameters
       ! nMGSteps (number of steps in the cycle strategy) and
       ! cycleStrategy (the cycle strategy itself).

       nMGSteps = 1
       allocate(cycleStrategy(1), stat=error)
       if(error /= 0) call terminate("extractMgInfo", &
            "allocation error for 1 integer")
       cycleStrategy(1) = 0

       ! Set the number of grid levels needed by the multigrid to 1.

       nMGLevels = 1

    else if(mgDescription(stringLen:stringLen) == "v") then

       ! Must be a v-cycle. The rest of the string should only contain
       ! digits. Check this.

       if(.not. digitsOnlyInString(mgDescription(:stringLen-1))) then
          write(errorMessage,*) "Invalid cycle strategy, ", &
               mgDescription(:stringLen), ", specified"
          if(myID == 0) call terminate("extractMgInfo", errorMessage)
       endif

       ! Read the number of levels in the cycle.

       read(mgDescription(:stringLen-1),*) nMGLevels

       ! Determine the number of steps in cycleStrategy and allocate
       ! the memory for it.

       nMGSteps = 4*nMGLevels - 4
       allocate(cycleStrategy(nMGSteps), stat=error)
       if(error /= 0) then
          write(errorMessage,*) "Allocation error for", nMGSteps, &
               "integers for the v-cycle ",       &
               mgDescription(:stringLen)
          call terminate("extractMgInfo", errorMessage)
       endif

       ! Set the values of cycleStrategy.

       ii = 1
       do i=1,(nMGLevels-1)
          cycleStrategy(ii)   = 0
          cycleStrategy(ii+1) = 1
          ii = ii+2
       enddo

       do i=1,(nMGLevels-1)
          cycleStrategy(ii)   =  0
          cycleStrategy(ii+1) = -1
          ii = ii+2
       enddo

    else if(mgDescription(stringLen:stringLen) == "w") then

       ! Must be a w-cycle. The rest of the string should only contain
       ! digits. Check this.

       if(.not. digitsOnlyInString(mgDescription(:stringLen-1))) then
          write(errorMessage,*) "Invalid cycle strategy, ", &
               mgDescription(:stringLen), ", specified"
          if(myID == 0) call terminate("extractMgInfo", errorMessage)
       endif

       ! Read the number of levels in the cycle.

       read(mgDescription(:stringLen-1),*) nMGLevels

       ! Determine the number of steps in cycleStrategy and allocate
       ! the memory for it.

       nMGSteps = computeNstepsWcycle(nMGLevels)
       allocate(cycleStrategy(nMGSteps), stat=error)
       if(error /= 0) then
          write(errorMessage,*) "Allocation error for", nMGSteps, &
               "integers for the w-cycle ",       &
               mgDescription(:stringLen)
          call terminate("extractMgInfo", errorMessage)
       endif

       ! Set the values of cycleStrategy.

       ii = 1
       call setEntriesWcycle(ii, nMGLevels)

    else

       ! The string must be a collection of 0's, -1's and 1's to
       ! describe the cycle strategy. Get rid of the internal spaces
       ! first and determine the amount of -'s.

       ii = 0
       nMinus = 0
       do i=1,stringLen
          if(mgDescription(i:i) /= " ") then
             ii = ii+1
             mgDescription(ii:ii) = mgDescription(i:i)
             if(mgDescription(ii:ii) == "-") nMinus = nMinus+1
          endif
       enddo
       stringLen = ii

       ! Determine the number of steps in the cycle strategy and
       ! allocate the memory for it.

       nMGSteps = ii - nMinus
       allocate(cycleStrategy(nMGSteps), stat=error)
       if(error /= 0) then
          write(errorMessage,*) "Allocation error for", nMGSteps, &
               "integers for the cycle strategy"
          call terminate("extractMgInfo", errorMessage)
       endif

       ! Determine the entries for cycleStrategy.

       i = 1
       nn = 1
       do
          ii = i
          if(mgDescription(i:i) == "-") i = i+1

          ! Determine the case we are having here.

          select case (mgDescription(ii:i))
          case ("0")
             cycleStrategy(nn) = 0
          case ("1")
             cycleStrategy(nn) = 1
          case ("-1")
             cycleStrategy(nn) = -1
          case default
             write(errorMessage, *) "Invalid character, ", &
                  mgDescription(ii:i), &
                  ", in the string describing &
                  &cycling strategy"
             if(myID == 0) call terminate("extractMgInfo", errorMessage)
          end select

          ! Update i and nn

          i  = i  + 1
          nn = nn + 1

          ! Exit the do loop in case i is larger than stringLen.

          if(i > stringLen) exit
       enddo

       ! Check if the string specified is valid and determine the
       ! maximum grid level needed in the cycle.

       nn = 0
       nMGLevels = 0
       do i=1,nMGSteps
          nn = nn + cycleStrategy(i)
          nMGLevels = max(nn, nMGLevels)
       enddo
       nMGLevels = nMGLevels + 1

       if(nn /= 0 .and. myID == 0) &
            call terminate("extractMgInfo", &
            "sum of coefficients in cycle strategy is not 0")
    endif

    ! Correct the value of mgStartlevel in case a nonpositive number
    ! has been specified. In that case it is set to -1, the default
    ! value.

    if(mgStartlevel <= 0) mgStartlevel = -1

    ! Determine the value of mgStartlevel. This parameter might be
    ! specified in the python script file and is checked here for
    ! consistency. If mgStartlevel has not been specified to any
    ! specific level it is set to the coarsest level in the mg cycle
    ! (starting from free stream) or to the finest level (restart).
    ! The restart is handled in python wrapper.

    if(mgStartlevel == -1) then

       ! Value has not been specified. Default value is set, see
       ! the comments above.
       mgStartlevel = nMGLevels

    endif

  end subroutine extractMGInfo

  !      ==================================================================

  logical function digitsOnlyInString(string)
    !
    !       digitsOnlyInString checks whether the given string contains
    !       digits only or if other character types are present. In the
    !       former case the function returns .True., otherwise .False.
    !
    implicit none
    !
    !      Subroutine argument                                              *
    !
    character (len=*), intent(in) :: string
    !
    !      Local variables
    !
    integer :: i, stringLen

    ! Initialize digitsOnlyInString to .True.

    digitsOnlyInString = .true.

    ! Determine the length of the string.

    stringLen = len_trim(string)

    ! Loop over the elements of the string and check if they are digits.

    do i=1,stringLen
       if(string(i:i) < "0" .or. string(i:i) > "9") &
            digitsOnlyInString = .false.
    enddo

  end function digitsOnlyInString

  !      ==================================================================

  recursive function computeNstepsWcycle(nLevels) result(nSteps)
    !
    !       computeNstepsWcycle is recursive function, which determines
    !       the number of entries of a w-cycle of a given level.
    !
    use constants
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Result variable
    !
    integer(kind=intType) :: nSteps
    !
    !      Function argument
    !
    integer(kind=intType), intent(in) :: nLevels
    !
    !      Local variables
    !
    character (len=maxStringLen) :: errorMessage

    ! Determine the case we are having here. For nLevels is less
    ! than 2 an error message is printed, in case nLevels is 2
    ! the recursion is broken and otherwise a recursive call is made.

    if(nLevels < 2) then
       write(errorMessage,*) "Wrong value of nLevels", nLevels
       if(myID == 0) call terminate("computeNstepsWcycle", errorMessage)
    else if(nLevels == 2) then
       nSteps = 4
    else
       nSteps = 4 + 2*computeNstepsWcycle(nLevels-1)
    endif

  end function computeNstepsWcycle

  !      ==================================================================

  recursive subroutine setEntriesWcycle(counter, nLevels)
    !
    !       setEntriesWcycle is a recursive subroutine, which actually
    !       fills the entries of cycleStrategy for a w-cycle.
    !
    use constants
    use inputIteration, only : cycleStrategy
    use communication, only : myID
    use utils, only : terminate
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(inout) :: counter
    integer(kind=intType), intent(in)    :: nLevels
    !
    !      Local variables
    !
    character (len=maxStringLen) :: errorMessage

    ! Determine the case we are having here. For nLevels is less
    ! than 2 an error message is printed, in case nLevels is 2
    ! the recursion is broken and otherwise a recursive call is made.

    if(nLevels < 2) then

       write(errorMessage,*) "Wrong value of nLevels", nLevels
       if(myID == 0) call terminate("setEntriesWcycle", errorMessage)

    else if(nLevels == 2) then

       cycleStrategy(counter)   =  0
       cycleStrategy(counter+1) =  1
       cycleStrategy(counter+2) =  0
       cycleStrategy(counter+3) = -1

       counter = counter + 4

    else

       cycleStrategy(counter)   =  0
       cycleStrategy(counter+1) =  1
       counter = counter + 2

       call setEntriesWcycle(counter, nLevels-1)
       call setEntriesWcycle(counter, nLevels-1)

       cycleStrategy(counter)   =  0
       cycleStrategy(counter+1) = -1
       counter = counter + 2

    endif

  end subroutine setEntriesWcycle
  subroutine isoVariables(variables)
    !
    !       isoVariables extracts from the given string the extra
    !       iso surface variables to be written to the solution file.
    !
    use constants
    use extraOutput
    use utils, only : convertToLowerCase, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(inout) :: variables
    !
    !      Local variables.
    !
    integer :: nVarSpecified, pos

    character(len=15) :: keyword
    character(len=maxStringLen) :: errorMessage

    ! Convert the string variables to lower case.

    call convertToLowerCase(variables)

    ! Initialize all the iso output variables to .False.
    isoWriteRho   = .false.
    isoWriteVx    = .false.
    isoWriteVy    = .false.
    isoWriteVz    = .false.
    isoWriteP     = .false.
    isoWriteTurb  = .false.

    isoWriteMx    = .false.
    isoWriteMy    = .false.
    isoWriteMz    = .false.
    isoWriteRhoe  = .false.
    isoWriteTemp  = .false.
    isoWriteVort  = .false.
    isoWriteVortx = .false.
    isoWriteVorty = .false.
    isoWriteVortz = .false.

    isoWriteCp       = .false.
    isoWriteMach     = .false.
    isoWriteMachTurb = .false.
    isoWritePtotloss = .false.

    isoWriteEddyVis      = .false.
    isoWriteRatioEddyVis = .false.
    isoWriteDist         = .false.

    isoWriteResRho  = .false.
    isoWriteResMom  = .false.
    isoWriteResRhoe = .false.
    isoWriteResTurb = .false.

    isoWriteShock = .false.
    isoWriteFilteredShock = .false.

    isoWriteBlank = .false.

    ! Initialize nVarSpecified to 0. This serves as a test
    ! later on.

    nVarSpecified = 0

    ! Loop to extract the info from the string variables.

    do
       ! Condition to exit the loop.

       if(len_trim(variables) == 0) exit

       ! Locate the first occurance of the _ in the string and
       ! determine the string keyword.

       pos = index(variables, "_")
       if(pos == 0) then
          keyword   = variables
          variables = ""
       else
          keyword   = variables(:pos-1)
          variables = variables(pos+1:)
       endif

       ! Check the keyword.

       select case (keyword)
       case ("")
          ! Multiple occurence of "_". Just ignore it.

       case("rho")
          isoWriteRho   = .true.
          nVarSpecified = nVarSpecified + 1

       case("vx")
          isoWriteVx   = .true.
          nVarSpecified = nVarSpecified + 1

       case("vy")
          isoWriteVy   = .true.
          nVarSpecified = nVarSpecified + 1

       case("vz")
          isoWriteVz   = .true.
          nVarSpecified = nVarSpecified + 1

       case("P")
          isoWriteP   = .true.
          nVarSpecified = nVarSpecified + 1

       case("turb")
          isoWriteTurb = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mx")
          isoWriteMx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("my")
          isoWriteMy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mz")
          isoWriteMz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvx")
          isoWriteRVx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvy")
          isoWriteRVy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvz")
          isoWriteRVz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rhoe")
          isoWriteRhoe = .true.
          nVarSpecified = nVarSpecified + 1

       case ("temp")
          isoWriteTemp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vort")
          isoWriteVort = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vortx")
          isoWriteVortx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vorty")
          isoWriteVorty = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vortz")
          isoWriteVortz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cp")
          isoWriteCp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mach")
          isoWriteMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rmach")
          isoWriteRMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("macht")
          isoWriteMachTurb = .true.
          nVarSpecified = nVarSpecified + 1

       case ("ptloss")
          isoWritePtotloss = .true.
          nVarSpecified = nVarSpecified + 1

       case ("eddy")
          isoWriteEddyVis = .true.
          nVarSpecified = nVarSpecified + 1

       case ("eddyratio")
          isoWriteRatioEddyVis = .true.
          nVarSpecified = nVarSpecified + 1

       case ("dist")
          isoWriteDist = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resrho")
          isoWriteResRho = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resmom")
          isoWriteResMom = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resrhoe")
          isoWriteResRhoe = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resturb")
          isoWriteResTurb = .true.
          nVarSpecified = nVarSpecified + 1

       case ("blank")
          isoWriteBlank = .true.
          nVarSpecified = nVarSpecified + 1

       case("shock")
          isoWriteShock = .true.
          nVarSpecified = nVarSpecified + 1

       case("filteredshock")
          isoWriteFilteredShock = .true.
          nVarSpecified = nVarSpecified + 1


       case default
          pos = len_trim(keyword)
          write(errorMessage,"(3a)" ) "Unknown extra iso output &
               &variable, ", trim(keyword), &
               ", specified"
          call terminate("isoVariables", errorMessage)

       end select

    enddo

    ! Set this to true regardless...it is possible no varibles were
    ! specified
    isoOutSpecified = .true.

  end subroutine isoVariables

  subroutine monitorVariables(variables)
    !
    !       monitorVariables extracts from the given string the variables
    !       to be monitored during the convergence.
    !
    use constants
    use cgnsNames
    use communication, only : myid, adflow_comm_world
    use monitor, only : monNames, nMOn , nMonMax, nMonSum, showCPU
    use utils, only : convertToLowerCase, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(inout) :: variables
    !
    !      Local parameter.
    !
    integer(kind=intType), parameter :: nVarMax = 21
    !
    !      Local variables.
    !
    integer :: pos, ierr

    character(len=15) :: keyword
    character(len=maxStringLen) :: errorMessage

    character(len=maxCGNSNameLen), dimension(nVarMax) :: tmpNames

    logical :: monDrho, monTotalR

    ! Check if the monitoring names have already been allocated.
    ! This happens when multiple lines for the monitoring variables
    ! are specified in the parameter file. If this happens the last
    ! value is taken and thus release the memory of previously
    ! specified names.

    if( allocated(monNames) ) then
       deallocate(monNames, stat=ierr)
       if(ierr /= 0) call terminate("monitorVariables", &
            "Deallocation error for monNames")
    endif

    ! Initialize monDrho, monDturb and showCPU to .false.

    monDrho  = .false.
    monTotalR = .false.
    monDturb = .false.
    showCPU  = .false.

    ! Initialize nMonSum, nMonMax and nMon to 0.

    nMonSum = 0
    nMonMax = 0
    nMon    = 0

    ! Convert the string variables to lower case.

    call convertToLowerCase(variables)

    ! Loop to extract the info from the string variables.

    do
       ! Condition to exit the loop.

       if(len_trim(variables) == 0) exit

       ! Locate the first occurance of the _ in the string and
       ! determine the string keyword.

       pos = index(variables, "_")
       if(pos == 0) then
          keyword   = variables
          variables = ""
       else
          keyword   = variables(:pos-1)
          variables = variables(pos+1:)
       endif

       ! Check the keyword.

       select case (keyword)
       case ("")
          ! Multiple occurence of "_". Just ignore it.

       case ("cpu")    ! only written to stdout.
          showCPU = .true.

       case ("resrho")
          monDrho = .true.
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsL2resRho

       case ("resmom")
          nMon = nMon + 3; nMonSum = nMonSum + 3
          tmpNames(nMon-2) = cgnsL2resMomx
          tmpNames(nMon-1) = cgnsL2resMomy
          tmpNames(nMon)   = cgnsL2resMomz

       case ("resrhoe")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsL2resRhoe

       case ("resturb")  ! special case, because the turbulence model
          ! Is not yet known. See checkMonitor.
          monDturb = .true.

       case ("cl")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCl

       case ("clp")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsClp

       case ("clv")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsClv

       case ("cd")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCd

       case ("cdp")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCdp

       case ("cdv")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCdv

       case ("cfx")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCfx

       case ("cfy")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCfy

       case ("cfz")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCfz

       case ("cmx")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCmx

       case ("cmy")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCmy

       case ("cmz")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCmz

       case ("hdiff")
          nMon = nMon + 1; nMonMax = nMonMax + 1
          tmpNames(nMon) = cgnsHdiffMax

       case ("mach")
          nMon = nMon + 1; nMonMax = nMonMax + 1
          tmpNames(nMon) = cgnsMachMax

       case ("yplus")
          nMon = nMon + 1; nMonMax = nMonMax + 1
          tmpNames(nMon) = cgnsYplusMax

       case ("eddyv")
          nMon = nMon + 1; nMonMax = nMonMax + 1
          tmpNames(nMon) = cgnsEddyMax

       case("totalr")
          monTotalR = .True.
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = 'totalR'

       case("sepsensor")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsSepSensor

       case("cavitation")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsCavitation

       case("axismoment")
          nMon = nMon + 1; nMonSum = nMonSum + 1
          tmpNames(nMon) = cgnsAxisMoment

       case default
          write(errorMessage,"(3a)") "Unknown monitoring variable, ", &
               trim(keyword), ", specified"
          if(myID == 0) &
               call terminate("monitorVariables", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)

       end select

    enddo

    ! If the density residual was not specified to be monitored,
    ! add it to tmpNames.

    if(.not. monDrho) then
       nMon = nMon + 1; nMonSum = nMonSum + 1
       tmpNames(nMon) = cgnsL2resRho
    endif

    ! If the total residual was not specified to be monitored, add
    ! it to tmpNames.

    if(.not. monTotalR) then
       nMon = nMon + 1; nMonSum = nMonSum + 1
       tmpNames(nMon) = "totalR"
    endif


    ! Allocate the memory for monNames. If the turbulent residuals
    ! must be monitored allocate some extra place.

    pos = nMon
    if( monDturb ) pos = nMon + 4
    allocate(monNames(pos), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("monitorVariables", &
         "Memory allocation failure for monNames")

    ! Copy the monitoring names into monNames.

    do pos=1,nMon
       monNames(pos) = tmpNames(pos)
    enddo

    ! Set monitorSpecified to .true. to indicate that monitoring
    ! variables have been specified.

    monitorSpecified = .true.

  end subroutine monitorVariables
  subroutine readCpTempCurveFits
    !
    !       readCpTempCurveFits reads the curve fits for the cp as a
    !       function of the temperature from the file cpFile.
    !
    use constants
    use communication, only : myid, adflow_comm_world
    use cpCurveFits, only : cpTrange, cpTempFit, cpEint, cpHint, cvn, cv0, &
         cpNParts
    use inputIO, only : cpFile
    use utils, only : terminate
    implicit none

    ! Local variables.

    integer, parameter :: readUnit = 32

    integer :: ios, ierr

    integer(kind=intType) :: nn, mm, kk, ii
    real(kind=realType)   :: T1, T2, e0

    character(len=2*maxStringLen) :: errorMessage
    character(len=512)            :: string

    ! Open the file for reading and check if it went okay. If the file
    ! is not found, processor 0 prints an error message.

    open(unit=readUnit, file=cpFile, status="old", &
         action="read", iostat=ios)

    if(ios /= 0) then

       write(errorMessage,*) "Cp curve fit file ", trim(cpFile), &
            " not found."
       if(myID == 0) &
            call terminate("readCpTempCurveFits", errorMessage)

       call mpi_barrier(ADflow_comm_world, ierr)
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
       call mpi_barrier(ADflow_comm_world, ierr)
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
             call mpi_barrier(ADflow_comm_world, ierr)
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
                call mpi_barrier(ADflow_comm_world, ierr)
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
                call mpi_barrier(ADflow_comm_world, ierr)
             endif
          endif
       enddo

    enddo nRanges

    ! Close the file

    close(unit=readUnit)
    !
    !       Compute the constants eint0, such that the internal energy is
    !       a continous function of the temperature.
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
    !       findNextInfoLine skips the comment lines in the given unit
    !       and finds the first line containing information.
    !
    use communication
    use utils, only : terminate
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

    ! Loop to skip the comment lines.

    do
       read(unit=readUnit, fmt="(a512)", iostat=ios) string

       ! Test if everything went okay.

       if(ios /= 0) then
          if(myID == 0)                        &
               call terminate("findNextInfoLine", &
               "Unexpected end of Cp curve fit file")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! Get rid of the leading and trailing spaces in string.

       string = adjustl(string)
       string = trim(string)

       ! Check if this is the correct line. If so, exit

       if((len_trim(string) > 0) .and. (string(:1) /= "#")) exit
    enddo

  end subroutine findNextInfoLine

  subroutine setEquationParameters
    !
    !       setEquationParameters sets the number of variables in the
    !       governing equations, the number of turbulent variables, etc.
    !
    use constants
    use paramTurb
    use turbCurveFits
    use flowVarRefState, only : nw, nwf, nt1, nt2, nwt, viscous, &
         eddyModel, kPresent
    use inputPhysics, only : equations, turbModel, wallFunctions, rvfN
    implicit none

    ! Set the number of flow variables to 5, nt1 to 6. This is valid
    ! for all governing equations. Furthermore initialize viscous,
    ! kPresent and eddyModel to .False., which indicates an inviscid
    ! computation. For ns and rans this will be corrected.

    nwf = 5
    nt1 = 6

    viscous    = .false.
    kPresent  = .false.
    eddyModel = .false.

    ! Determine the set of governing equations to solve for and set
    ! the parameters accordingly.

    select case (equations)
    case (EulerEquations)
       nw  = 5
       nt2 = 5

       !===============================================================

    case (NSEquations)
       nw  = 5
       nt2 = 5

       viscous = .true.

       !===============================================================

    case (RANSEquations)

       viscous = .true.

       select case(turbModel)

       case (spalartAllmaras)
          nw  = 6
          nt2 = 6

          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataSa

          !===========================================================

       case (spalartAllmarasEdwards)
          nw  = 6
          nt2 = 6

          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataSae

          !===========================================================

       case (komegaWilcox)
          nw  = 7
          nt2 = 7

          kPresent  = .true.
          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataKw

          !===========================================================

       case (komegaModified)
          nw  = 7
          nt2 = 7

          kPresent  = .true.
          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataKwMod

          !===========================================================

       case (menterSST)
          nw  = 7
          nt2 = 7

          kPresent  = .true.
          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataSST

          !===========================================================

       case (ktau)
          nw  = 7
          nt2 = 7

          kPresent  = .true.
          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataKtau

          !===========================================================

       case (v2f)
          nw  = 9
          nt2 = 9

          rvfLimitK = 1.e-25_realType
          rvfLimitE = 1.e-25_realType

          if(rvfN == 6) then
             rvfCmu   = rvfN6Cmu
             rvfCl    = rvfN6Cl
          else
             rvfCmu   = rvfN1Cmu
             rvfCl    = rvfN1Cl
          endif

          kPresent  = .true.
          eddyModel = .true.
          if( wallFunctions ) call initCurveFitDataVf

       end select

    end select

    ! Determine the number of turbulent variables.

    nwt = nw - nwf

  end subroutine setEquationParameters
  subroutine setStageCoeffExplicitRK
    !
    !       setStageCoeffExplicitRK determines the coefficients of the
    !       stages for the explicit Runge Kutta time integration schemes
    !       for unsteady problems.
    !
    use constants
    use inputUnsteady, only : timeAccuracy, betaRKUnsteady, &
         gammaRKUnsteady, nRKStagesUnsteady
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! Determine the number of Runge Kutta stages as a function of
    ! the accuracy.

    select case (timeAccuracy)
    case (firstOrder)
       nRKStagesUnsteady = 1

    case (secondOrder)
       nRKStagesUnsteady = 2

    case (thirdOrder)
       nRKStagesUnsteady = 3

    case default
       call terminate("setStageCoeffExplicitRK", &
            "No higher order stuff yet")
    end select

    ! Allocate and determine betaRKUnsteady and gammaRKUnsteady.

    allocate(betaRKUnsteady(nRKStagesUnsteady,nRKStagesUnsteady), &
         gammaRKUnsteady(nRKStagesUnsteady), stat=ierr)
    if(ierr /= 0)                               &
         call terminate("setStageCoeffExplicitRK", &
         "Memory allocation failure for betaRKUnsteady &
         &and gammaRKUnsteady.")

    betaRKUnsteady = zero

    select case (timeAccuracy)
    case (firstOrder)

       ! Just the forward Euler time integration scheme.

       betaRKUnsteady(1,1) = 1.0_realType
       gammaRKUnsteady(1)  = 0.0_realType

       !==============================================================

    case (secondOrder)

       ! The TVD Runge Kutta scheme which allows for the maximum
       ! CFL number (1.0).

       betaRKUnsteady(1,1) =  1.0_realType
       betaRKUnsteady(2,1) = -0.5_realType
       betaRKUnsteady(2,2) =  0.5_realType

       gammaRKUnsteady(1)  = 0.0_realType
       gammaRKUnsteady(2)  = 1.0_realType

       !==============================================================

    case (thirdOrder)

       ! Low storage (although not exploited in this implemetation)
       ! 3 stage scheme of Le and Moin.

       betaRKUnsteady(1,1) =   8.0_realType/15.0_realType
       betaRKUnsteady(2,1) = -17.0_realType/60.0_realType
       betaRKUnsteady(2,2) =   5.0_realType/12.0_realType
       betaRKUnsteady(3,2) =  -5.0_realType/12.0_realType
       betaRKUnsteady(3,3) =   3.0_realType/ 4.0_realType

       gammaRKUnsteady(1)  = 0.0_realType
       gammaRKUnsteady(2)  = 8.0_realType/15.0_realType
       gammaRKUnsteady(3)  = 2.0_realType/ 3.0_realType

       ! The TVD Runge Kutta scheme which allows for the maximum
       ! CFL number (1.0).

       ! betaRKUnsteady(1,1) =  1.0_realType
       ! betaRKUnsteady(2,1) = -3.0_realType/ 4.0_realType
       ! betaRKUnsteady(2,2) =  1.0_realType/ 4.0_realType
       ! betaRKUnsteady(3,1) = -1.0_realType/12.0_realType
       ! betaRKUnsteady(3,2) = -1.0_realType/12.0_realType
       ! betaRKUnsteady(3,3) =  2.0_realType/ 3.0_realType

       ! gammaRKUnsteady(1)  = 0.0_realType
       ! gammaRKUnsteady(2)  = 1.0_realType
       ! gammaRKUnsteady(3)  = 0.5_realType

       !==============================================================

    case default
       call terminate("setStageCoeffExplicitRK", &
            "No higher order stuff yet")
    end select

  end subroutine setStageCoeffExplicitRK
  subroutine surfaceVariables(variables)
    !
    !       surfaceVariables extracts from the given string the surface
    !       variables to be written to the solution file.
    !
    use constants
    use extraOutput
    use utils, only : convertToLowerCase, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(inout) :: variables
    !
    !      Local variables.
    !
    integer :: nVarSpecified, pos

    character(len=15) :: keyword
    character(len=maxStringLen) :: errorMessage

    ! Convert the string variables to lower case.

    call convertToLowerCase(variables)

    ! Initialize all the surface output variables to .false.

    surfWriteRho  = .false.
    surfWriteP    = .false.
    surfWriteTemp = .false.
    surfWriteVx   = .false.
    surfWriteVy   = .false.
    surfWriteVz   = .false.
    surfWriteRVx   = .false.
    surfWriteRVy   = .false.
    surfWriteRVz   = .false.

    surfWriteCp       = .false.
    surfWritePtotloss = .false.
    surfWriteMach     = .false.
    surfWriteRMach     = .false.

    surfWriteCf    = .false.
    surfWriteCh    = .false.
    surfWriteYplus = .false.
    surfWriteCfx   = .false.
    surfWriteCfy   = .false.
    surfWriteCfz   = .false.

    surfWriteBlank = .false.
    surfWriteSepSensor = .false.
    surfWriteCavitation = .false.
    surfWriteAxisMoment = .false.
    surfWriteGC = .false.

    ! Initialize nVarSpecified to 0. This serves as a test
    ! later on.

    nVarSpecified = 0

    ! Loop to extract the info from the string variables.

    do
       ! Condition to exit the loop.

       if(len_trim(variables) == 0) exit

       ! Locate the first occurance of the _ in the string and
       ! determine the string keyword.

       pos = index(variables, "_")
       if(pos == 0) then
          keyword   = variables
          variables = ""
       else
          keyword   = variables(:pos-1)
          variables = variables(pos+1:)
       endif

       ! Check the keyword.

       select case (keyword)
       case ("")
          ! Multiple occurence of "_". Just ignore it.

       case ("rho")
          surfWriteRho = .true.
          nVarSpecified = nVarSpecified + 1

       case ("p")
          surfWriteP = .true.
          nVarSpecified = nVarSpecified + 1

       case ("temp")
          surfWriteTemp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vx")
          surfWriteVx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vy")
          surfWriteVy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vz")
          surfWriteVz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvx")
          surfWriteRVx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvy")
          surfWriteRVy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvz")
          surfWriteRVz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cp")
          surfWriteCp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("ptloss")
          surfWritePtotloss = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mach")
          surfWriteMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rmach")
          surfWriteRMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cf")
          surfWriteCf = .true.
          nVarSpecified = nVarSpecified + 1

       case ("ch")
          surfWriteCh = .true.
          nVarSpecified = nVarSpecified + 1

       case ("yplus")
          surfWriteYplus = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cfx")
          surfWriteCfx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cfy")
          surfWriteCfy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cfz")
          surfWriteCfz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("blank")
          surfWriteBlank = .true.
          nVarSpecified = nVarSpecified + 1

       case ("sepsensor")
          surfWriteSepSensor = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cavitation")
          surfWriteCavitation = .true.
          nVarSpecified = nVarSpecified + 1

       case ("axismoment")
          surfWriteAxisMoment = .true.
          nVarSpecified = nVarSpecified + 1

       case ("gc")
          surfWriteGC = .True.
          nVarSpecified = nVarSpecified + 1

       case default
          pos = len_trim(keyword)
          write(errorMessage,"(3a)") "Unknown surface output &
               &variable, ", trim(keyword), &
               ", specified"
          call terminate("surfaceVariables", errorMessage)

       end select

    enddo

    ! Set surfaceOutSpecified to .true. if variables were specified.
    ! If not, later on the defaults will be set.

    if(nVarSpecified > 0) surfaceOutSpecified = .true.

  end subroutine surfaceVariables

  subroutine volumeVariables(variables)
    !
    !       volumeVariables extracts from the given string the extra
    !       volume variables to be written to the solution file.
    !
    use constants
    use extraOutput
    use utils, only : convertToLowerCase, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(inout) :: variables
    !
    !      Local variables.
    !
    integer :: nVarSpecified, pos

    character(len=15) :: keyword
    character(len=maxStringLen) :: errorMessage

    ! Convert the string variables to lower case.

    call convertToLowerCase(variables)

    ! Initialize all the volume output variables to .False.

    volWriteMx    = .false.
    volWriteMy    = .false.
    volWriteMz    = .false.
    volWriteRhoe  = .false.
    volWriteTemp  = .false.
    volWriteVort  = .false.
    volWriteVortx = .false.
    volWriteVorty = .false.
    volWriteVortz = .false.

    volWriteCp       = .false.
    volWriteMach     = .false.
    volWriteMachTurb = .false.
    volWritePtotloss = .false.

    volWriteEddyVis      = .false.
    volWriteRatioEddyVis = .false.
    volWriteDist         = .false.

    volWriteResRho  = .false.
    volWriteResMom  = .false.
    volWriteResRhoe = .false.
    volWriteResTurb = .false.

    volWriteShock   = .false.
    volWriteFilteredShock = .false.

    volWriteBlank = .false.
    volWriteGC = .false.
    volWriteStatus = .false.
    volWriteIntermittency = .false.


    ! Initialize nVarSpecified to 0. This serves as a test
    ! later on.

    nVarSpecified = 0

    ! Loop to extract the info from the string variables.

    do
       ! Condition to exit the loop.

       if(len_trim(variables) == 0) exit

       ! Locate the first occurance of the _ in the string and
       ! determine the string keyword.

       pos = index(variables, "_")
       if(pos == 0) then
          keyword   = variables
          variables = ""
       else
          keyword   = variables(:pos-1)
          variables = variables(pos+1:)
       endif

       ! Check the keyword.

       select case (keyword)
       case ("")
          ! Multiple occurence of "_". Just ignore it.

       case ("mx")
          volWriteMx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("my")
          volWriteMy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mz")
          volWriteMz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvx")
          volWriteRVx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvy")
          volWriteRVy = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rvz")
          volWriteRVz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rhoe")
          volWriteRhoe = .true.
          nVarSpecified = nVarSpecified + 1

       case ("temp")
          volWriteTemp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vort")
          volWriteVort = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vortx")
          volWriteVortx = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vorty")
          volWriteVorty = .true.
          nVarSpecified = nVarSpecified + 1

       case ("vortz")
          volWriteVortz = .true.
          nVarSpecified = nVarSpecified + 1

       case ("cp")
          volWriteCp = .true.
          nVarSpecified = nVarSpecified + 1

       case ("mach")
          volWriteMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("rmach")
          volWriteRMach = .true.
          nVarSpecified = nVarSpecified + 1

       case ("macht")
          volWriteMachTurb = .true.
          nVarSpecified = nVarSpecified + 1

       case ("ptloss")
          volWritePtotloss = .true.
          nVarSpecified = nVarSpecified + 1

       case ("eddy")
          volWriteEddyVis = .true.
          nVarSpecified = nVarSpecified + 1

       case ("eddyratio")
          volWriteRatioEddyVis = .true.
          nVarSpecified = nVarSpecified + 1

       case ("dist")
          volWriteDist = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resrho")
          volWriteResRho = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resmom")
          volWriteResMom = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resrhoe")
          volWriteResRhoe = .true.
          nVarSpecified = nVarSpecified + 1

       case ("resturb")
          volWriteResTurb = .true.
          nVarSpecified = nVarSpecified + 1

       case ("shock")
          volWriteShock = .true.
          nVarSpecified = nVarSpecified + 1

       case ("filteredshock")
          volWriteFilteredShock = .true.
          nVarSpecified = nVarSpecified + 1

       case ("blank")
          volWriteBlank = .true.
          nVarSpecified = nVarSpecified + 1

       case ("gc")
          volWriteGC = .true.
          nVarSpecified = nVarSpecified + 1

       case ("status")
          volWriteStatus = .true.
          nVarSpecified = nVarSpecified + 1

       case("intermittency")
          volWriteIntermittency = .true.
          nVarSpecified = nVarSpecified + 1

       case default
          pos = len_trim(keyword)
          write(errorMessage,"(3a)" ) "Unknown extra volume output &
               &variable, ", trim(keyword), &
               ", specified"
          call terminate("volumeVariables", errorMessage)

       end select

    enddo

    ! Set volumeOutSpecified to .true. if variables were specified.
    ! If not, later on the defaults will be set.

    if(nVarSpecified > 0) volumeOutSpecified = .true.

  end subroutine volumeVariables

  subroutine checkInputParam
    !
    !       checkInputParam checks if all necessary data has been
    !       specified. If some key data is missing an error message will
    !       be printed and the program will exit. Key data depends on the
    !       case to be solved. E.g. for the Navier Stokes equations it is
    !       necessary to specify the Reynolds number, but for Euler this
    !       can be omitted.
    !       Furthermore warnings are printed in case parameters have been
    !       specified that are ignored, e.g. Mach number for internal flow
    !       computations.
    !       Note that only processor 0 prints warning and error messages,
    !       such that the output does not become messy.
    !
    use constants
    !  --------- Bare imports...too many to list -------
    use inputDiscretization
    use inputIO
    use inputIteration
    use inputMotion
    use inputOverset
    use inputParallel
    use inputPhysics
    use inputTimeSpectral
    use inputUnsteady
    use inputADjoint
    use inputTSStabDeriv
    !  ------------------------------------------------
    use communication, only : myid, adflow_comm_world
    use iteration, only : coefTime, coefTimeALE, coefMeshALE, &
         oldSolWritten, nALEMeshes, nALESteps, nOldLevels
    use monitor, only : nTimeStepsRestart
    use utils, only : terminate
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    integer(kind=intType) :: nn, oldSolWrittenSize

    real(kind=realType) :: vecLength, dot

    logical :: gridPrecisionWarning, solPrecisionWarning

    !       Discretization parameters. Check if the key parameters have
    !       been specified and set some coarse grid parameters in case
    !       these have not been specified.
    !
    if(spaceDiscr == none) then
       if(myID == 0)                       &
            call terminate("checkInputParam", &
            "Discretization scheme not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(spaceDiscrCoarse == none) spaceDiscrCoarse = spaceDiscr

    if(riemannCoarse == none) riemannCoarse = riemann

    ! Set dirScaling to .false. if a scheme other than scalar
    ! dissipation is used.

    if(spaceDiscr /= dissScalar) dirScaling = .false.

    ! Determine whether or not the spectral radIi are needed for
    ! the flux computations.

    radiiNeededFine = .false.
    if(spaceDiscr == dissScalar) radiiNeededFine = .true.

    radiiNeededCoarse = .false.
    if(spaceDiscrCoarse == dissScalar) radiiNeededCoarse = .true.
    !
    !       IO parameters. Check if the grid file has been specified
    !       Possibly correct the
    !       value of restart. Note that restart got the default value of
    !       .true. in case no restart file has been specified it is now
    !       set to false. Set the names of the solution files if not
    !       specified and check if a cp curve fit file has been specified
    !       if curve fits must be used.
    !       If the code has been compiled without cgns check that the file
    !       format is not cgns.
    !       Overwrite storeConvInnerIter to .true. if this is not an
    !       unsteady computation.
    !
    if(gridFile == "") then
       if(myID == 0) &
            call terminate("checkInputParam", "Grid file not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(newGridFile == "") then
       newGridFile = "NewGrid.cgns"
    endif

    if(solFile == "") then
       solFile = "SolADflow.cgns"
    endif

    if(surfaceSolFile == "") &
         surfaceSolFile = trim(solfile)//"Surface"

    if(cpModel == cpTempCurveFits .and. cpFile == "") then
       if(myID == 0)                        &
            call terminate("checkInputParam", &
            "Cp curve fit file not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

#ifdef USE_NO_CGNS

    if(fileFormatRead  == cgnsFormat .or. &
         fileFormatWrite == cgnsFormat) then
       if(myID == 0)                       &
            call terminate("checkInputParam", &
            "cgns support disabled during compile time")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

#endif

    if(equationMode == unsteady) then
       if(timeIntegrationScheme == explicitRK) &
            storeConvInnerIter = .false.
    else
       storeConvInnerIter = .true.
    endif
    !
    !       Iteration parameters. Check if the key parameters have specified
    !       been and set some coarse grid parameters in case these
    !       have not been specified.
    !
    if(equationMode          == unsteady .and. &
         timeIntegrationScheme == explicitRK) then
       smoother = none
    else
       if(smoother == none) then
          if(myID == 0) &
               call terminate("checkInputParam", "Smoother not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(ncycles < 0) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Number of multigrid cycles not or wrongly &
               &specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(cfl < zero) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "cfl number not or wrongly specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(l2Conv <= zero .or. L2Conv >= one) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Relative L2 norm for convergence must be a &
               & number between 0 and 1.")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(l2ConvCoarse <= zero .or. L2ConvCoarse >= one) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Relative L2 norm for convergence coarse grid &
               &must be a number between 0 and 1.")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif
    endif
    !
    !       Grid motion parameters. These can only be specified for an
    !       external flow problem.
    !
    if(flowType == internalFlow .and. gridMotionSpecified) then
       if(myID == 0) &
            call terminate("checkInputParam", &
            "Grid motion specified for an internal flow; &
            &this is not possible")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif
    !
    !       Physics parameters. Check if the key parameters have been
    !       specified and set the unit vector for the free-stream velocity.
    !
    if(equations == none) then
       if(myID == 0) &
            call terminate("checkInputParam", "Equations not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(equationMode == none) then
       if(myID == 0) &
            call terminate("checkInputParam", "Mode not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(flowType == none) then
       if(myID == 0) &
            call terminate("checkInputParam", "Flow type not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(Mach < zero .and. flowType == externalFlow) then
       if(myID == 0)                        &
            call terminate("checkInputParam", &
            "Mach not or wrongly specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    if(equations == RANSEquations .and. turbModel == none) then
       if(myID == 0)                        &
            call terminate("checkInputParam", &
            "Turbulence model not specified")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Create a unit vector for the free stream velocity. It is checked
    ! if the vector specified is a valid one. If not processor 0 prints
    ! an error message. Only for external flows.

    if(flowType == externalFlow) then
       vecLength = sqrt(velDirFreestream(1)*velDirFreestream(1) &
            +      velDirFreestream(2)*velDirFreestream(2) &
            +      velDirFreestream(3)*velDirFreestream(3))
       if(vecLength < eps) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Free stream velocity direction wrongly &
               &specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       vecLength = one/vecLength
       velDirFreestream(1) = velDirFreestream(1)*vecLength
       velDirFreestream(2) = velDirFreestream(2)*vecLength
       velDirFreestream(3) = velDirFreestream(3)*vecLength
    else
       ! Internal flow; simply reset the velocity direction. The value
       ! will be determined later from the inflow boundary conditions.

       velDirFreestream(1) = one
       velDirFreestream(2) = zero
       velDirFreestream(3) = zero
    endif

    ! Set the drag direction to the velocity direction.

    dragDirection = velDirFreestream

    ! Check the lift direction if it was specified for an external
    ! flow. Otherwise set the default direction.

    if(liftDirSpecified .and. flowType == externalFlow) then

       ! Create a unit vector. Perform the same check as for
       ! for the free stream velocity direction.

       vecLength = sqrt(liftDirection(1)*liftDirection(1) &
            +      liftDirection(2)*liftDirection(2) &
            +      liftDirection(3)*liftDirection(3))
       if(vecLength < eps) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Lift direction wrongly specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       vecLength = one/vecLength
       liftDirection(1) = liftDirection(1)*vecLength
       liftDirection(2) = liftDirection(2)*vecLength
       liftDirection(3) = liftDirection(3)*vecLength

       ! Check the orthogonality with the drag direction.

       dot = liftDirection(1)*dragDirection(1) &
            + liftDirection(2)*dragDirection(2) &
            + liftDirection(3)*dragDirection(3)

       if(abs(dot) > 1.e-3_realType) then
          if(myID == 0)                       &
               call terminate("checkInputParam", &
               "Lift direction not orthogonal to &
               &free-stream")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

    else

       ! Lift direction not specified. Set the default direction.
       ! It will have a zero component in the y-direction and a positive
       ! one in the z-direction.

       liftDirection(1) = -dragDirection(3)
       liftDirection(2) =  zero
       liftDirection(3) =  dragDirection(1)

       if(liftDirection(3) < zero) then
          liftDirection(1) = -liftDirection(1)
          liftDirection(3) = -liftDirection(3)
       endif
    endif

    ! Set the Mach number for the coefficients equal to the Mach
    ! number if it was not specified. For internal flow field this
    ! will again be changed in initFlo.

    if(MachCoef < zero) MachCoef = Mach
    !
    !       Time spectral parameters. They only need to be specified for a
    !       time spectral computation.
    !
    testSpectral: if(equationMode == timeSpectral) then

       ! Check if the number of time intervals was specified.

       if(nTimeIntervalsSpectral < 0) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Number time intervals spectral not or &
               &wrongly specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! If an unsteady restart solution file must be written, check
       ! if the corresponding time step has been specified.

       if( writeUnsteadyRestartSpectral ) then
          if(dtUnsteadyRestartSpectral <= zero) then
             if(myID == 0)                        &
                  call terminate("checkInputParam", &
                  "Time step (in sec) for unsteady restart &
                  &not or wrongly specified.")
             call mpi_barrier(ADflow_comm_world, ierr)
          endif
       endif

       ! If solution files (for postprocessing) must be written,
       ! check if the number has been specified.

       if( writeUnsteadyVolSpectral .or. &
            writeUnsteadySurfSpectral) then
          if(nunsteadySolSpectral <= 0) then
             if(myID == 0)                        &
                  call terminate("checkInputParam", &
                  "Number of unsteady solution files &
                  &not or wrongly specified.")
             call mpi_barrier(ADflow_comm_world, ierr)
          endif
       endif

    else testSpectral

       ! No spectral method. Set nTimeIntervalsSpectral to 1.

       nTimeIntervalsSpectral = 1

    endif testSpectral
    !
    !       Unsteady parameters. They only need to be specified for an
    !       unsteady computation.
    !
    testUnsteady: if(equationMode == unsteady) then

       ! Physical time step parameters.

       if(nTimeStepsFine < 0) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Number of unsteady time steps fine grid &
               &not or wrongly specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(nTimeStepsCoarse < 0) nTimeStepsCoarse = nTimeStepsFine

       if(deltaT < 0) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Unsteady time step (in sec) &
               &not or wrongly specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! Check if the rigid body rotation parameters are consistent.
       ! The polynomial rotation coefficients.

       if(degreePolXRot >= 0 .and. &
            .not. allocated(coefPolXRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Polynomial coefficients x-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreePolYRot >= 0 .and. &
            .not. allocated(coefPolYRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Polynomial coefficients y-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreePolZRot >= 0 .and. &
            .not. allocated(coefPolZRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Polynomial coefficients z-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! The fourier rotation coefficients.

       if(degreeFourXRot >= 0 .and. &
            .not. allocated(cosCoefFourXRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier cosine coefficients x-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreeFourXRot >= 1 .and. &
            .not. allocated(sinCoefFourXRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier sine coefficients x-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreeFourYRot >= 0 .and. &
            .not. allocated(cosCoefFourYRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier cosine coefficients y-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreeFourYRot >= 1 .and. &
            .not. allocated(sinCoefFourYRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier sine coefficients y-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreeFourZRot >= 0 .and. &
            .not. allocated(cosCoefFourZRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier cosine coefficients z-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       if(degreeFourZRot >= 1 .and. &
            .not. allocated(sinCoefFourZRot)) then
          if(myID == 0)                        &
               call terminate("checkInputParam", &
               "Fourier sine coefficients z-rotation &
               &not specified")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

    endif testUnsteady
    !
    !                             Warning messages.
    !
    ! Check for an invisid problem if the Reynolds number is specified.
    ! If so, print a Warning that this info is ignored.

    if(myID == 0 .and. equations == EulerEquations .and. &
         Reynolds > zero) then

       print "(a)", "#"
       print "(a)", "#                      Warning"
       print "(a)", "# Reynolds number specified for the Euler &
            &equations."
       print "(a)", "# This information is ignored."
       print "(a)", "#"

    endif

    ! Check if the Mach and Reynolds number are specified for an
    ! internal flow problem. If so, print a Warning message that this
    ! info is ignored.

    if(flowType == internalFlow) then

       ! Check whether a viscous or an inviscid problem is to be solved.
       ! For an inviscid problem you do not want to mention that the
       ! Reynolds number is ignored, because this has already been
       ! taken care of.

       if((equations == NSEquations .or.     &
            equations == RANSEquations) .and. &
            Mach > zero .and. Reynolds > zero) then

          ! Viscous problem, where both the Mach and Reynolds were
          ! specified. Processor 0 prints the Warning.

          if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print "(a)", "# Mach and Reynolds number specified &
                  &for an internal flow problem."
             print "(a)", "# This information is ignored."
             print "(a)", "#"
          endif

       else if(Mach > zero) then

          ! The Mach number has been specified. Processor 0 prints
          ! a Warning.

          if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print "(a)", "# Mach number specified for an internal &
                  &flow problem."
             print "(a)", "# This information is ignored."
             print "(a)", "#"
          endif

       endif

    endif

    ! For a steady computation possible specified rigid body
    ! rotation info is ignored. Processor 0 will print the Warning.

    if(degreePolXRot  >= 0 .or. degreePolYRot  >= 0 .or. &
         degreePolZRot  >= 0 .or. degreeFourXRot >= 0 .or. &
         degreeFourYRot >= 0 .or. degreeFourZRot >= 0) then

       if(equationMode == steady .and. myID == 0) then
          print "(a)", "#"
          print "(a)", "#                      Warning"
          print "(a)", "# Rigid body rotation info specified for &
               &a steady computation."
          print "(a)", "# This information is ignored."
          print "(a)", "#"
       endif
    endif

    ! Print warning messages if the precision to be written
    ! is larger than the precision used in the computation.

    gridPrecisionWarning = .false.
    solPrecisionWarning  = .false.

#ifdef USE_SINGLE_PRECISION
    if(precisionGrid == precisionDouble) gridPrecisionWarning = .true.
    if(precisionSol  == precisionDouble) solPrecisionWarning = .true.
#endif

    if(gridPrecisionWarning .and. myID == 0) then
       print "(a)", "#"
       print "(a)", "#                      Warning"
       print "(a)", "# Precision of the grid file to write is &
            &bigger than used in the computation."
       print "(a)", "# This does not make sense and is a waste &
            &of disk space"
       print "(a)", "#"
    endif

    if(solPrecisionWarning .and. myID == 0) then
       print "(a)", "#"
       print "(a)", "#                      Warning"
       print "(a)", "# Precision of the solution file to write is &
            &bigger than used in the computation."
       print "(a)", "# This does not make sense and is a waste &
            &of disk space"
       print "(a)", "#"
    endif
    !
    !       Wall functions can only be used if the RANS equations are to
    !       be solved. If no wall functions are used the wall offset is
    !       set to zero.
    !
    if(equations /= RANSEquations) wallFunctions = .false.
    if(.not. wallFunctions) wallOffset = zero
    !
    !       Check whether or not the wall distance is needed for the
    !       turbulence model.
    !
    if(equations == RANSEquations) then

       ! RANS simulation. Determine if the turbulence model is
       ! wall distance free. Note that updateWallDistanceUnsteady is
       ! NOT overruled, because this is just the case for which this
       ! parameter was intended.

       select case (turbModel)
       case (komegaWilcox, komegaModified, ktau)

          ! Wall distance free turbulence models.

          wallDistanceNeeded = .false.

          !=============================================================

       case default

          ! The turbulence model needs the wall distance

          wallDistanceNeeded = .true.

       end select

    else

       ! Laminar or inviscid computation. Simply initialize the
       ! logicals for the wall distance to .false.

       wallDistanceNeeded         = .false.
       updateWallDistanceUnsteady = .false.

    endif
    !
    !       Parallelization parameters. Set the minimum load imbalance to
    !       3 percent to avoid any problems.
    !
    loadImbalance = max(loadImbalance, 0.03_realType)
    !
    !       Some default parameters, which depend on other parameters.
    !       Only if these have not been specified of course.
    !
    if(nsgStartup < 0)    nsgStartup    = 0
    if(ncyclesCoarse < 0) nCyclesCoarse = nCycles
    if(cflCoarse < zero)  cflCoarse     = cfl
    if(betaTurb  < zero)  betaTurb      = alfaTurb

    if(turbRelax == turbRelaxNotDefined) then
       turbRelax = turbRelaxImplicit
       if(turbModel == v2f) turbRelax = turbRelaxExplicit
    endif

    ! V2f should only be solved with explicit underrelaxation.

    if(equations == RANSEquations .and. turbModel == v2f .and. &
         turbRelax == turbRelaxImplicit) then

       turbRelax = turbRelaxExplicit

       if(myID == 0) then
          print "(a)", "#"
          print "(a)", "#                      Warning"
          print "(a)", "# Implicit underrelaxation specified for &
               &the v2f model."
          print "(a)", "# This is overwritten to explicit &
               &underrelaxation."
          print "(a)", "#"
       endif

    endif

    if(nsaveVolume <= 0) then
       select case (equationMode)
       case (steady, timeSpectral)
          nSaveVolume = nCycles + nCyclesCoarse + nsgStartup + 1

       case (unsteady)
          nSaveVolume = nTimeStepsFine + nTimeStepsCoarse &
               + nTimeStepsRestart + 1
       end select
    endif

    if(nsaveSurface <= 0) nSaveSurface  = nSaveVolume

    if(eddyVisInfRatio < zero) then

       ! Default value depends on the turbulence model.

       select case (turbModel)

       case (spalartAllmaras, spalartAllmarasEdwards)
          eddyVisInfRatio = 0.009_realType

       case default
          eddyVisInfRatio = 0.1_realType

       end select
    endif
    !
    !       Determine the number of old grid levels needed for the BDF
    !       time integration of unsteady problems and allocate the memory
    !       for the coefficients. The actual values are not yet set,
    !       because in the first (and possibly second) time step a reduced
    !       order must be used, because the older states are not available
    !       yet. Also allocate the memory for the logicals to indicate
    !       whether or not old solutions have been written.
    !       If a Runge Kutta scheme must be used for the time integration,
    !       either explicit or implicit, a separate routine is called to
    !       set all the necessary variables.
    !
    select case (timeIntegrationScheme)
    case (BDF, MD)

       ! First check if the accuracy is okay.

       if(timeAccuracy > thirdOrder) then
          if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print "(a)", "# Maximum third order possible for BDF."
             print "(a)", "# Order has been reduced to third."
             print "(a)", "#"
          endif

          timeAccuracy = thirdOrder
       endif

       ! Determine the accuracy and set nOldLevels accordingly.

       select case (timeAccuracy)
       case (firstOrder)
          nOldLevels = 1

       case (secondOrder)
          nOldLevels = 2

       case (thirdOrder)
          nOldLevels = 3
       end select

       ! Allocate the memory for coefTime.
       if( allocated(coefTime)) deallocate(coefTime)
       allocate(coefTime(0:nOldLevels), stat=ierr)
       if(ierr /= 0)                       &
            call terminate("checkInputParam", &
            "Memory allocation error for coefTime")

       ! Determine the accuracy and set ALE parameters accordingly.
       if (useALE) then
          select case (timeAccuracy)
          case (firstOrder)
             nALEMeshes = 1
             nALESteps  = 2

          case (secondOrder)
             nALEMeshes = 2
             nALESteps  = 4

          case (thirdOrder)
             call terminate("checkInputParam", &
                  "ALE can only use 1st and 2nd order time accuracy")
          end select

          if( allocated(coefTimeALE)) deallocate(coefTimeALE)
          allocate(coefTimeALE(1:nALEsteps), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("checkInputParam", &
               "Memory allocation error for coefTimeALE")

          if( allocated(coefMeshALE)) deallocate(coefMeshALE)
          allocate(coefMeshALE(1:nALEMeshes,2), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("checkInputParam", &
               "Memory allocation error for coefMeshALE")
       end if

       !===============================================================

    case (explicitRK)
       nOldLevels = 1
       call setStageCoeffExplicitRK


    end select

    ! Set the logicals whether or not the old solutions have been
    ! written. Note that this is only used for the second and
    ! higher order BDF schemes. However it is allocated with a
    ! minimum size of 1 to avoid problems.

    oldSolWrittenSize = max(nOldLevels-1_intType, 1_intType)

    !check allocations for multipile succesive calls
    if (allocated(oldSolWritten)) deallocate(oldSolWritten)
    allocate(oldSolWritten(oldSolWrittenSize), stat=ierr)
    if(ierr /= 0)                       &
         call terminate("checkInputParam", &
         "Memory allocation error for oldSolWritten")

    do nn=1,oldSolWrittenSize
       oldSolWritten(nn) = .false.
    enddo
    !
    !       Determine the values of the runge kutta parameters, depending
    !       on the number of stages specified.
    !
    ! Limit the number of stages between 1 and 6 and allocate the
    ! memory.

    nRKStages = min(6_intType,max(1_intType,nRKStages))

    !check allocations for multipile succesive calls
    if (allocated(etaRk)) deallocate(etaRk)
    if (allocated(cdisRK)) deallocate(cdisRK)

    allocate(etaRk(nRKStages), cdisRK(nRKStages), stat=ierr)
    if(ierr /= 0) &
         call terminate("checkInputParam", &
         "Memory allocation error for etaRK and cdisRK")

    ! Determine the case we are having here.

    select case (nRKStages)
    case (1_intType)
       etaRK(1) = one

       cdisRK(1) = one

    case (2_intType)
       etaRK(1) = 0.2222_realType
       etaRK(2) = one

       cdisRK(1) = one
       cdisRK(2) = one

    case (3_intType)
       etaRK(1) = 0.2846_realType
       etaRK(2) = 0.6067_realType
       etaRK(3) = one

       cdisRK(1) = one
       cdisRK(2) = one
       cdisRK(3) = one

    case (4_intType)
       etaRK(1) = 0.33333333_realType
       etaRK(2) = 0.26666667_realType
       etaRK(3) = 0.55555555_realType
       etaRK(4) = one

       cdisRK(1) = one
       cdisRK(2) = half
       cdisRK(3) = zero
       cdisRK(4) = zero

    case (5_intType)
       etaRK(1) = fourth
       etaRK(2) = 0.16666667_realType !1/6
       etaRK(3) = 0.37500000_realType !3/8
       etaRK(4) = half
       etaRK(5) = one

       cdisRK(1) = one
       cdisRK(2) = zero
       cdisRK(3) = 0.56_realType
       cdisRK(4) = zero
       cdisRK(5) = 0.44_realType

    case (6_intType)
       etaRK(1) = 0.0722_realType
       etaRK(2) = 0.1421_realType
       etaRK(3) = 0.2268_realType
       etaRK(4) = 0.3425_realType
       etaRK(5) = 0.5349_realType
       etaRK(6) = one

       cdisRK(1) = one
       cdisRK(2) = one
       cdisRK(3) = one
       cdisRK(4) = one
       cdisRK(5) = one
       cdisRK(6) = one
    end select
    !
    !       To avoid any problems later on, allocate the memory for the
    !       rigid body motion parameters if these values were not present
    !       in the parameter file.
    !
    if(.not. allocated(coefPolXRot) ) then
       allocate(coefPolXRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for coefPolXRot")
       coefPolXRot = zero
    endif

    if(.not. allocated(coefPolYRot) ) then
       allocate(coefPolYRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for coefPolYRot")
       coefPolYRot = zero
    endif

    if(.not. allocated(coefPolZRot) ) then
       allocate(coefPolZRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for coefPolZRot")
       coefPolZRot = zero
    endif

    if(.not. allocated(cosCoefFourXRot) ) then
       allocate(cosCoefFourXRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &cosCoefFourXRot")
       cosCoefFourXRot = zero
    endif

    if(.not. allocated(sinCoefFourXRot) ) then
       allocate(sinCoefFourXRot(1), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &sinCoefFourXRot")
       sinCoefFourXRot = zero
    endif

    if(.not. allocated(cosCoefFourYRot) ) then
       allocate(cosCoefFourYRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &cosCoefFourYRot")
       cosCoefFourYRot = zero
    endif

    if(.not. allocated(sinCoefFourYRot) ) then
       allocate(sinCoefFourYRot(1), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &sinCoefFourYRot")
       sinCoefFourYRot = zero
    endif

    if(.not. allocated(cosCoefFourZRot) ) then
       allocate(cosCoefFourZRot(0:0), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &cosCoefFourZRot")
       cosCoefFourZRot = zero
    endif

    if(.not. allocated(sinCoefFourZRot) ) then
       allocate(sinCoefFourZRot(1), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("checkInputParam", &
            "Memory allocation failure for &
            &sinCoefFourZRot")
       sinCoefFourZRot = zero
    endif

  end subroutine checkInputParam
  subroutine setDefaultValues
    !
    !       setDefaultValues sets the default values for the input
    !       parameters where-ever possible. The parameters that must be
    !       set by the user are initialized such a check can be performed
    !       later.
    !
    use constants

    !  --------- Bare imports...too many to list -------
    use inputDiscretization
    use inputIO
    use inputIteration
    use inputMotion
    use inputOverset
    use inputParallel
    use inputPhysics
    use inputTimeSpectral
    use inputUnsteady
    use inputADjoint
    use inputTSStabDeriv
    !  ------------------------------------------------
    use flowVarRefState, only : Lref, lRefSpecified, pRef, rhoRef, &
         TinfDim, Tref
    use iteration, only : nOldSolAvail, timeSpectralGridsNotWritten
    use monitor, only : monMassSliding, nTimeStepsRestart, timeUnsteadyRestart
    use killSignals, only : fatalFail, routineFailed
    use ADjointPETSc, only : adjointPETScVarsAllocated
    use inputCostFunctions
    implicit none

    ! Initialize monitoring the turbulent residuals as well as the
    ! monitoring of mass flow of the sliding interfaces to .false.

    monDturb       = .false.
    monMassSliding = .false.

    ! Initialize the logicals to check whether or not monitoring,
    ! surface output and volume output variables were specified to
    ! .false.

    monitorSpecified    = .false.
    surfaceOutSpecified = .false.
    volumeOutSpecified  = .false.
    isoOutSpecified     = .false.
    !
    !       Set the default values for the discretization parameters.
    !
    spaceDiscr = none                  ! Serves as a check later on.
    orderTurb  = firstOrder            ! First order discretization.
    ! Of turbulent advective terms.
    riemann     = Roe
    limiter     = noLimiter            ! No limiter in upwind schemes.
    precond     = noPrecond            ! No preconditioning.

    eulerWallBCTreatment = normalMomentum   ! Normal momentum equation is
    ! Used to determine ghost
    ! cell pressure.

    viscWallBCTreatment = constantPressure  ! Normal momentum equation is
    ! Used to determine ghost
    ! cell pressure.

    outflowTreatment = constantExtrapol ! Constant extrapolation at
    ! outflow boundaries.

    spaceDiscrCoarse = none             ! Serves as a check. If nothing
    riemannCoarse    = none             ! is specified the fine grid
    ! parameter is taken.

    nonMatchTreatment = NonConservative ! Non conservative treatment
    ! of non-matching block to
    ! block boundaries.

    vortexCorr = .false.                ! No vortex correction is
    ! applied.

    vis2       = half
    vis4       = one/64.0_realType
    vis2Coarse = half

    dirScaling = .true.                 ! Apply isotropic directional
    adis       = two*third              ! scaling in the artificial
    ! dissipation schemes.

    hScalingInlet = .false.             ! No total enthalpy scaling.

    kappaCoef = third
    !
    !       Set the default values for the IO-parameters.

    gridFile       = ""          ! Serves as a check later on.

    checkRestartSol = .true.     ! Restart solution is checked for
    ! correct nonDimensionalization.

    newGridFile = ""             ! This will be corrected later on
    solFile     = ""             ! if nothing is specified. The
    ! default names depend on the
    ! format used

    surfaceSolFile = ""          ! This will be corrected later if no
    ! surface solution file is specified.

    storeRindLayer = .True.     ! No halo cells in solution files.

    autoParameterUpdate = .true. ! Update the input parameter file
    ! when a restart file is written.
    writeCoorMeter = .false.     ! Use original coordinate units
    ! when writing solution files.

    cpFile = ""                  ! Serves as a check later on.

    storeConvInnerIter = .false. ! Do not store the convergence of
    ! the inner iterations in unsteady
    ! mode.

#ifdef USE_SINGLE_PRECISION
    precisionGrid = precisionSingle   ! Default IO precision depends
    precisionSol  = precisionSingle   ! on the default floating
                                      ! point type used. Note that
#else
    precisionGrid = precisionDouble   ! for quadrupole precision the
    precisionSol  = precisionDouble   ! IO takes place in double
                                      ! precision.
#endif

    ! Surface solution defaults to single precision
    precisionSurfGrid = precisionSingle
    precisionSurfSol  = precisionSingle

    !
    !       Set the default values for the iteration parameters.
    !
    nCycles       = -1    ! Serves as a check later on.
    nsgStartup    =  0    ! No single grid startup iterations.
    nSubIterTurb  =  0    ! No additional turbulent subiterations.
    nUpdateBleeds = 50    ! Update the bleeds every 50 iterations.

    nSaveVolume  = 1      ! Only save at the end of the computation.
    nSaveSurface = 1

    smoother  = none
    nRKStages = 5
    nSubiterations = 1

    !resAveraging =  noResAveraging ! No residual averaging.
    resAveraging =  noResAveraging
    smoop        = 1.5_realType

    turbTreatment     = decoupled      ! Decoupled solver for the
    ! turbulent equations
    turbSmoother      = adi            ! solved using an adi scheme.
    freezeTurbSource = .true.          ! Freeze the coarse grid source
    ! terms for a coupled solver.
    turbRelax = turbRelaxNotDefined    ! Will be set later, depending
    ! on the turbulence model.

    cfl = -one                         ! Serves as a check later on.

    relaxBleeds = 0.1_realType         ! Relaxation factor for the
    ! bleed boundary conditions.

    alfaTurb =  0.8_realType
    betaTurb = -one                    ! Serves as a check later on.

    L2Conv        = 1.e-6_realType     ! Six orders of magnitude for
    ! convergence.
    L2ConvCoarse = 1.e-2_realType      ! Only two on coarse grids in
    ! full mg.

    maxL2DeviationFactor = 1_realType
    nCyclesCoarse = -1             ! If these parameters are not
    cflCoarse     = -one           ! specified the corresponding fine
    ! grid values are taken.

    fcoll = one       ! No relaxation when restricting the residuals.

    mgBoundCorr = bcDirichlet0 ! Zero out the boundary halo's for
    ! the multigrid corrections.

    mgStartlevel = -1    ! Start at the coarsest grid of the mg cycle
    ! when no restart is performed.
    mgDescription = "sg" ! Single grid computation.
    !
    !       Set the default values for the motion parameters,
    !       i.e. no motion.
    !
    ! Translation data.


    ! Rotation data.

    rotPoint = zero

    degreePolXRot = -1      ! -1, because the start index is 0.
    degreePolYRot = -1
    degreePolZRot = -1

    degreeFourXRot = -1     ! -1, because the start index is 0,
    ! at least of the cosine part.
    degreeFourYRot = -1
    degreeFourZRot = -1

    omegaFourXRot = zero
    omegaFourYRot = zero
    omegaFourZRot = zero

    ! The logical to determine whether or not a motion is specified.
    ! Initialize it to .false.

    gridMotionSpecified = .false.
    !
    !       Set the default values for the parallel parameters.
    !
    loadImbalance = 0.1_realType  ! Allow 10 percent load imbalance.
    splitBlocks   = .true.        ! Allow the splitting of blocks to
    ! obtain a better load balancing.
    loadbalanceiter = 2           ! Do two iterations
    !
    !       Set the default values for the physics parameters.
    !
    equations     = none       ! These are parameters that must be
    equationMode = none        ! specified. If not, the program
    flowType     = none        ! exits.
    turbModel    = none

    cpModel = cpConstant       ! Constant cp.

    turbProd = strain          ! Strain is used in the production
    ! term of transport turbulence models.

    wallFunctions = .false.    ! No wall functions used.

    Mach     = -one            ! Both parameters must be specified
    Reynolds = -one            ! for external flows. The -1. serves
    ! as a check later on.

    MachCoef = -one            ! If not specified MachCoef will
    ! be set to Mach.

    velDirFreestream(1) = one  ! Free stream velocity
    velDirFreestream(2) = zero ! is specified in the
    velDirFreestream(3) = zero ! x-axis direction.

    liftDirSpecified = .false. ! Lift direction not specified.

    ReynoldsLength = one
    TinfDim = 288.15_realType
    gammaConstant  = 1.4_realType
    RGasDim        = 287.87_realType

    prandtl     = 0.72_realType
    prandtlTurb = 0.90_realType
    pklim       = 20.0_realType
    wallOffset  = zero

    SSuthDim    = 110.55_realType
    muSuthDim   = 1.716e-5_realType
    TSuthDim    = 273.15_realType

    rvfN = 1                             ! Version 1 of the v2f
    ! model is used.
    rvfB = .true.                        ! An upper bound is used
    ! in the v2f scales.
    eddyVisInfRatio = -one               ! Default value depends on
    ! the turbulence model.
    turbIntensityInf = 0.001_realType

    surfaceRef = one
    lengthRef  = one

    pointRef(1) = zero
    pointRef(2) = zero
    pointRef(3) = zero

    momentAxis(1,1) = zero
    momentAxis(1,2) = one
    momentAxis(2,1) = zero
    momentAxis(2,2) = zero
    momentAxis(3,1) = zero
    momentAxis(3,2) = zero

    !
    !       Set the default values for the time spectral parameters.
    !
    nTimeIntervalsSpectral   = -1   ! Serves as a check later on.

    nUnsteadySolSpectral = -1       ! Serves as a check later on.

    writeUnsteadyVolSpectral  = .false. ! No writing of the files
    writeUnsteadySurfSpectral = .false. ! for postprocessing.

    writeUnsteadyRestartSpectral = .false. ! No writing of an unsteady
    ! mode restart file.

    dtUnsteadyRestartSpectral    = -one    ! Is checked later on.
    !
    !       Set the default values for the unsteady parameters.
    !
    timeAccuracy = secondOrder  ! Second order time accuracy.

    nTimeStepsCoarse = -1       ! Serves as a check later on.
    nTimeStepsFine   = -1       ! Serves as a check later on.

    deltaT = -one               ! Serves as a check later on.

    useALE = .True.             ! Use the ALE scheme by default.

    updateWallDistanceUnsteady = .true.  ! This default value is
    ! overruled for models that
    ! are wall distance free.
    !
    !       The reference state variables. Set them to -1, such that they
    !       can be checked later on.
    !
    pRef   = -one
    rhoRef = -one
    TRef   = -one
    !
    !       The conversion factor of the grid units to meters. Default 1.
    !
    LRef           = one
    LRefSpecified = .false.
    !
    !       Initialization of some unsteady restart parameters. These will
    !       be overwritten when an actual unsteady restart is performed.
    !
    nOldSolAvail        = 1
    nTimeStepsRestart   = 0
    timeUnsteadyRestart = zero
    !
    !       Variables needed for the writing of grid and solution files.
    !
    timeSpectralGridsNotWritten = .true.

    ! Additional Paramters Requiring Defaults
    printIterations = .True.
    routineFailed = .False.
    fatalFail     = .False.
    lumpedDiss    = .False.
    approxSA      = .False.
    useApproxWallDistance = .False.
    cflLimit = 3.0
    adjointPETScVarsAllocated = .False.
    usematrixfreedrdw = .False.
    sepSensorOffset = zero
    sepSensorSharpness = 10_realType
  end subroutine setDefaultValues

  subroutine initializeIsoSurfaceVariables(values, nValues)
    !
    !       isoVariables extracts from the given string the extra
    !       iso surface variables to be written to the solution file.
    !
    use constants
    use extraOutput, only : isoValues, isoSurfaceNames, nIsoSurface
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), dimension(nValues), intent(in) :: values
    integer(kind=intType), intent(in) :: nValues

    ! Basically just copy into module
    if (allocated(isoValues)) then
       deallocate(isoValues)
    end if

    if (allocated(isoSurfaceNames)) then
       deallocate(isoSurfaceNames)
    end if

    nIsoSurface = nValues
    allocate(isoValues(nIsoSurface))
    allocate(isoSurfaceNames(nIsoSurface))

    isoValues = values

  end subroutine initializeIsoSurfaceVariables

  subroutine setIsoSurfaceVariable(variable, iVar)

    ! Set variable to iVar. initializeIsoSurfaceVariables MUST be called
    ! first with the desired number of values to set.

    use constants
    use cgnsNames
    use extraOutput
    use communication, only : myID
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(in):: variable
    integer(kind=intType) :: iVar

    select case (variable)
    case("rho")
       isoSurfaceNames(iVar) = cgnsDensity
    case("vx")
       isoSurfaceNames(iVar) = cgnsVelX
    case("vy")
       isoSurfaceNames(iVar) = cgnsVelY
    case("vz")
       isoSurfaceNames(iVar) = cgnsVelZ
    case("P")
       isoSurfaceNames(iVar) = cgnsPressure
    case ("mx")
       isoSurfaceNames(iVar) = cgnsMomX
    case ("my")
       isoSurfaceNames(iVar) = cgnsMomY
    case ("mz")
       isoSurfaceNames(iVar) = cgnsMomZ
    case ("rvx")
       isoSurfaceNames(iVar) = cgnsRelVelX
    case ("rvy")
       isoSurfaceNames(iVar) = cgnsRelVelY
    case ("rvz")
       isoSurfaceNames(iVar) = cgnsRelVelZ
    case ("rhoe")
       isoSurfaceNames(iVar) = cgnsEnergy
    case ("temp")
       isoSurfaceNames(iVar) = cgnsTemp
    case ("vort")
       isoSurfaceNames(iVar) = cgnsVortMagn
    case ("vortx")
       isoSurfaceNames(iVar) = cgnsVortX
    case ("vorty")
       isoSurfaceNames(iVar) = cgnsVortY
    case ("vortz")
       isoSurfaceNames(iVar) = cgnsVortZ
    case ("cp")
       isoSurfaceNames(iVar) = cgnsCp
    case ("mach")
       isoSurfaceNames(iVar) = cgnsMach
    case ("rmach")
       isoSurfaceNames(iVar) = cgnsRelMach
    case ("macht")
       isoSurfaceNames(iVar) = cgnsMachTurb
    case ("ptloss")
       isoSurfaceNames(iVar) = cgnsPTotLoss
    case ("eddy")
       isoSurfaceNames(iVar) = cgnsEddy
    case ("eddyratio")
       isoSurfaceNames(iVar) = cgnsEddyRatio
    case ("dist")
       isoSurfaceNames(iVar) = cgnsWallDist
    case ("resrho")
       isoSurfaceNames(iVar) = cgnsResRho
    case("shock")
       isoSurfaceNames(iVar) = cgnsShock
    case("filteredShock")
       isoSurfaceNames(iVar) = cgnsFilteredShock
    case default

       if(myID == 0) Then
          print *,'Error: ', variable, 'cannot be used as an isoSurface'
       end if
       call EChk(-99, __FILE__, __LINE__)
    end select
  end subroutine setIsoSurfaceVariable

end module inputParamRoutines
