module surfaceIntegrations

contains

  subroutine getCostFunctions(globalVals, funcValues)

    use constants
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : pRef, rhoRef, tRef, LRef, gammaInf, pInf, uRef, uInf
    use inputPhysics, only : liftDirection, dragDirection, surfaceRef, machCoef, lengthRef, alpha, beta, liftIndex
    use inputTSStabDeriv, only : TSstability
    use utils, only : computeTSDerivatives
    use flowUtils, only : getDirVector
    implicit none

    ! Input/Output
    real(kind=realType), intent(in), dimension(:, :) :: globalVals
    real(Kind=realType), intent(out), dimension(:) :: funcValues

    ! Working
    real(kind=realType) :: fact, factMoment, ovrNTS
    real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, forceP, forceV, forceM, &
         moment, cForce, cForceP, cForceV, cForceM, cMoment
    real(kind=realType), dimension(3) :: VcoordRef, VFreestreamRef
    real(kind=realType) ::  mAvgPtot, mAvgTtot, mAvgRho, mAvgPs, mFlow, mAvgMn, mAvga, &
                            mAvgVx, mAvgVy, mAvgVz, gArea

    real(kind=realType) ::  vdotn, mag, u, v, w
    integer(kind=intType) :: sps
    real(kind=realType), dimension(8):: dcdq, dcdqdot
    real(kind=realType), dimension(8):: dcdalpha,dcdalphadot
    real(kind=realType), dimension(8):: Coef0

    ! Factor used for time-averaged quantities.
    ovrNTS = one/nTimeIntervalsSpectral

    ! Sum pressure and viscous contributions
    force = globalvals(iFp:iFp+2, :) + globalvals(iFv:iFv+2, :) + globalvals(iFlowFm:iFlowFm+2, :)
    forceP = globalvals(iFp:iFp+2, :)
    forceV = globalvals(iFv:iFv+2, :)
    forceM = globalvals(iFlowFm:iFlowFm+2, :)

    Moment = globalvals(iMp:iMp+2, :) + globalvals(iMv:iMv+2, :) + globalvals(iFlowMm:iFlowMm+2, :)

    fact = two/(gammaInf*MachCoef*MachCoef &
         *surfaceRef*LRef*LRef*pRef)
    cForce = fact*force
    cForceP = fact*forceP
    cForceV = fact*forceV
    cForceM = fact*forceM

    ! Moment factor has an extra lengthRef
    fact = fact/(lengthRef*LRef)
    cMoment = fact*Moment

    ! Zero values since we are summing.
    funcValues = zero

    ! Here we finally assign the final function values
    !$AD II-LOOP
    do sps=1, nTimeIntervalsSpectral
       funcValues(costFuncForceX) = funcValues(costFuncForceX) + ovrNTS*force(1, sps)
       funcValues(costFuncForceY) = funcValues(costFuncForceY) + ovrNTS*force(2, sps)
       funcValues(costFuncForceZ) = funcValues(costFuncForceZ) + ovrNTS*force(3, sps)

       funcValues(costFuncForceXPressure) = funcValues(costFuncForceXPressure) + ovrNTS*forceP(1, sps)
       funcValues(costFuncForceYPressure) = funcValues(costFuncForceYPressure) + ovrNTS*forceP(2, sps)
       funcValues(costFuncForceZPressure) = funcValues(costFuncForceZPressure) + ovrNTS*forceP(3, sps)

       funcValues(costFuncForceXViscous) = funcValues(costFuncForceXViscous) + ovrNTS*forceV(1, sps)
       funcValues(costFuncForceYViscous) = funcValues(costFuncForceYViscous) + ovrNTS*forceV(2, sps)
       funcValues(costFuncForceZViscous) = funcValues(costFuncForceZViscous) + ovrNTS*forceV(3, sps)

       funcValues(costFuncForceXMomentum) = funcValues(costFuncForceXMomentum) + ovrNTS*forceM(1, sps)
       funcValues(costFuncForceYMomentum) = funcValues(costFuncForceYMomentum) + ovrNTS*forceM(2, sps)
       funcValues(costFuncForceZMomentum) = funcValues(costFuncForceZMomentum) + ovrNTS*forceM(3, sps)

       ! ------------

       funcValues(costFuncForceXCoef) = funcValues(costFuncForceXCoef) + ovrNTS*cForce(1, sps)
       funcValues(costFuncForceYCoef) = funcValues(costFuncForceYCoef) + ovrNTS*cForce(2, sps)
       funcValues(costFuncForceZCoef) = funcValues(costFuncForceZCoef) + ovrNTS*cForce(3, sps)

       funcValues(costFuncForceXCoefPressure) = funcValues(costFuncForceXCoefPressure) + ovrNTS*cForceP(1, sps)
       funcValues(costFuncForceYCoefPressure) = funcValues(costFuncForceYCoefPressure) + ovrNTS*cForceP(2, sps)
       funcValues(costFuncForceZCoefPressure) = funcValues(costFuncForceZCoefPressure) + ovrNTS*cForceP(3, sps)

       funcValues(costFuncForceXCoefViscous) = funcValues(costFuncForceXCoefViscous) + ovrNTS*cForceV(1, sps)
       funcValues(costFuncForceYCoefViscous) = funcValues(costFuncForceYCoefViscous) + ovrNTS*cForceV(2, sps)
       funcValues(costFuncForceZCoefViscous) = funcValues(costFuncForceZCoefViscous) + ovrNTS*cForceV(3, sps)

       funcValues(costFuncForceXCoefMomentum) = funcValues(costFuncForceXCoefMomentum) + ovrNTS*cForceM(1, sps)
       funcValues(costFuncForceYCoefMomentum) = funcValues(costFuncForceYCoefMomentum) + ovrNTS*cForceM(2, sps)
       funcValues(costFuncForceZCoefMomentum) = funcValues(costFuncForceZCoefMomentum) + ovrNTS*cForceM(3, sps)

       ! ------------

       funcValues(costFuncMomX) = funcValues(costFuncMomX) + ovrNTS*moment(1, sps)
       funcValues(costFuncMomY) = funcValues(costFuncMomY) + ovrNTS*moment(2, sps)
       funcValues(costFuncMomZ) = funcValues(costFuncMomZ) + ovrNTS*moment(3, sps)

       funcValues(costFuncMomXCoef) = funcValues(costFuncMomXCoef) + ovrNTS*cMoment(1, sps)
       funcValues(costFuncMomYCoef) = funcValues(costFuncMomYCoef) + ovrNTS*cMoment(2, sps)
       funcValues(costFuncMomZCoef) = funcValues(costFuncMomZCoef) + ovrNTS*cMoment(3, sps)

       funcValues(costFuncSepSensor) = funcValues(costFuncSepSensor) + ovrNTS*globalVals(iSepSensor, sps)
       funcValues(costFuncCavitation) = funcValues(costFuncCavitation) + ovrNTS*globalVals(iCavitation, sps)
       funcValues(costFuncAxisMoment) = funcValues(costFuncAxisMoment) + ovrNTS*globalVals(iAxisMoment, sps)
       funcValues(costFuncSepSensorAvgX) = funcValues(costFuncSepSensorAvgX) + ovrNTS*globalVals(iSepAvg  , sps)
       funcValues(costFuncSepSensorAvgY) = funcValues(costFuncSepSensorAvgY) + ovrNTS*globalVals(iSepAvg+1, sps)
       funcValues(costFuncSepSensorAvgZ) = funcValues(costFuncSepSensorAvgZ) + ovrNTS*globalVals(iSepAvg+2, sps)
       funcValues(costFuncArea)    = funcValues(costFuncArea) + ovrNTS*globalVals(iArea, sps)
       funcValues(costFuncFlowPower) = funcValues(costFuncFlowPower) + ovrNTS*globalVals(iPower, sps)

       funcValues(costFuncCpError2) = funcValues(costFuncCpError2) + ovrNTS*globalVals(iCpError2,sps)

       ! Mass flow like objective
       mFlow = globalVals(iMassFlow, sps)
       if (mFlow /= zero) then
          mAvgPtot = globalVals(iMassPtot, sps)/mFlow
          mAvgTtot = globalVals(iMassTtot, sps)/mFlow
          mAvgrho   = globalVals(iMassRho, sps)/mFlow
          mAvgPs   = globalVals(iMassPs, sps)/mFlow
          mAvgMn   = globalVals(iMassMn, sps)/mFlow
          mAvga   = globalVals(iMassa, sps)/mFlow

          mAvgVx   = globalVals(iMassVx, sps)/mFlow
          mAvgVy   = globalVals(iMassVy, sps)/mFlow
          mAvgVz   = globalVals(iMassVz, sps)/mFlow

          mag = sqrt(globalVals(iMassnx, sps)**2 + &
                    globalVals(iMassny, sps)**2 + &
                    globalVals(iMassnz, sps)**2)

       else
          mAvgPtot = zero
          mAvgTtot = zero
          mAvgrho = zero
          mAvgPs = zero
          mAvgMn = zero
          mAvga = zero
          mAvgVx   = zero
          mAvgVy   = zero
          mAvgVz   = zero

       end if

       ! area averaged objectives
       gArea = globalVals(iArea, sps)
       if (gArea /= zero) then
          ! area averaged pressure
          funcValues(costFuncAAvgPTot) = funcValues(costFuncAAvgPTot) + ovrNTS*globalVals(iAreaPTot, sps) / gArea
          funcValues(costFuncAAvgPs) = funcValues(costFuncAAvgPs) + ovrNTS*globalVals(iAreaPs, sps) / gArea
       end if

       funcValues(costFuncMdot)      = funcValues(costFuncMdot) + ovrNTS*mFlow
       funcValues(costFuncMavgPtot ) = funcValues(costFuncMavgPtot) + ovrNTS*mAvgPtot
       funcValues(costFuncMavgTtot)  = funcValues(costFuncMavgTtot) + ovrNTS*mAvgTtot
       funcValues(costFuncMavgRho)  = funcValues(costFuncMavgRho) + ovrNTS*mAvgRho
       funcValues(costFuncMavgPs)    = funcValues(costFuncMAvgPs) + ovrNTS*mAvgPs
       funcValues(costFuncMavgMn)    = funcValues(costFuncMAvgMn) + ovrNTS*mAvgMn
       funcValues(costFuncMavga)    = funcValues(costFuncMAvga) + ovrNTS*mAvga
       funcValues(costfuncmavgvx)    = funcValues(costfuncmavgvx) + ovrNTS*mAvgVx
       funcValues(costfuncmavgvy)    = funcValues(costfuncmavgvy) + ovrNTS*mAvgVy
       funcValues(costfuncmavgvz)    = funcValues(costfuncmavgvz) + ovrNTS*mAvgVz
       ! Bending moment calc - also broken.
       ! call computeRootBendingMoment(cForce, cMoment, liftIndex, bendingMoment)
       ! funcValues(costFuncBendingCoef) = funcValues(costFuncBendingCoef) + ovrNTS*bendingMoment

    end do

    ! Lift and Drag (coefficients): Dot product with the lift/drag direction.
    funcValues(costFuncLift) = &
         funcValues(costFuncForceX)*liftDirection(1) + &
         funcValues(costFuncForceY)*liftDirection(2) + &
         funcValues(costFuncForceZ)*liftDirection(3)

    funcValues(costFuncLiftPressure) = &
         funcValues(costFuncForceXPressure)*liftDirection(1) + &
         funcValues(costFuncForceYPressure)*liftDirection(2) + &
         funcValues(costFuncForceZPressure)*liftDirection(3)

    funcValues(costFuncLiftViscous) = &
         funcValues(costFuncForceXViscous)*liftDirection(1) + &
         funcValues(costFuncForceYViscous)*liftDirection(2) + &
         funcValues(costFuncForceZViscous)*liftDirection(3)

    funcValues(costFuncLiftMomentum) = &
         funcValues(costFuncForceXMomentum)*liftDirection(1) + &
         funcValues(costFuncForceYMomentum)*liftDirection(2) + &
         funcValues(costFuncForceZMomentum)*liftDirection(3)

    !-----

    funcValues(costFuncDrag) = &
         funcValues(costFuncForceX)*dragDirection(1) + &
         funcValues(costFuncForceY)*dragDirection(2) + &
         funcValues(costFuncForceZ)*dragDirection(3)

    funcValues(costFuncDragPressure) = &
         funcValues(costFuncForceXPressure)*dragDirection(1) + &
         funcValues(costFuncForceYPressure)*dragDirection(2) + &
         funcValues(costFuncForceZPressure)*dragDirection(3)

    funcValues(costFuncDragViscous) = &
         funcValues(costFuncForceXViscous)*dragDirection(1) + &
         funcValues(costFuncForceYViscous)*dragDirection(2) + &
         funcValues(costFuncForceZViscous)*dragDirection(3)

    funcValues(costFuncDragMomentum) = &
         funcValues(costFuncForceXMomentum)*dragDirection(1) + &
         funcValues(costFuncForceYMomentum)*dragDirection(2) + &
         funcValues(costFuncForceZMomentum)*dragDirection(3)

    !-----

    funcValues(costFuncLiftCoef) = &
         funcValues(costFuncForceXCoef)*liftDirection(1) + &
         funcValues(costFuncForceYCoef)*liftDirection(2) + &
         funcValues(costFuncForceZCoef)*liftDirection(3)

    funcValues(costFuncLiftCoefPressure) = &
         funcValues(costFuncForceXCoefPressure)*liftDirection(1) + &
         funcValues(costFuncForceYCoefPressure)*liftDirection(2) + &
         funcValues(costFuncForceZCoefPressure)*liftDirection(3)

    funcValues(costFuncLiftCoefViscous) = &
         funcValues(costFuncForceXCoefViscous)*liftDirection(1) + &
         funcValues(costFuncForceYCoefViscous)*liftDirection(2) + &
         funcValues(costFuncForceZCoefViscous)*liftDirection(3)

    funcValues(costFuncLiftCoefMomentum) = &
         funcValues(costFuncForceXCoefMomentum)*liftDirection(1) + &
         funcValues(costFuncForceYCoefMomentum)*liftDirection(2) + &
         funcValues(costFuncForceZCoefMomentum)*liftDirection(3)

    !-----

    funcValues(costFuncDragCoef) = &
         funcValues(costFuncForceXCoef)*dragDirection(1) + &
         funcValues(costFuncForceYCoef)*dragDirection(2) + &
         funcValues(costFuncForceZCoef)*dragDirection(3)

    funcValues(costFuncDragCoefPressure) = &
         funcValues(costFuncForceXCoefPressure)*dragDirection(1) + &
         funcValues(costFuncForceYCoefPressure)*dragDirection(2) + &
         funcValues(costFuncForceZCoefPressure)*dragDirection(3)

    funcValues(costFuncDragCoefViscous) = &
         funcValues(costFuncForceXCoefViscous)*dragDirection(1) + &
         funcValues(costFuncForceYCoefViscous)*dragDirection(2) + &
         funcValues(costFuncForceZCoefViscous)*dragDirection(3)

    funcValues(costFuncDragCoefMomentum) = &
         funcValues(costFuncForceXCoefMomentum)*dragDirection(1) + &
         funcValues(costFuncForceYCoefMomentum)*dragDirection(2) + &
         funcValues(costFuncForceZCoefMomentum)*dragDirection(3)


    ! -------------------- Time Spectral Objectives ------------------

    if (TSSTability) then
       print *,'Error: TSStabilityDerivatives are *BROKEN*. They need to be '&
            &'completely verifed from scratch'
       stop

       call computeTSDerivatives(force, moment, coef0, dcdalpha, &
            dcdalphadot, dcdq, dcdqdot)

       funcValues( costFuncCl0  )       = coef0(1)
       funcValues( costFuncCd0 )        = coef0(2)
       funcValues( costFuncCFy0 )       = coef0(4)
       funcValues( costFuncCm0 )        = coef0(8)

       funcValues( costFuncClAlpha)     = dcdalpha(1)
       funcValues( costFuncCdAlpha)     = dcdalpha(2)
       funcValues( costFuncCFyAlpha)    = dcdalpha(4)
       funcValues( costFuncCmzAlpha)    = dcdalpha(8)

       funcValues( costFuncClAlphaDot)     = dcdalphadot(1)
       funcValues( costFuncCdAlphaDot)     = dcdalphadot(2)
       funcValues( costFuncCFyAlphaDot)    = dcdalphadot(4)
       funcValues( costFuncCmzAlphaDot)    = dcdalphadot(8)

       funcValues( costFuncClq)         = dcdq(1)
       funcValues( costFuncCdq)         = dcdq(2)
       funcValues( costFuncCfyq)        = dcdq(4)
       funcValues( costFuncCmzq)        = dcdq(8)

       funcValues( costFuncClqDot)         = dcdqdot(1)
       funcValues( costFuncCdqDot)         = dcdqdot(2)
       funcValues( costFuncCfyqDot)        = dcdqdot(4)
       funcValues( costFuncCmzqDot)        = dcdqdot(8)
    end if

  end subroutine getCostFunctions

  subroutine wallIntegrationFace(localValues, mm)
    !
    !       wallIntegrations computes the contribution of the block
    !       given by the pointers in blockPointers to the force and
    !       moment of the geometry. A distinction is made
    !       between the inviscid and viscous parts. In case the maximum
    !       yplus value must be monitored (only possible for rans), this
    !       value is also computed. The separation sensor and the cavita-
    !       tion sensor is also computed
    !       here.
    !
    use constants
    use communication
    use blockPointers
    use flowVarRefState
    use inputCostFunctions
    use inputPhysics, only : MachCoef, pointRef, velDirFreeStream, equations, momentAxis, cavitationnumber
    use BCPointers
    implicit none

    ! Input/output variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType) :: mm

    ! Local variables.
    real(kind=realType), dimension(3)  :: Fp, Fv, Mp, Mv
    real(kind=realType)  :: yplusMax, sepSensor, sepSensorAvg(3), Cavitation
    integer(kind=intType) :: i, j, ii, blk

    real(kind=realType) :: pm1, fx, fy, fz, fn
    real(kind=realType) :: xc, yc, zc, qf(3), r(3), n(3), L
    real(kind=realType) :: fact, rho, mul, yplus, dwall
    real(kind=realType) :: V(3), sensor, sensor1, Cp, tmp, plocal
    real(kind=realType) :: tauXx, tauYy, tauZz
    real(kind=realType) :: tauXy, tauXz, tauYz

    real(kind=realType), dimension(3) :: refPoint
    real(kind=realType), dimension(3,2) :: axisPoints
    real(kind=realType) :: mx, my, mz, cellArea, m0x, m0y, m0z, Mvaxis, Mpaxis
    real(kind=realType) :: CpError, CpError2

    select case (BCFaceID(mm))
    case (iMin, jMin, kMin)
       fact = -one
    case (iMax, jMax, kMax)
       fact = one
    end select

    ! Determine the reference point for the moment computation in
    ! meters.

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    ! Determine the points defining the axis about which to compute a moment
    axisPoints(1,1) = LRef*momentAxis(1,1)
    axisPoints(1,2) = LRef*momentAxis(1,2)
    axisPoints(2,1) = LRef*momentAxis(2,1)
    axisPoints(2,2) = LRef*momentAxis(2,2)
    axisPoints(3,1) = LRef*momentAxis(3,1)
    axisPoints(3,2) = LRef*momentAxis(3,2)

    ! Initialize the force and moment coefficients to 0 as well as
    ! yplusMax.

    Fp = zero; Fv = zero;
    Mp = zero; Mv = zero;
    yplusMax = zero
    sepSensor = zero
    Cavitation = zero
    sepSensorAvg = zero
    Mpaxis = zero; Mvaxis = zero;
    CpError2 = zero;

    !
    !         Integrate the inviscid contribution over the solid walls,
    !         either inviscid or viscous. The integration is done with
    !         cp. For closed contours this is equal to the integration
    !         of p; for open contours this is not the case anymore.
    !         Question is whether a force for an open contour is
    !         meaningful anyway.
    !


    ! Loop over the quadrilateral faces of the subface. Note that
    ! the nodal range of BCData must be used and not the cell
    ! range, because the latter may include the halo's in i and
    ! j-direction. The offset +1 is there, because inBeg and jnBeg
    ! refer to nodal ranges and not to cell ranges. The loop
    ! (without the AD stuff) would look like:
    !
    ! do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
    !    do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

    !$AD II-LOOP
    do ii=0,(BCData(mm)%jnEnd - bcData(mm)%jnBeg)*(bcData(mm)%inEnd - bcData(mm)%inBeg) -1
       i = mod(ii, (bcData(mm)%inEnd-bcData(mm)%inBeg)) + bcData(mm)%inBeg + 1
       j = ii/(bcData(mm)%inEnd-bcData(mm)%inBeg) + bcData(mm)%jnBeg + 1

       ! Compute the average pressure minus 1 and the coordinates
       ! of the centroid of the face relative from from the
       ! moment reference point. Due to the usage of pointers for
       ! the coordinates, whose original array starts at 0, an
       ! offset of 1 must be used. The pressure is multipled by
       ! fact to account for the possibility of an inward or
       ! outward pointing normal.

       pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)*pRef

       tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
       cp = (half*(pp2(i,j) + pp1(i,j)) - pInf)*tmp
       CpError = (cp - BCData(mm)%CpTarget(i,j))
       CPError2 = CpError2 + CpError*CpError

       xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
            +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
       yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
            +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
       zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
            +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

       ! Compute the force components.
       blk = max(BCData(mm)%iblank(i,j), 0)
       fx = pm1*ssi(i,j,1)
       fy = pm1*ssi(i,j,2)
       fz = pm1*ssi(i,j,3)

       ! Update the inviscid force and moment coefficients. Iblank as we sum
       Fp(1) = Fp(1) + fx*blk
       Fp(2) = Fp(2) + fy*blk
       Fp(3) = Fp(3) + fz*blk

       mx = yc*fz - zc*fy
       my = zc*fx - xc*fz
       mz = xc*fy - yc*fx

       Mp(1) = Mp(1) + mx*blk
       Mp(2) = Mp(2) + my*blk
       Mp(3) = Mp(3) + mz*blk

       ! Compute the r and n vectores for the moment around an
       ! axis computation where r is the distance from the
       ! force to the first point on the axis and n is a unit
       ! normal in the direction of the axis
       r(1) = fourth*(xx(i,j,   1) + xx(i+1,j,   1)&
            + xx(i,j+1,1) + xx(i+1,j+1,1)) - axisPoints(1,1)
       r(2) = fourth*(xx(i,j,   2) + xx(i+1,j,   2)&
            + xx(i,j+1,2) + xx(i+1,j+1,2)) - axisPoints(2,1)
       r(3) = fourth*(xx(i,j,   3) + xx(i+1,j,   3)&
            + xx(i,j+1,3) + xx(i+1,j+1,3)) - axisPoints(3,1)

       L = sqrt((axisPoints(1,2) - axisPoints(1,1)) **2 &
            + (axisPoints(2,2) - axisPoints(2,1)) **2 &
            + (axisPoints(3,2) - axisPoints(3,1)) **2)

       n(1) = (axisPoints(1,2) - axisPoints(1,1)) / L
       n(2) = (axisPoints(2,2) - axisPoints(2,1)) / L
       n(3) = (axisPoints(3,2) - axisPoints(3,1)) / L

       ! Compute the moment of the force about the first point
       ! used to define the axis, and the project that axis in
       ! the n direction
       m0x = r(2)*fz - r(3)*fy
       m0y = r(3)*fx - r(1)*fz
       m0z = r(1)*fy - r(2)*fx
       Mpaxis = Mpaxis +(m0x*n(1) + m0y*n(2) + m0z*n(3))*blk

       ! Save the face-based forces and area
       bcData(mm)%Fp(i, j, 1) = fx
       bcData(mm)%Fp(i, j, 2) = fy
       bcData(mm)%Fp(i, j, 3) = fz
       cellArea = sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2)

       bcData(mm)%area(i, j) = cellArea

       ! Get normalized surface velocity:
       v(1) = ww2(i, j, ivx)
       v(2) = ww2(i, j, ivy)
       v(3) = ww2(i, j, ivz)
       v = v / (sqrt(v(1)**2 + v(2)**2 + v(3)**2) + 1e-16)

       ! Dot product with free stream
       sensor = -(v(1)*velDirFreeStream(1) + v(2)*velDirFreeStream(2) + &
            v(3)*velDirFreeStream(3))

       !Now run through a smooth heaviside function:
       sensor = one/(one + exp(-2*sepSensorSharpness*(sensor-sepSensorOffset)))

       ! And integrate over the area of this cell and save, blanking as we go.
       sensor = sensor * cellArea * blk
       sepSensor = sepSensor + sensor

       ! Also accumulate into the sepSensorAvg
       xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
            +         xx(i,j+1,1) + xx(i+1,j+1,1))
       yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
            +         xx(i,j+1,2) + xx(i+1,j+1,2))
       zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
            +         xx(i,j+1,3) + xx(i+1,j+1,3))

       sepSensorAvg(1) = sepSensorAvg(1)  + sensor * xc
       sepSensorAvg(2) = sepSensorAvg(2)  + sensor * yc
       sepSensorAvg(3) = sepSensorAvg(3)  + sensor * zc

       if (computeCavitation) then
          plocal = pp2(i,j)
          tmp = two/(gammaInf*MachCoef*MachCoef)
          Cp = tmp*(plocal-pinf)
          Sensor1 = -Cp - cavitationnumber
          Sensor1 = one/(one+exp(-2*10*Sensor1))
          Sensor1 = Sensor1 * cellArea * blk
          Cavitation = Cavitation + Sensor1
       end if
    enddo

    !
    ! Integration of the viscous forces.
    ! Only for viscous boundaries.
    !
    visForce: if( BCType(mm) == NSWallAdiabatic .or. &
         BCType(mm) == NSWallIsoThermal) then

       ! Initialize dwall for the laminar case and set the pointer
       ! for the unit normals.

       dwall = zero

       ! Loop over the quadrilateral faces of the subface and
       ! compute the viscous contribution to the force and
       ! moment and update the maximum value of y+.

       !$AD II-LOOP
       do ii=0,(BCData(mm)%jnEnd - bcData(mm)%jnBeg)*(bcData(mm)%inEnd - bcData(mm)%inBeg) -1
          i = mod(ii, (bcData(mm)%inEnd-bcData(mm)%inBeg)) + bcData(mm)%inBeg + 1
          j = ii/(bcData(mm)%inEnd-bcData(mm)%inBeg) + bcData(mm)%jnBeg + 1

          ! Store the viscous stress tensor a bit easier.
          blk = max(BCData(mm)%iblank(i,j), 0)

          tauXx = viscSubface(mm)%tau(i,j,1)
          tauYy = viscSubface(mm)%tau(i,j,2)
          tauZz = viscSubface(mm)%tau(i,j,3)
          tauXy = viscSubface(mm)%tau(i,j,4)
          tauXz = viscSubface(mm)%tau(i,j,5)
          tauYz = viscSubface(mm)%tau(i,j,6)

          ! Compute the viscous force on the face. A minus sign
          ! is now present, due to the definition of this force.

          fx = -fact*(tauXx*ssi(i,j,1) + tauXy*ssi(i,j,2) &
               +        tauXz*ssi(i,j,3))*pRef
          fy = -fact*(tauXy*ssi(i,j,1) + tauYy*ssi(i,j,2) &
               +        tauYz*ssi(i,j,3))*pRef
          fz = -fact*(tauXz*ssi(i,j,1) + tauYz*ssi(i,j,2) &
               +        tauZz*ssi(i,j,3))*pRef

          ! Compute the coordinates of the centroid of the face
          ! relative from the moment reference point. Due to the
          ! usage of pointers for xx and offset of 1 is present,
          ! because x originally starts at 0.

          xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
               +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
          yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
               +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
          zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
               +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

          ! Update the viscous force and moment coefficients, blanking as we go.

          Fv(1) = Fv(1) + fx * blk
          Fv(2) = Fv(2) + fy * blk
          Fv(3) = Fv(3) + fz * blk

          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx

          Mv(1) = Mv(1) + mx * blk
          Mv(2) = Mv(2) + my * blk
          Mv(3) = Mv(3) + mz * blk

          ! Compute the r and n vectors for the moment around an
          ! axis computation where r is the distance from the
          ! force to the first point on the axis and n is a unit
          ! normal in the direction of the axis
          r(1) = fourth*(xx(i,j,   1) + xx(i+1,j,   1)&
               + xx(i,j+1,1) + xx(i+1,j+1,1)) - axisPoints(1,1)
          r(2) = fourth*(xx(i,j,   2) + xx(i+1,j,   2)&
               + xx(i,j+1,2) + xx(i+1,j+1,2)) - axisPoints(2,1)
          r(3) = fourth*(xx(i,j,   3) + xx(i+1,j,   3)&
               + xx(i,j+1,3) + xx(i+1,j+1,3)) - axisPoints(3,1)

          L = sqrt((axisPoints(1,2) - axisPoints(1,1)) ** 2 &
               + (axisPoints(2,2) - axisPoints(2,1)) ** 2 &
               + (axisPoints(3,2) - axisPoints(3,1)) ** 2 )

          n(1) = (axisPoints(1,2) - axisPoints(1,1)) / L
          n(2) = (axisPoints(2,2) - axisPoints(2,1)) / L
          n(3) = (axisPoints(3,2) - axisPoints(3,1)) / L

          ! Compute the moment of the force about the first point
          ! used to define the axis, and then project that axis in
          ! the n direction
          m0x = r(2)*fz - r(3)*fy
          m0y = r(3)*fx - r(1)*fz
          m0z = r(1)*fy - r(2)*fx
          Mvaxis = Mvaxis + (m0x*n(1) + m0y*n(2) + m0z*n(3))*blk

          ! Save the face based forces for the slice operations
          bcData(mm)%Fv(i, j, 1) = fx
          bcData(mm)%Fv(i, j, 2) = fy
          bcData(mm)%Fv(i, j, 3) = fz

          ! Compute the tangential component of the stress tensor,
          ! which is needed to monitor y+. The result is stored
          ! in fx, fy, fz, although it is not really a force.
          ! As later on only the magnitude of the tangential
          ! component is important, there is no need to take the
          ! sign into account (it should be a minus sign).

          fx = tauXx*BCData(mm)%norm(i,j,1) + tauXy*BCData(mm)%norm(i,j,2) &
               + tauXz*BCData(mm)%norm(i,j,3)
          fy = tauXy*BCData(mm)%norm(i,j,1) + tauYy*BCData(mm)%norm(i,j,2) &
               + tauYz*BCData(mm)%norm(i,j,3)
          fz = tauXz*BCData(mm)%norm(i,j,1) + tauYz*BCData(mm)%norm(i,j,2) &
               + tauZz*BCData(mm)%norm(i,j,3)

          fn = fx*BCData(mm)%norm(i,j,1) + fy*BCData(mm)%norm(i,j,2) + fz*BCData(mm)%norm(i,j,3)

          fx = fx - fn*BCData(mm)%norm(i,j,1)
          fy = fy - fn*BCData(mm)%norm(i,j,2)
          fz = fz - fn*BCData(mm)%norm(i,j,3)

          ! Compute the local value of y+. Due to the usage
          ! of pointers there is on offset of -1 in dd2Wall..
#ifndef USE_TAPENADE
          if(equations == RANSEquations) then
             dwall = dd2Wall(i-1,j-1)
             rho   = half*(ww2(i,j,irho) + ww1(i,j,irho))
             mul   = half*(rlv2(i,j) + rlv1(i,j))
             yplus = sqrt(rho*sqrt(fx*fx + fy*fy + fz*fz))*dwall/mul

             ! Store this value if this value is larger than the
             ! currently stored value. Blank non-active cells.

             yplusMax = max(yplusMax, yplus*blk)
          end if
#endif
       enddo
    else
       ! If we had no viscous force, set the viscous component to zero
       bcData(mm)%Fv = zero
    end if visForce

    ! Increment the local values array with the values we computed here.
    localValues(iFp:iFp+2) = localValues(iFp:iFp+2) + Fp
    localValues(iFv:iFv+2) = localValues(iFv:iFv+2) + Fv
    localValues(iMp:iMp+2) = localValues(iMp:iMp+2) + Mp
    localValues(iMv:iMv+2) = localValues(iMv:iMv+2) + Mv
    localValues(iSepSensor) = localValues(iSepSensor) + sepSensor
    localValues(iCavitation) = localValues(iCavitation) + cavitation
    localValues(iSepAvg:iSepAvg+2) = localValues(iSepAvg:iSepAvg+2) + sepSensorAvg
    localValues(iAxisMoment) = localValues(iAxisMoment) + Mpaxis + Mvaxis
    localValues(iCpError2) = localValues(iCpError2) + CpError2

#ifndef USE_TAPENADE
    localValues(iyPlus) = max(localValues(iyPlus), yplusMax)
#endif
  end subroutine wallIntegrationFace

  subroutine flowIntegrationFace(isInflow, localValues, mm)

    use constants
    use blockPointers, only : BCType, BCFaceID, BCData, addGridVelocities
    use flowVarRefState, only : pRef, pInf, rhoRef, timeRef, LRef, TRef, RGas, uRef, uInf, rhoInf
    use inputPhysics, only : pointRef, flowType
    use flowUtils, only : computePtot, computeTtot
    use BCPointers, only : ssi, sFace, ww1, ww2, pp1, pp2, xx, gamma1, gamma2
    use utils, only : mynorm2
    implicit none

    ! Input/Output variables
    logical, intent(in) :: isInflow
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), intent(in) :: mm

    ! Local variables
    real(kind=realType) ::  massFlowRate, mass_Ptot, mass_Ttot, mass_Ps, mass_MN, mass_a, mass_rho, &
                            mass_Vx, mass_Vy, mass_Vz, mass_nx, mass_ny, mass_nz
    real(kind=realType) ::  area_Ptot, area_Ps
    real(kind=realType) ::  mReDim
    integer(kind=intType) :: i, j, ii, blk
    real(kind=realType) :: internalFlowFact, inFlowFact, fact, xc, yc, zc, mx, my, mz
    real(kind=realType) :: sF, vmag, vnm, vnmFreeStreamRef, vxm, vym, vzm, Fx, Fy, Fz, u, v, w
    real(kind=realType) :: pm, Ptot, Ttot, rhom, gammam, am
    real(kind=realType) :: area, cellArea, overCellArea
    real(kind=realType), dimension(3) :: Fp, Mp, FMom, MMom, refPoint, sFaceCoordRef
    real(kind=realType) :: MNm, massFlowRateLocal

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    ! Note that these are *opposite* of force integrations. The reason
    ! is that we want positive mass flow into the domain and negative
    ! mass flow out of the domain. Since the low faces have ssi
    ! vectors pointining into the domain, this is correct. The high
    ! end faces need to flip this.
    select case (BCFaceID(mm))
    case (iMin, jMin, kMin)
       fact = one
    case (iMax, jMax, kMax)
       fact = -one
    end select

    ! the sign of momentum forces are flipped for internal flows
    internalFlowFact = one
    if (flowType == internalFlow) then
      internalFlowFact = -one
    end if

    inFlowFact = one
    if (isInflow) then
      inflowFact= -one
    end if

    ! Loop over the quadrilateral faces of the subface. Note that
    ! the nodal range of BCData must be used and not the cell
    ! range, because the latter may include the halo's in i and
    ! j-direction. The offset +1 is there, because inBeg and jnBeg
    ! refer to nodal ranges and not to cell ranges. The loop
    ! (without the AD stuff) would look like:
    !
    ! do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
    !    do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

    mReDim = sqrt(pRef*rhoRef)
    Fp = zero
    Mp = zero
    FMom = zero
    MMom = zero

    massFlowRate = zero
    area = zero
    mass_Ptot = zero
    mass_Ttot = zero
    mass_Ps = zero
    mass_MN = zero
    mass_a = zero
    mass_rho = zero

    mass_Vx = zero
    mass_Vy = zero
    mass_Vz = zero
    mass_nx = zero
    mass_ny = zero
    mass_nz = zero

    area_Ptot = zero
    area_Ps   = zero

    !$AD II-LOOP
    do ii=0,(BCData(mm)%jnEnd - bcData(mm)%jnBeg)*(bcData(mm)%inEnd - bcData(mm)%inBeg) -1
      i = mod(ii, (bcData(mm)%inEnd-bcData(mm)%inBeg)) + bcData(mm)%inBeg + 1
      j = ii/(bcData(mm)%inEnd-bcData(mm)%inBeg) + bcData(mm)%jnBeg + 1

      if( addGridVelocities ) then
        sF = sFace(i,j)
      else
        sF = zero
      end if

      ! Compute the force components.
      blk = max(BCData(mm)%iblank(i,j), 0) ! iBlank forces for overset stuff

      vxm = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
      vym = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
      vzm = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
      rhom = half*(ww1(i,j,irho) + ww2(i,j,irho))
      pm = half*(pp1(i,j)+ pp2(i,j))
      gammam = half*(gamma1(i,j) + gamma2(i,j))

      vnm = vxm*ssi(i,j,1) + vym*ssi(i,j,2) + vzm*ssi(i,j,3)  - sF
      vmag = sqrt((vxm**2 + vym**2 + vzm**2)) - sF
      am = sqrt(gammam*pm/rhom)
      MNm = vmag/am

      cellArea = sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2)
      area = area + cellArea*blk
      overCellArea = 1/cellArea

      call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
      call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)

      massFlowRateLocal = rhom*vnm*blk*fact*mReDim



      massFlowRate = massFlowRate + massFlowRateLocal

      ! re-dimentionalize quantities
      pm = pm*pRef

      mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
      mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
      mass_rho  = mass_rho + rhom * massFlowRateLocal * rhoRef
      mass_a  = mass_a + am * massFlowRateLocal * uRef

      mass_Ps = mass_Ps + pm*massFlowRateLocal
      mass_MN = mass_MN + MNm*massFlowRateLocal

      area_pTot = area_pTot + Ptot * Pref * cellArea * blk
      area_Ps = area_Ps + pm * cellArea * blk

      sFaceCoordRef(1) = sF * ssi(i,j,1)*overCellArea
      sFaceCoordRef(2) = sF * ssi(i,j,2)*overCellArea
      sFaceCoordRef(3) = sF * ssi(i,j,3)*overCellArea

      mass_Vx = mass_Vx + (vxm*uRef - sFaceCoordRef(1)) *massFlowRateLocal
      mass_Vy = mass_Vy + (vym*uRef - sFaceCoordRef(2)) *massFlowRateLocal
      mass_Vz = mass_Vz + (vzm*uRef - sFaceCoordRef(3)) *massFlowRateLocal

      mass_nx = mass_nx + ssi(i,j,1)*overCellArea * massFlowRateLocal
      mass_ny = mass_ny + ssi(i,j,2)*overCellArea * massFlowRateLocal
      mass_nz = mass_nz + ssi(i,j,3)*overCellArea * massFlowRateLocal

      xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
          +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
      yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
          +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
      zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
          +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

      ! Pressure forces. Note that these need a *negative* and to subtract
      ! the reference pressure sign to be consistent with the force
      ! computation on the walls.
      pm = -(pm-pInf*pRef)*fact*blk

      fx = pm*ssi(i,j,1)
      fy = pm*ssi(i,j,2)
      fz = pm*ssi(i,j,3)

      ! Update the pressure force and moment coefficients.
      Fp(1) = Fp(1) + fx
      Fp(2) = Fp(2) + fy
      Fp(3) = Fp(3) + fz

      mx = yc*fz - zc*fy
      my = zc*fx - xc*fz
      mz = xc*fy - yc*fx

      Mp(1) = Mp(1) + mx
      Mp(2) = Mp(2) + my
      Mp(3) = Mp(3) + mz

      ! Momentum forces are a little tricky.  We negate because
      ! have to re-apply fact to massFlowRateLocal to undoo it, because
      ! we need the signed behavior of ssi to get the momentum forces correct.
      ! Also, the sign is flipped between inflow and outflow types

      massFlowRateLocal = massFlowRateLocal*fact/timeRef*blk/cellArea*internalFlowFact*inFlowFact

      fx = massFlowRateLocal * ssi(i,j,1)*vxm
      fy = massFlowRateLocal * ssi(i,j,2)*vym
      fz = massFlowRateLocal * ssi(i,j,3)*vzm

      FMom(1) = FMom(1) + fx
      FMom(2) = FMom(2) + fy
      FMom(3) = FMom(3) + fz

      mx = yc*fz - zc*fy
      my = zc*fx - xc*fz
      mz = xc*fy - yc*fx

      MMom(1) = MMom(1) + mx
      MMom(2) = MMom(2) + my
      MMom(3) = MMom(3) + mz

    enddo

    ! Increment the local values array with what we computed here
    localValues(iMassFlow) = localValues(iMassFlow) + massFlowRate
    localValues(iArea) = localValues(iArea) + area
    localValues(iMassRho) = localValues(iMassRho) + mass_rho
    localValues(iMassa) = localValues(iMassa) + mass_a
    localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
    localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
    localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps
    localValues(iMassMN)   = localValues(iMassMN)   + mass_MN
    localValues(iFp:iFp+2)   = localValues(iFp:iFp+2) + Fp
    localValues(iFlowFm:iFlowFm+2)   = localValues(iFlowFm:iFlowFm+2) + FMom
    localValues(iFlowMp:iFlowMp+2)   = localValues(iFlowMp:iFlowMp+2) + Mp
    localValues(iFlowMm:iFlowMm+2)   = localValues(iFlowMm:iFlowMm+2) + MMom

    localValues(iAreaPTot) = localValues(iAreaPTot) + area_pTot
    localValues(iAreaPs) = localValues(iAreaPs) + area_Ps

    localValues(iMassVx)   = localValues(iMassVx)   + mass_Vx
    localValues(iMassVy)   = localValues(iMassVy)   + mass_Vy
    localValues(iMassVz)   = localValues(iMassVz)   + mass_Vz
    localValues(iMassnx)   = localValues(iMassnx)   + mass_nx
    localValues(iMassny)   = localValues(iMassny)   + mass_ny
    localValues(iMassnz)   = localValues(iMassnz)   + mass_nz

  end subroutine flowIntegrationFace


  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE

  subroutine getSolutionWrap(famLists, funcValues, nCost, nGroups, nFamMax)

    use constants
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    implicit none
    ! Input/output Variables
    integer(kind=intType), dimension(nGroups, nFamMax) :: famLists
    real(kind=realType), dimension(nCost, nGroups), intent(out)  :: funcValues
    integer(kind=intType) :: nGroups, nCost, nFamMax

    ! Local variable

    call getSolution(famLists, funcValues)
  end subroutine getSolutionWrap

  subroutine getSolution(famLists, funcValues, globalValues)
    !--------------------------------------------------------------
    ! Manual Differentiation Warning: Modifying this routine requires
    ! modifying the hand-written forward and reverse routines.
    ! --------------------------------------------------------------

    use constants
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    use communication, only : adflow_comm_world
    use blockPointers, only : nDom
    use utils, only : setPointers, EChk
    use zipperIntegrations, only :integrateZippers
    use userSurfaceIntegrations, only : integrateUserSurfaces
    use actuatorRegion, only : integrateActuatorRegions
    implicit none

    ! Input/Output Variables
    integer(kind=intType), dimension(:, :), intent(in), target :: famLists
    real(kind=realType), dimension(:, :), intent(out) :: funcValues
    real(kind=realType), optional, intent(out), dimension(:, :, :) :: globalValues

    ! Working
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal
    integer(kind=intType) :: nn, sps, ierr, iGroup, nFam
    integer(kind=intType), dimension(:), pointer :: famList
    ! Master loop over the each of the groups we have

    groupLoop: do iGroup=1, size(famLists, 1)

       ! Extract the current family list
       nFam = famLists(iGroup, 1)
       famList => famLists(iGroup, 2:2+nFam-1)
       funcValues = zero
       localVal = zero

       do sps=1, nTimeIntervalsSpectral
          ! Integrate the normal block surfaces.
          do nn=1, nDom
             call setPointers(nn, 1, sps)
             call integrateSurfaces(localval(:, sps), famList)
          end do

          ! Integrate any zippers we have
          call integrateZippers(localVal(:, sps), famList, sps)

          ! Integrate any user-supplied surfaces as have as well.
          call integrateUserSurfaces(localVal(:, sps), famList, sps)

          ! Integrate any actuator regions we have
          call integrateActuatorRegions(localVal(:, sps), famList, sps)
       end do

       ! Now we need to reduce all the cost functions
       call mpi_allreduce(localval, globalVal, nLocalValues*nTimeIntervalsSpectral, adflow_real, &
            MPI_SUM, adflow_comm_world, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Call the final routine that will comptue all of our functions of
       ! interest.
       call getCostFunctions(globalVal, funcValues(:, iGroup))

       if (present(globalValues)) then
          globalValues(:, :, iGroup) = globalVal
       end if
    end do groupLoop
  end subroutine getSolution

  subroutine integrateSurfaces(localValues, famList)
    !--------------------------------------------------------------
    ! Manual Differentiation Warning: Modifying this routine requires
    ! modifying the hand-written forward and reverse routines.
    ! --------------------------------------------------------------
    !
    ! This is a shell routine that calls the specific surface
    ! integration routines. Currently we have have the forceAndMoment
    ! routine as well as the flow properties routine. This routine
    ! takes care of setting pointers, while the actual computational
    ! routine just acts on a specific fast pointed to by pointers.

    use constants
    use blockPointers, only : nBocos, BCData, BCType, sk, sj, si, x, rlv, &
         sfacei, sfacej, sfacek, gamma, rev, p, viscSubface
    use utils, only : setBCPointers, isWallType
    use sorting, only : famInList
    ! Tapenade needs to see these modules that the callees use.
    use BCPointers
    use flowVarRefState
    use inputPhysics

    implicit none

    ! Input/output Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList

    ! Working variables
    integer(kind=intType) :: mm

    ! Loop over all possible boundary conditions
    bocos: do mm=1, nBocos

       ! Determine if this boundary condition is to be incldued in the
       ! currently active group
       famInclude: if (famInList(BCData(mm)%famID, famList)) then

          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible.
          call setBCPointers(mm, .True.)

          ! no post gathered integrations currently
          isWall: if( isWallType(BCType(mm)) ) then
             call wallIntegrationFace(localvalues, mm)
          end if isWall

          isInflowOutflow: if (BCType(mm) == SubsonicInflow .or. &
               BCType(mm) == SupersonicInflow) then
             call flowIntegrationFace(.true., localValues, mm)
          else if (BCType(mm) == SubsonicOutflow .or. &
               BCType(mm) == SupersonicOutflow) then
             call flowIntegrationFace(.false., localValues, mm)
          end if isInflowOutflow

       end if famInclude
    end do bocos

  end subroutine integrateSurfaces

#ifndef USE_COMPLEX
  subroutine integrateSurfaces_d(localValues, localValuesd, famList)
    !------------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    ! Forward mode linearization of integrateSurfaces
    use constants
    use blockPointers, only : nBocos, BCData, BCType
    use utils, only : setBCPointers_d, isWallType
    use sorting, only : famInList
    use surfaceIntegrations_d, only : wallIntegrationFace_d, flowIntegrationFace_d
    implicit none

    ! Input/output Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
    integer(kind=intType), dimension(:), intent(in) :: famList

    ! Working variables
    integer(kind=intType) :: mm

    ! Loop over all possible boundary conditions
    do mm=1, nBocos
       ! Determine if this boundary condition is to be incldued in the
       ! currently active group
       famInclude: if (famInList(BCData(mm)%famID, famList)) then

          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible.
          call setBCPointers_d(mm, .True.)

          ! not post gathered integrations currently
          isWall: if( isWallType(BCType(mm)) ) then
             call wallIntegrationFace_d(localValues, localValuesd, mm)
          end if isWall

          isInflowOutflow: if (BCType(mm) == SubsonicInflow .or. &
               BCType(mm) == SupersonicInflow) then
             call flowIntegrationFace_d(.true., localValues, localValuesd, mm)
            else if (BCType(mm) == SubsonicOutflow .or. &
               BCType(mm) == SupersonicOutflow) then

               call flowIntegrationFace_d(.false., localValues, localValuesd, mm)

          end if isInflowOutflow
       end if famInclude
    end do
  end subroutine integrateSurfaces_d

  subroutine integrateSurfaces_b(localValues, localValuesd, famList)
    !------------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    ! Reverse mode linearization of integrateSurfaces
    use constants
    use blockPointers, only : nBocos, BCData, BCType, bcDatad
    use utils, only : setBCPointers_d, isWallType
    use sorting, only : famInList
    use surfaceIntegrations_b, only : wallIntegrationFace_b, flowIntegrationFace_b
    implicit none

    ! Input/output Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
    integer(kind=intType), dimension(:), intent(in) :: famList
    ! Working variables
    integer(kind=intType) :: mm

    ! Call the individual integration routines.
    do mm=1, nBocos
       ! Determine if this boundary condition is to be incldued in the
       ! currently active group
       famInclude: if (famInList(BCData(mm)%famID, famList)) then

          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible.
          call setBCPointers_d(mm, .True.)

          ! not post gathered integrations currently
          isWall: if( isWallType(BCType(mm)) ) then
             call wallIntegrationFace_b(localValues, localValuesd, mm)
          end if isWall

          isInflowOutflow: if (BCType(mm) == SubsonicInflow .or. &
               BCType(mm) == SupersonicInflow) then
             call flowIntegrationFace_b(.true., localValues, localValuesd, mm)
          else if (BCType(mm) == SubsonicOutflow .or. &
               BCType(mm) == SupersonicOutflow) then
             call flowIntegrationFace_b(.false., localValues, localValuesd, mm)
          end if isInflowOutflow
       end if famInclude
    end do
  end subroutine integrateSurfaces_b

  subroutine getSolution_d(famLists, funcValues, funcValuesd)
    !------------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    use constants
    use inputTSStabDeriv, only : TSSTability
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    use communication, only : adflow_comm_world
    use blockPointers, only : nDom
    use utils, only : setPointers_d, EChk
    use surfaceIntegrations_d, only : getCostFunctions_d
    use zipperIntegrations, only :integrateZippers_d
    use userSurfaceIntegrations, only : integrateUserSurfaces_d
    use actuatorRegion, only : integrateActuatorRegions_d
    implicit none

    ! Input/Output Variables
    integer(kind=intType), dimension(:, :), target, intent(in) :: famLists
    real(kind=realType), dimension(:, :), intent(out) :: funcValues, funcValuesd

    ! Working
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVald, globalVald
    integer(kind=intType) :: nn, sps, ierr, iGroup, nFam
    integer(kind=intType), dimension(:), pointer :: famList

    groupLoop: do iGroup=1, size(famLists, 1)

       ! Extract the current family list
       nFam = famLists(iGroup, 1)
       famList => famLists(iGroup, 2:2+nFam-1)

       localVal = zero
       localVald = zero
       do sps=1, nTimeIntervalsSpectral
          ! Integrate the normal block surfaces.
          do nn=1, nDom
             call setPointers_d(nn, 1, sps)
             call integrateSurfaces_d(localval(:, sps), localvald(:, sps), famList)
          end do

          ! Integrate any zippers we have
          call integrateZippers_d(localVal(:, sps), localVald(:, sps), famList, sps)

          ! Integrate any user-supplied surface as have as well.
          call integrateUserSurfaces_d(localVal(:, sps), localVald(:, sps), famList, sps)

          ! Integrate any actuator regions we have
          call integrateActuatorRegions_d(localVal(:, sps), localVald(:, sps), famList, sps)
       end do

       ! Now we need to reduce all the cost functions
       call mpi_allreduce(localval, globalVal, nLocalValues*nTimeIntervalsSpectral, adflow_real, &
            MPI_SUM, adflow_comm_world, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Now we need to reduce all the cost functions
       call mpi_allreduce(localvald, globalVald, nLocalValues*nTimeIntervalsSpectral, adflow_real, &
            MPI_SUM, adflow_comm_world, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Call the final routine that will comptue all of our functions of
       ! interest.

       call getCostFunctions_d(globalVal, globalVald, funcValues(:, iGroup), funcValuesd(:, iGroup))

       ! if (present(globalValues)) then
       !    globalValues = globalVal
       !    globalValuesd = globalVald
       ! end if
    end do groupLoop

  end subroutine getSolution_d

  subroutine getSolution_b(famLists, funcValues, funcValuesd)
    ! -----------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    use constants
    use communication, only : myid
    use inputTSStabDeriv, only : TSSTability
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    use communication, only : adflow_comm_world
    use blockPointers, only : nDom
    use utils, only : setPointers_b, EChk, setPointers
    use surfaceIntegrations_b, only : getCostFunctions_b
    use zipperIntegrations, only :integrateZippers_b
    use userSurfaceIntegrations, only : integrateUserSurfaces_b
    use actuatorRegion, only : integrateActuatorRegions_b
    implicit none

    ! Input/Output Variables
    integer(kind=intType), dimension(:, :), target, intent(in) :: famLists
    real(kind=realType), dimension(:, :) :: funcValues, funcValuesd

    ! Working
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVald, globalVald
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral, size(famLists, 1)) :: globalValues
    integer(kind=intType) :: nn, sps, ierr, iGroup, nFam
    integer(kind=intType), dimension(:), pointer :: famList


    call getSolution(famLists, funcValues, globalValues)

    groupLoop: do iGroup=1, size(famLists, 1)

       ! Extract the current family list
       nFam = famLists(iGroup, 1)
       famList => famLists(iGroup, 2:2+nFam-1)

       localVal = zero
       localVald = zero

       ! Retrive the forward pass values from getSolution
       globalVal = globalValues(:, :, iGroup)

       if (myid == 0) then
          call getCostFunctions_b(globalVal, globalVald, funcValues(:, iGroup), funcValuesd(:, iGroup))
          localVald = globalVald
       end if

       ! Now we need to bcast out the localValues to all procs.
       call mpi_bcast(localVald, nLocalValues*nTimeIntervalsSpectral, &
            adflow_real, 0, adflow_comm_world, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       do sps=1, nTimeIntervalsSpectral

          ! Integrate any actuator regions we have:
          call integrateActuatorRegions_b(localVal(:, sps), localVald(:, sps), famList, sps)

          ! Integrate any user-supplied planes as have as well.
          call integrateUserSurfaces_b(localVal(:, sps), localVald(:, sps), famList, sps)

          ! Integrate any zippers we have
          call integrateZippers_b(localVal(:, sps), localVald(:, sps), famList, sps)

          ! Integrate the normal block surfaces.
          do nn=1, nDom
             call setPointers_b(nn, 1, sps)
             call integrateSurfaces_b(localval(:, sps), localVald(:, sps), famList)
          end do

       end do
    end do groupLoop
  end subroutine getSolution_b
#endif
#endif
end module surfaceIntegrations
