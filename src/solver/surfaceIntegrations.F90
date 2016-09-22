module surfaceIntegrations

contains

  subroutine integrateSurfaces(localValues)
    ! This is a shell routine that calls the specific surface
    ! integration routines. Currently we have have the forceAndMoment
    ! routine as well as the flow properties routine. This routine
    ! takes care of setting pointers, while the actual computational
    ! routine just acts on a specific fast pointed to by pointers. 

    use constants
    use blockPointers, only : nBocos, BCData, BCType, sk, sj, si, x, rlv, &
         sfacei, sfacej, sfacek, gamma, rev, p, viscSubface
    use surfaceFamilies, only : famGroups
    use utils, only : setBCPointers, isWallType
    use sorting, only : bsearchIntegers
    use costFunctions, only : nLocalValues
    ! Tapenade needs to see these modules that the callees use.
    use BCPointers 
    use flowVarRefState
    use inputPhysics

    implicit none

    ! Input/output Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues

    ! Working variables
    integer(kind=intType) :: mm

    ! Loop over all possible boundary conditions
    bocos: do mm=1, nBocos
       
       ! Determine if this boundary condition is to be incldued in the
       ! currently active group
       famInclude: if (bsearchIntegers(BCdata(mm)%famID, &
            famGroups, size(famGroups)) > 0) then
          
          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible. 
          call setBCPointers(mm, .True.)

          isWall: if( isWallType(BCType(mm))) then 
             call wallIntegrationFace(localvalues, mm)
          end if isWall

#ifndef USE_TAPENADE
          isInflowOutflow: if (BCType(mm) == SubsonicInflow .or. &
               BCType(mm) == SubsonicOutflow .or. &
               BCType(mm) == SupersonicInflow .or. &
               BCType(mm) == SupersonicOutflow) then 
             call flowIntegrationFace(localValues, mm)
          end if isInflowOutflow
#endif

       end if famInclude
    end do bocos

  end subroutine integrateSurfaces

  subroutine flowIntegrationFace(localValues, mm)

    use constants
    use costFunctions
    use blockPointers, only : BCFaceID, BCData, addGridVelocities
    use costFunctions, onlY : nLocalValues, iMassFlow, iMassPtot, iMassTtot, iMassPs
    use sorting, only : bsearchIntegers
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef
    use inputPhysics, only : pointREf
    use flowUtils, only : computePtot, computeTtot
    use BCPointers, only : ssi, sFace, ww1, ww2, pp1, pp2, xx
    implicit none

    ! Input/Output variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), intent(in) :: mm

    ! Local variables
    real(kind=realType) ::  massFlowRate, mass_Ptot, mass_Ttot, mass_Ps
    integer(kind=intType) :: i, j, ii
    real(kind=realType) :: fact, xc, yc, zc, cellArea, mx, my, mz
    real(kind=realType) :: sF, vnm, vxm, vym, vzm, mReDim, Fx, Fy, Fz
    real(kind=realType) :: pm, Ptot, Ttot, rhom, massFlowRateLocal
    real(kind=realType), dimension(3) :: Fp, Mp, FMom, MMom, refPoint

    massFlowRate = zero
    mass_Ptot = zero
    mass_Ttot = zero
    mass_Ps = zero

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    select case (BCFaceID(mm))
    case (iMin, jMin, kMin)
       fact = -one
    case (iMax, jMax, kMax)
       fact = one
    end select

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

    !$AD II-LOOP
    do ii=0,(BCData(mm)%jnEnd - bcData(mm)%jnBeg)*(bcData(mm)%inEnd - bcData(mm)%inBeg) -1
       i = mod(ii, (bcData(mm)%inEnd-bcData(mm)%inBeg)) + bcData(mm)%inBeg + 1
       j = ii/(bcData(mm)%inEnd-bcData(mm)%inBeg) + bcData(mm)%jnBeg + 1
       
       if( addGridVelocities ) then 
          sF = sFace(i,j)
       else
          sF = zero
       end if

       vxm = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
       vym = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
       vzm = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
       rhom = half*(ww1(i,j,irho) + ww2(i,j,irho))
       pm = half*(pp1(i,j)+ pp2(i,j))

       vnm = vxm*ssi(i,j,1) + vym*ssi(i,j,2) + vzm*ssi(i,j,3)  - sF
       
       call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
       call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)

       pm = pm*pRef
       
       massFlowRateLocal = rhom*vnm*fact*mReDim
       massFlowRate = massFlowRate + massFlowRateLocal

       mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
       mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
       mass_Ps = mass_Ps + pm*massFlowRateLocal

       xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
            +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
       yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
            +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
       zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
            +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

       ! Compute the force components.
       ! blk = max(BCData(mm)%iblank(i,j), 0) ! iBlank forces for overset stuff
       
       fx = pm*ssi(i,j,1)
       fy = pm*ssi(i,j,2)
       fz = pm*ssi(i,j,3)

       ! Pressure forces
       ! fx = fx*blk
       ! fy = fy*blk
       ! fz = fz*blk

       ! Update the pressure force and moment coefficients.
       Fp(1) = Fp(1) + fx*fact
       Fp(2) = Fp(2) + fy*fact
       Fp(3) = Fp(3) + fz*fact
                 
       mx = yc*fz - zc*fy
       my = zc*fx - xc*fz
       mz = xc*fy - yc*fx
       
       Mp(1) = Mp(1) + mx
       Mp(2) = Mp(2) + my
       Mp(3) = Mp(3) + mz
       
       ! Momentum forces
       cellArea = sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2)

       fx = massFlowRateLocal*BCData(mm)%norm(i,j,1)*vxm/timeRef
       fy = massFlowRateLocal*BCData(mm)%norm(i,j,2)*vym/timeRef
       fz = massFlowRateLocal*BCData(mm)%norm(i,j,3)*vzm/timeRef

       ! fx = fx*blk
       ! fy = fy*blk
       ! fz = fz*block

       ! Note: momentum forces have opposite sign to pressure forces
       FMom(1) = FMom(1) - fx*fact
       FMom(2) = FMom(2) - fy*fact
       FMom(3) = FMom(3) - fz*fact

       mx = yc*fz - zc*fy
       my = zc*fx - xc*fz
       mz = xc*fy - yc*fx

       MMom(1) = MMom(1) + mx
       MMom(2) = MMom(2) + my
       MMom(3) = MMom(3) + mz
       
    enddo

   ! Increment the local values array with what we computed here
    localValues(iMassFlow) = localValues(iMassFlow) + massFlowRate
    localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
    localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
    localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps
    localValues(iFlowFp:iFlowFp+2)   = localValues(iFlowFp:iFlowFp+2) + Fp
    localValues(iFlowFm:iFlowFm+2)   = localValues(iFlowFm:iFlowFm+2) + FMom
    localValues(iFlowMp:iFlowMp+2)   = localValues(iFlowMp:iFlowMp+2) + Mp
    localValues(iFlowMm:iFlowMm+2)   = localValues(iFlowMm:iFlowMm+2) + MMom

  end subroutine flowIntegrationFace

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
    use costFunctions
    use communication
    use blockPointers
    use flowVarRefState
    use inputPhysics, only : MachCoef, pointRef, velDirFreeStream, equations
    use costFunctions, only : nLocalValues, iFp, iFv, iMp, iMv, iSepSensor, &
         iSepAvg, iCavitation, sepSensorSharpness, sepSensorOffset, iYplus
    use sorting, only :bsearchIntegers
    use BCPointers
    implicit none
  
    ! Input/output variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType) :: mm

    ! Local variables.
    real(kind=realType), dimension(3)  :: Fp, Fv, Mp, Mv
    real(kind=realType)  :: yplusMax, sepSensor, sepSensorAvg(3), Cavitation
    integer(kind=intType) :: i, j, ii, blk

    real(kind=realType) :: pm1, fx, fy, fz, fn, sigma
    real(kind=realType) :: xc, yc, zc, qf(3)
    real(kind=realType) :: fact, rho, mul, yplus, dwall
    real(kind=realType) :: V(3), sensor, sensor1, Cp, tmp, plocal
    real(kind=realType) :: tauXx, tauYy, tauZz
    real(kind=realType) :: tauXy, tauXz, tauYz

    real(kind=realType), dimension(3) :: refPoint
    real(kind=realType) :: mx, my, mz, cellArea

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

    ! Initialize the force and moment coefficients to 0 as well as
    ! yplusMax.

    Fp = zero; Fv = zero;
    Mp = zero; Mv = zero;
    yplusMax = zero
    sepSensor = zero
    Cavitation = zero
    sepSensorAvg = zero

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
       
       ! iBlank forces
       fx = fx*blk
       fy = fy*blk
       fz = fz*blk

       ! Update the inviscid force and moment coefficients.
       Fp(1) = Fp(1) + fx
       Fp(2) = Fp(2) + fy
       Fp(3) = Fp(3) + fz
       
       mx = yc*fz - zc*fy
       my = zc*fx - xc*fz
       mz = xc*fy - yc*fx

       Mp(1) = Mp(1) + mx
       Mp(2) = Mp(2) + my
       Mp(3) = Mp(3) + mz
       
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
       
       ! And integrate over the area of this cell and save:
       sensor = sensor * cellArea
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

       plocal = pp2(i,j)
       tmp = two/(gammaInf*MachCoef*MachCoef)
       Cp = tmp*(plocal-pinf)
       Sigma = 1.4
       Sensor1 = -Cp - Sigma
       Sensor1 = one/(one+exp(-2*10*Sensor1))
       Sensor1 = Sensor1 * cellArea
       Cavitation = Cavitation + Sensor1
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
          
          ! iBlank forces after saving for zipper mesh
          tauXx = tauXx*blk
          tauYy = tauYy*blk
          tauZz = tauZz*blk
          tauXy = tauXy*blk
          tauXz = tauXz*blk
          tauYz = tauYz*blk
          
          fx = fx*blk
          fy = fy*blk
          fz = fz*blk

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
          
          ! Update the viscous force and moment coefficients.

          Fv(1) = Fv(1) + fx
          Fv(2) = Fv(2) + fy
          Fv(3) = Fv(3) + fz
          
          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx

          Mv(1) = Mv(1) + mx
          Mv(2) = Mv(2) + my
          Mv(3) = Mv(3) + mz

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
             ! currently stored value.
             
             yplusMax = max(yplusMax, yplus)
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
#ifndef USE_TAPENADE
    localValues(iyPlus) = max(localValues(iyPlus), yplusMax)
#endif
  end subroutine wallIntegrationFace

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------


#ifndef USE_TAPENADE
  subroutine computeAeroCoef(globalCFVals,sps)
    !
    !      Compute the aerodynamic coefficients from the force and moment
    !      produced by the pressure and shear stresses on the body walls:
    !
    use constants
    use blockPointers
    use communication, only : adflow_comm_world, myid
    use inputPhysics
    use iteration
    use costFunctions
    use inputTimeSpectral
    use flowVarRefState
    use costFunctions, only : nLocalValues
    use overset, only : oversetPresent
    use utils, only : EChk, setPointers
    implicit none

    ! Input/Ouput Variables
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), intent(out), dimension(nCostFunction)::globalCFVals

    !      Local variables.
    integer(kind=intType) :: nn, ierr
    real(kind=realType), dimension(3) :: force, moment, cForce, cMoment
    real(kind=realType) :: fact, cd, cl, lift, drag
    real(kind=realType) ::localValues(nLocalValues)
    real(kind=realType), dimension(nCostFunction)::localCFVals

    localValues = zero
    domains: do nn=1,nDom
       call setPointers(nn,1_intType,sps)
       call integrateSurfaces(localValues)

    end do domains

    if (oversetPresent) then
       call wallIntegrationsZipper(localValues, sps)
    end if

    ! Sum pressure and viscous contributions from the walls, and
    ! pressure and momentum contributions from the inflow/outflow
    ! conditions. 
    Force = localValues(iFp:iFp+2) + localValues(iFv:iFv+2) + &
            localValues(iFlowFp:iFlowFp+2) + localValues(iFlowFm:iFlowFm+2)
    Moment = localValues(iMp:iMp+2) + localValues(iMv:iMv+2) + &
             localValues(iFlowMp:iFlowMp+2) + localValues(iFlowMm:iFlowMm+2)

    fact = two/(gammaInf*MachCoef*MachCoef &
         *surfaceRef*LRef*LRef*pRef)
    cForce = fact*force

    ! Moment factor has an extra lengthRef
    fact = fact/(lengthRef*LRef)
    cMoment = fact*Moment

    ! Get Lift coef and Drag coef
    CD =  cForce(1)*dragDirection(1) &
         + cForce(2)*dragDirection(2) &
         + cForce(3)*dragDirection(3)

    CL =  cForce(1)*liftDirection(1) &
         + cForce(2)*liftDirection(2) &
         + cForce(3)*liftDirection(3)

    Drag = Force(1)*dragDirection(1) &
         + Force(2)*dragDirection(2) &
         + Force(3)*dragDirection(3)

    Lift=  Force(1)*liftDirection(1) &
         + Force(2)*liftDirection(2) &
         + Force(3)*liftDirection(3)

    localCFVals(costFuncLift) = Lift
    localCFVals(costFuncDrag) = Drag
    localCFVals(costFuncLiftCoef) = Cl
    localCFVals(costFuncDragCoef) = Cd
    localCFVals(costFuncForceX) = Force(1)
    localCFVals(costFuncForceY) = Force(2)
    localCFVals(costFuncForceZ) = Force(3)
    localCFVals(costFuncForceXCoef) = cForce(1)
    localCFVals(costFuncForceYCoef) = cForce(2)
    localCFVals(costFuncForceZCoef) = cForce(3)
    localCFVals(costFuncMomX) = moment(1)
    localCFVals(costFuncMomY) = moment(2)
    localCFVals(costFuncMomZ) = moment(3)
    localCFVals(costFuncMomXCoef) = cmoment(1)
    localCFVals(costFuncMomYCoef) = cmoment(2)
    localCFVals(costFuncMomZCoef) = cmoment(3)
    localCFVals(costFuncSepSensor) = localValues(iSepSensor)
    localCFVals(costFuncSepSensorAvgX) = localValues(iSepAvg+0)
    localCFVals(costFuncSepSensorAvgY) = localValues(iSepAvg+1)
    localCFVals(costFuncSepSensorAvgZ) = localValues(iSepAvg+2)
    localCFVals(costFuncCavitation) = localValues(iCavitation)
    localCFVals(costFuncMdot) = localValues(iMassFlow)
    localCFVals(costFuncMavgPtot) = localValues(iMassPtot)
    localCFVals(costFuncMavgTtot) = localValues(iMassTtot)
    localCFVals(costFuncMavgPs) = localValues(iMassPs)

    ! Now we will mpi_allReduce them into globalCFVals
    call mpi_allreduce(localCFVals, globalCFVals, nCostFunction, adflow_real, &
         mpi_sum, ADflow_comm_world, ierr)

    globalCFVals(costFuncMavgPtot) = globalCFVals(costFuncMavgPtot)/globalCFVals(costFuncMdot)
    globalCFVals(costFuncMavgTtot) = globalCFVals(costFuncMavgTtot)/globalCFVals(costFuncMdot)
    globalCFVals(costFuncMavgPs) = globalCFVals(costFuncMavgPs)/globalCFVals(costFuncMdot)
    globalCFVals(costFuncMdot) = globalCFVals(costFuncMdot)

  end subroutine computeAeroCoef

  subroutine getSolution(sps)
    use constants
    use costFunctions
    use inputTSStabDeriv, only : TSSTability
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    use communication, only : adflow_comm_world
    use blockPointers, only : nDom
    use utils, only : computeTSDerivatives, computeRootBendingMoment, getDirAngle
    use inputPhysics, only : velDirFreeStream, liftDirection, dragDirection
    implicit none

    ! Input Variables
    integer(kind=intType) :: sps

    !  Local variables.
    real(kind=realType)   :: alpha, beta
    real(kind=realType), dimension(8) ::  dcdq, dcdqdot
    real(kind=realType), dimension(8) :: dcdalpha, dcdalphadot
    real(kind=realType), dimension(8) :: Coef0
    real(kind=realType), dimension(nCostFunction)::globalCFVals
    real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment
    real(kind=realType)::bendingMoment,bendingSum, cf(3), cm(3)
    integer(kind=intType) :: i

    funcValues(:) = 0.0

    bendingSum = 0.0
    do i =1,nTimeIntervalsSpectral
       call computeAeroCoef(globalCFVals,i)

       force(1, i) = globalCFVals(costFuncForceX)
       force(2, i) = globalCFVals(costFuncForceY)
       force(3, i) = globalCFVals(costFuncForceZ)
       moment(1, i) = globalCFVals(costFuncMomX)
       moment(2, i) = globalCFVals(costFuncMomY)
       moment(3, i) = globalCFVals(costFuncMomZ)

       cf = (/globalCFVals(costFuncForceXCoef), &
            globalCFVals(costFuncForceYCoef), &
            globalCFVals(costFuncForceZCoef)/)

       cm = (/globalCFVals(costFuncMomXCoef), &
            globalCFVals(costFuncMomYCoef), &
            globalCFVals(costFuncMomZCoef)/)
       call computeRootBendingMoment(cf, cm, bendingMoment)
       bendingsum = bendingsum+bendingMoment
    end do

    call computeAeroCoef(globalCFVals, sps)
    funcValues(costFuncBendingCoef)=bendingSum/nTimeIntervalsSpectral

    funcValues(costFuncLift) = globalCFVals(costFuncLift)
    funcValues(costFuncDrag) = globalCFVals(costFuncDrag)
    funcValues(costFuncLiftCoef) = globalCFVals(costFuncLiftCoef)
    funcValues(costFuncDragCoef) = globalCFVals(costFuncDragCoef)
    funcValues(costFuncForceX) = globalCFVals(costFuncForceX)
    funcValues(costFuncForceY) = globalCFVals(costFuncForceY)
    funcValues(costFuncForceZ) = globalCFVals(costFuncForceZ)
    funcValues(costFuncForceXCoef) = globalCFVals(costFuncForceXCoef)
    funcValues(costFuncForceYCoef) = globalCFVals(costFuncForceYCoef)
    funcValues(costFuncForceZCoef) = globalCFVals(costFuncForceZCoef)
    funcValues(costFuncMomX) = globalCFVals(costFuncMomX)
    funcValues(costFuncMomY) = globalCFVals(costFuncMomY)
    funcValues(costFuncMomZ) = globalCFVals(costFuncMomZ)
    funcValues(costFuncMomXCoef) = globalCFVals(costFuncMomXCoef)
    funcValues(costFuncMomYCoef) = globalCFVals(costFuncMomYCoef)
    funcValues(costFuncMomZCoef) = globalCFVals(costFuncMomZCoef)
    funcValues(costFuncSepSensor) = globalCFVals(costFuncSepSensor)
    funcValues(costFuncSepSensorAvgX) = globalCFVals(costFuncSepSensorAvgX)
    funcValues(costFuncSepSensorAvgY) = globalCFVals(costFuncSepSensorAvgY)
    funcValues(costFuncSepSensorAvgZ) = globalCFVals(costFuncSepSensorAvgZ)

    funcValues(costFuncCavitation) = globalCFVals(costFuncCavitation)
    funcValues(costFuncMdot) = globalCFVals(costFuncMdot)
    funcValues(costFuncMavgPtot) = globalCFVals(costFuncMavgPtot)
    funcValues(costFuncMavgTtot) = globalCFVals(costFuncMavgTtot)
    funcValues(costFuncMavgPs) = globalCFVals(costFuncMavgPs)

    if(TSStability)then

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
  end subroutine getSolution

  subroutine wallIntegrationsZipper(localValues, sps)

    use communication
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use costFunctions
    use overset, only : nodeZipperScatter, globalNodalVec, localZipperNodes, localZipperTp, localZipperTv
    use inputTimeSpectral
    use inputIteration
    use costFunctions, only : nLocalValues
    use utils, only : EChk, setPointers, myNorm2
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    ! Input/Output
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), intent(inout) :: localValues(nLocalValues)
    ! Working
    real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv

    integer(kind=intType) :: i, j, k, l, ierr, nn, mm, gind, ind, iVar
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, rowStart, rowEnd
    real(kind=realType), dimension(:, :, :), pointer :: xx
    real(kind=realType), dimension(:), pointer :: xPtr, pPtr, vPtr, localPtr
    integer(kind=intType), dimension(:, :), pointer :: gnp

    real(kind=realType), dimension(3) :: x1, x2, x3, xc, ss, norm, refPoint
    real(kind=realType), dimension(3) :: pp1, pp2, pp3, presForce
    real(kind=realType), dimension(3) :: vv1, vv2, vv3, viscForce
    real(kind=realType), dimension(3) ::  MTmp
    real(kind=realType) :: fact, triArea

    ! PETsc
    Vec, pointer :: localVec

    ! Determine the reference point for the moment computation in
    ! meters.
    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    ! Communicate the nodes, pressure tractions and viscous tractions to
    ! the root proc for the triangles. We do a generic loop

    call VecGetOwnershipRange(globalNodalVec, rowStart, rowEnd, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    do iVar=1, 3

       call VecGetArrayF90(globalNodalVec, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       domainLoop: do nn=1, nDom
          call setPointers(nn, 1, sps)

          bocoLoop: do mm=1, nBocos

             if(BCType(mm) == EulerWall .or. &
                  BCType(mm) == NSWallAdiabatic .or. &
                  BCType(mm) == NSWallIsothermal) then

                select case (BCFaceID(mm))
                case (iMin)
                   xx => x(1, :, :, :)
                   gnp => globalNode(1, :, :)
                case (iMax)
                   xx => x(il, :, :, :)
                   gnp => globalNode(il, :, :)
                case (jMin)
                   xx => x(:, 1, :, :)
                   gnp => globalNode(:, 1, :)
                case (jMax)
                   xx => x(:, jl, :, :)
                   gnp => globalNode(:, jl, :)
                case (kMin)
                   xx => x(:, :, 1, :)
                   gnp => globalNode(:, :, 1)
                case (kMax)
                   xx => x(:, :, kl, :)
                   gnp => globalNode(:, :, kl)
                end select

                jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
                iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd

                ! Loop over the nodes of the subface:
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      gInd = gnp(i+1, j+1)
                      ind = (gInd-rowStart/3)*3+1
                      if (ind < 0) then
                         print *,'something wrong:', myid, gind, rowstart, ind
                         stop
                      end if

                      select case(iVar)
                      case (1) ! Nodes
                         localPtr(ind:ind+2) = xx(i+1, j+1, :)
                      case (2)
                         localPtr(ind:ind+2) = bcData(mm)%Tp(i, j, :)
                      case (3)
                         localPtr(ind:ind+2) = bcData(mm)%Tv(i, j, :)
                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       call vecRestoreArrayF90(globalNodalVec, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       select case(iVar)
       case(1)
          localVec => localZipperNodes
       case(2)
          localVec => localZipperTp
       case(3)
          localVec => localZipperTv
       end select

       ! Perform the scatter from the global x vector to zipperNodes
       call VecScatterBegin(nodeZipperScatter, globalNodalVec, &
            localVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(nodeZipperScatter, globalNodalVec, &
            localVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

    end do


    if (myid == 0) then

       Fp = zero
       Fv = zero
       Mp = zero
       Mv = zero

       ! Puck out pointers for the nodes and tractions
       call VecGetArrayF90(localZipperTp, pPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecGetArrayF90(localZipperTv, vPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecGetArrayF90(localZipperNodes, xPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Number of triangles is the length of localZipperNodes /9

       do i=1, size(xPtr)/9

          ! Nodes for the triangles
          x1 = xPtr((i-1)*9+1:(i-1)*9+3)
          x2 = xPtr((i-1)*9+4:(i-1)*9+6)
          x3 = xPtr((i-1)*9+7:(i-1)*9+9)

          ! Nodal pressure tractions
          pp1 = pPtr((i-1)*9+1:(i-1)*9+3)
          pp2 = pPtr((i-1)*9+4:(i-1)*9+6)
          pp3 = pPtr((i-1)*9+7:(i-1)*9+9)

          ! Nodal viscous tractions
          vv1 = vPtr((i-1)*9+1:(i-1)*9+3)
          vv2 = vPtr((i-1)*9+4:(i-1)*9+6)
          vv3 = vPtr((i-1)*9+7:(i-1)*9+9)

          ! Compute area
          call cross_prod(x2-x1, x3-x1, norm)
          ss = half * norm
          triArea = mynorm2(ss)

          ! This is the actual integration
          presForce = third*(pp1 + pp2 + pp3) * triArea
          viscForce = third*(vv1 + vv2 + vv3) * triArea

          ! Add to Fp and Fv
          Fp = Fp + presForce
          Fv = Fv + viscForce

          ! Compute cell center
          xc(:) = third*(x1 + x2 + x3) - refPoint(:)

          ! Add to Mp and Mv
          call cross_prod(xc, presForce, Mtmp)
          Mp = Mp + Mtmp

          call cross_prod(xc, viscForce, Mtmp)
          Mv = Mv + Mtmp

       end do

       ! Return the array pointers
       call VecRestoreArrayF90(localZipperNodes, xPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecRestoreArrayF90(localZipperTp, pPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecRestoreArrayF90(localZipperTv, vPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Increment into the local vector
       localValues(iFp:iFp+2) = localValues(iFp:iFp+2) + Fp
       localValues(iFv:iFv+2) = localValues(iFv:iFv+2) + Fv
       localValues(iMp:iMp+2) = localValues(iMp:iMp+2) + Mp
       localValues(iMv:iMv+2) = localValues(iMv:iMv+2) + Mv


    end if

  end subroutine wallIntegrationsZipper
#endif
end module surfaceIntegrations
