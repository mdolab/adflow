module surfaceIntegrations

  use constants
  use communication, only : commType, internalCommType

  type userSurfCommType
     ! Data required on each proc:

     ! nDonor: The number of donor points the proc will provide
     ! frac (3, nDonor) : The uvw coordinates of the interpolation point
     ! donorInfo(4, nDonor) : Donor information. 1 is the local block ID and 2-4 is the 
     !    starting i,j,k indices for the interpolation. 
     ! procSizes(0:nProc-1) : The number of donors on each proc
     ! procDisps(0:nProc) : Cumulative form of procSizes

     ! inv(nConn) : Array allocated only on root processor used to
     ! reorder the nodes or elements back to the original order. 

     integer(kind=intType) :: nDonor
     real(kind=realType), dimension(:,:), allocatable :: frac
     integer(kind=intType), dimension(:, :), allocatable :: donorInfo
     integer(kind=intTYpe), dimension(:), allocatable :: procSizes, procDisps
     integer(kind=intTYpe), dimension(:), allocatable :: inv
     logical, dimension(:), allocatable :: valid

  end type userSurfCommType

  type userIntSurf

     character(len=maxStringLen) :: famName
     integer(Kind=intType) :: famID
     real(kind=realType), dimension(:, :), allocatable :: pts
     integer(kind=intType), dimension(:, :), allocatable :: conn
     
     ! Two separate commes: One for the nodes (based on the primal
     ! mesh) and one for the variables (based on the dual mesh)
     type(userSurfCommType) :: nodeComm, faceComm

  end type userIntSurf

  integer(kind=intType), parameter :: nUserIntSurfsMax=25
  type(userIntSurf), dimension(nUserIntSurfsMax), target :: userIntSurfs
  integer(kind=intTYpe) :: nUserIntSurfs=0

contains

  subroutine integrateSurfaces(localValues, famList)
    ! This is a shell routine that calls the specific surface
    ! integration routines. Currently we have have the forceAndMoment
    ! routine as well as the flow properties routine. This routine
    ! takes care of setting pointers, while the actual computational
    ! routine just acts on a specific fast pointed to by pointers. 

    use constants
    use blockPointers, only : nBocos, BCData, BCType, sk, sj, si, x, rlv, &
         sfacei, sfacej, sfacek, gamma, rev, p, viscSubface
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
    integer(kind=intType), dimension(:), intent(in) :: famList
    ! Working variables
    integer(kind=intType) :: mm

    ! Loop over all possible boundary conditions
    bocos: do mm=1, nBocos
       
       ! Determine if this boundary condition is to be incldued in the
       ! currently active group
       famInclude: if (bsearchIntegers(BCdata(mm)%famID, famList) > 0) then
          
          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible. 
          call setBCPointers(mm, .True.)

          isWall: if( isWallType(BCType(mm))) then 
             call wallIntegrationFace(localvalues, mm)
          end if isWall

          isInflowOutflow: if (BCType(mm) == SubsonicInflow .or. &
               BCType(mm) == SubsonicOutflow .or. &
               BCType(mm) == SupersonicInflow .or. &
               BCType(mm) == SupersonicOutflow) then 
             call flowIntegrationFace(localValues, mm)
          end if isInflowOutflow
       end if famInclude
    end do bocos

  end subroutine integrateSurfaces

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

  subroutine flowIntegrationFace(localValues, mm)

    use constants
    use costFunctions
    use blockPointers, only : BCFaceID, BCData, addGridVelocities
    use costFunctions, onlY : nLocalValues, iMassFlow, iMassPtot, iMassTtot, iMassPs
    use sorting, only : bsearchIntegers
    use flowVarRefState, only : pRef, rhoRef, timeRef, LRef, TRef, RGas
    use inputPhysics, only : pointRef
    use flowUtils, only : computePtot, computeTtot
    use BCPointers, only : ssi, sFace, ww1, ww2, pp1, pp2, xx
    implicit none

    ! Input/Output variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), intent(in) :: mm

    ! Local variables
    real(kind=realType) ::  massFlowRate, mass_Ptot, mass_Ttot, mass_Ps
    integer(kind=intType) :: i, j, ii, blk
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

       vnm = vxm*ssi(i,j,1) + vym*ssi(i,j,2) + vzm*ssi(i,j,3)  - sF
       
       call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
       call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)

       pm = pm*pRef
       
       massFlowRateLocal = rhom*vnm*fact*mReDim*blk
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

       ! Pressure forces. Note that these need a *negative* sign to be
       ! consistent with the force computation. 
       fx = -pm*ssi(i,j,1)*blk*fact
       fy = -pm*ssi(i,j,2)*blk*fact
       fz = -pm*ssi(i,j,3)*blk*fact

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
       
       ! Momentum forces. These gets a little more complex: TheB
       ! BCData(mm)%norm already is setup such that it points *out* of
       ! the domain. The use of fact here *reverses* the factor that
       ! has already been taken into account on massFLowRateLocal.
       fx = -massFlowRateLocal * BCData(mm)%norm(i, j, 1)*vxm/timeRef*fact*blk
       fy = -massFlowRateLocal * BCData(mm)%norm(i, j, 2)*vym/timeRef*fact*blk
       fz = -massFlowRateLocal * BCData(mm)%norm(i, j, 3)*vzm/timeRef*fact*blk

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
    localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
    localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
    localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps
    localValues(iFlowFp:iFlowFp+2)   = localValues(iFlowFp:iFlowFp+2) + Fp
    localValues(iFlowFm:iFlowFm+2)   = localValues(iFlowFm:iFlowFm+2) + FMom
    localValues(iFlowMp:iFlowMp+2)   = localValues(iFlowMp:iFlowMp+2) + Mp
    localValues(iFlowMm:iFlowMm+2)   = localValues(iFlowMm:iFlowMm+2) + MMom

  end subroutine flowIntegrationFace
  subroutine flowIntegrationZipper(vars, localValues, famList, sps)

    ! Integrate over the trianges for the inflow/outflow conditions. 

    use constants
    use costFunctions, only : nLocalValues, iMassFlow, iMassPtot, iMassTtot, iMassPs, &
         iFlowMm, iFlowMp, iFlowFm, iFlowFp
    use sorting, only : bsearchIntegers
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef, rGas
    use inputPhysics, only : pointRef
    use flowUtils, only : computePtot, computeTtot
    use overset, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : mynorm2, cross_prod
    implicit none

    ! Input/output Variables
    real(kind=realType), dimension(:, :), intent(in) :: vars
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps

    ! Working variables
    integer(kind=intType) :: i, j
    real(kind=realType) :: sF, vnm, vxm, vym, vzm, mReDim, Fx, Fy, Fz
    real(kind=realType), dimension(3) :: Fp, Mp, FMom, MMom, refPoint, ss, x1, x2, x3, norm
    real(kind=realType) :: pm, Ptot, Ttot, rhom, massFlowRateLocal
    real(kind=realType) ::  massFlowRate, mass_Ptot, mass_Ttot, mass_Ps
    real(kind=realType) :: fact, xc, yc, zc, cellArea, mx, my, mz

    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    zipper => zipperMeshes(iBCGroupInflowOutFlow)

    massFlowRate = zero
    mass_Ptot = zero
    mass_Ttot = zero
    mass_Ps = zero

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    mReDim = sqrt(pRef*rhoRef)
    Fp = zero
    Mp = zero
    FMom = zero
    MMom = zero

    !$AD II-LOOP
    do i=1, size(zipper%conn, 2)
       if (bsearchIntegers(zipper%fam(i), famList) > 0) then 
          ! Compute the averaged values for this trianlge
          vxm = zero; vym = zero; vzm = zero; rhom = zero; pm = zero;
          sF = zero
          do j=1,3
             rhom = rhom + vars(zipper%conn(j, i), iRho)
             vxm = vxm + vars(zipper%conn(j, i), iVx)
             vym = vym + vars(zipper%conn(j, i), iVy)
             vzm = vzm + vars(zipper%conn(j, i), iVz)
             pm = pm + vars(zipper%conn(j, i), iRhoE)
             sF = sF + vars(zipper%conn(j, i), 6)
          end do

          ! Divide by 3 due to the summation above:
          rhom = third*rhom
          vxm = third*vxm
          vym = third*vym
          vzm = third*vzm
          pm = third*pm
          sF = third*sF

          ! Get the nodes of triangle.
          x1 = vars(zipper%conn(1, i), 7:9)
          x2 = vars(zipper%conn(2, i), 7:9)
          x3 = vars(zipper%conn(3, i), 7:9)
          call cross_prod(x2-x1, x3-x1, norm)
          ss = -half*norm

          call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
          call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)
          
          pm = pm*pRef
          vnm = vxm*ss(1) + vym*ss(2) + vzm*ss(3)  - sF
          
          massFlowRateLocal = rhom*vnm*mReDim
          massFlowRate = massFlowRate + massFlowRateLocal
          
          mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
          mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
          mass_Ps = mass_Ps + pm*massFlowRateLocal
          
          ! Compute the average cell center. 
          xc = zero
          yc = zero
          zc = zero
          do j=1,3
             xc = xc + (vars(zipper%conn(1, i), 7)) 
             yc = yc + (vars(zipper%conn(2, i), 8)) 
             zc = zc + (vars(zipper%conn(3, i), 9)) 
          end do

          ! Finish average for cell center
          xc = third*xc
          yc = third*yc
          zc = third*zc

          xc = xc - refPoint(1)
          yc = yc - refPoint(2)
          zc = zc - refPoint(3)

          fx = pm*ss(1)
          fy = pm*ss(2)
          fz = pm*ss(3)

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

          ! Momentum forces

          ! Get unit normal vector. 
          ss = ss/mynorm2(ss)

          fx = massFlowRateLocal*ss(1) * vxm/timeRef
          fy = massFlowRateLocal*ss(2) * vym/timeRef
          fz = massFlowRateLocal*ss(3) * vzm/timeRef

          ! Note: momentum forces have opposite sign to pressure forces
          FMom(1) = FMom(1) - fx
          FMom(2) = FMom(2) - fy
          FMom(3) = FMom(3) - fz

          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx
          
          MMom(1) = MMom(1) + mx
          MMom(2) = MMom(2) + my
          MMom(3) = MMom(3) + mz
       end if
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

  end subroutine flowIntegrationZipper

  subroutine wallIntegrationZipper(vars, localValues, famList, sps)

    use constants
    use costFunctions

    use sorting, only : bsearchIntegers
    use flowVarRefState, only : LRef
    use inputPhysics, only : pointRef
    use overset, only : zipperMeshes, zipperMesh
    use utils, only : mynorm2, cross_prod
    implicit none

    ! Input/Output
    real(kind=realType), intent(in), dimension(:, :) :: vars
    real(kind=realType), intent(inout) :: localValues(nLocalValues)
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps

    ! Working
    real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv
    type(zipperMesh), pointer :: zipper

    integer(kind=intType) :: i, j
    real(kind=realType), dimension(3) :: ss, norm, refPoint
    real(kind=realType), dimension(3) :: p1, p2, p3, v1, v2, v3, x1, x2, x3
    real(kind=realType) :: fact, triArea, fx, fy, fz, mx, my, mz, xc, yc, zc

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    zipper => zipperMeshes(iBCGroupWalls)

    ! Determine the reference point for the moment computation in
    ! meters.
    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)
    Fp = zero
    Fv = zero
    Mp = zero
    Mv = zero

    !$AD II-LOOP
    do i=1, size(zipper%conn, 2)
       if (bsearchIntegers(zipper%fam(i), famList) > 0) then 

          ! Get the nodes of triangle. The *3 is becuase of the
          ! blanket third above. 
          x1 = vars(zipper%conn(1, i), 7:9)
          x2 = vars(zipper%conn(2, i), 7:9)
          x3 = vars(zipper%conn(3, i), 7:9)
          call cross_prod(x2-x1, x3-x1, norm)
          ss = half*norm
          ! The third here is to account for the summation of P1, p2
          ! and P3
          triArea = mynorm2(ss)*third

          ! Compute the average cell center. 
          xc = zero
          yc = zero
          zc = zero
          do j=1,3
             xc = xc + (vars(zipper%conn(1, i), 7)) 
             yc = yc + (vars(zipper%conn(2, i), 8)) 
             zc = zc + (vars(zipper%conn(3, i), 9)) 
          end do

          ! Finish average for cell center          
          xc = third*xc
          yc = third*yc
          zc = third*zc

          xc = xc - refPoint(1)
          yc = yc - refPoint(2)
          zc = zc - refPoint(3)

          ! Update the pressure force and moment coefficients.
          p1 = vars(zipper%conn(1, i), 1:3) 
          p2 = vars(zipper%conn(2, i), 1:3)
          p3 = vars(zipper%conn(3, i), 1:3)

          fx = (p1(1) + p2(1) + p3(1))*triArea
          fy = (p1(2) + p2(2) + p3(2))*triArea
          fz = (p1(3) + p2(3) + p3(3))*triArea

          Fp(1) = Fp(1) + fx
          Fp(2) = Fp(2) + fy
          Fp(3) = Fp(3) + fz
                 
          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx
       
          Mp(1) = Mp(1) + mx
          Mp(2) = Mp(2) + my
          Mp(3) = Mp(3) + mz

          ! Update the viscous force and moment coefficients
          v1 = vars(zipper%conn(1, i), 4:6) 
          v2 = vars(zipper%conn(2, i), 4:6)
          v3 = vars(zipper%conn(3, i), 4:6)

          fx = (v1(1) + v2(1) + v3(1))*triArea
          fy = (v1(2) + v2(2) + v3(2))*triArea
          fz = (v1(3) + v2(3) + v3(3))*triArea

          ! Note: momentum forces have opposite sign to pressure forces
          Fv(1) = Fv(1) + fx
          Fv(2) = Fv(2) + fy
          Fv(3) = Fv(3) + fz

          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx
          
          Mv(1) = Mv(1) + mx
          Mv(2) = Mv(2) + my
          Mv(3) = Mv(3) + mz
       end if
    enddo

    ! Increment into the local vector
    localValues(iFp:iFp+2) = localValues(iFp:iFp+2) + Fp
    localValues(iFv:iFv+2) = localValues(iFv:iFv+2) + Fv
    localValues(iMp:iMp+2) = localValues(iMp:iMp+2) + Mp
    localValues(iMv:iMv+2) = localValues(iMv:iMv+2) + Mv
    
  end subroutine wallIntegrationZipper

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE

  subroutine integrateZippers(localValues, famList, sps)

    ! Integrate over the triangles formed by the zipper mesh. This
    ! will perform both all necesasry zipper integrations. Currently
    ! this includes the wall force integrations as well as the
    ! flow-though surface integration. 

    use constants
    use costFunctions, only : nLocalValues
    use overset, only : zipperMeshes
    use haloExchange, only : wallIntegrationZipperComm, flowIntegrationZipperComm
    implicit none

    ! Input Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), dimension(:, :), allocatable :: vars

    ! Determine if we have a wall Zipper:
    if (zipperMeshes(iBCGroupWalls)%allocated) then 
       
       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call wallIntegrationZipperComm(vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call wallIntegrationZipper(vars, localValues, famList, sps)

       ! Cleanup vars
       deallocate(vars)
    end if

    ! Determine if we have a flowthrough Zipper:
    if (zipperMeshes(iBCGroupInflowOutflow)%allocated) then 
       
       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call flowIntegrationZipperComm(vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call flowIntegrationZipper(vars, localValues, famList, sps)

       ! Cleanup vars
       deallocate(vars)
    end if
  end subroutine integrateZippers

  ! subroutine integrateZippers_d(localValues, localValuesd, famList, sps)

  !   ! Forward mode linearization of the zipper integration. 

  !   use constants
  !   use costFunctions, only : nLocalValues
  !   use overset, only : zipperMeshes
  !   use haloExchange, only : wallIntegrationZipperComm_d, flowIntegrationZipperComm_d
  !   implicit none

  !   ! Input Variables
  !   real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
  !   integer(kind=intType), dimension(:), intent(in) :: famList
  !   integer(kind=intType), intent(in) :: sps
  !   real(kind=realType), dimension(:, :), allocatable :: vars, varsd

  !   ! Determine if we have a wall Zipper:
  !   if (zipperMeshes(iBCGroupWalls)%allocated) then 
       
  !      ! Gather up the required variables in "vars" on the root
  !      ! proc. This routine is differientated by hand. 
  !      call wallIntegrationZipperComm_d(vars, varsd, sps)

  !      ! Perform actual integration. Tapenade ADs this routine.
  !      call wallIntegrationZipper_d(vars, varsd, localValues, localValuesd, famList)

  !      ! Cleanup vars
  !      deallocate(vars, varsd)
  !   end if

  !   ! Determine if we have a flowthrough Zipper:
  !   if (zipperMeshes(iBCGroupInflowOutflow)%allocated) then 
       
  !      ! Gather up the required variables in "vars" on the root
  !      ! proc. This routine is differientated by hand. 
  !      call flowIntegrationZipperComm_d(vars, varsd, sps)

  !      ! Perform actual integration. Tapenade ADs this routine.
  !      call flowIntegrationZipper_d(vars, varsd, localValues, localValuesd, famList)

  !      ! Cleanup vars
  !      deallocate(vars, varsd)
  !   end if
  ! end subroutine integrateZippers_d

  ! subroutine integrateZippers_b(localValues, localValuesd, famList, sps)

  !   ! Forward mode linearization of the zipper integration. 

  !   use constants
  !   use costFunctions, only : nLocalValues
  !   use overset, only : zipperMeshes
  !   !use haloExchange, only : wallIntegrationZipperComm_b, flowIntegrationZipperComm_b
    
  !   implicit none

  !   ! Input Variables
  !   real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
  !   integer(kind=intType), dimension(:), intent(in) :: famList
  !   integer(kind=intType), intent(in) :: sps
  !   real(kind=realType), dimension(:, :), allocatable :: vars, varsd

  !   ! Determine if we have a wall Zipper:
  !   if (zipperMeshes(iBCGroupWalls)%allocated) then 

  !      ! Perform actual integration. Tapenade ADs this routine.
  !      call wallIntegrationZipper_b(vars, varsd, localValues, localValuesd, famList)
       
  !      ! Scatter (becuase we are reverse) the values from the root
  !      ! back out to all necessary procs.  This routine is
  !      ! differientated by hand.
  !      !call wallIntegrationZipperComm_b(vars, varsd, sps)

  !      ! Cleanup vars
  !      deallocate(vars, varsd)
  !   end if

  !   ! Determine if we have a flowthrough Zipper:
  !   if (zipperMeshes(iBCGroupInflowOutflow)%allocated) then 

  !      ! Perform actual integration. Tapenade ADs this routine.
  !      call flowIntegrationZipper_b(vars, varsd, localValues, localValuesd, famList)
       
  !      ! Scatter (becuase we are reverse) the values from the root
  !      ! back out to all necessary procs.  This routine is
  !      ! differientated by hand.
  !      !call flowIntegrationZipperComm_b(vars, varsd, sps)

  !      ! Cleanup vars
  !      deallocate(vars, varsd)
  !   end if
  ! end subroutine integrateZippers_b


  subroutine computeAeroCoef(globalCFVals,sps, famList)
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
    integer(kind=intType), dimension(:), intent(in) :: famList

    !      Local variables.
    integer(kind=intType) :: nn, ierr
    real(kind=realType), dimension(3) :: force, moment, cForce, cMoment
    real(kind=realType) :: fact, cd, cl, lift, drag
    real(kind=realType) ::localValues(nLocalValues)
    real(kind=realType), dimension(nCostFunction)::localCFVals

    localValues = zero
    domains: do nn=1,nDom
       call setPointers(nn,1_intType,sps)
       call integrateSurfaces(localValues, famList)

    end do domains

    if (oversetPresent) then
       call integrateZippers(localValues, famList, sps)
    end if

    ! Integrate any user-supplied planes as well:
    call evalUserIntegrationSurfaces(localvalues, famList, sps)

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

  subroutine getSolution(sps, famList)
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
    integer(kind=intType), dimension(:), intent(in) :: famList

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
       call computeAeroCoef(globalCFVals,i, famList)

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

    call computeAeroCoef(globalCFVals, sps, famList)
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


  subroutine addIntegrationSurface(pts, conn, famName, famID, nPts, nConn)
    ! Add a user-supplied integration surface. 

    use communication, only : myID
    use constants

    implicit none

    ! Input variables
    real(kind=realType), dimension(3, nPts), intent(in) :: pts
    integer(kind=intType), dimension(4, nConn), intent(in) :: conn
    integer(kind=intType), intent(in) :: nPts, nConn, famID
    character(len=*) :: famName
    type(userIntSurf), pointer :: surf

    ! Not really much to do here...we just have to save the data
    ! into the data structure untilly we actual have to do the
    ! search.
    nUserIntSurfs = nUserIntSurfs + 1
    if (nUserIntSurfs > nUserIntSurfsMax) then 
       print *,"Error: Exceeded the maximum number of user-supplied "&
            &"integration slices. Increase nUserIntSurfsMax"
       stop
    end if
    
    if (myid == 0) then 
       surf => userIntSurfs(nUserIntSurfs)

       allocate(surf%pts(3, nPts), surf%conn(4, nConn))
       surf%pts = pts
       surf%conn = conn
       surf%famName = famName
       surf%famID = famID
    end if

  end subroutine addIntegrationSurface

  subroutine buildVolumeADTs(oBlocks, useDual)

    ! This builds volume ADTs for the the owned blocks. It will build
    ! either the dual mesh or the primal mesh depening on the flag
    ! useDual. 

    use constants
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl,  vol, ib, jb, kb
    use adtBuild, only : buildSerialHex
    use utils, only : setPointers, EChk
    implicit none

    ! Input/Output Parameters
    type(oversetBlock), dimension(:), target :: oBlocks
    logical :: useDual

    ! Working Parameters
    integer(kind=intType) :: nInterpol, nn, i, j, k, iii, jjj, kkk
    integer(kind=intType) :: ii, jj, kk, mm, nADT, nHexa, planeOffset
    type(oversetBlock), pointer :: oBlock

    nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

    domainLoop: do nn=1, nDom

       call setPointers(nn, 1, 1)
       oBlock => oBlocks(nn)

       primalOrDual: if (useDual) then 

          ! Now setup the data for the ADT
          nHexa = il * jl * kl
          nADT = ie * je * ke
          oBlock%il = il
          oBlock%jl = jl
          oBlock%kl = kl
          
          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))
          ! Fill up the xADT using cell centers (dual mesh)
          mm = 0
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = eighth*(&
                        x(i-1, j-1, k-1, :) + &
                        x(i  , j-1, k-1, :) + &
                        x(i-1, j  , k-1, :) + &
                        x(i  , j  , k-1, :) + &
                        x(i-1, j-1, k  , :) + &
                        x(i  , j-1, k  , :) + &
                        x(i-1, j  , k  , :) + &
                        x(i  , j  , k  , :))
                   oBlock%qualDonor(1, mm) = vol(i, j, k)
                end do
             end do
          end do
          
          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ie * je
          do k=2, ke
             do j=2, je
                do i=2, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 
                   
                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       else
          ! Note that we will be including the halo primal cells. This
          ! should slightly increase robusness for viscous off-wall
          ! spacing. This means the primal mesh has 1 MORE node/cell
          ! in each direction. 

          ! Now setup the data for the ADT
          nHexa = ie * je * ke
          nADT = ib * jb * kb
          oBlock%il = ie
          oBlock%jl = je
          oBlock%kl = ke

          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))

          oBlock%qualDonor = zero
          ! Fill up the xADT using the primal nodes
          mm = 0
          do k=0, ke
             do j=0, je
                do i=0, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = x(i, j, k, :)

                   ! Since we don't have all 8 volumes surrounding the
                   ! halo nodes, clip the volumes to be between 0 and ib etc.
                   do iii=0,1
                      do jjj=0,1
                         do kkk=0,1
                            ii = min(max(0, iii+i), ib)
                            jj = min(max(0, jjj+j), jb)
                            kk = min(max(0, kkk+k), kb)
                            
                            oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) + &
                                 vol(ii, jj, kk)
                         end do
                      end do
                   end do

                   ! Dividing by 8 isn't strictly necessary but we'll
                   ! do it anyway.
                   oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) * eighth
                end do
             end do
          end do
          
          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ib * jb
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-1)*planeOffset + (j-1)*ib + (i-1) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ib
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 
                   
                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       end if primalOrDual

       ! Call the custom build routine -- Serial only, only Hexa volumes,
       ! we supply our own ADT Type

       call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)
    end do domainLoop

  end subroutine buildVolumeADTs

  subroutine performInterpolation(pts, oBlocks, useDual, comm)!, oSurfs, comm)

    ! This routine performs the actual searches for the slices. It is
    ! generic in the sense that it will search an arbtitrary set of
    ! points on either the primal or dual meshes. The final required
    ! communication data is then written into the supplied comm. 
    
    use constants
    use block, only : fringeType
    use communication, only : adflow_comm_world, myid, nProc
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl, x, iBlank, vol
    use adtLocalSearch, only :  mindistancetreesearchsinglepoint, &
         containmenttreesearchsinglepoint
    use adtUtils, only : stack
    use adtData, only : adtBBoxTargetType
    use utils, only : setPointers, mynorm2, EChk
    use inputOverset, only : oversetProjTol
    use oversetUtilities, only : fracToWeights2, qsortFringeType, getCumulativeForm
    
    implicit none

    ! Input parameters:
    real(kind=realType), dimension(:, :), intent(in) :: pts
    type(userIntSurf) :: surf
    type(oversetBlock), dimension(:), target, intent(in) :: oBlocks
    logical, intent(in) :: useDual
    !type(oversetSurf), dimension(:), intent(in) :: oSurfs
    type(userSurfCommType) :: comm

    ! Working parameters
    type(oversetBlock), pointer :: oBlock
    type(fringeType), dimension(:), allocatable :: surfFringes

    integer(Kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, nn, mm 
    integer(kind=intType) :: iSurf, ierr, nInterpol, iProc
    integer(kind=intType) :: nHexa, nAdt, planeOffset, elemID, nPts
    real(kind=realType) :: xc(4), weight(8)
    integer status(MPI_STATUS_SIZE) 

    real(kind=realType) :: uvw(5), uvw2(5), donorQual, xcheck(3)
    integer(kind=intType) :: intInfo(3), intInfo2(3)
    logical :: failed, invalid

    integer(kind=intType), dimension(:, :), allocatable :: donorInfo, intSend
    integer(kind=intType), dimension(:), allocatable :: procSizes
    real(kind=realType), dimension(:, :), allocatable :: donorFrac, realSend

    ! Variables we have to pass the ADT search routine
    integer(kind=intType), dimension(:), pointer :: BB
    type(adtBBoxTargetType), dimension(:), pointer :: BB2
    integer(kind=intType), dimension(:), pointer :: frontLeaves
    integer(kind=intType), dimension(:), pointer :: frontLeavesNew
   
    ! Data for the search
    allocate(BB(20), frontLeaves(25), frontLeavesNew(25), stack(100))

    nPts = size(pts, 2)
              
    ! Allocate donor information arrays
    allocate(donorFrac(4, nPts), donorInfo(5, nPts))
    donorInfo = -1
    donorFrac(4, :) = large
    nInterpol = 1
    domainSearch: do nn=1, nDom
       oBlock => oBlocks(nn)
       call setPointers(nn, 1, 1)

       ! Search the supplied pts one at a time
       elemLoop: do i=1, nPts

          xc(1:3) = pts(:, i)
             
          ! Call the standard tree search
          call containmentTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
               oBlock%qualDonor, nInterpol, BB, frontLeaves, frontLeavesNew, failed)
          
          ! Make sure this point is not garbage.
          if (intInfo(1) >= 0) then 
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do
             
             if (mynorm2(xcheck - xc(1:3)) > oversetProjTol) then 
                failed = .True.
             end if
          end if
             
          if (intInfo(1) >= 0 .and. failed) then 
             ! we "found" a point but it is garbage. Do the failsafe search
             xc(4) = large
             call minDistanceTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
                  oBlock%qualDonor, nInterpol, BB2, frontLeaves, frontLeavesNew)
             
             ! Check this one:
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do
             
             ! Since this is the last line of defence, relax the tolerance a bit
             if (mynorm2(xcheck - xc(1:3)) > 100*oversetProjTol) then 
                ! This fringe has not found a donor
                intInfo(1) = -1
             else
                ! This one has now passed.
                
                ! Important! uvw(4) is the distance squared for this search
                ! not
                uvw(4) = uvw(5)
             end if
          end if
             
          elemFound: if (intInfo(1) >= 0) then 

             ! Donor and block and index information for this donor. 
             donorQual = uvw(4)
             elemID = intInfo(3) - 1 ! Make it zero based

             ! The dual mesh needs an offset of 1 becuase it only used
             ! 1:ie values. This is not necessary for the primal. 
             if (useDual) then 
                ii = mod(elemID, oBlock%il) + 1
                jj = mod(elemID/oBlock%il, oBlock%jl) + 1
                kk = elemID/(oBlock%il*oBlock%jl) + 1
             else
                ii = mod(elemID, oBlock%il) 
                jj = mod(elemID/oBlock%il, oBlock%jl) 
                kk = elemID/(oBlock%il*oBlock%jl) 
             end if
             ! Rememebr donorFrac(4, i) is the current best quality 
             if ( donorQual < donorFrac(4, i)) then 
                
                invalid = .False.

                ! For the dual mesh search, we have to make sure the
                ! potential donors are valid. Such a check is not
                ! necessary for the primal search since all nodes are
                ! considered valid. 

                if (useDual) then 
                   ! Check if the point is invalid. We can do this
                   ! with i-blank array. We can only accept compute
                   ! cells (iblank=1) or interpolated
                   ! cells(iblank=-1).
                   do kkk=0,1
                      do jjj=0,1
                         do iii=0,1
                            if (.not. (iblank(ii+iii, jj+jjj, kk+kkk) == 1 .or. &
                                 iblank(ii+iii, jj+jjj, kk+kkk) == -1)) then 
                               invalid = .True.
                            end if
                         end do
                      end do
                   end do
                end if
                
                if (.not. invalid) then 
                   
                   ! Set the quality of the donor to the one we
                   ! just found. Save the rest of the necessary
                   ! information. 
                   donorInfo(1, i) = myid
                   donorInfo(2, i) = nn
                   donorInfo(3, i) = ii
                   donorInfo(4, i) = jj
                   donorInfo(5, i) = kk
                   donorFrac(1:3, i) = uvw(1:3)
                   donorFrac(4, i) = donorQual
                end if
             end if
          end if elemFound
       end do elemLoop
    end do domainSearch
      
    ! Next count up the number of valid donors we've found and compact
    ! the info back to that length.
    if (myid /=0) then 
       j = 0
       do i=1, nPts
          if (donorInfo(1, i) /= -1) then 
             j = j + 1
          end if
       end do
       allocate(intSend(6, j), realSend(4, i))
       if (j > 0) then 
          j = 0
          do i=1, nPts
             if (donorInfo(1, i) /= -1) then 
                j = j + 1
                intSend(1:5, j) = donorInfo(:, i)
                intSend(6, j) = i
                realSend(:, j) = donorFrac(:, i)
             end if
          end do
       end if
    else
       ! On the root proc, use intSend and realSend as the receiver
       ! buffer. These can be at most nPts sized.
       allocate(intSend(6, nPts), realSend(4, nPts))
    end if
       
    ! Gather up the sizes (j) to the root processor so he know who to
    ! expect data from.
    allocate(procSizes(0:nProc-1))
    
    call mpi_gather(j, 1, adflow_integer, procSizes, 1, &
         adflow_integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    
    ! Next all the procs need to send all the information back to the
    ! root processor where we will determine the proper donors for
    ! each of the cells
    
    ! All procs except root fire off their data.
    if (myid >= 1) then 
       if (j > 0) then 
          call mpi_send(intSend, j*6, adflow_integer, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          
          call mpi_send(realSend, j*4, adflow_real, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end if
    end if
    
    ! And the root processor recieves it...
    if (myid == 0) then 
       do iProc=1, nProc-1
          ! Determine if this proc has sent anything:
          if (procSizes(iProc) /= 0) then 
             
             call MPI_recv(intSend, 6*nPts, adflow_integer, iProc, iProc,&
                  adflow_comm_world, status, ierr)
             call EChk(ierr,__FILE__,__LINE__)
             
             call MPI_recv(realSend, 4*nPts, adflow_real, iProc, iProc,&
                  adflow_comm_world, status, ierr)
             call EChk(ierr,__FILE__,__LINE__)
             
             ! Now process the data (intSend and realSend) that we
             ! just received. We don't need to check the status for
             ! the sizes becuase we already know the sizes from the
             ! initial gather we did. 
                
             do i=1, procSizes(iProc)
                ii = intSend(6, i)
                
                if (realSend(4, i) < donorFrac(4, ii)) then 
                   ! The incoming quality is better. Accept it. 
                   donorInfo(1:5, ii) = intSend(1:5, i)
                   donorFrac(:, ii) = realSend(:, i)
                end if
             end do
          end if
       end do
          
       ! To make this easier, convert the information we have to a
       ! 'fringeType' array so we can use the pre-existing sorting
       ! routine.
       allocate(surfFringes(nPts))
       do i=1, nPts
          surfFringes(i)%donorProc = donorInfo(1, i)
          surfFringes(i)%donorBlock = donorInfo(2, i)
          surfFringes(i)%dI = donorInfo(3, i)
          surfFringes(i)%dJ = donorInfo(4, i)
          surfFringes(i)%dK = donorInfo(5, i)
          surfFringes(i)%donorFrac = donorFrac(1:3, i)
          ! Use the myBlock attribute to keep track of the original
          ! index. When we sort the fringes, they will no longer be
          ! in the same order
          surfFringes(i)%myBlock = i
       end do
       
       ! Perform the actual sort. 
       call qsortFringeType(surfFringes, nPts)
       
       ! We will reuse-proc sizes to now mean the number of elements
       ! that the processor *actually* has to send. We will include
       ! the root proc itself in the calc becuase that will tell us
       ! the size of the internal comm structure. 
       
       procSizes = 0
       allocate(comm%valid(nPts))
       comm%valid = .True.
       do i=1, nPts
          if (surfFringes(i)%donorProc < 0) then 
             ! We dont have a donor. Flag this point as invalid
             comm%valid(i) = .False.
          end if

          ! Dump the points without donors on the root proc by making
          ! sure j is at least 0 for the root proc. These will just
          ! simply be ignored during the comm. 
          j = max(surfFringes(i)%donorProc, 0)
          procSizes(j) = procSizes(j) + 1
       end do
    end if
       
    ! Simply broadcast out the the proc sizes back to everyone so all
    ! processors know if they are to receive anything back. 
    call mpi_bcast(procSizes, nProc, adflow_Integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
       
    ! We can now save some of the final comm data required on the
    ! root proc for this surf. 
    allocate(comm%procSizes(0:nProc-1), comm%procDisps(0:nProc))
    
    ! Copy over procSizes and generate the cumulative form of
    ! the size array, procDisps
    comm%procSizes = procSizes
    call getCumulativeForm(comm%procSizes, nProc, comm%procDisps)

    ! Record the elemInverse which is necessary to index into
    ! the original conn array. 
    if (myid == 0) then 
       allocate(comm%Inv(nPts))
       do i=1, nPts
          comm%Inv(i) = surfFringes(i)%myBlock
       end do
    end if

    ! Now we can send out the final donor information to the
    ! processors that must supply it. 
    comm%nDonor = procSizes(myID)
    allocate(comm%frac(3, comm%nDonor), comm%donorInfo(4, comm%nDonor))
    
    if (myid >= 1) then 
       if (comm%nDonor > 0) then 
          ! We are responible for at least 1 donor. We have to make
          ! use of the intSend and realSend buffers again (which
          ! are guaranteed to be big enough). The reason we can't
          ! dump the data in directlyis that intSend and realSend
          ! have a different leading index than we need on the
          ! final data structure. 
          
          call MPI_recv(intSend, 6*comm%nDonor, adflow_integer, 0, myid, &
               adflow_comm_world, status, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          
          call MPI_recv(realSend, 4*comm%nDonor, adflow_real, 0, myID, &
               adflow_comm_world, status, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy into final structure
          do i=1, comm%nDonor
             comm%donorInfo(:, i) = intSend(1:4, i)
             comm%frac(:, i) = realSend(1:3, i)
          end do
       end if
    else
       ! We are the root processor.
       if (comm%nDonor > 0) then 
          ! We need to copy out our donor info on the root proc if we have any
          do i= comm%procDisps(myID)+1, comm%procDisps(myID+1)
             comm%donorInfo(1, i) = surfFringes(i)%donorBlock
             comm%donorInfo(2, i) = surfFringes(i)%dI
             comm%donorInfo(3, i) = surfFringes(i)%dJ
             comm%donorInfo(4, i) = surfFringes(i)%dK
             comm%frac(1:3, i)    = surfFringes(i)%donorFrac
          end do
       end if
          
       ! Now loop over the rest of the procs and send out the info we
       ! need. We have to temporarily copy the data back out of
       ! fringes to the intSend and realSend arrays
       do iProc=1, nProc-1
          
          if (comm%procSizes(iProc) > 0) then 
             ! Have something to send here:
             j = 0
             do i=comm%procDisps(iProc)+1, comm%procDisps(iProc+1)
                j = j + 1
                
                intSend(1, j) = surfFringes(i)%donorBlock
                intSend(2, j) = surfFringes(i)%dI      
                intSend(3, j) = surfFringes(i)%dJ
                intSend(4, j) = surfFringes(i)%dK
                realSend(1:3, j) = surfFringes(i)%donorFrac
             end do
             
             call mpi_send(intSend, j*6, adflow_integer, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)
             
             call mpi_send(realSend, j*4, adflow_real, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end if
       end do
       
       ! Deallocate data allocatd only on root proc
       deallocate(surfFringes)
    end if
       
    ! Nuke rest of allocated on all procs
    deallocate(intSend, realSend, procSizes, donorInfo, donorFrac)
 
  end subroutine performInterpolation

  subroutine interpolateIntegrationSurfaces

    ! This routine performs the actual searches for the slices. We
    ! reuse much of the same machinery as is used in the overset code.

    use constants
    use communication, only : adflow_comm_world, myid
    use overset, only : oversetBlock
    use blockPointers, only : nDom, ie, je, ke, il, jl, kl
    use adtBuild, only : buildSerialHex, destroySerialHex
    use utils, only : setPointers, EChk

    implicit none

    ! Working parameters
    type(oversetBlock), dimension(nDom), target :: oBlocks
    type(userIntSurf), pointer :: surf

    integer(Kind=intType) :: iSurf, ii, i, nn, nPts, ierr
    real(kind=realType), dimension(:, :), allocatable :: pts
    logical :: useDual

    primalDualLoop: do ii=1, 2
       if (ii==1) then 
          useDual = .True.
       else
          useDual = .False.
       end if

       call buildVolumeADTs(oBlocks, useDual)
       
       masterLoop: do iSurf=1, nUserIntSurfs

          surf => userIntSurfs(iSurf)

          if (myid == 0) then 

             if (ii==1) then 
                nPts = size(surf%conn, 2)
                allocate(pts(3, nPts))
                ! Use the dual (cell center values)
                elemLoop: do i=1, nPts

                   ! Compute the center of the cell. 
                   pts(:, i) = fourth*( & 
                        surf%pts(:, surf%conn(1, i)) + &
                        surf%pts(:, surf%conn(2, i)) + &
                        surf%pts(:, surf%conn(3, i)) + &
                        surf%pts(:, surf%conn(4, i)))
                end do elemLoop
             else if (ii==2) then 
                nPts = size(surf%pts, 2)
                allocate(pts(3, nPts))

                ! Use the primal nodal values
                nodeLoop: do i=1, nPts
                   ! Compute the center of the cell. 
                   pts(:, i) = surf%pts(:, i)
                end do nodeLoop
             end if
          end if

          ! Send the number of points back to all procs:
          call mpi_bcast(nPts, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
           
          ! All other procs except the root allocate space and receive
          ! the pt array. 
          if (myid /= 0) then 
             allocate(pts(3, nPts))
          end if
          
          call mpi_bcast(pts, 3*nPts, adflow_real, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          
          ! Call the actual interpolation routine
          if (ii==1) then 
             call performInterpolation(pts, oBlocks, .True., surf%faceComm)
          else
             call performInterpolation(pts, oBlocks, .False., surf%nodeComm)
          end if
          
          deallocate(pts)
       end do masterLoop

       ! Destroy the ADT Data and allocated values
       do nn=1, nDom
          call destroySerialHex(oBlocks(nn)%ADT)
          deallocate(oBlocks(nn)%xADT, oBlocks(nn)%hexaConn, oBlocks(nn)%qualDonor)
       end do
    end do primalDualLoop
  end subroutine interpolateIntegrationSurfaces

  subroutine commUsetIntegrationSurfaceVars(recvBuffer, varStart, varEnd, comm)

    use constants
    use block, onlY : flowDoms, nDom
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk
    use oversetUtilities, only :fracToWeights
    
    implicit none

    ! Input/Output
    real(kind=realType), dimension(:) :: recvBuffer
    integer(kind=intType), intent(in) :: varStart, varEnd
    type(userSurfCommType) :: comm

    ! Working
    real(kind=realType), dimension(:), allocatable :: sendBuffer
    integer(Kind=intType) :: d1, i1, j1, k1, jj, k, nvar, i, ierr
    real(kind=realType), dimension(8) :: weight

    ! The number of variables we are transferring:
    nVar = varEnd - varStart + 1
    
    ! We assume that the pointers to the realCommVars have already been set. 
    
    allocate(sendBuffer(nVar*comm%nDonor))

    ! First generate the interpolated data necessary
    jj = 0
    donorLoop: do i=1, comm%nDonor
       ! Convert the frac to weights
       call fracToWeights(comm%frac(:, i), weight)
       
       ! Block and indices for easier reading. The +1 is due to the
       ! pointer offset on realCommVars.

       d1 = comm%donorInfo(1, i) ! Block Index
       i1 = comm%donorInfo(2, i)+1 ! donor I index
       j1 = comm%donorInfo(3, i)+1 ! donor J index
       k1 = comm%donorInfo(4, i)+1 ! donor K index
       
       ! We are interpolating nVar variables
       do k=varStart, varEnd
          jj = jj + 1
          if (d1 > 0) then ! Is this pt valid?
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1+1)
          end if
       end do
    end do donorLoop
    
    ! Now we can do an mpi_gatherv to the root proc:
    call mpi_gatherv(sendBuffer, nVar*comm%nDonor, adflow_real, recvBuffer, &
         nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    deallocate(sendBuffer)

  end subroutine commUsetIntegrationSurfaceVars

  subroutine evalUserIntegrationSurfaces(localValues, famList, sps)

    use constants
    use block, onlY : flowDoms, nDom
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk, mynorm2
    use costFunctions, only : nLocalValues, iMassFlow, iMassPs, iMassPTot, iMassTtot
    use flowUtils, only : computePtot, computeTtot
    use sorting, only : bsearchIntegers
    implicit none

    ! Input Parameters
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps

    ! Working parameters
    integer(kind=intType) :: iSurf, i, j, k, jj, ierr, nn, iDim
    real(kind=realType), dimension(:), allocatable :: recvBuffer1, recvBuffer2
    type(userIntSurf), pointer :: surf
    real(kind=realType) :: mReDim, rho, vx, vy, vz, PP, massFLowRate, vn
    real(Kind=realType) :: pTot, Ttot, mass_Ptot, mass_Ttot, mass_ps, massFlowRateLocal
    real(kind=realType), dimension(3) :: v1, v2, x1, x2, x3, x4, n
    logical :: valid
    logical, dimension(:), allocatable :: ptValid
    ! Set the pointers for the required communication variables: rho,
    ! vx, vy, vz, P and the x-coordinates
    
    domainLoop:do nn=1, nDom
       flowDoms(nn, 1, 1)%realCommVars(1)%var => flowDoms(nn, 1, 1)%w(:, :, :, 1)
       flowDoms(nn, 1, 1)%realCommVars(2)%var => flowDoms(nn, 1, 1)%w(:, :, :, 2)
       flowDoms(nn, 1, 1)%realCommVars(3)%var => flowDoms(nn, 1, 1)%w(:, :, :, 3)
       flowDoms(nn, 1, 1)%realCommVars(4)%var => flowDoms(nn, 1, 1)%w(:, :, :, 4)
       flowDoms(nn, 1, 1)%realCommVars(5)%var => flowDoms(nn, 1, 1)%P(:, :, :)
       flowDoms(nn, 1, 1)%realCommVars(6)%var => flowDoms(nn, 1, 1)%x(:, :, :, 1)
       flowDoms(nn, 1, 1)%realCommVars(7)%var => flowDoms(nn, 1, 1)%x(:, :, :, 2)
       flowDoms(nn, 1, 1)%realCommVars(8)%var => flowDoms(nn, 1, 1)%x(:, :, :, 3)
    end do domainLoop

    masterLoop: do iSurf=1, nUserIntSurfs

       ! Pointer for easier reading
       surf => userIntSurfs(iSurf)
       
       ! Do we need to include this surface?
       famInclude: if (bsearchIntegers(surf%famID, famList) > 0) then
          
          ! Communicate the face values and the nodal values
          if (myid == 0) then 
             allocate(recvBuffer1(5*size(surf%conn, 2)))
             allocate(recvBuffer2(3*size(surf%pts, 2)))
          else
             allocate(recvBuffer1(0))
             allocate(recvBuffer2(0))
          end if
          
          call commUsetIntegrationSurfaceVars(recvBuffer1, 1, 5, surf%faceComm)
          call commUsetIntegrationSurfaceVars(recvBuffer2, 6, 8, surf%nodeComm)
          
          ! *Finally* we can do the actual integrations
          if (myid == 0)  then 
             allocate(ptValid(size(surf%pts, 2)))
             
             ! Before we do the integration, update the pts array from
             ! the values in recvBuffer2
             do i=1, size(surf%pts, 2)
                j = surf%nodeComm%inv(i)
                surf%pts(:, j) = recvBuffer2(3*i-2:3*i)
                ptValid(j) = surf%nodeComm%valid(i)
             end do
             
             massFlowRate = zero
             mass_Ptot = zero
             mass_Ttot = zero
             mass_Ps = zero
             
             mReDim = sqrt(pRef*rhoRef)
             
             elemLoop: do i=1, size(surf%conn, 2)
                
                ! First check if this cell is valid to integrate. Both
                ! the cell center *and* all nodes must be valid
                
                ! j is now the element index into the original conn array. 
                j = surf%faceComm%Inv(i)
                
                valid = surf%faceComm%valid(i) .and. &
                     ptValid(surf%conn(1, j)) .and. ptValid(surf%conn(2, j)) .and. &
                     ptValid(surf%conn(3, j)) .and. ptValid(surf%conn(4, j))
                if (valid) then 
                   ! Extract out the interpolated quantities
                   rho = recvBuffer1(5*i-4)
                   vx  = recvBuffer1(5*i-3)
                   vy  = recvBuffer1(5*i-2)
                   vz  = recvBuffer1(5*i-1)
                   PP  = recvBuffer1(5*i  )
                   
                   call computePtot(rho, vx, vy, vz, pp, Ptot)
                   call computeTtot(rho, vx, vy, vz, pp, Ttot)
                   
                   ! Get the coordinates of the quad
                   x1 = surf%pts(:, surf%conn(1, j))
                   x2 = surf%pts(:, surf%conn(2, j))
                   x3 = surf%pts(:, surf%conn(3, j))
                   x4 = surf%pts(:, surf%conn(4, j))
                   
                   ! Diagonal vectors
                   v1 = x3 - x1
                   v2 = x4 - x2
                   
                   ! Take the cross product of the two diagaonal vectors.
                   n(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
                   n(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
                   n(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))
                
                   ! Normal face velocity
                   vn = vx*n(1) + vy*n(2) + vz*n(3)
                   
                   ! Finally the mass flow rate
                   massFlowRateLocal = rho*vn*mReDim
                   massFlowRate = massFlowRate + massFlowRateLocal
                   mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
                   mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
                   mass_Ps = mass_Ps + pp*massFlowRateLocal*Pref
                end if
                
             end do elemLoop
             deallocate(ptValid)
          end if

          ! Accumulate the final values
          localValues(iMassFlow) = localValues(iMassFlow) + massFlowRate
          localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
          localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
          localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps
                    
          deallocate(recvBuffer1, recvBuffer2)


       end if famInclude
    end do masterLoop

  end subroutine evalUserIntegrationSurfaces


#endif
end module surfaceIntegrations
