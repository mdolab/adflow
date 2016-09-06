module surfaceIntegrations

contains


  subroutine flowProperties(massFlowRate, mass_Ptot, mass_Ttot, mass_Ps)

    use constants
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use bcroutines
    use costFunctions
    use surfaceFamilies
    use sorting, only : bsearchIntegers
    use utils, only : setBCPointers, resetBCPointers
    use flowUtils, only : computePtot, computeTtot
    use BCPointers
    implicit none

    !
    !      Subroutine arguments
    !
    real(kind=realType), intent(out) :: massFlowRate, mass_Ptot, mass_Ttot, mass_Ps
    integer(kind=intType) :: nn, i, j, ii
    real(kind=realType) :: fact
    real(kind=realType) :: sF, vnm, vxm, vym, vzm
    real(kind=realType) :: pm, Ptot, Ttot, rhom, massFlowRateLocal, tmp

    massFlowRate = zero
    mass_Ptot = zero
    mass_Ttot = zero
    mass_Ps = zero
    sF = zero

    bocos: do nn=1,nBocos
       famInclude: if (bsearchIntegers(BCdata(nn)%famID, famGroups, size(famGroups)) > 0) then
          call setBCPointers(nn, .True.)

          select case (BCFaceID(nn))
          case (iMin)
             fact = -one
          case (iMax)
             fact = one
          case (jMin)
             fact = -one
          case (jMax)
             fact = one
          case (kMin)
             fact = -one
          case (kMax)
             fact = one
          end select


          ! Loop over the quadrilateral faces of the subface. Note that
          ! the nodal range of BCData must be used and not the cell
          ! range, because the latter may include the halo's in i and
          ! j-direction. The offset +1 is there, because inBeg and jnBeg
          ! refer to nodal ranges and not to cell ranges. The loop
          ! (without the AD stuff) would look like:
          !
          ! do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
          !    do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

          !$AD II-LOOP
          do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1
             i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
             j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

             if( addGridVelocities ) sF = sFace(i,j)
             vxm = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
             vym = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
             vzm = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
             rhom = half*(ww1(i,j,irho) + ww2(i,j,irho))
             pm = half*(pp1(i,j)+ pp2(i,j))

             vnm = vxm*ssi(i,j,1) + vym*ssi(i,j,2) + vzm*ssi(i,j,3)  - sF

             massFlowRateLocal = rhom*vnm

             !  vn1 = ww1(i,j,ivx)*ssi(i,j,1) + ww1(i,j,ivy)*ssi(i,j,2) &
             ! + ww1(i,j,ivz)*ssi(i,j,3) - sF
             !  vn2 = ww2(i,j,ivx)*ssi(i,j,1) + ww2(i,j,ivy)*ssi(i,j,2) &
             ! + ww2(i,j,ivz)*ssi(i,j,3) - sF

             ! massFlowRateLocal = half*(ww1(i,j,irho)*vn1 &
             !                     + ww2(i,j,irho)*vn2)

             massFlowRate = massFlowRate + massFlowRateLocal

             call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
             call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)

             mass_Ptot = mass_pTot + Ptot * massFlowRateLocal
             mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal
             mass_Ps = mass_Ps + pm*massFlowRateLocal
             print *, nn, i, j, pTot, massFlowRateLocal

          enddo

          massFlowRate = massFlowRate*fact
          mass_Ptot = mass_pTot*fact
          mass_Ttot = mass_Ttot*fact
          mass_Ps = mass_Ps*fact

       end if famInclude
    end do bocos
  end subroutine flowProperties

  subroutine forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, &
       sepSensorAvg, Cavitation)
    !
    !       forcesAndMoments computes the contribution of the block
    !       given by the pointers in blockPointers to the force and
    !       moment coefficients of the geometry. A distinction is made
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
    use inputPhysics
    use bcroutines
    use costFunctions
    use surfaceFamilies
    use sorting, only :bsearchIntegers
    use utils, only : setBCPointers, resetBCPointers
    use BCPointers
    implicit none
    !
    !      Subroutine arguments
    !
    real(kind=realType), dimension(3), intent(out) :: cFp, cFv
    real(kind=realType), dimension(3), intent(out) :: cMp, cMv

    real(kind=realType), intent(out) :: yplusMax, sepSensor
    real(kind=realType), intent(out) :: sepSensorAvg(3), Cavitation
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, i, j, ii, blk

    real(kind=realType) :: pm1, fx, fy, fz, fn, sigma
    real(kind=realType) :: xc, yc, zc, qf(3)
    real(kind=realType) :: fact, rho, mul, yplus, dwall
    real(kind=realType) :: scaleDim, V(3), sensor, sensor1, Cp, tmp, plocal
    real(kind=realType) :: tauXx, tauYy, tauZz
    real(kind=realType) :: tauXy, tauXz, tauYz

    real(kind=realType), dimension(3) :: refPoint
    real(kind=realType) :: mx, my, mz, cellArea
    logical :: viscousSubface

    ! Set the actual scaling factor such that ACTUAL forces are computed
    scaleDim = pRef/pInf

    ! Determine the reference point for the moment computation in
    ! meters.

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    ! Initialize the force and moment coefficients to 0 as well as
    ! yplusMax.

    cFp(1) = zero; cFp(2) = zero; cFp(3) = zero
    cFv(1) = zero; cFv(2) = zero; cFv(3) = zero
    cMp(1) = zero; cMp(2) = zero; cMp(3) = zero
    cMv(1) = zero; cMv(2) = zero; cMv(3) = zero

    yplusMax = zero
    sepSensor = zero
    Cavitation = zero
    sepSensorAvg = zero

    ! Loop over the boundary subfaces of this block.

    bocos: do nn=1,nBocos
       !
       !         Integrate the inviscid contribution over the solid walls,
       !         either inviscid or viscous. The integration is done with
       !         cp. For closed contours this is equal to the integration
       !         of p; for open contours this is not the case anymore.
       !         Question is whether a force for an open contour is
       !         meaningful anyway.
       !

       famInclude: if (bsearchIntegers(BCdata(nn)%famID, &
            famGroups, size(famGroups)) > 0) then

          invForce: if(BCType(nn) == EulerWall .or. &
               BCType(nn) == NSWallAdiabatic .or. &
               BCType(nn) == NSWallIsothermal) then
             ! Subface is a wall. Check if it is a viscous wall.

             viscousSubface = .true.
             if(BCType(nn) == EulerWall) viscousSubface = .false.

             ! Set a bunch of pointers depending on the face id to make
             ! a generic treatment possible. The routine setBcPointers
             ! is not used, because quite a few other ones are needed.
             call setBCPointers(nn, .True.)

             select case (BCFaceID(nn))
             case (iMin)
                fact = -one
             case (iMax)
                fact = one
             case (jMin)
                fact = -one
             case (jMax)
                fact = one
             case (kMin)
                fact = -one
             case (kMax)
                fact = one
             end select

             ! Loop over the quadrilateral faces of the subface. Note that
             ! the nodal range of BCData must be used and not the cell
             ! range, because the latter may include the halo's in i and
             ! j-direction. The offset +1 is there, because inBeg and jnBeg
             ! refer to nodal ranges and not to cell ranges. The loop
             ! (without the AD stuff) would look like:
             !
             ! do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             !    do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

             !$AD II-LOOP
             do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1
                i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
                j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

                ! Compute the average pressure minus 1 and the coordinates
                ! of the centroid of the face relative from from the
                ! moment reference point. Due to the usage of pointers for
                ! the coordinates, whose original array starts at 0, an
                ! offset of 1 must be used. The pressure is multipled by
                ! fact to account for the possibility of an inward or
                ! outward pointing normal.

                pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)*scaleDim

                xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
                     +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
                yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
                     +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
                zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
                     +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

                ! Compute the force components.
                blk = max(BCData(nn)%iblank(i,j), 0)
                fx = pm1*ssi(i,j,1)
                fy = pm1*ssi(i,j,2)
                fz = pm1*ssi(i,j,3)

                ! iBlank forces
                fx = fx*blk
                fy = fy*blk
                fz = fz*blk

                ! Update the inviscid force and moment coefficients.
                cFp(1) = cFp(1) + fx
                cFp(2) = cFp(2) + fy
                cFp(3) = cFp(3) + fz

                mx = yc*fz - zc*fy
                my = zc*fx - xc*fz
                mz = xc*fy - yc*fx

                cMp(1) = cMp(1) + mx
                cMp(2) = cMp(2) + my
                cMp(3) = cMp(3) + mz

                ! Save the face-based forces and area
                bcData(nn)%Fp(i, j, 1) = fx
                bcData(nn)%Fp(i, j, 2) = fy
                bcData(nn)%Fp(i, j, 3) = fz
                cellArea = sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2)
                bcData(nn)%area(i, j) = cellArea

                ! Get normalized surface velocity:
                v(1) = ww2(i, j, ivx)
                v(2) = ww2(i, j, ivy)
                v(3) = ww2(i, j, ivz)
                v = v / (sqrt(v(1)**2 + v(2)**2 + v(3)**2) + 1e-16)

                ! Dot product with free stream
                sensor = -(v(1)*velDirFreeStream(1) + &
                     v(2)*velDirFreeStream(2) + &
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
                tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
                Cp = tmp*(plocal-pinf)
                Sigma = 1.4
                Sensor1 = -Cp - Sigma
                Sensor1 = one/(one+exp(-2*10*Sensor1))
                Sensor1 = Sensor1 * cellArea
                Cavitation = Cavitation + Sensor1
             enddo
             !
             !           Integration of the viscous forces.
             !           Only for viscous boundaries.
             !
             visForce: if( viscousSubface ) then

                ! Initialize dwall for the laminar case and set the pointer
                ! for the unit normals.

                dwall = zero
                ! Replace norm with BCData norm - Peter Lyu
                !norm => BCData(nn)%norm

                ! Loop over the quadrilateral faces of the subface and
                ! compute the viscous contribution to the force and
                ! moment and update the maximum value of y+.

                do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1
                   i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
                   j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

                   ! Store the viscous stress tensor a bit easier.
                   blk = max(BCData(nn)%iblank(i,j), 0)

                   tauXx = viscSubface(nn)%tau(i,j,1)
                   tauYy = viscSubface(nn)%tau(i,j,2)
                   tauZz = viscSubface(nn)%tau(i,j,3)
                   tauXy = viscSubface(nn)%tau(i,j,4)
                   tauXz = viscSubface(nn)%tau(i,j,5)
                   tauYz = viscSubface(nn)%tau(i,j,6)

                   ! Compute the viscous force on the face. A minus sign
                   ! is now present, due to the definition of this force.

                   fx = -fact*(tauXx*ssi(i,j,1) + tauXy*ssi(i,j,2) &
                        +        tauXz*ssi(i,j,3))*scaleDim
                   fy = -fact*(tauXy*ssi(i,j,1) + tauYy*ssi(i,j,2) &
                        +        tauYz*ssi(i,j,3))*scaleDim
                   fz = -fact*(tauXz*ssi(i,j,1) + tauYz*ssi(i,j,2) &
                        +        tauZz*ssi(i,j,3))*scaleDim

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

                   cFv(1) = cFv(1) + fx
                   cFv(2) = cFv(2) + fy
                   cFv(3) = cFv(3) + fz

                   mx = yc*fz - zc*fy
                   my = zc*fx - xc*fz
                   mz = xc*fy - yc*fx

                   cMv(1) = cMv(1) + mx
                   cMv(2) = cMv(2) + my
                   cMv(3) = cMv(3) + mz

                   ! Save the face based forces for the slice operations
                   bcData(nn)%Fv(i, j, 1) = fx
                   bcData(nn)%Fv(i, j, 2) = fy
                   bcData(nn)%Fv(i, j, 3) = fz

                   ! Compute the tangential component of the stress tensor,
                   ! which is needed to monitor y+. The result is stored
                   ! in fx, fy, fz, although it is not really a force.
                   ! As later on only the magnitude of the tangential
                   ! component is important, there is no need to take the
                   ! sign into account (it should be a minus sign).

                   fx = tauXx*BCData(nn)%norm(i,j,1) + tauXy*BCData(nn)%norm(i,j,2) &
                        + tauXz*BCData(nn)%norm(i,j,3)
                   fy = tauXy*BCData(nn)%norm(i,j,1) + tauYy*BCData(nn)%norm(i,j,2) &
                        + tauYz*BCData(nn)%norm(i,j,3)
                   fz = tauXz*BCData(nn)%norm(i,j,1) + tauYz*BCData(nn)%norm(i,j,2) &
                        + tauZz*BCData(nn)%norm(i,j,3)

                   fn = fx*BCData(nn)%norm(i,j,1) + fy*BCData(nn)%norm(i,j,2) + fz*BCData(nn)%norm(i,j,3)

                   fx = fx - fn*BCData(nn)%norm(i,j,1)
                   fy = fy - fn*BCData(nn)%norm(i,j,2)
                   fz = fz - fn*BCData(nn)%norm(i,j,3)

                   ! Compute the local value of y+. Due to the usage
                   ! of pointers there is on offset of -1 in dd2Wall..
#ifndef TAPENADE_REVERSE
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
                bcData(nn)%Fv = zero
             end if visForce

             call resetBCPointers(nn, .True.)

          end if invForce
       else
          ! If it wasn't included, but still a wall...zero
          if(BCType(nn) == EulerWall .or. &
               BCType(nn) == NSWallAdiabatic .or. &
               BCType(nn) == NSWallIsothermal) then
             bcData(nn)%area = zero
             bcData(nn)%Fp = zero
             bcData(nn)%Fv = zero
          end if
       end if famInclude
    end do bocos

    ! Currently the coefficients only contain the surface integral
    ! of the pressure tensor. These values must be scaled to
    ! obtain the correct coefficients.

    fact = two/(gammaInf*pInf*MachCoef*MachCoef &
         *surfaceRef*LRef*LRef*scaleDim)

    cFp(1) = cFp(1)*fact; cFp(2) = cFp(2)*fact; cFp(3) = cFp(3)*fact
    cFv(1) = cFv(1)*fact; cFv(2) = cFv(2)*fact; cFv(3) = cFv(3)*fact

    fact = fact/(lengthRef*LRef)
    cMp(1) = cMp(1)*fact; cMp(2) = cMp(2)*fact; cMp(3) = cMp(3)*fact
    cMv(1) = cMv(1)*fact; cMv(2) = cMv(2)*fact; cMv(3) = cMv(3)*fact

  end subroutine forcesAndMoments

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
    use communication, only : sumb_comm_world, myid
    use inputPhysics
    use iteration
    use costFunctions
    use inputTimeSpectral
    use flowVarRefState
    use overset, only : oversetPresent
    use utils, only : EChk, setPointers
    implicit none

    ! Input/Ouput Variables
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), intent(out), dimension(nCostFunction)::globalCFVals

    !      Local variables.
    integer(kind=intType) :: nn, ierr
    real(kind=realType) :: force(3), cforce(3), Lift, Drag, CL, CD
    real(kind=realType) :: Moment(3),cMoment(3), fact, scaleDim
    real(kind=realType), dimension(3) :: cFpLocal, cFvLocal, cMpLocal, cMvLocal, sepSensorAvgLocal
    real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yPlusMax, sepSensorLocal, cavitationLocal
    real(Kind=realType) :: sepSensor, sepSensorAvg(3), Cavitation
    real(Kind=realType) :: massFlowRate, mass_Ptot, mass_Ttot, mass_Ps
    real(kind=realType), dimension(nCostFunction)::localCFVals


    cFpLocal = zero
    cFvLocal = zero
    cMpLocal = zero
    cMvLocal = zero
    sepSensorLocal = zero
    sepSensorAvgLocal = zero
    cavitationLocal = zero
    domains: do nn=1,nDom
       call setPointers(nn,1_intType,sps)

       call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, &
            sepSensorAvg, Cavitation)

       cFpLocal = cFpLocal + cFp
       cFvLocal = cFvLocal + cFv
       cMpLocal = cMpLocal + cMp
       cMvLocal = cMvLocal + cMv
       sepSensorLocal = sepSensorLocal + sepSensor
       sepSensorAvgLocal = sepSensorAvgLocal + sepSensorAvg
       cavitationLocal = cavitationLocal + cavitation

       call flowProperties(massFlowRate, mass_Ptot, mass_Ttot, mass_Ps)
       ! Compute the dimensional mass flow
    end do domains

    if (oversetPresent) then
       call forcesAndMomentsZipper(cfp, cfv, cmp, cmv, sps)

       ! Add the extra zipper component on the root proc
       if (myid == 0) then
          cFpLocal = cFpLocal + cFp
          cFvLocal = cFvLocal + cFv
          cMpLocal = cMpLocal + cMp
          cMvLocal = cMvLocal + cMv
       end if
    end if

    ! Now compute the other variables
    cFp = cFpLocal
    cFv = cFvLocal
    cMp = cMpLocal
    cMv = cMvLocal
    sepSensor = sepSensorLocal
    cavitation = cavitationLocal
    sepSensorAvg = sepSensorAvgLocal


    scaleDim = pRef/pInf

    ! Sum pressure and viscous contributions
    cForce = cFp + cFv
    cMoment = cMp + cMv

    ! Get Lift coef and Drag coef
    CD =  cForce(1)*dragDirection(1) &
         + cForce(2)*dragDirection(2) &
         + cForce(3)*dragDirection(3)

    CL =  cForce(1)*liftDirection(1) &
         + cForce(2)*liftDirection(2) &
         + cForce(3)*liftDirection(3)

    ! Divide by fact to get the forces, Lift and Drag back
    fact = two/(gammaInf*pInf*MachCoef*MachCoef &
         *surfaceRef*LRef*LRef*scaleDim)
    Force = cForce / fact
    Lift  = CL / fact
    Drag  = CD / fact

    ! Moment factor has an extra lengthRef
    fact = fact/(lengthRef*LRef)
    Moment = cMoment / fact

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
    localCFVals(costFuncSepSensor) = sepSensor
    localCFVals(costFuncSepSensorAvgX) = sepSensorAvg(1)
    localCFVals(costFuncSepSensorAvgY) = sepSensorAvg(2)
    localCFVals(costFuncSepSensorAvgZ) = sepSensorAvg(3)
    localCFVals(costFuncCavitation) = Cavitation
    localCFVals(costFuncMdot) = localCFVals(costFuncMdot) + massFlowRate
    localCFVals(costFuncMavgPtot) = localCFVals(costFuncMavgPtot) + mass_Ptot
    localCFVals(costFuncMavgTtot) = localCFVals(costFuncMavgTtot) + mass_Ttot
    localCFVals(costFuncMavgPs) = localCFVals(costFuncMavgPs) + mass_Ps

    ! Now we will mpi_allReduce them into globalCFVals
    call mpi_allreduce(localCFVals, globalCFVals, nCostFunction, sumb_real, &
         mpi_sum, SUmb_comm_world, ierr)

    globalCFVals(costFuncMavgPtot) = globalCFVals(costFuncMavgPtot)/globalCFVals(costFuncMdot)*pRef
    globalCFVals(costFuncMavgTtot) = globalCFVals(costFuncMavgTtot)/globalCFVals(costFuncMdot)*tRef
    globalCFVals(costFuncMavgPs) = globalCFVals(costFuncMavgPs)/globalCFVals(costFuncMdot)*pRef
    globalCFVals(costFuncMdot) = globalCFVals(costFuncMdot)*sqrt(pRef*rhoRef)

  end subroutine computeAeroCoef

  subroutine getSolution(sps)
    use constants
    use costFunctions
    use inputTSStabDeriv, only : TSSTability
    use inputTimeSpectral , only : nTimeIntervalsSpectral
    use communication, only : sumb_comm_world
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
    integer(kind=intType) :: i, liftIndex


    call getDirAngle(velDirFreestream, LiftDirection,&
         liftIndex, alpha, beta)

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
       call computeRootBendingMoment(cf, cm, liftIndex, bendingMoment)
       bendingsum = bendingsum+bendingMoment
    end do

    call computeAeroCoef(globalCFVals,sps)
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

       call computeTSDerivatives(force, moment, liftIndex, coef0, dcdalpha, &
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

  subroutine forcesAndMomentsZipper(cFp, cFv, cMp, cMv, sps)

    use communication
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use costFunctions
    use overset, only : nodeZipperScatter, globalNodalVec, localZipperNodes, localZipperTp, localZipperTv
    use inputTimeSpectral
    use inputIteration
    use utils, only : EChk, setPointers, myNorm2
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    ! Input/Output
    real(kind=realType), dimension(3), intent(out) :: cFp, cFv
    real(kind=realType), dimension(3), intent(out) :: cMp, cMv
    integer(kind=intType), intent(in) :: sps

    ! Working
    integer(kind=intType) :: i, j, k, l, ierr, nn, mm, gind, ind, iVar
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, rowStart, rowEnd
    real(kind=realType), dimension(:, :, :), pointer :: xx
    real(kind=realType), dimension(:), pointer :: xPtr, pPtr, vPtr, localPtr
    integer(kind=intType), dimension(:, :), pointer :: gnp

    real(kind=realType), dimension(3) :: x1, x2, x3, xc, ss, norm, refPoint
    real(kind=realType), dimension(3) :: pp1, pp2, pp3, presForce
    real(kind=realType), dimension(3) :: vv1, vv2, vv3, viscForce
    real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv, MTmp
    real(kind=realType) :: fact, scaleDim, triArea

    ! PETsc
    Vec, pointer :: localVec

    ! Set the actual scaling factor such that ACTUAL forces are computed
    scaleDim = pRef/pInf

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

    ! Puck out pointers for the nodes and tractions
    cFp = zero
    cFv = zero
    cMp = zero
    cMv = zero
    if (myid == 0) then

       Fp = zero
       Fv = zero
       Mp = zero
       Mv = zero

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

       fact = two/(gammaInf*pInf*MachCoef*MachCoef &
            *surfaceRef*LRef*LRef*scaleDim)

       ! Convert the values to coefficients
       cFp(:) =  fact*Fp
       cFv(:) =  fact*Fv

       fact = fact/(lengthRef*LRef)
       cMp(:) = Mp(:)*fact
       cMv(:) = Mv(:)*fact

    end if

  end subroutine forcesAndMomentsZipper
#endif
end module surfaceIntegrations
