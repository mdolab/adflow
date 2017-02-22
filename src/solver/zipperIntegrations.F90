module zipperIntegrations
contains

  subroutine flowIntegrationZipper(isInflow, conn, fams, vars, localValues, famList, sps, &
       withGathered, funcValues, ptValid)

    ! Integrate over the trianges for the inflow/outflow conditions. 

    use constants
    use blockPointers, only : BCType
    use sorting, only : famInList
    use flowVarRefState, only : pRef, pInf, rhoRef, pRef, timeRef, LRef, TRef, rGas, uRef, uInf
    use inputPhysics, only : pointRef, flowType, velDirFreeStream, alpha, beta, liftIndex
    use flowUtils, only : computePtot, computeTtot, getDirVector
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : mynorm2, cross_prod
    implicit none

    ! Input/output Variables
    logical, intent(in) :: isInflow
    integer(kind=intType), intent(in), dimension(:, :) :: conn
    integer(kind=intType), intent(in), dimension(:) :: fams
    real(kind=realType), dimension(:, :), intent(in) :: vars
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    logical, intent(in) :: withGathered
    real(kind=realType), optional, dimension(:), intent(in) :: funcValues
    logical(kind=intType), dimension(:), optional, intent(in) :: ptValid
    ! Working variables
    integer(kind=intType) :: i, j
    real(kind=realType) :: sF, vmag, vnm, vxm, vym, vzm, Fx, Fy, Fz, cellArea, u, v, w, vnmFreeStreamRef
    real(kind=realType), dimension(3) :: Fp, Mp, FMom, MMom, refPoint, ss, x1, x2, x3, & 
      norm, VcoordRef, VFreestreamRef, sFaceFreestreamRef, normFreeStreamRef
    real(kind=realType) :: pm, Ptot, Ttot, rhom, gammam, MNm, massFlowRateLocal
    real(kind=realType) ::  massFlowRate, mass_Ptot, mass_Ttot, mass_Ps, mass_MN
    real(kind=realType) ::  mReDim, pk, sigma_MN, sigma_Ptot
    real(kind=realType) :: internalFlowFact, inflowFact, xc, yc, zc, mx, my, mz
    real(kind=realType) :: edotA, edotV, edotP

    logical :: triIsValid

    mReDim = sqrt(pRef*rhoRef)
    Fp = zero
    Mp = zero
    FMom = zero
    MMom = zero

    massFlowRate = zero
    mass_Ptot = zero
    mass_Ttot = zero
    mass_Ps = zero

    edotA = zero
    edotV = zero
    edotP = zero

    sigma_Mn = zero
    sigma_Ptot = zero
    pk = zero

    refPoint(1) = LRef*pointRef(1)
    refPoint(2) = LRef*pointRef(2)
    refPoint(3) = LRef*pointRef(3)

    internalFlowFact = one
    if (flowType == internalFlow) then 
       internalFlowFact = -one
    end if

    inFlowFact = one
    if (isInflow) then 
       inflowFact=-one
    end if

    !$AD II-LOOP
    do i=1, size(conn, 2)
       if (famInList(fams(i), famList)) then 

          ! If the ptValid list is given, check if we should integrate
          ! this triangle.
          triIsValid = .True. 
          if (present(ptValid)) then 
             ! Check if each of the three nodes are valid
             if ((ptValid(conn(1, i)) .eqv. .False.) .or. & 
                 (ptValid(conn(2, i)) .eqv. .False.) .or. &
                 (ptValid(conn(3, i)) .eqv. .False.)) then 
                triIsValid = .False.
             end if
          end if

          validTrianlge: if (triIsValid) then 

             ! Compute the averaged values for this triangle
             vxm = zero; vym = zero; vzm = zero; rhom = zero; pm = zero; MNm = zero; gammam = zero;
             sF = zero
             do j=1,3
                rhom = rhom + vars(conn(j, i), iRho)
                vxm = vxm + vars(conn(j, i), iVx)
                vym = vym + vars(conn(j, i), iVy)
                vzm = vzm + vars(conn(j, i), iVz)
                pm = pm + vars(conn(j, i), iRhoE)
                gammam = gammam + vars(conn(j, i), iZippFlowGamma)
                sF = sF + vars(conn(j, i), iZippFlowSFace)
             end do
             
             ! Divide by 3 due to the summation above:
             rhom = third*rhom
             vxm = third*vxm
             vym = third*vym
             vzm = third*vzm
             pm = third*pm
             gammam = third*gammam
             sF = third*sF
             
             ! Get the nodes of triangle.
             x1 = vars(conn(1, i), iZippFLowX:iZippFlowZ)
             x2 = vars(conn(2, i), iZippFlowX:iZippFlowZ)
             x3 = vars(conn(3, i), iZippFlowX:iZippFlowZ)
             call cross_prod(x2-x1, x3-x1, norm)
             ss = half*norm

             call computePtot(rhom, vxm, vym, vzm, pm, Ptot)
             call computeTtot(rhom, vxm, vym, vzm, pm, Ttot)
             
             vnm = vxm*ss(1) + vym*ss(2) + vzm*ss(3)  - sF
             
             vmag = sqrt((vxm**2 + vym**2 + vzm**2)) - sF
             ! a = sqrt(gamma*p/rho); sqrt(v**2/a**2)
             MNm = vmag/sqrt(gammam*pm/rhom)

             cellArea = sqrt(ss(1)**2 + ss(2)**2 + ss(3)**2)

             massFlowRateLocal = rhom*vnm*mReDim
             
             if (withGathered) then 
                
                sigma_Mn = sigma_Mn  + massFlowRateLocal*(MNm - funcValues(costFuncMavgMN))**2
                Ptot = Ptot * pRef
                sigma_Ptot = sigma_Ptot + massFlowRateLocal*(Ptot - funcValues(costFuncMavgPtot))**2
                
             else 
                massFlowRate = massFlowRate + massFlowRateLocal
                
                pk = pk + ((pm-pInf) + half*rhom*(vmag**2 - uInf**2)) * vnm * pRef * uRef * internalFlowFact

                ! computes the normalized vector maped into the freestream direction, so we multiply by the magnitude after
                VcoordRef(1) = vxm
                VcoordRef(2) = vym
                VcoordRef(3) = vzm

                call getDirVector(VcoordRef, -alpha, -beta, VFreestreamRef, liftIndex)
                VFreestreamRef = VFreestreamRef * vmag

                !project the face normal into the freestream velocity and scale by the face
                call getDirVector(ss, -alpha, -beta, normFreeStreamRef, liftIndex)
                ! note that normFreeStreamRef is of magnitude 1 now
                sFaceFreestreamRef = normFreeStreamRef * sF

                ! compute the pertubations of the flow from the free-stream velocity
                u = VFreestreamRef(1) - sFaceFreestreamRef(1) - uInf
                v = VFreestreamRef(2) - sFaceFreestreamRef(2)
                w = VFreestreamRef(3) - sFaceFreestreamRef(3)

                vnmFreeStreamRef =  (u+uInf)*normFreeStreamRef(1) + v*normFreeStreamRef(2) + w*normFreeStreamRef(3)
                vnmFreeStreamRef = vnmFreeStreamRef * cellArea


                edotA = edotA + half * rhom*u**2 * vnmFreeStreamRef * pref*uRef * internalFlowFact
                edotV = edotV + half * rhom*(v**2+w**2) * vnmFreeStreamRef * pref*uRef * internalFlowFact
                edotP = edotP + (pm-pInf) * (vnm - uInf*normFreeStreamRef(1)*cellArea) * pref*uRef * internalFlowFact
                
                pm = pm*pRef
                
                mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
                mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
                mass_Ps = mass_Ps + pm*massFlowRateLocal
                mass_MN = mass_MN + MNm*massFlowRateLocal
                
                ! Compute the average cell center. 
                xc = zero
                yc = zero
                zc = zero
                do j=1,3
                   xc = xc + (vars(conn(1, i), iZippFlowX)) 
                   yc = yc + (vars(conn(2, i), iZippFlowY)) 
                   zc = zc + (vars(conn(3, i), iZippFlowZ)) 
                end do
                
                ! Finish average for cell center
                xc = third*xc
                yc = third*yc
                zc = third*zc
                
                xc = xc - refPoint(1)
                yc = yc - refPoint(2)
                zc = zc - refPoint(3)
                
                pm = -(pm-pInf*pRef)
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
                ss = ss/cellArea
                massFlowRateLocal = massFlowRateLocal/timeRef*internalFlowFact*inflowFact
                
                fx = massFlowRateLocal*ss(1) * vxm/timeRef
                fy = massFlowRateLocal*ss(2) * vym/timeRef
                fz = massFlowRateLocal*ss(3) * vzm/timeRef
                
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
          end if validTrianlge
       end if
    enddo

    if (withGathered) then 
        localValues(isigmaMN) = localValues(isigmaMN) + sigma_Mn
        localValues(isigmaPtot) = localValues(isigmaPtot) + sigma_Ptot
     else
       ! Increment the local values array with what we computed here
       localValues(iMassFlow) = localValues(iMassFlow) + massFlowRate
       localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
       localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
       localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps
       localValues(iMassMN)   = localValues(iMassMN)   + mass_MN
       localValues(iPk)   = localValues(iPk)   + pk
       localValues(iEdotA) = localValues(iEdotA) + edotA
       localValues(iEdotV) = localValues(iEdotV) + edotV
       localValues(iEdotP) = localValues(iEdotP) + edotP
       localValues(iEdot)   = localValues(iEdot)   + edotA + edotV + edotP
       localValues(iFp:iFp+2)   = localValues(iFp:iFp+2) + Fp
       localValues(iFlowFm:iFlowFm+2)   = localValues(iFlowFm:iFlowFm+2) + FMom
       localValues(iFlowMp:iFlowMp+2)   = localValues(iFlowMp:iFlowMp+2) + Mp
       localValues(iFlowMm:iFlowMm+2)   = localValues(iFlowMm:iFlowMm+2) + MMom
    end if

  end subroutine flowIntegrationZipper

  subroutine wallIntegrationZipper(conn, fams, vars, localValues, famList, sps)

    use constants
    use sorting, only : famInList
    use flowVarRefState, only : LRef
    use inputPhysics, only : pointRef
    use utils, only : mynorm2, cross_prod
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in), dimension(:, :) :: conn
    integer(kind=intType), intent(in), dimension(:) :: fams
    real(kind=realType), intent(in), dimension(:, :) :: vars
    real(kind=realType), intent(inout) :: localValues(nLocalValues)
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps

    ! Working
    real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv


    integer(kind=intType) :: i, j
    real(kind=realType), dimension(3) :: ss, norm, refPoint
    real(kind=realType), dimension(3) :: p1, p2, p3, v1, v2, v3, x1, x2, x3
    real(kind=realType) :: fact, triArea, fx, fy, fz, mx, my, mz, xc, yc, zc
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
    do i=1, size(conn, 2)
       if (famInList(fams(i), famList)) then 

          ! Get the nodes of triangle. 
          x1 = vars(conn(1, i), iZippWallX:iZIppWallZ)
          x2 = vars(conn(2, i), iZippWallX:iZIppWallZ)
          x3 = vars(conn(3, i), iZippWallX:iZIppWallZ)
          call cross_prod(x2-x1, x3-x1, norm)
          ss = half*norm
          ! The third here is to account for the summation of P1, p2
          ! and P3
          triArea = mynorm2(ss)*third
       
          ! Compute the average cell center. 
          xc = third*(x1(1)+x2(1)+x3(1))
          yc = third*(x1(2)+x2(2)+x3(2))
          zc = third*(x1(3)+x2(3)+x3(3))

          xc = xc - refPoint(1)
          yc = yc - refPoint(2)
          zc = zc - refPoint(3)

          ! Update the pressure force and moment coefficients.
          p1 = vars(conn(1, i), iZippWallTpx:iZippWallTpz)
          p2 = vars(conn(2, i), iZippWallTpx:iZippWallTpz)
          p3 = vars(conn(3, i), iZippWallTpx:iZippWallTpz)

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
          v1 = vars(conn(1, i), iZippWallTvx:iZippWallTvz)
          v2 = vars(conn(2, i), iZippWallTvx:iZippWallTvz)
          v3 = vars(conn(3, i), iZippWallTvx:iZippWallTvz)

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

#ifndef USE_TAPENADE

  subroutine integrateZippers(localValues, famList, sps, withGathered, funcValues)
    !--------------------------------------------------------------
    ! Manual Differentiation Warning: Modifying this routine requires
    ! modifying the hand-written forward and reverse routines. 
    ! --------------------------------------------------------------

    ! Integrate over the triangles formed by the zipper mesh. This
    ! will perform both all necesasry zipper integrations. Currently
    ! this includes the wall force integrations as well as the
    ! flow-though surface integration. 

    use constants
    use overset, only : zipperMeshes, zipperMesh
    use haloExchange, only : wallIntegrationZipperComm, flowIntegrationZipperComm
    implicit none

    ! Input Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), dimension(:, :), allocatable :: vars
    type(zipperMesh), pointer :: zipper
    logical, intent(in) :: withGathered
    real(kind=realType), intent(in), dimension(:) :: funcValues

    ! Determine if we have a wall Zipper:
    zipper => zipperMeshes(iBCGroupWalls)
    if (zipper%allocated .and. (withGathered .eqv. .False.)) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippWallComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call wallIntegrationZipperComm(vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call wallIntegrationZipper(zipper%conn, zipper%fam, vars, localValues, famList, sps)

       ! Cleanup vars
       deallocate(vars)
    end if

    zipper => zipperMeshes(iBCGroupInflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call flowIntegrationZipperComm(.true., vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call flowIntegrationZipper(.true., zipper%conn, zipper%fam, vars, &
            localValues, famList, sps, withGathered, funcValues)

       ! Cleanup vars
       deallocate(vars)
    end if

    zipper => zipperMeshes(iBCGroupOutflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call flowIntegrationZipperComm(.false., vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call flowIntegrationZipper(.false., zipper%conn, zipper%fam, vars, &
            localValues, famList, sps, withGathered, funcValues)

       ! Cleanup vars
       deallocate(vars)
    end if
  end subroutine integrateZippers


#ifndef USE_COMPLEX
  subroutine integrateZippers_d(localValues, localValuesd, famList, sps, withGathered, &
       funcValues, funcValuesd)
    !------------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    ! Forward mode linearization of the zipper integration. 

    use constants
    use overset, only : zipperMeshes, zipperMesh
    use haloExchange, only : wallIntegrationZipperComm_d, flowIntegrationZipperComm_d
    use zipperIntegrations_d, only : flowIntegrationZipper_d, wallIntegrationZipper_d
    implicit none

    ! Input Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), dimension(:, :), allocatable :: vars, varsd
    type(zipperMesh), pointer :: zipper
    logical, optional, intent(in) :: withGathered
    real(kind=realType), intent(in), dimension(:) :: funcValues, funcValuesd

    zipper => zipperMeshes(iBCGroupWalls)
    if (zipper%allocated .and. (withGathered .eqv. .False.)) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippWallComm), varsd(size(zipper%indices), nZippWallComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call wallIntegrationZipperComm_d(vars, varsd, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call wallIntegrationZipper_d(zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if

    zipper => zipperMeshes(iBCGroupInflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm), varsd(size(zipper%indices), nZippFlowComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call flowIntegrationZipperComm_d(.true., vars, varsd, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call flowIntegrationZipper_d(.true., zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps, withGathered, funcValues, funcValuesd)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if

    zipper => zipperMeshes(iBCGroupOutflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm), varsd(size(zipper%indices), nZippFlowComm))

       ! Gather up the required variables in "vars" on the root
       ! proc. This routine is differientated by hand. 
       call flowIntegrationZipperComm_d(.false., vars, varsd, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       call flowIntegrationZipper_d(.false., zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps, withGathered, funcValues, funcValuesd)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if
  end subroutine integrateZippers_d

  subroutine integrateZippers_b(localValues, localValuesd, famList, sps, withGathered, &
       funcValues, funcValuesd)
    !------------------------------------------------------------------------
    ! Manual Differentiation Warning: This routine is differentiated by hand.
    ! -----------------------------------------------------------------------

    ! Reverse mode linearization of the zipper integration. 

    use constants
    use overset, only : zipperMeshes, zipperMesh
    use haloExchange, only : wallIntegrationZipperComm_b, flowIntegrationZipperComm_b, &
         wallIntegrationZipperComm, flowIntegrationZipperComm
    use zipperIntegrations_b, only : wallIntegrationZipper_b, flowIntegrationZipper_b
    implicit none

    ! Input Variables
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), dimension(:, :), allocatable :: vars, varsd
    type(zipperMesh), pointer :: zipper
    logical,  intent(in) :: withGathered
    real(kind=realType), dimension(:) :: funcValues, funcValuesd

    ! Determine if we have a wall Zipper:
    zipper => zipperMeshes(iBCGroupWalls)
    if (zipper%allocated .and. (withGathered .eqv. .False.)) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippWallComm), varsd(size(zipper%indices), nZippWallComm))

       ! Set the forward variables
       call wallIntegrationZipperComm(vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       varsd = zero
       call wallIntegrationZipper_b(zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps)

       ! Scatter (becuase we are reverse) the values from the root
       ! back out to all necessary procs.  This routine is
       ! differientated by hand.
       call wallIntegrationZipperComm_b(vars, varsd, sps)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if

    zipper => zipperMeshes(iBCGroupInflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm), varsd(size(zipper%indices), nZippFlowComm))

       ! Set the forward variables
       call flowIntegrationZipperComm(.true., vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       varsd = zero
       call flowIntegrationZipper_b(.true., zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps, withGathered, funcValues, funcValuesd)

       ! Scatter (becuase we are reverse) the values from the root
       ! back out to all necessary procs.  This routine is
       ! differientated by hand.
       call flowIntegrationZipperComm_b(.true., vars, varsd, sps)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if

    zipper => zipperMeshes(iBCGroupOutflow)
    ! Determine if we have a flowthrough Zipper:
    if (zipper%allocated) then 

       ! Allocate space necessary to store variables. Only non-zero on
       ! root proc. 
       allocate(vars(size(zipper%indices), nZippFlowComm), varsd(size(zipper%indices), nZippFlowComm))

       ! Set the forward variables
       call flowIntegrationZipperComm(.false., vars, sps)

       ! Perform actual integration. Tapenade ADs this routine.
       varsd = zero
       call flowIntegrationZipper_b(.false., zipper%conn, zipper%fam, vars, varsd, &
            localValues, localValuesd, famList, sps, withGathered, funcValues, funcValuesd)

       ! Scatter (becuase we are reverse) the values from the root
       ! back out to all necessary procs.  This routine is
       ! differientated by hand.
       call flowIntegrationZipperComm_b(.false., vars, varsd, sps)

       ! Cleanup vars
       deallocate(vars, varsd)
    end if
  end subroutine integrateZippers_b
#endif
#endif

end module zipperIntegrations
