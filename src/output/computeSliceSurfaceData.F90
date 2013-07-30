subroutine computeSliceSurfaceData(sps, nFields)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * The purpose of this routine is to compute all the cell         *
  !      * centered data on all wall surfaces that will be communicated   *
  !      * and used for the user-defined slices. All the variables        *
  !      * that are used for the surface solution file are included in    *
  !      * addition to the three components of the pressure and viscous   *
  !      * forces. On exit, the localValues() array is filled and can     *
  !      * be communicated.                                               *
  !      *                                                                *
  !      ******************************************************************
  !

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  use liftDistributionData
  use cgnsNames
  implicit none

  ! Input Params
  integer(kind=intType), intent(in) :: sps

  ! Output Params
  integer(kind=intType), intent(out) :: nFields

  ! Working parameters
  integer(kind=intType) :: mm, nn, i, j, jj
  integer(kind=intType) :: iFace, nSolVar, iFaceStart, iSolVar, ivar
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yplusmax, qf(3)
  character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames

  integer(kind=intType), dimension(:,:), pointer :: viscPointer
  integer(kind=intType), dimension(:,:), pointer :: iblank2
       
  real(kind=realType) :: fact, gm1, ptotInf, ptot, psurf, rsurf
  real(kind=realType) :: usurf, vsurf, wsurf, m2surf, musurf
  real(kind=realType) :: fx, fy, fz, fn, a2Tot, a2, qw
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz
  real(kind=realType) :: pm1, scaleDim, a
  real(kind=realType), dimension(3) :: norm

  real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
  real(kind=realType), dimension(:,:,:), pointer :: ss1, ss2, ss
  real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
  
  real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
  real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall
  
  ! We will do the pressure and viscous forces first:
  iFace = 0
  domains1: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax)

     ! Loop over the number of boundary subfaces of this block.
     bocos1: do mm=1,nBocos

        ! Check for wall
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! Face Loop:
           do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                 iFace = iFace + 1
                 localData(ifxp, iFace) = bcData(mm)%fP(i, j, 1)
                 localData(ifyp, iFace) = bcData(mm)%fP(i, j, 2)
                 localData(ifzp, iFace) = bcData(mm)%fP(i, j, 3)

                 localData(ifxv, iFace) = bcData(mm)%fV(i, j, 1)
                 localData(ifyv, iFace) = bcData(mm)%fV(i, j, 2)
                 localData(ifzv, iFace) = bcData(mm)%fV(i, j, 3)
              end do
           end do
        end if
     end do bocos1
  end do domains1

  ! Get the number and names of surface solution variables to be
  ! included on the 
  call numberOfSurfSolVariables(nSolVar)
  allocate(solNames(nSolVar))
  call surfSolNames(solNames)
  nFields = 6 + nSolVar

  ! Next we loop back over the surfaces, but this time we
  ! incrementally add any of the surface solution data to the 
  iFace = 0
  domains2: do nn=1,nDom
     call setPointers(nn,1_intType,sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos2: do mm=1,nBocos

        ! Check for wall
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! Set pointers for the correct block face
           select case (bcFaceID(mm))
           case (iMin)
              ww1    => w(1,1:,1:,:);   ww2    => w(2,1:,1:,:)
              pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)
              ss => si(1,:,:,:) ; fact = -one
              gamma1 => gamma(1,1:,1:); gamma2 => gamma(2,1:,1:)

              if( blockIsMoving)then
                 ss1    => s(1,1:,1:,:);   ss2    => s(2,1:,1:,:)
              end if

              iblank2 => iblank(2,1:,1:)
              viscPointer => viscIminPointer

              if( viscous ) then
                 rlv1 => rlv(1,1:,1:); rlv2 => rlv(2,1:,1:)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)

              !===============================================================

           case (iMax)
              ww1    => w(ie,1:,1:,:);   ww2    => w(il,1:,1:,:)
              ss => si(il,:,:,:) ; fact = one
              pp1    => p(ie,1:,1:);     pp2    => p(il,1:,1:)
              gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)

              if( blockIsMoving)then
                 ss1    => s(ie-1,1:,1:,:);   ss2    => s(ie,1:,1:,:)
              end if

              iblank2 => iblank(il,1:,1:)
              viscPointer => viscImaxPointer

              if( viscous ) then
                 rlv1 => rlv(ie,1:,1:); rlv2 => rlv(il,1:,1:)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)

              !===============================================================

           case (jMin)
              ww1    => w(1:,1,1:,:);   ww2    => w(1:,2,1:,:)
              ss => sj(:,1,:,:) ; fact = -one
              pp1    => p(1:,1,1:);     pp2    => p(1:,2,1:)
              gamma1 => gamma(1:,1,1:); gamma2 => gamma(1:,2,1:)

              if( blockIsMoving)then
                 ss1    => s(1:,1,1:,:);   ss2    => s(1:,2,1:,:)
              end if

              iblank2 => iblank(1:,2,1:)
              viscPointer => viscJminPointer

              if( viscous ) then
                 rlv1 => rlv(1:,1,1:); rlv2 => rlv(1:,2,1:)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)

              !===============================================================

           case (jMax)
              ww1    => w(1:,je,1:,:);   ww2    => w(1:,jl,1:,:)
              ss => sj(:,jl,:,:); fact = one
              pp1    => p(1:,je,1:);     pp2    => p(1:,jl,1:)
              gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)

              if( blockIsMoving)then
                 ss1    => s(1:,je-1,1:,:);   ss2    => s(1:,je,1:,:)
              end if

              iblank2 => iblank(1:,jl,1:)
              viscPointer => viscJmaxPointer

              if( viscous ) then
                 rlv1 => rlv(1:,je,1:); rlv2 => rlv(1:,jl,1:)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)

              !===============================================================

           case (kMin)
              ww1    => w(1:,1:,1,:);   ww2    => w(1:,1:,2,:)
              ss => sk(:,:,1,:);  fact = -one
              pp1    => p(1:,1:,1);     pp2    => p(1:,1:,2)
              gamma1 => gamma(1:,1:,1); gamma2 => gamma(1:,1:,2)

              if( blockIsMoving)then
                 ss1    => s(1:,1:,1,:);   ss2    => s(1:,1:,2,:)
              end if

              iblank2 => iblank(1:,1:,2)
              viscPointer => viscKminPointer

              if( viscous ) then
                 rlv1 => rlv(1:,1:,1); rlv2 => rlv(1:,1:,2)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)

              !===============================================================

           case (kMax)
              ww1    => w(1:,1:,ke,:);   ww2    => w(1:,1:,kl,:)
              ss => sk(:,:,kl,:);  fact = one
              pp1    => p(1:,1:,ke);     pp2    => p(1:,1:,kl)
              gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)

              if( blockIsMoving)then
                 ss1    => s(1:,1:,ke-1,:);   ss2    => s(1:,1:,ke,:)
              end if

              iblank2 => iblank(1:,1:,kl)
              viscPointer => viscKmaxPointer

              if( viscous ) then
                 rlv1 => rlv(1:,1:,ke); rlv2 => rlv(1:,1:,kl)
              endif

              if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)

           end select

           ! Now we have the pointers set. We can loop through and
           ! determine the variable to be written 

           ! Store the index of the current face since we will most
           ! likely going though multiple storage variables
           iFaceStart = iFace

           do iSolVar=1,nSolVar
              ! Reset to current starting face
              iFace = iFaceStart

              ! Index is 6 plus the current index
              iVar = iSolVar + ifzv
              varName: select case (solNames(iSolVar))

              case (cgnsDensity)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,irho) + ww2(i,j,irho))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsPressure)

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(pp1(i,j) + pp2(i,j))
                    enddo
                 enddo
                 !===============================================================

              case (cgnsTemp)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = (pp1(i,j) + pp2(i,j)) &
                            / (RGas*(ww1(i,j,irho) + ww2(i,j,irho)))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsVelx)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsVely)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsVelz)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsRelVelx)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsRelVely)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
                    enddo
                 enddo

                 !===============================================================

              case (cgnsRelVelz)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
                    enddo
                 enddo

                 !================================================================

              case (cgnsCp)

                 ! Factor multiplying p-pInf

                 fact = two/(gammaInf*pInf*MachCoef*MachCoef)
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = fact*(half*(pp1(i,j) + pp2(i,j)) - pInf)
                    enddo
                 enddo

                 !===============================================================

              case (cgnsPtotloss)

                 ! First compute the total pressure of the free stream.

                 call computePtot(rhoInf, uInf, zero, zero, &
                      pInf, ptotInf, 1_intType)
                 ptotInf = one/ptotInf

                 ! Loop over the faces and compute the total pressure loss.
                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                       psurf = half*(pp1(i,j) + pp2(i,j))
                       rsurf = half*(ww1(i,j,irho) + ww2(i,j,irho))
                       usurf = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
                       vsurf = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
                       wsurf = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))

                       call computePtot(rsurf, usurf, vsurf, wsurf, &
                            psurf, ptot, 1_intType)

                       iFace = iFace + 1
                       localData(iVar, iFace) = one - ptot*ptotInf
                    enddo
                 enddo

                 !===============================================================

              case (cgnsMach)

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                       psurf  = half*(pp1(i,j) + pp2(i,j))
                       rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
                       usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
                       vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
                       wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))
                       m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                            / (half*(gamma1(i,j) + gamma2(i,j))*psurf)

                       iFace = iFace + 1
                       localData(iVar, iFace) = sqrt(m2surf)
                    enddo
                 enddo

                 !===============================================================

              case (cgnsRelMach)

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                       psurf  = half*(pp1(i,j) + pp2(i,j))
                       rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
                       usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
                       vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
                       wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
                       m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                            / (half*(gamma1(i,j) + gamma2(i,j))*psurf)

                       iFace = iFace + 1
                       localData(iVar, iFace) = sqrt(m2surf)
                    enddo
                 enddo

                 !        ================================================================

              case (cgnsSkinFmag, cgnsYplus, &
                   cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

                 ! To avoid a lot of code duplication these 5 variables are
                 ! treated together.

                 ! Multiplication factor to obtain the skin friction from
                 ! the wall shear stress.

                 fact = two/(gammaInf*pInf*MachCoef*MachCoef)

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                       ! Determine the viscous subface on which this
                       ! face is located.
                       
                       jj = viscPointer(i,j)

                       ! Store the 6 components of the viscous stress tensor
                       ! a bit easier.

                       tauxx = viscSubface(jj)%tau(i,j,1)
                       tauyy = viscSubface(jj)%tau(i,j,2)
                       tauzz = viscSubface(jj)%tau(i,j,3)
                       tauxy = viscSubface(jj)%tau(i,j,4)
                       tauxz = viscSubface(jj)%tau(i,j,5)
                       tauyz = viscSubface(jj)%tau(i,j,6)

                       ! Compute the "unit" force on this face. The unit normal
                       ! is outward pointing per definition. A minus sign is
                       ! present, because of the definition of the viscous
                       ! stress tensor. 

                       norm(1) = BCData(jj)%norm(i,j,1)
                       norm(2) = BCData(jj)%norm(i,j,2)
                       norm(3) = BCData(jj)%norm(i,j,3)

                       fx = -(tauxx*norm(1) + tauxy*norm(2) + tauxz*norm(3))
                       fy = -(tauxy*norm(1) + tauyy*norm(2) + tauyz*norm(3))
                       fz = -(tauxz*norm(1) + tauyz*norm(2) + tauzz*norm(3))

                       fn = fx*norm(1) + fy*norm(2) + fz*norm(3)

                       fx = fx - fn*norm(1)
                       fy = fy - fn*norm(2)
                       fz = fz - fn*norm(3)

                       ! Determine the variable to be stored and compute it.
                       ! Note that an offset of -1 must be used in dd2Wall,
                       ! because the original array, d2Wall, starts at 2.
                       ! First update the counter nn.

                       iFace = iFace + 1

                       select case (solNames(iSolVar))
                       case (cgnsSkinFmag)
                          localData(iVar, iFace) = fact*sqrt(fx*fx + fy*fy + fz*fz)

                       case (cgnsSkinFx)
                          localData(iVar, iFace) = fact*fx

                       case (cgnsSkinFy)
                          localData(iVar, iFace) = fact*fy

                       case (cgnsSkinFz)
                          localData(iVar, iFace) = fact*fz

                       case (cgnsYplus)
                          rsurf      = half*(ww1(i,j,irho) + ww2(i,j,irho))
                          musurf     = half*(rlv1(i,j)     + rlv2(i,j))
                          localData(iVar, iFace) = sqrt(rsurf*sqrt(fx*fx + fy*fy + fz*fz)) &
                               * dd2Wall(i-1,j-1)/musurf
                       end select

                    enddo
                 enddo
                 
                 !        ================================================================

              case (cgnsStanton)

                 ! Some constants needed to compute the stanton number.

                 gm1   = gammaInf - one
                 a2Tot = gammaInf*pInf*(one + half*gm1*MachCoef*MachCoef) &
                      / rhoInf
                 fact   = MachCoef*sqrt(gammaInf*pInf*rhoInf)/gm1

                 ! Loop over the given range of faces. As the viscous data is
                 ! only present in the owned faces, the values of the halo's
                 ! are set equal to the nearest physical face. Therefore the
                 ! working indices are ii and jj.

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                       ! Determine the viscous subface on which this
                       ! face is located.

                       jj = viscPointer(i,j)

                       ! Compute the heat flux. Multipy with the sign of the
                       ! normal to obtain the correct value.

                       qw = viscSubface(jj)%q(i,j,1)*BCData(jj)%norm(i,j,1) &
                            + viscSubface(jj)%q(i,j,2)*BCData(jj)%norm(i,j,2) &
                            + viscSubface(jj)%q(i,j,3)*BCData(jj)%norm(i,j,3)

                       ! Compute the speed of sound squared at the wall and
                       ! the stanton number, which is stored in buffer.

                       a2 = half*(gamma1(i,j)   + gamma2(i,j)) &
                            *      (pp1(i,j)      + pp2(i,j))    &
                            /      (ww1(i,j,irho) + ww2(i,j,irho))

                       iFace = iFace + 1
                       localData(ivar, iFace) = qw/(fact*(a2Tot-a2))

                    enddo
                 enddo

                 !        ================================================================

              case (cgnsBlank)

                 ! Loop over the given range of faces. Since iblanks are set
                 ! to 2 for boundary conditions and >= 10 for the boundary,
                 ! take the minimum of the value and 1, so that cells with
                 ! valid data always have an iblank of 1.

                 do j= BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                    do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                       iFace = iFace + 1
                       localData(iVar, iFace) = real(min(iblank2(i,j), 1_intType), realType)
                    enddo
                 enddo

              end select varName
           end do ! Solution Variable loop
        end if ! Wall Type
     end do bocos2 ! Boco loop
  end do domains2 ! Domain loop

  ! Destroy the surface variable names
  deallocate(solNames)
end subroutine computeSliceSurfaceData
