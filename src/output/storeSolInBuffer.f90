!
!      ******************************************************************
!      *                                                                *
!      * File:          storeSolInBuffer.f90                            *
!      * Author:        Edwin van der Weide, Georgi Kalitzin,           *
!      *                Steve Repsher                                   *
!      * Starting date: 04-14-2003                                      *
!      * Last modified: 10-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine storeSolInBuffer(buffer, copyInBuffer, solName, &
                                   iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
!
!      ******************************************************************
!      *                                                                *
!      * StoreSolInBuffer stores the given range of the variable        *
!      * indicated by solName in IOVar and copies it into buffer if     *
!      * desired. It is assumed that the variables in blockPointers     *
!      * already point to the correct block.                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use cgnsNames
       use flowVarRefState
       use inputPhysics
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType), intent(in) :: kBeg, kEnd

       real(kind=realType), dimension(*), intent(out) :: buffer
       character(len=*), intent(in)                   :: solName

       logical, intent(in) :: copyInBuffer
!
!      Local parameters
!
       real(kind=realType), parameter :: plim   = 0.001_realType
       real(kind=realType), parameter :: rholim = 0.001_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, kk, nn

       real(kind=realType) :: uuy, uuz, vvx, vvz, wwx, wwy, tmp
       real(kind=realType) :: vortx, vorty, vortz, a2, ptotInf, ptot
       real(kind=realType) :: a, UovA(3), gradP(3)

       real(kind=realType), dimension(:,:,:,:), pointer :: wIO
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointer to the correct entry of IOVar. I'm cheating a
       ! bit here, because I know that only memory has been allocated
       ! for the first solution ID of IOVar.

       wIO => IOVar(nbkLocal,1)%w

       ! Determine the variable to be stored, compute it and store
       ! it in the 1D array buffer.

       select case(solName)

         case (cgnsDensity)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,irho)
               enddo
             enddo
           enddo

         case (cgnsMomx)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivx)
               enddo
             enddo
           enddo

         case (cgnsMomy)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivy)
               enddo
             enddo
           enddo

         case (cgnsMomz)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivz)
               enddo
             enddo
           enddo

         case (cgnsEnergy)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,irhoE)
               enddo
             enddo
           enddo

         case (cgnsTurbSaNu,cgnsTurbK)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,itu1)
               enddo
             enddo
           enddo

         case (cgnsTurbOmega,cgnsTurbTau,cgnsTurbEpsilon)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,itu2)
               enddo
             enddo
           enddo

         case (cgnsTurbV2)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,itu3)
               enddo
             enddo
           enddo

         case (cgnsTurbF)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,itu4)
               enddo
             enddo
           enddo

         case (cgnsVelx)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivx)
               enddo
             enddo
           enddo

         case (cgnsVely)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivy)
               enddo
             enddo
           enddo

         case (cgnsVelz)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivz)
               enddo
             enddo
           enddo

         case (cgnsRelVelx)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivx)-s(i,j,k,1)
               enddo
             enddo
           enddo

         case (cgnsRelVely)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivy)-s(i,j,k,2)
               enddo
             enddo
           enddo

         case (cgnsRelVelz)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = w(i,j,k,ivz)-s(i,j,k,3)
               enddo
             enddo
           enddo

         case (cgnsPressure)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = p(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsTemp)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = p(i,j,k)/(RGas*w(i,j,k,irho))
               enddo
             enddo
           enddo

         case (cgnsCp)
           tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = tmp*(p(i,j,k) - pInf)
               enddo
             enddo
           enddo

         case (cgnsMach)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                     / max(w(i,j,k,irho),rholim)
                 tmp = (w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                     +  w(i,j,k,ivz)**2)/a2
                 wIO(i,j,k,1) = sqrt(max(zero,tmp))
               enddo
             enddo
           enddo

         case (cgnsRelMach)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                     / max(w(i,j,k,irho),rholim)
                 tmp = ((w(i,j,k,ivx)-s(i,j,k,1))**2 +&
                      (w(i,j,k,ivy)-s(i,j,k,2))**2 &
                      +(w(i,j,k,ivz)-s(i,j,k,3))**2)/a2
                 wIO(i,j,k,1) = sqrt(max(zero,tmp))
               enddo
             enddo
           enddo


         case (cgnsMachTurb)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 tmp = w(i,j,k,irho)*w(i,j,k,itu1) &
                     / (gamma(i,j,k)*max(p(i,j,k),plim))
                 wIO(i,j,k,1) = sqrt(max(zero,tmp))
               enddo
             enddo
           enddo

         case (cgnsEddy)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = rev(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsEddyRatio)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = rev(i,j,k)/rlv(i,j,k)
               enddo
             enddo
           enddo

         case (cgNSWallDist)
           do k=kBeg,kEnd
             kk = max(2_intType,k); kk = min(kl,kk)
             do j=jBeg,jEnd
               jj = max(2_intType,j); jj = min(jl,jj)
               do i=iBeg,iEnd
                 ii = max(2_intType,i); ii = min(il,ii)
                 wIO(i,j,k,1) = d2Wall(ii,jj,kk)
               enddo
             enddo
           enddo

         case (cgnsVortMagn)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 tmp = half/vol(i,j,k)
                 uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                    - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                    + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                    - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                    + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                    - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                 uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                    - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                    + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                    - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                    + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                    - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                 vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                    - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                    + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                    - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                    + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                    - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                 vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                    - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                    + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                    - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                    + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                    - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                 wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                    - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                    + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                    - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                    + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                    - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                 wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                    - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                    + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                    - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                    + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                    - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                 vortx = wwy - vvz; vorty = uuz - wwx; vortz = vvx - uuy

                 wIO(i,j,k,1) = tmp*sqrt(vortx**2 + vorty**2 + vortz**2)
               enddo
             enddo
           enddo

         case (cgnsVortx)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 tmp = half/vol(i,j,k)
                 vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                      - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                      + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                      - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                      + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                      - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                 wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                    - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                    + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                    - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                    + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                    - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                 wIO(i,j,k,1) = tmp*(wwy - vvz)
               enddo
             enddo
           enddo

         case (cgnsVorty)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 tmp = half/vol(i,j,k)
                 uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                    - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                    + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                    - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                    + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                    - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                 wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                    - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                    + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                    - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                    + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                    - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                 wIO(i,j,k,1) = tmp*(uuz - wwx)
               enddo
             enddo
           enddo

         case (cgnsVortz)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 tmp = half/vol(i,j,k)
                 uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                    - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                    + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                    - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                    + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                    - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                 vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                    - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                    + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                    - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                    + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                    - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                 wIO(i,j,k,1) = tmp*(vvx - uuy)
               enddo
             enddo
           enddo

         case (cgnsPtotloss)

           ! Compute the free stream total pressure.

           call computePtot(rhoInf, uInf, zero, zero, &
                            pInf, ptotInf, 1_intType)
           ptotInf = one/ptotInf

           ! Loop over the cell centers and compute the
           ! total pressure loss.

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 call computePtot(w(i,j,k,irho), w(i,j,k,ivx), &
                                  w(i,j,k,ivy),  w(i,j,k,ivz), &
                                  p(i,j,k),      ptot, 1_intType)

                 wIO(i,j,k,1) = one - ptot*ptotInf
               enddo
             enddo
           enddo

         case (cgnsResRho)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,irho)
                 wIO(i,j,k,1) = dw(i,j,k,irho)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResMomx)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,imx)
                 wIO(i,j,k,1) = dw(i,j,k,imx)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResMomy)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,imy)
                 wIO(i,j,k,1) = dw(i,j,k,imy)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResMomz)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,imz)
                 wIO(i,j,k,1) = dw(i,j,k,imz)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResRhoE)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,irhoE)
                 wIO(i,j,k,1) = dw(i,j,k,irhoE)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResNu,cgnsResK)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,itu1)
                 wIO(i,j,k,1) = dw(i,j,k,itu1)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResOmega,cgnsResTau,cgnsResEpsilon)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,itu2)
                 wIO(i,j,k,1) = dw(i,j,k,itu2)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResV2)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,itu3)
                 wIO(i,j,k,1) = dw(i,j,k,itu3)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsResF)

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
             !   wIO(i,j,k,1) = dw(i,j,k,itu4)
                 wIO(i,j,k,1) = dw(i,j,k,itu4)/vol(i,j,k)
               enddo
             enddo
           enddo

         case (cgnsBlank)
           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 wIO(i,j,k,1) = real(iblank(i,j,k),realType)
               enddo
             enddo
           enddo

         case (cgnsGC)
            do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                  wIO(i,j,k,1) = real(globalcell(i,j,k),realType)
               enddo
             enddo
           enddo

         case (cgnsstatus)
            do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                  wIO(i,j,k,1) = real(fringes(i,j,k)%status)
               enddo
             enddo
           enddo

        case (cgnsShock)
           
           do k=kBeg,kEnd
              do j=jBeg,jEnd
                 do i=iBeg,iEnd

                    ! Here we compute U/a <dot> grad P / ||grad P||
                    ! Whre U is the velocity vector, a is the speed of
                    ! sound and P is the pressure. 

                    ! U / a
                    a  = sqrt(gamma(i,j,k)*max(p(i,j,k),plim) &
                         / max(w(i,j,k,irho),rholim))
                    UovA = (/w(i,j,k,ivx), w(i,j,k,ivy), w(i,j,k,ivz)/)/a

                    ! grad P / ||grad P||

                    gradP(1) = si(i,  j,k,1)*P(i+1,j,k) &
                         - si(i-1,j,k,1)*P(i-1,j,k) &
                         + sj(i,j,  k,1)*P(i,j+1,k) &
                         - sj(i,j-1,k,1)*P(i,j-1,k) &
                         + sk(i,j,k,  1)*P(i,j,k+1) &
                         - sk(i,j,k-1,1)*P(i,j,k-1)

                    gradP(2) = si(i,  j,k,2)*P(i+1,j,k) &
                         - si(i-1,j,k,2)*P(i-1,j,k) &
                         + sj(i,j,  k,2)*P(i,j+1,k) &
                         - sj(i,j-1,k,2)*P(i,j-1,k) &
                         + sk(i,j,k,  2)*P(i,j,k+1) &
                         - sk(i,j,k-1,2)*P(i,j,k-1)

                    gradP(3) = si(i,  j,k,3)*P(i+1,j,k) &
                         - si(i-1,j,k,3)*P(i-1,j,k) &
                         + sj(i,j,  k,3)*P(i,j+1,k) &
                         - sj(i,j-1,k,3)*P(i,j-1,k) &
                         + sk(i,j,k,  3)*P(i,j,k+1) &
                         - sk(i,j,k-1,3)*P(i,j,k-1)

                    gradP = gradP / sqrt(gradP(1)**2 + gradP(2)**2 + gradP(3)**2 )
                    ! Dot product
                    wIO(i,j,k,1) = UovA(1)*gradP(1) + UovA(2)*gradP(2) + UovA(3)*gradP(3)
                 end do
              end do
           end do
           

         case default
           call returnFail("storeSolInBuffer", &
                          "This should not happen")

       end select

       ! Copy the data in the 1D buffer, if desired.

       if( copyInBuffer ) then
         nn = 0
         do k=kBeg,kEnd
           do j=jBeg,jEnd
             do i=iBeg,iEnd
               nn = nn + 1
               buffer(nn) = wIO(i,j,k,1)
             enddo
           enddo
         enddo
       endif

       end subroutine storeSolInBuffer
