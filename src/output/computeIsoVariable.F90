!
!      ******************************************************************
!      *                                                                *
!      * File:          computeIsoVariable.F90                          *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 07-21-2013                                      *
!      * Last modified: 07-21-2013                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine computeIsoVariable(solName, sps, isoVal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * computeIsoVar computes NODE centered values for the given      *
  !      * solName variable. It is essentially equilivent to              *
  !      * sotreSolInBuffer. It is assumed blockPointers are already      *
  !      * set to the correct block.                                      *
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
  character(len=*), intent(in)                   :: solName
  integer(kind=intType), intent(in) :: sps
  real(kind=realType), intent(in) :: isoVal
  !
  !      Local parameters
  !
  real(kind=realType), parameter :: plim   = 0.001_realType
  real(kind=realType), parameter :: rholim = 0.001_realType
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, ii, jj, kk, nn

  real(kind=realType) :: uy, uz, vx, vz, wx, wy, tmp
  real(kind=realType) :: vortx, vorty, vortz, a2, ptotInf, ptot, uova(3), gradP(3), a
  real(kind=realType), dimension(:, :, :), pointer :: fc, fn
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  do nn=1,nDom
     call setPointers(nn, 1, sps)
     fc => flowDoms(nn, 1, sps)%fc
     fn => flowDoms(nn, 1, sps)%fn

     select case(solName)

     case (cgnsDensity)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,irho)
              enddo
           enddo
        enddo

     case (cgnsMomx)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivx)
              enddo
           enddo
        enddo

     case (cgnsMomy)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivy)
              enddo
           enddo
        enddo

     case (cgnsMomz)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivz)
              enddo
           enddo
        enddo

     case (cgnsEnergy)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,irhoE)
              enddo
           enddo
        enddo

     case (cgnsTurbSaNu,cgnsTurbK)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,itu1)
              enddo
           enddo
        enddo

     case (cgnsTurbOmega,cgnsTurbTau,cgnsTurbEpsilon)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,itu2)
              enddo
           enddo
        enddo

     case (cgnsTurbV2)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,itu3)
              enddo
           enddo
        enddo

     case (cgnsTurbF)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,itu4)
              enddo
           enddo
        enddo

     case (cgnsVelx)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivx)
              enddo
           enddo
        enddo

     case (cgnsVely)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivy)
              enddo
           enddo
        enddo

     case (cgnsVelz)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivz)
              enddo
           enddo
        enddo

     case (cgnsRelVelx)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivx)-s(i,j,k,1)
              enddo
           enddo
        enddo

     case (cgnsRelVely)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivy)-s(i,j,k,2)
              enddo
           enddo
        enddo

     case (cgnsRelVelz)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = w(i,j,k,ivz)-s(i,j,k,3)
              enddo
           enddo
        enddo

     case (cgnsPressure)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = p(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsTemp)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = p(i,j,k)/(RGas*w(i,j,k,irho))
              enddo
           enddo
        enddo

     case (cgnsCp)
        tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = tmp*(p(i,j,k) - pInf)
              enddo
           enddo
        enddo

     case (cgnsMach)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                      / max(w(i,j,k,irho),rholim)
                 tmp = (w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                      +  w(i,j,k,ivz)**2)/a2
                 fc(i,j,k) = sqrt(max(zero,tmp))
              enddo
           enddo
        enddo

     case (cgnsRelMach)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                      / max(w(i,j,k,irho),rholim)
                 tmp = ((w(i,j,k,ivx)-s(i,j,k,1))**2 +&
                      (w(i,j,k,ivy)-s(i,j,k,2))**2 &
                      +(w(i,j,k,ivz)-s(i,j,k,3))**2)/a2
                 fc(i,j,k) = sqrt(max(zero,tmp))
              enddo
           enddo
        enddo


     case (cgnsMachTurb)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 tmp = w(i,j,k,irho)*w(i,j,k,itu1) &
                      / (gamma(i,j,k)*max(p(i,j,k),plim))
                 fc(i,j,k) = sqrt(max(zero,tmp))
              enddo
           enddo
        enddo

     case (cgnsEddy)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = rev(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsEddyRatio)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = rev(i,j,k)/rlv(i,j,k)
              enddo
           enddo
        enddo

     case (cgNSWallDist)
        do k=1,ke
           kk = max(2_intType,k); kk = min(kl,kk)
           do j=1,je
              jj = max(2_intType,j); jj = min(jl,jj)
              do i=1,ie
                 ii = max(2_intType,i); ii = min(il,ii)
                 fc(i,j,k) = d2Wall(ii,jj,kk)
              enddo
           enddo
        enddo

     case (cgnsVortMagn)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 tmp = half/vol(i,j,k)
                 uy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                      - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                      + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                      - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                      + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                      - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                 uz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                      - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                      + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                      - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                      + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                      - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                 vx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                      - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                      + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                      - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                      + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                      - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                 vz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                      - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                      + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                      - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                      + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                      - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                 wx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                      - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                      + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                      - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                      + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                      - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                 wy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                      - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                      + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                      - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                      + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                      - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                 vortx = wy - vz; vorty = uz - wx; vortz = vx - uy

                 fc(i,j,k) = tmp*sqrt(vortx**2 + vorty**2 + vortz**2)
              enddo
           enddo
        enddo

     case (cgnsVortx)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 tmp = half/vol(i,j,k)
                 vz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                      - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                      + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                      - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                      + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                      - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                 wy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                      - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                      + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                      - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                      + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                      - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                 fc(i,j,k) = tmp*(wy - vz)
              enddo
           enddo
        enddo

     case (cgnsVorty)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 tmp = half/vol(i,j,k)
                 uz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                      - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                      + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                      - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                      + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                      - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                 wx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                      - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                      + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                      - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                      + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                      - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                 fc(i,j,k) = tmp*(uz - wx)
              enddo
           enddo
        enddo

     case (cgnsVortz)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 tmp = half/vol(i,j,k)
                 uy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                      - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                      + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                      - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                      + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                      - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                 vx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                      - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                      + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                      - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                      + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                      - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                 fc(i,j,k) = tmp*(vx - uy)
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

        do k=1,ke
           do j=1,je
              do i=1,ie
                 call computePtot(w(i,j,k,irho), w(i,j,k,ivx), &
                      w(i,j,k,ivy),  w(i,j,k,ivz), &
                      p(i,j,k),      ptot, 1_intType)

                 fc(i,j,k) = one - ptot*ptotInf
              enddo
           enddo
        enddo

     case (cgnsResRho)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,irho)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResMomx)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,imx)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResMomy)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,imy)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResMomz)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,imz)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResRhoE)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,irhoE)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResNu,cgnsResK)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,itu1)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResOmega,cgnsResTau,cgnsResEpsilon)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,itu2)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResV2)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,itu3)/vol(i,j,k)
              enddo
           enddo
        enddo

     case (cgnsResF)

        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = dw(i,j,k,itu4)/vol(i,j,k)
              enddo
           enddo
        enddo

        case (cgnsShock)
           
           do k=1,ke
              do j=1,je
                 do i=1,ie

                    ! Here we compute U/a <dot> grad P / ||grad P||
                    ! Whre U is the velocity vector, a is the speed of
                    ! sound and P is the pressure. 

                    ! U / a
                    a  = sqrt(gamma(i,j,k)*max(p(i,j,k),plim) &
                         / max(w(i,j,k,irho),rholim))
                    
                    if (addGridVelocities) then
                       UovA = (/w(i,j,k,ivx)-s(i,j,k,1), &
                            w(i,j,k,ivy)-s(i,j,k,2), &
                            w(i,j,k,ivz)-s(i,j,k,3)/)/a
                    else
                       UovA = (/w(i,j,k,ivx),w(i,j,k,ivy), w(i,j,k,ivz)/)/a
                    end if
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

                    ! Protect against divide by zero
                    gradP = gradP / sqrt(gradP(1)**2 + gradP(2)**2 + gradP(3)**2 + 1e-12)

                    ! Dot product
                    fc(i,j,k) = UovA(1)*gradP(1) + UovA(2)*gradP(2) + UovA(3)*gradP(3)
                 end do
              end do
           end do

     case (cgnsBlank)
        do k=1,ke
           do j=1,je
              do i=1,ie
                 fc(i,j,k) = real(min(iblank(i,j,k),1_intType),realType)
              enddo
           enddo
        enddo

     case default
        call terminate("computeIsoVariable", &
             "This should not happen")

     end select

     ! We now create nodal values from the cell centered
     ! values. This was the reason for going from 1 to ie
     ! etc.

     do k=1,kl
        do j=1,jl 
           do i=1,il
              fn(i,j,k) = eighth*( &
                   fc(i  , j  , k  ) + &
                   fc(i+1, j  , k  ) + &
                   fc(i  , j+1, k  ) + &
                   fc(i+1, j+1, k  ) + &
                   fc(i  , j  , k+1) + &
                   fc(i+1, j  , k+1) + &
                   fc(i  , j+1, k+1) + &
                   fc(i+1, j+1, k+1)) - isoVal
           end do
        end do
     end do
  end do
end subroutine computeIsoVariable
