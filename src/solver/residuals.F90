module residuals
contains

  subroutine residual_block
    !
    !       residual computes the residual of the mean flow equations on   
    !       the current MG level.                                          
    !
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputIteration
    use inputDiscretization
    use inputTimeSpectral
    use inputUnsteady ! Added by HDN
    use iteration
    use inputAdjoint
    use flowUtils, only : computeSpeedOfSoundSquared, allNodalGradients
    use fluxes
    use ALEUtils, only : interpLevelALE_block, recoverLevelALE_block
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: discr
    integer(kind=intType) :: i, j, k, l
    integer(kind=intType) :: iale, jale, kale, lale, male ! For loops of ALE
    real(kind=realType), parameter :: K1 = 1.05_realType
    real(kind=realType), parameter :: K2 = 0.6_realType ! Random given number
    real(kind=realType), parameter :: M0 = 0.2_realType ! Mach number preconditioner activation
    real(kind=realType), parameter :: alpha = 0_realType
    real(kind=realType), parameter :: delta = 0_realType
    !real(kind=realType), parameter :: hinf = 2_realType ! Test phase 
    real(kind=realType), parameter :: Cpres = 4.18_realType ! Test phase
    real(kind=realType), parameter :: TEMP= 297.15_realType

    !
    !     Local variables
    !
    real(kind=realType) :: K3, h, velXrho, velYrho, velZrho,SoS,hinf
    real(kind=realType) :: resM,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35
    real(kind=realType) :: A41,A42,A43,A44,A45,A51,A52,A53,A54,A55,B11,B12,B13,B14,B15
    real(kind=realType) :: B21,B22,B23,B24,B25,B31,B32,B33,B34,B35
    real(kind=realType) :: B41,B42,B43,B44,B45,B51,B52,B53,B54,B55
    real(kind=realType) :: rhoHdash, betaMr2
    real(kind=realType) :: G, q
    real(kind=realType) :: b1, b2, b3, b4, b5
    real(kind=realType) :: dwo(nwf)
    logical :: fineGrid


    ! Set the value of rFil, which controls the fraction of the old
    ! dissipation residual to be used. This is only for the runge-kutta
    ! schemes; for other smoothers rFil is simply set to 1.0.
    ! Note the index rkStage+1 for cdisRK. The reason is that the
    ! residual computation is performed before rkStage is incremented.

    if(smoother == RungeKutta) then
       rFil = cdisRK(rkStage+1)
    else
       rFil = one
    endif

    ! Set the value of the discretization, depending on the grid level,
    ! and the logical fineGrid, which indicates whether or not this
    ! is the finest grid level of the current mg cycle.

    discr = spaceDiscrCoarse
    if(currentLevel == 1) discr = spaceDiscr

    fineGrid = .false.
    if(currentLevel == groundLevel) fineGrid = .true.


    ! ===========================================================
    !
    ! Assuming ALE has nothing to do with MG
    ! The geometric data will be interpolated if in MD mode
    !
    ! ===========================================================
#ifndef USE_TAPENADE
    call interpLevelALE_block
#endif
    ! ===========================================================
    !
    ! The fluxes are calculated as usual
    !
    ! ===========================================================

    call inviscidCentralFlux

    select case (discr)

    case (dissScalar) ! Standard scalar dissipation scheme.

       if( fineGrid) then 
          if (.not. lumpedDiss) then
             call inviscidDissFluxScalar
          else
             call inviscidDissFluxScalarApprox
          end if
       else
#ifndef USE_TAPENADE
          call inviscidDissFluxScalarCoarse
#endif
       endif

       !===========================================================

    case (dissMatrix) ! Matrix dissipation scheme.

       if( fineGrid ) then
          if (.not. lumpedDiss) then 
             call inviscidDissFluxMatrix
          else
             call inviscidDissFluxMatrixApprox
          end if
       else
#ifndef USE_TAPENADE
          call inviscidDissFluxMatrixCoarse
#endif
       endif

       !===========================================================

    case (upwind) ! Dissipation via an upwind scheme.

       call inviscidUpwindFlux(fineGrid)

    end select

    !-------------------------------------------------------
    ! Lastly, recover the old s[I,J,K], sFace[I,J,K]
    ! This shall be done before difussive and source terms
    ! are computed.
    !-------------------------------------------------------
#ifndef USE_TAPENADE
    call recoverLevelALE_block
#endif


    if( viscous ) then 
       ! Only compute viscous fluxes if rFil > 0
       if(abs(rFil) > thresholdReal) then 
          ! not lumpedDiss means it isn't the PC...call the vicousFlux
          if (.not. lumpedDiss) then 
             call computeSpeedOfSoundSquared
             call allNodalGradients
             call viscousFlux
          else
             ! This is a PC calc...only include viscous fluxes if viscPC
             ! is used
             call computeSpeedOfSoundSquared
             if (viscPC) then
                call allNodalGradients
                call viscousFlux 
             else
                call viscousFluxApprox
             end if
          end if
       end if
    end if

    !===========================================================

    ! Add the dissipative and possibly viscous fluxes to the
    ! Euler fluxes. Loop over the owned cells and add fw to dw.
    ! Also multiply by iblank so that no updates occur in holes
    if ( lowspeedpreconditioner ) then
       do k=2,kl
          do j=2,jl
             do i=2,il
                !    Compute speed of sound
                SoS = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))

                ! Coompute velocities without rho from state vector
                velXrho = w(i,j,k,ivx)
                velYrho = w(i,j,k,ivy)
                velZrho = w(i,j,k,ivz)

                q = (velXrho**2 + velYrho**2 + velZrho**2)
                resM=sqrt(q)/SoS

                !
                !    Compute K3
                K3 = K1 * ( 1 + ((1-K1*M0**2)*resM**2)/(K1*M0**4) )
                !    Compute BetaMr2
                betaMr2 = min( max( K3*(velXrho**2 + velYrho**2  &
                     + velZrho**2), ((K2)*(wInf(ivx)**2 &
                     + wInf(ivy)**2 + wInf(ivz)**2))) , SoS**2 )

                A11=  (betaMr2)*(1/SoS**4)
                A12 = zero
                A13 = zero
                A14 = zero
                A15 = (-betaMr2)/SoS**4 

                A21 =one*velXrho/SoS**2 
                A22 = one*w(i,j,k,irho)
                A23 = zero
                A24 = zero
                A25 = one*(-velXrho)/SoS**2

                A31 = one*velYrho/SoS**2
                A32 = zero
                A33 = one*w(i,j,k,irho)
                A34 = zero 
                A35 = one*(-velYrho)/SoS**2

                A41 = one*velZrho/SoS**2 
                A42 = zero
                A43 = zero
                A44 = one*w(i,j,k,irho) 
                A45 = zero + one *(-velZrho)/SoS**2 

                A51=  one*((1/(gamma(i,j,k)-1))+(resM**2)/2)
                A52 = one * w(i,j,k,irho)*velXrho
                A53 = one * w(i,j,k,irho)*velYrho
                A54 = one * w(i,j,k,irho)*velzrho
                A55 = one * ((-(resM**2))/2)

                B11 = A11*(gamma(i,j,k)-1)*q/2+A12*(-velXrho)&
                     /w(i,j,k,irho)+A13*(-velYrho)/w(i,j,k,irho)+A14*(-velZrho)/w(i,j,k,irho)&
                     + A15* (((gamma(i,j,k)-1)*q/2)-SoS**2)
                B12 = A11*(1-gamma(i,j,k))*velXrho+A12*1/w(i,j,k,irho)&
                     +A15*(1-gamma(i,j,k))*velXrho
                B13 = A11*(1-gamma(i,j,k))*velYrho+A13&
                     /w(i,j,k,irho)+A15*(1-gamma(i,j,k))*velYrho
                B14 = A11*(1-gamma(i,j,k))*velZrho&
                     +A14/w(i,j,k,irho)+A15*(1-gamma(i,j,k))*velZrho
                B15 = A11*(gamma(i,j,k)-1)+A15*(gamma(i,j,k)-1)

                B21 = A21*(gamma(i,j,k)-1)*q/2+A22*(-velXrho)&
                     /w(i,j,k,irho)+A23*(-velYrho)/w(i,j,k,irho)+A24*(-velZrho)&
                     /w(i,j,k,irho)+ A25* (((gamma(i,j,k)-1)*q/2)-SoS**2)
                B22 = A21*(1-gamma(i,j,k))*velXrho+A22&
                     /w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velXrho
                B23 = A21*(1-gamma(i,j,k))*velYrho&
                     +A23*1/w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velYrho
                B24 = A21*(1-gamma(i,j,k))*velZrho&
                     +A24*1/w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velZrho
                B25 = A21*(gamma(i,j,k)-1)+A25*(gamma(i,j,k)-1)

                B31 = A31*(gamma(i,j,k)-1)*q/2+A32*(-velXrho)&
                     /w(i,j,k,irho)+A33*(-velYrho)/w(i,j,k,irho)+A34*(-velZrho)/w(i,j,k,irho)&
                     + A35* (((gamma(i,j,k)-1)*q/2)-SoS**2)
                B32 = A31*(1-gamma(i,j,k))*velXrho+A32&
                     /w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velXrho
                B33 = A31*(1-gamma(i,j,k))*velYrho&
                     +A33*1/w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velYrho
                B34 = A31*(1-gamma(i,j,k))*velZrho&
                     +A34*1/w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velZrho
                B35 = A31*(gamma(i,j,k)-1)+A35*(gamma(i,j,k)-1)

                B41 = A41*(gamma(i,j,k)-1)*q/2+A42*(-velXrho)&
                     /w(i,j,k,irho)+A43*(-velYrho)/w(i,j,k,irho)+A44*(-velZrho) &
                     /w(i,j,k,irho)+ A45* (((gamma(i,j,k)-1)*q/2)-SoS**2)
                B42 = A41*(1-gamma(i,j,k))*velXrho+A42&
                     /w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velXrho
                B43 = A41*(1-gamma(i,j,k))*velYrho&
                     +A43*1/w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velYrho
                B44 = A41*(1-gamma(i,j,k))*velZrho&
                     +A44*1/w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velZrho
                B45 = A41*(gamma(i,j,k)-1)+A45*(gamma(i,j,k)-1)

                B51 = A51*(gamma(i,j,k)-1)*q/2+A52*(-velXrho)&
                     /w(i,j,k,irho)+A53*(-velYrho)/w(i,j,k,irho)+A54*(-velZrho) &
                     /w(i,j,k,irho)+ A55* (((gamma(i,j,k)-1)*q/2)-SoS**2)
                B52 = A51*(1-gamma(i,j,k))*velXrho+A52&
                     /w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velXrho
                B53 = A51*(1-gamma(i,j,k))*velYrho&
                     +A53*1/w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velYrho
                B54 = A51*(1-gamma(i,j,k))*velZrho&
                     +A54*1/w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velZrho
                B55 = A51*(gamma(i,j,k)-1)+A55*(gamma(i,j,k)-1)

                ! dwo is the orginal redisual
                do l=1,nwf
                   dwo(l) = (dw(i,j,k,l) + fw(i,j,k,l))* max(real(iblank(i,j,k), realType), zero)
                end do

                dw(i,j,k,1)=B11*dwo(1) + B12*dwo(2)+ B13*dwo(3) + B14*dwo(4) + B15*dwo(5)
                dw(i,j,k,2)=B21*dwo(1) + B22*dwo(2)+ B23*dwo(3) + B24*dwo(4) + B25*dwo(5)
                dw(i,j,k,3)=B31*dwo(1) + B32*dwo(2)+ B33*dwo(3) + B34*dwo(4) + B35*dwo(5)
                dw(i,j,k,4)=B41*dwo(1) + B42*dwo(2)+ B43*dwo(3) + B44*dwo(4) + B45*dwo(5)
                dw(i,j,k,5)=B51*dwo(1) + B52*dwo(2)+ B53*dwo(3) + B54*dwo(4) + B55*dwo(5)

             enddo
          enddo
       enddo
    else
       do l=1,nwf
          do k=2,kl
             do j=2,jl
                do i=2,il
                   dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
                        * max(real(iblank(i,j,k), realType), zero)
                enddo
             enddo
          enddo
       enddo
    endif

  end subroutine residual_block

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE
  subroutine initres(varStart, varEnd)
    !
    ! Shell function to call initRes_block on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use iteration
    use section
    use utils, only : setPointers
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(in) :: varStart, varEnd
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, nn

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, currentLevel, sps)

          call initres_block(varStart, varEnd, nn, sps)

       end do domains

    end do spectralLoop

  end subroutine initRes

  subroutine initres_block(varStart, varEnd, nn, sps)
    !
    !       initres initializes the given range of the residual. Either to 
    !       zero, steady computation, or to an unsteady term for the time  
    !       spectral and unsteady modes. For the coarser grid levels the   
    !       residual forcing term is taken into account.                   
    !
    use blockPointers
    use flowVarRefState
    use inputIteration
    use inputPhysics
    use inputTimeSpectral
    use inputUnsteady
    use iteration

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: varStart, varEnd, nn, sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm, ll, ii, jj, i, j, k, l, m
    real(kind=realType)   :: oneOverDt, tmp

    real(kind=realType), dimension(:,:,:,:), pointer :: ww, wsp, wsp1
    real(kind=realType), dimension(:,:,:),   pointer :: volsp

    ! Return immediately of no variables are in the range.

    if(varEnd < varStart) return

    ! Determine the equation mode and act accordingly.

    select case (equationMode)
    case (steady)

       ! Steady state computation.
       ! Determine the currently active multigrid level.

       steadyLevelTest: if(currentLevel == groundLevel) then

          ! Ground level of the multigrid cycle. Initialize the
          ! owned residuals to zero.

          do l=varStart,varEnd
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      dw(i,j,k,l) = zero
                   enddo
                enddo
             enddo
          enddo
       else steadyLevelTest

          ! Coarse grid level. Initialize the owned cells to the
          ! residual forcing terms.

          do l=varStart,varEnd
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      dw(i,j,k,l) = wr(i,j,k,l)
                   enddo
                enddo
             enddo
          enddo
       endif steadyLevelTest

       !===========================================================

    case (unsteady)

       ! Unsteady computation.
       ! A further distinction must be made.

       select case(timeIntegrationScheme)

       case (explicitRK)

          ! We are always on the finest grid.
          ! Initialize the residual to zero.

          do l=varStart,varEnd
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      dw(i,j,k,l) = zero
                   enddo
                enddo
             enddo
          enddo

          !=======================================================

       case (MD,BDF) ! Modified by HDN

          ! Store the inverse of the physical nonDimensional
          ! time step a bit easier.

          oneOverDt = timeRef/deltaT

          ! Store the pointer for the variable to be used to compute
          ! the unsteady source term. For a runge-kutta smoother this
          ! is the solution of the zeroth runge-kutta stage. As for
          ! rkStage == 0 this variable is not yet set w is used.
          ! For other smoothers w is to be used as well.

          if(smoother == RungeKutta .and. rkStage > 0) then
             ww => wn
          else
             ww => w
          endif

          ! Determine the currently active multigrid level.

          unsteadyLevelTest: if(currentLevel == groundLevel) then

             ! Ground level of the multigrid cycle. Initialize the
             ! owned cells to the unsteady source term. First the
             ! term for the current time level. Note that in w the
             ! velocities are stored and not the momentum variables.
             ! Therefore the if-statement is present to correct this.

             do l=varStart,varEnd

                if(l == ivx .or. l == ivy .or. l == ivz) then

                   ! Momentum variables.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                                 * ww(i,j,k,l)*ww(i,j,k,irho)
                         enddo
                      enddo
                   enddo

                else

                   ! Non-momentum variables, for which the variable
                   ! to be solved is stored; for the flow equations this
                   ! is the conservative variable, for the turbulent
                   ! equations the primitive variable.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                                 * ww(i,j,k,l)
                         enddo
                      enddo
                   enddo

                endif

             enddo

             ! The terms from the older time levels. Here the
             ! conservative variables are stored. In case of a
             ! deforming mesh, also the old volumes must be taken.

             deformingTest: if( deforming_Grid ) then

                ! Mesh is deforming and thus the volumes can change.
                ! Use the old volumes as well.

                do m=1,nOldLevels
                   do l=varStart,varEnd
                      do k=2,kl
                         do j=2,jl
                            do i=2,il
                               dw(i,j,k,l) = dw(i,j,k,l)                 &
                                    + coefTime(m)*volOld(m,i,j,k) &
                                    * wOld(m,i,j,k,l)
                            enddo
                         enddo
                      enddo
                   enddo
                enddo

             else deformingTest

                ! Rigid mesh. The volumes remain constant.

                do m=1,nOldLevels
                   do l=varStart,varEnd
                      do k=2,kl
                         do j=2,jl
                            do i=2,il
                               dw(i,j,k,l) = dw(i,j,k,l)            &
                                    + coefTime(m)*vol(i,j,k) &
                                    * wOld(m,i,j,k,l)
                            enddo
                         enddo
                      enddo
                   enddo
                enddo

             endif deformingTest

             ! Multiply the time derivative by the inverse of the
             ! time step to obtain the true time derivative.
             ! This is done after the summation has been done, because
             ! otherwise you run into finite accuracy problems for
             ! very small time steps.

             do l=varStart,varEnd
                do k=2,kl
                   do j=2,jl
                      do i=2,il
                         dw(i,j,k,l) = oneOverDt*dw(i,j,k,l)
                      enddo
                   enddo
                enddo
             enddo

          else unsteadyLevelTest

             ! Coarse grid level. Initialize the owned cells to the
             ! residual forcing term plus a correction for the
             ! multigrid treatment of the time derivative term.
             ! As the velocities are stored instead of the momentum,
             ! these terms must be multiplied by the density.

             tmp = oneOverDt*coefTime(0)

             do l=varStart,varEnd

                if(l == ivx .or. l == ivy .or. l == ivz) then

                   ! Momentum variables.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = tmp*vol(i,j,k)               &
                                 * (ww(i,j,k,l)*ww(i,j,k,irho)  &
                                 -  w1(i,j,k,l)*w1(i,j,k,irho))
                            dw(i,j,k,l) = dw(i,j,k,l) + wr(i,j,k,l)
                         enddo
                      enddo
                   enddo

                else

                   ! Non-momentum variables.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = tmp*vol(i,j,k)             &
                                 * (ww(i,j,k,l) - w1(i,j,k,l))
                            dw(i,j,k,l) =  dw(i,j,k,l) + wr(i,j,k,l)
                         enddo
                      enddo
                   enddo

                endif

             enddo

          endif unsteadyLevelTest

       end select

       !===========================================================

    case (timeSpectral)

       ! Time spectral computation. The time derivative of the
       ! current solution is given by a linear combination of
       ! all other solutions, i.e. a matrix vector product.

       ! First store the section to which this block belongs
       ! in jj.

       jj = sectionID

       ! Determine the currently active multigrid level.

       spectralLevelTest: if(currentLevel == groundLevel) then

          ! Finest multigrid level. The residual must be
          ! initialized to the time derivative.

          ! Initialize it to zero.

          do l=varStart,varEnd
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      dw(i,j,k,l) = zero
                   enddo
                enddo
             enddo
          enddo

          ! Loop over the number of terms which contribute
          ! to the time derivative.

          timeLoopFine: do mm=1,nTimeIntervalsSpectral

             ! Store the pointer for the variable to be used to
             ! compute the unsteady source term and the volume.
             ! Also store in ii the offset needed for vector
             ! quantities.

             wsp   => flowDoms(nn,currentLevel,mm)%w
             volsp => flowDoms(nn,currentLevel,mm)%vol
             ii    =  3*(mm-1)

             ! Loop over the number of variables to be set.

             varLoopFine: do l=varStart,varEnd

                ! Test for a momentum variable.

                if(l == ivx .or. l == ivy .or. l == ivz) then

                   ! Momentum variable. A special treatment is
                   ! needed because it is a vector and the velocities
                   ! are stored instead of the momentum. Set the
                   ! coefficient ll, which defines the row of the
                   ! matrix used later on.

                   if(l == ivx) ll = 3*sps - 2
                   if(l == ivy) ll = 3*sps - 1
                   if(l == ivz) ll = 3*sps

                   ! Loop over the owned cell centers to add the
                   ! contribution from wsp.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il

                            ! Store the matrix vector product with the
                            ! velocity in tmp.

                            tmp = dvector(jj,ll,ii+1)*wsp(i,j,k,ivx) &
                                 + dvector(jj,ll,ii+2)*wsp(i,j,k,ivy) &
                                 + dvector(jj,ll,ii+3)*wsp(i,j,k,ivz)

                            ! Update the residual. Note the
                            ! multiplication with the density to obtain
                            ! the correct time derivative for the
                            ! momentum variable.

                            dw(i,j,k,l) = dw(i,j,k,l) &
                                 + tmp*volsp(i,j,k)*wsp(i,j,k,irho)

                         enddo
                      enddo
                   enddo

                else

                   ! Scalar variable.  Loop over the owned cells to
                   ! add the contribution of wsp to the time
                   ! derivative.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = dw(i,j,k,l)        &
                                 + dscalar(jj,sps,mm) &
                                 * volsp(i,j,k)*wsp(i,j,k,l)

                         enddo
                      enddo
                   enddo

                endif

             enddo varLoopFine

          enddo timeLoopFine
       else spectralLevelTest

          ! Coarse grid level. Initialize the owned cells to the
          ! residual forcing term plus a correction for the
          ! multigrid treatment of the time derivative term.

          ! Initialization to the residual forcing term.

          do l=varStart,varEnd
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      dw(i,j,k,l) = wr(i,j,k,l)
                   enddo
                enddo
             enddo
          enddo

          ! Loop over the number of terms which contribute
          ! to the time derivative.

          timeLoopCoarse: do mm=1,nTimeIntervalsSpectral

             ! Store the pointer for the variable to be used to
             ! compute the unsteady source term and the pointers
             ! for wsp1, the solution when entering this MG level
             ! and for the volume.
             ! Furthermore store in ii the offset needed for
             ! vector quantities.

             wsp   => flowDoms(nn,currentLevel,mm)%w
             wsp1  => flowDoms(nn,currentLevel,mm)%w1
             volsp => flowDoms(nn,currentLevel,mm)%vol
             ii    =  3*(mm-1)

             ! Loop over the number of variables to be set.

             varLoopCoarse: do l=varStart,varEnd

                ! Test for a momentum variable.

                if(l == ivx .or. l == ivy .or. l == ivz) then

                   ! Momentum variable. A special treatment is
                   ! needed because it is a vector and the velocities
                   ! are stored instead of the momentum. Set the
                   ! coefficient ll, which defines the row of the
                   ! matrix used later on.

                   if(l == ivx) ll = 3*sps - 2
                   if(l == ivy) ll = 3*sps - 1
                   if(l == ivz) ll = 3*sps

                   ! Add the contribution of wps to the correction
                   ! of the time derivative. The difference between
                   ! the current time derivative and the one when
                   ! entering this grid level must be added, because
                   ! the residual forcing term only takes the spatial
                   ! part of the coarse grid residual into account.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il

                            ! Store the matrix vector product with the
                            ! momentum in tmp.

                            tmp = dvector(jj,ll,ii+1)                &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivx)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivx)) &
                                 +  dvector(jj,ll,ii+2)               &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivy)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivy)) &
                                 + dvector(jj,ll,ii+3)                &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivz)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivz))

                            ! Add tmp to the residual. Multiply it by
                            ! the volume to obtain the finite volume
                            ! formulation of the  derivative of the
                            ! momentum.

                            dw(i,j,k,l) = dw(i,j,k,l) + tmp*volsp(i,j,k)

                         enddo
                      enddo
                   enddo

                else

                   ! Scalar variable. Loop over the owned cells
                   ! to add the contribution of wsp to the correction
                   ! of the time derivative.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dw(i,j,k,l) = dw(i,j,k,l)        &
                                 + dscalar(jj,sps,mm) &
                                 * volsp(i,j,k)       &
                                 * (wsp(i,j,k,l) - wsp1(i,j,k,l))
                         enddo
                      enddo
                   enddo

                endif

             enddo varLoopCoarse

          enddo timeLoopCoarse
       endif spectralLevelTest

    end select

    ! Set the residual in the halo cells to zero. This is just
    ! to avoid possible problems. Their values do not matter.

    do l=varStart,varEnd
       do k=0,kb
          do j=0,jb
             dw(0,j,k,l)  = zero
             dw(1,j,k,l)  = zero
             dw(ie,j,k,l) = zero
             dw(ib,j,k,l) = zero
          enddo
       enddo

       do k=0,kb
          do i=2,il
             dw(i,0,k,l)  = zero
             dw(i,1,k,l)  = zero
             dw(i,je,k,l) = zero
             dw(i,jb,k,l) = zero
          enddo
       enddo

       do j=2,jl
          do i=2,il
             dw(i,j,0,l)  = zero
             dw(i,j,1,l)  = zero
             dw(i,j,ke,l) = zero
             dw(i,j,kb,l) = zero
          enddo
       enddo
    enddo

  end subroutine initres_block

  subroutine residual
    !
    ! Shell function to call residual_block on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use Iteration
    use utils, only : setPointers
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, nn

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, currentLevel, sps)

          call residual_block

       end do domains

    end do spectralLoop
  end subroutine residual

  subroutine computedwDADI
    !
    !       executeRkStage executes one runge kutta stage. The stage       
    !       number, rkStage, is defined in the local module iteration.     
    !
    use blockPointers
    use constants
    use flowVarRefState
    use inputIteration
    use inputPhysics
    use inputTimeSpectral
    use inputUnsteady
    use iteration
    implicit none
    !
    !      Local parameter.
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, n
    real(kind=realType) :: gm1, epsval, fac, eps2
    real(kind=realType) :: uvel, vvel, wvel, cijk, cijkinv, c2inv, uvw
    real(kind=realType) :: ri1, ri2, ri3, rj1, rj2, rj3, rk1, rk2, rk3, uu
    real(kind=realType) :: ri, rj, rk, qsi, qsj, qsk, currentCfl
    real(kind=realType) :: xfact,  cInf, cInf2
    real(kind=realType) :: dw1, dw2, dw3, dw4, dw5
    real(kind=realType) :: a1, a2, a3, a4, a5, a6, a7, mut, ge
    real(kind=realType) :: viscTerm1, viscTerm2, viscTerm3
    real(kind=realType) :: metterm
    real(kind=realType) :: unsteadyImpl, mult
    real(kind=realType) :: sqrt2, sqrt2inv, alph, alphinv
    real(kind=realType) :: volhalf, volhalfrho, volfact
    real(kind=realType) :: mutpI, mutmI, mutpJ, mutmJ, mutpK, mutmK, &
         mettermp,mettermm, &
         viscTermI,viscTermJ,viscTermK

    real(kind=realType), dimension(:,:,:), pointer:: qq_i,qq_j,qq_k, &
         cc_i,cc_j,cc_k,spectral_i,spectral_j,spectral_k,dual_dt
    real(kind=realType), dimension(ie,5) :: bbi,cci,ddi,ffi
    real(kind=realType), dimension(je,5) :: bbj,ccj,ddj,ffj
    real(kind=realType), dimension(ke,5) :: bbk,cck,ddk,ffk
    real(kind=realType), dimension(ie) :: mettermi
    real(kind=realType), dimension(je) :: mettermj
    real(kind=realType), dimension(ke) :: mettermk
    real(kind=realType), dimension(5) :: diagPlus,diagMinus


    ! Set the value of the current cfl number,
    ! depending on the situation. On the finest grid in the mg cycle
    ! the second halo is computed, otherwise not.

    currentCfl = cflCoarse
    if (currentLevel == 1) then
       currentCfl = cfl
    end if
    qq_i => scratch(:,:,:,1)
    qq_j => scratch(:,:,:,2)
    qq_k => scratch(:,:,:,3)
    cc_i => scratch(:,:,:,4)
    cc_j => scratch(:,:,:,5)
    cc_k => scratch(:,:,:,6)
    spectral_i => scratch(:,:,:,7)
    spectral_j => scratch(:,:,:,8)
    spectral_k => scratch(:,:,:,9)
    dual_dt => scratch(:,:,:,10)

    !  havent thought about iblank

    !   these are factors for robustness. 

    epsval = 0.08
    fac    = 1.05
    mut    = zero
    cInf2  = gammaInf*pInf/rhoInf
    cInf   = sqrt(cInf2)


    qsi = zero
    qsj = zero
    qsk = zero
    viscTermi=zero
    viscTermj=zero
    viscTermk=zero
    sqrt2=sqrt(two)
    sqrt2inv=one/sqrt2


    if(equationMode.eq.steady) then
       do k=2,kl
          do j=2,jl
             do i=2,il
                dual_dt(i,j,k)= currentCfl*dtl(i,j,k)*vol(i,j,k)
             enddo
          enddo
       enddo
    else  
       do k=2,kl
          do j=2,jl
             do i=2,il
                unsteadyImpl = coefTime(0)*timeRef/deltaT
                mult = currentCfl*dtl(i,j,k)
                mult = mult/(mult*unsteadyImpl*vol(i,j,k) + one)
                dual_dt(i,j,k)   = mult*vol(i,j,k)
             enddo
          enddo
       enddo
    endif

    !     Set up some arrays
    do k=2,kl
       do j=2,jl
          do i=2,il
             volhalf    = half/vol(i,j,k)
             volhalfrho = volhalf/w(i,j,k,irho)
             uvel  = w(i,j,k,ivx)
             vvel  = w(i,j,k,ivy)
             wvel  = w(i,j,k,ivz)

             cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))

             ri1  = volhalf*(si(i,j,k,1) + si(i-1,j,k,1))
             ri2  = volhalf*(si(i,j,k,2) + si(i-1,j,k,2))
             ri3  = volhalf*(si(i,j,k,3) + si(i-1,j,k,3))

             rj1  = volhalf*(sj(i,j,k,1) + sj(i,j-1,k,1))
             rj2  = volhalf*(sj(i,j,k,2) + sj(i,j-1,k,2))
             rj3  = volhalf*(sj(i,j,k,3) + sj(i,j-1,k,3))

             rk1  = volhalf*(sk(i,j,k,1) + sk(i,j,k-1,1))
             rk2  = volhalf*(sk(i,j,k,2) + sk(i,j,k-1,2))
             rk3  = volhalf*(sk(i,j,k,3) + sk(i,j,k-1,3))

             if( addGridVelocities ) then
                qsi = (sFaceI(i-1,j,k) + sFaceI(i,j,k))*volhalf
                qsj = (sFaceJ(i,j-1,k) + sFaceJ(i,j,k))*volhalf
                qsk = (sFaceK(i,j,k-1) + sFaceK(i,j,k))*volhalf
             endif


             qq_i(i,j,k) = ri1*w(i,j,k,ivx) + ri2*w(i,j,k,ivy) +  &
                  ri3*w(i,j,k,ivz) - qsi
             qq_j(i,j,k) = rj1*w(i,j,k,ivx) + rj2*w(i,j,k,ivy) +  &
                  rj3*w(i,j,k,ivz) - qsj
             qq_k(i,j,k) = rk1*w(i,j,k,ivx) + rk2*w(i,j,k,ivy) +  &
                  rk3*w(i,j,k,ivz) - qsk

             ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
             rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
             rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)

             cc_i(i,j,k) = cijk*ri
             cc_j(i,j,k) = cijk*rj
             cc_k(i,j,k) = cijk*rk
             if( viscous )  then
                mutpI = rlv(i,j,k)+rlv(i+1,j,k)
                mutmI = rlv(i,j,k)+rlv(i-1,j,k)
                mutpJ = rlv(i,j,k)+rlv(i,j+1,k)
                mutmJ = rlv(i,j,k)+rlv(i,j-1,k)
                mutpK = rlv(i,j,k)+rlv(i,j,k+1)
                mutmK = rlv(i,j,k)+rlv(i,j,k-1)
                if( eddyModel ) then
                   mutpI = mutpI+rev(i,j,k)+rev(i+1,j,k)
                   mutmI = mutpI+rev(i,j,k)+rev(i-1,j,k)
                   mutpJ = mutpJ+rev(i,j,k)+rev(i,j+1,k)
                   mutmJ = mutpJ+rev(i,j,k)+rev(i,j-1,k)
                   mutpK = mutpK+rev(i,j,k)+rev(i,j,k+1)
                   mutmK = mutpK+rev(i,j,k)+rev(i,j,k-1)
                endif

                volfact  =   two/(vol(i,j,k)+vol(i+1,j,k))
                mettermp =   (si(i,j,k,1)*si(i,j,k,1) &
                     +si(i,j,k,2)*si(i,j,k,2) &
                     +si(i,j,k,3)*si(i,j,k,3))*mutpI*volfact
                volfact  =   two/(vol(i,j,k)+vol(i-1,j,k))
                mettermm =   (si(i-1,j,k,1)*si(i-1,j,k,1) &
                     +si(i-1,j,k,2)*si(i-1,j,k,2) &
                     +si(i-1,j,k,3)*si(i-1,j,k,3))*mutmI*volfact
                viscTermi=   (mettermp+mettermm)*volhalfrho

                volfact  =   two/(vol(i,j,k)+vol(i,j+1,k))
                mettermp =   (sj(i,j,k,1)*sj(i,j,k,1) &
                     +sj(i,j,k,2)*sj(i,j,k,2) &
                     +sj(i,j,k,3)*sj(i,j,k,3))*mutpJ*volfact
                volfact  =   two/(vol(i,j,k)+vol(i,j-1,k))
                mettermm =   (sj(i,j-1,k,1)*sj(i,j-1,k,1) &
                     +sj(i,j-1,k,2)*sj(i,j-1,k,2) &
                     +sj(i,j-1,k,3)*sj(i,j-1,k,3))*mutmJ*volfact
                viscTermj=   (mettermp+mettermm)*volhalfrho

                volfact  =   two/(vol(i,j,k)+vol(i,j,k+1))
                mettermp =   (sk(i,j,k,1)*sk(i,j,k,1) &
                     +sk(i,j,k,2)*sk(i,j,k,2) &
                     +sk(i,j,k,3)*sk(i,j,k,3))*mutpK*volfact
                volfact  =   two/(vol(i,j,k)+vol(i,j,k-1))
                mettermm =   (sk(i,j,k-1,1)*sk(i,j,k-1,1) &
                     +sk(i,j,k-1,2)*sk(i,j,k-1,2) &
                     +sk(i,j,k-1,3)*sk(i,j,k-1,3))*mutmK*volfact
                viscTermk=   (mettermp+mettermm)*volhalfrho

             endif

             eps2 = ri*cInf*epsval
             spectral_i(i,j,k) = (fac*sqrt((abs(qq_i(i,j,k))+cc_i(i,j,k))**2 &
                  + eps2**2)+viscTermi)*dual_dt(i,j,k)
             eps2 = rj*cInf*epsval
             spectral_j(i,j,k) = (fac*sqrt((abs(qq_j(i,j,k))+cc_j(i,j,k))**2 &
                  + eps2**2)+viscTermj)*dual_dt(i,j,k)
             eps2 = rk*cInf*epsval
             spectral_k(i,j,k) = (fac*sqrt((abs(qq_k(i,j,k))+cc_k(i,j,k))**2 & 
                  + eps2**2)+viscTermk)*dual_dt(i,j,k)
             spectral_i(i,j,k)=spectral_i(i,j,k)*zero
             spectral_j(i,j,k)=spectral_j(i,j,k)*zero
             spectral_k(i,j,k)=spectral_k(i,j,k)*zero
          enddo
       enddo
    enddo

    !       Multiply by T_eta^inv
    do k=2,kl
       do j=2,jl
          do i=2,il

             gm1   = gamma(i,j,k)-one
             cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
             c2inv = one/(cijk*cijk)
             xfact = two*cijk
             alphinv  = sqrt2*cijk/w(i,j,k,irho)

             uvel  = w(i,j,k,ivx)
             vvel  = w(i,j,k,ivy)
             wvel  = w(i,j,k,ivz)
             uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)

             rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
             rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
             rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
             rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
             uu   = uvel*rj1+vvel*rj2+wvel*rj3
             rj1  = rj1/rj
             rj2  = rj2/rj
             rj3  = rj3/rj

             dw1   = dw(i,j,k,1)
             dw2   = dw(i,j,k,2)
             dw3   = dw(i,j,k,3)
             dw4   = dw(i,j,k,4)
             dw5   = dw(i,j,k,5)

             a1    =  dw2*uvel+dw3*vvel+dw4*wvel-dw5
             a1    =  a1*gm1*c2inv+dw1*(one-uvw*gm1*c2inv)

             a2    =  (rj2*wvel-rj3*vvel)*dw1+rj3*dw3-rj2*dw4
             a3    =  (rj3*uvel-rj1*wvel)*dw1+rj1*dw4-rj3*dw2
             a4    =  (rj1*vvel-rj2*uvel)*dw1+rj2*dw2-rj1*dw3

             a5    =  uvw*dw1-uvel*dw2-vvel*dw3-wvel*dw4+dw5
             a5    =  a5*gm1*c2inv

             a6    =  uu*dw1/rj-rj1*dw2-rj2*dw3-rj3*dw4

             dw(i,j,k,1)  = a1*rj1+a2/w(i,j,k,irho)
             dw(i,j,k,2)  = a1*rj2+a3/w(i,j,k,irho)
             dw(i,j,k,3)  = a1*rj3+a4/w(i,j,k,irho)
             dw(i,j,k,4)  = (half*a5-a6/xfact)*alphinv
             dw(i,j,k,5)  = (half*a5+a6/xfact)*alphinv

          enddo
       enddo
    enddo

    if(jl.gt.2) then

       !  Inversion in j

       do k=2,kl
          do i=2,il
             do j=1,jl
                if( viscous )   mut = rlv(i,j,k)+rlv(i,j+1,k)
                if( eddyModel ) mut = mut+rev(i,j,k)+rev(i,j+1,k)

                volfact=one/(vol(i,j,k)+vol(i,j+1,k))
                metterm =   (sj(i,j,k,1)*sj(i,j,k,1) &
                     +sj(i,j,k,2)*sj(i,j,k,2) &
                     +sj(i,j,k,3)*sj(i,j,k,3))
                mettermj(j)= metterm*mut*volfact
             enddo

             do j=2,jl

                viscTerm1 = mettermj(j)  /vol(i,j,k)/w(i,j,k,irho)
                viscTerm3 = mettermj(j-1)/vol(i,j,k)/w(i,j,k,irho)
                viscTerm2 = viscTerm1+viscTerm3

                volhalf = half/vol(i,j,k)
                rj1  = volhalf*(sj(i,j,k,1) + sj(i,j-1,k,1))
                rj2  = volhalf*(sj(i,j,k,2) + sj(i,j-1,k,2))
                rj3  = volhalf*(sj(i,j,k,3) + sj(i,j-1,k,3))
                metterm     = rj1*rj1+rj2*rj2+rj3*rj3
                eps2        = epsval*epsval*cInf2*metterm
                diagPlus(1) =half*(qq_j(i,j,k)+fac*sqrt(qq_j(i,j,k)**2+eps2))
                diagPlus(2) =diagPlus(1)
                diagPlus(3) =diagPlus(1)
                diagPlus(4) =half*(qq_j(i,j,k)+cc_j(i,j,k)   &
                     +fac*sqrt((qq_j(i,j,k)+cc_j(i,j,k))**2+eps2))
                diagPlus(5) =half*(qq_j(i,j,k)-cc_j(i,j,k)   &
                     +fac*sqrt((qq_j(i,j,k)-cc_j(i,j,k))**2+eps2))
                diagMinus(1) =half*(qq_j(i,j,k)-fac*sqrt(qq_j(i,j,k)**2+eps2))
                diagMinus(2) =diagMinus(1)
                diagMinus(3) =diagMinus(1)
                diagMinus(4) =half*(qq_j(i,j,k)+cc_j(i,j,k)   &
                     -fac*sqrt((qq_j(i,j,k)+cc_j(i,j,k))**2+eps2))
                diagMinus(5) =half*(qq_j(i,j,k)-cc_j(i,j,k)   &
                     -fac*sqrt((qq_j(i,j,k)-cc_j(i,j,k))**2+eps2))

                do n=1,5
                   bbj(j+1,n)= -viscTerm1-diagPlus(n)
                   ddj(j-1,n)= -viscTerm3+diagMinus(n)
                   ccj(j  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n) 
                enddo
             enddo

             do n=1,5
                bbj(je,n)=zero
                ddj(1 ,n)=zero
                do j=2,jl
                   bbj(j,n)=bbj(j,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   ddj(j,n)=ddj(j,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   ccj(j,n)=one+ccj(j,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)+spectral_i(i,j,k)+spectral_k(i,j,k)
                   ffj(j,n)=dw(i,j,k,n)
                enddo
             enddo

             call tridiagsolve(bbj,ccj,ddj,ffj,jl)

             do n=1,5
                do j=2,jl
                   dw(i,j,k,n)=ffj(j,n)
                enddo
             enddo

          enddo
       enddo

    endif
    !       Multiply by T_xi^inv T_eta


    do k=2,kl 
       do j=2,jl
          do i=2,il
             ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
             ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
             ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
             ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
             ri1  = ri1/ri
             ri2  = ri2/ri
             ri3  = ri3/ri

             rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
             rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
             rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
             rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
             rj1  = rj1/rj
             rj2  = rj2/rj
             rj3  = rj3/rj

             dw1   = dw(i,j,k,1)
             dw2   = dw(i,j,k,2)
             dw3   = dw(i,j,k,3)
             dw4   = dw(i,j,k,4)
             dw5   = dw(i,j,k,5)

             a1 = ( ri1*rj1 + ri2*rj2 + ri3*rj3 )
             a2 = ( ri1*rj2 - rj1*ri2 )
             a3 = ( ri3*rj2 - rj3*ri2 )
             a4 = ( ri1*rj3 - rj1*ri3 )
             a5 = (dw4-dw5)*sqrt2inv
             a6 = (dw4+dw5)*half
             a7 = (a3*dw1+a4*dw2-a2*dw3-a5*a1)*sqrt2inv

             dw(i,j,k,1)=  a1*dw1 + a2*dw2 + a4*dw3 + a5*a3
             dw(i,j,k,2)=- a2*dw1 + a1*dw2 - a3*dw3 + a5*a4
             dw(i,j,k,3)=- a4*dw1 + a3*dw2 + a1*dw3 - a5*a2
             dw(i,j,k,4)= -a7+a6
             dw(i,j,k,5)=  a7+a6

          enddo
       enddo
    enddo

    !  Multiply by diagonal

    do k=2,kl
       do j=2,jl
          do i=2,il
             xfact = one+spectral_i(i,j,k)+spectral_j(i,j,k) &
                  +spectral_k(i,j,k)	
             dw(i,j,k,1)=dw(i,j,k,1)*xfact
             dw(i,j,k,2)=dw(i,j,k,2)*xfact
             dw(i,j,k,3)=dw(i,j,k,3)*xfact
             dw(i,j,k,4)=dw(i,j,k,4)*xfact
             dw(i,j,k,5)=dw(i,j,k,5)*xfact
          enddo
       enddo
    enddo



    if(il.gt.2) then

       !  Inversion in i

       do k=2,kl
          do j=2,jl
             do i=1,il
                if( viscous )   mut = rlv(i,j,k)+rlv(i+1,j,k)
                if( eddyModel ) mut = mut+rev(i,j,k)+rev(i+1,j,k)

                volfact=one/(vol(i,j,k)+vol(i+1,j,k))
                metterm    =(si(i,j,k,1)*si(i,j,k,1) &
                     +si(i,j,k,2)*si(i,j,k,2) &
                     +si(i,j,k,3)*si(i,j,k,3))
                mettermi(i)=metterm*mut*volfact
             enddo

             do i=2,il

                viscTerm1 = mettermi(i)  /vol(i,j,k)/w(i,j,k,irho)
                viscTerm3 = mettermi(i-1)/vol(i,j,k)/w(i,j,k,irho)
                viscTerm2 = viscTerm1+viscTerm3

                volhalf = half/vol(i,j,k)
                ri1  = volhalf*(si(i,j,k,1) + si(i-1,j,k,1))
                ri2  = volhalf*(si(i,j,k,2) + si(i-1,j,k,2))
                ri3  = volhalf*(si(i,j,k,3) + si(i-1,j,k,3))
                metterm     = ri1*ri1+ri2*ri2+ri3*ri3
                eps2        = epsval*epsval*cInf2*metterm
                diagPlus(1) =half*(qq_i(i,j,k)+fac*sqrt(qq_i(i,j,k)**2+eps2))
                diagPlus(2) =diagPlus(1)
                diagPlus(3) =diagPlus(1)
                diagPlus(4) =half*(qq_i(i,j,k)+cc_i(i,j,k)   &
                     +fac*sqrt((qq_i(i,j,k)+cc_i(i,j,k))**2+eps2))
                diagPlus(5) =half*(qq_i(i,j,k)-cc_i(i,j,k)   &
                     +fac*sqrt((qq_i(i,j,k)-cc_i(i,j,k))**2+eps2))
                diagMinus(1) =half*(qq_i(i,j,k)-fac*sqrt(qq_i(i,j,k)**2+eps2))
                diagMinus(2) =diagMinus(1)
                diagMinus(3) =diagMinus(1)
                diagMinus(4) =half*(qq_i(i,j,k)+cc_i(i,j,k)   &
                     -fac*sqrt((qq_i(i,j,k)+cc_i(i,j,k))**2+eps2))
                diagMinus(5) =half*(qq_i(i,j,k)-cc_i(i,j,k)   &
                     -fac*sqrt((qq_i(i,j,k)-cc_i(i,j,k))**2+eps2))
                do n=1,5
                   bbi(i+1,n)= -viscTerm1-diagPlus(n)
                   ddi(i-1,n)= -viscTerm3+diagMinus(n)
                   cci(i  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n)
                enddo
             enddo

             do n=1,5
                bbi(ie ,n)=zero
                ddi(1 ,n)=zero
                do i=2,il
                   bbi(i,n)=bbi(i,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   ddi(i,n)=ddi(i,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   cci(i,n)=one+cci(i,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)+spectral_j(i,j,k)+spectral_k(i,j,k)
                   ffi(i,n)=dw(i,j,k,n)
                enddo
             enddo


             call tridiagsolve(bbi,cci,ddi,ffi,il)

             do n=1,5
                do i=2,il
                   dw(i,j,k,n)=ffi(i,n)
                enddo
             enddo

          enddo
       enddo

    endif
    !       Multiply by T_zeta^inv T_xi


    do k=2,kl 
       do j=2,jl
          do i=2,il
             ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
             ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
             ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
             ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
             ri1  = ri1/ri
             ri2  = ri2/ri
             ri3  = ri3/ri

             rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
             rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
             rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))
             rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)
             rk1  = rk1/rk
             rk2  = rk2/rk
             rk3  = rk3/rk

             dw1   = dw(i,j,k,1)
             dw2   = dw(i,j,k,2)
             dw3   = dw(i,j,k,3)
             dw4   = dw(i,j,k,4)
             dw5   = dw(i,j,k,5)

             a1 =  ri1*rk1 + ri2*rk2 + ri3*rk3 
             a2 =  rk1*ri2 - ri1*rk2 
             a3 =  rk3*ri2 - ri3*rk2 
             a4 =  rk1*ri3 - ri1*rk3 
             a5 = (dw4-dw5)*sqrt2inv
             a6 = (dw4+dw5)*half
             a7 = (a3*dw1+a4*dw2-a2*dw3-a5*a1)*sqrt2inv

             dw(i,j,k,1)=  a1*dw1 + a2*dw2 + a4*dw3 + a5*a3
             dw(i,j,k,2)=- a2*dw1 + a1*dw2 - a3*dw3 + a5*a4
             dw(i,j,k,3)=- a4*dw1 + a3*dw2 + a1*dw3 - a5*a2
             dw(i,j,k,4)= -a7+a6
             dw(i,j,k,5)=  a7+a6
          enddo
       enddo
    enddo

    !  Multiply by diagonal

    do k=2,kl
       do j=2,jl
          do i=2,il
             xfact = one+spectral_i(i,j,k)+spectral_j(i,j,k) &
                  +spectral_k(i,j,k)	
             dw(i,j,k,1)=dw(i,j,k,1)*xfact
             dw(i,j,k,2)=dw(i,j,k,2)*xfact
             dw(i,j,k,3)=dw(i,j,k,3)*xfact
             dw(i,j,k,4)=dw(i,j,k,4)*xfact
             dw(i,j,k,5)=dw(i,j,k,5)*xfact
          enddo
       enddo
    enddo

    if(kl.gt.2) then

       !  Inversion in k

       do j=2,jl
          do i=2,il
             do k=1,kl
                if( viscous )   mut = rlv(i,j,k)+rlv(i,j,k+1)
                if( eddyModel ) mut = mut+rev(i,j,k)+rev(i,j,k+1)

                volfact=1./(vol(i,j,k)+vol(i,j,k+1))
                metterm    = sk(i,j,k,1)*sk(i,j,k,1)  &
                     +sk(i,j,k,2)*sk(i,j,k,2)  &
                     +sk(i,j,k,3)*sk(i,j,k,3)
                mettermk(k)= metterm*mut*volfact
             enddo

             do k=2,kl

                viscTerm1 = mettermk(k)  /vol(i,j,k)/w(i,j,k,irho)
                viscTerm3 = mettermk(k-1)/vol(i,j,k)/w(i,j,k,irho)
                viscTerm2 = viscTerm1+viscTerm3

                volhalf = half/vol(i,j,k)
                rk1  = volhalf*(sk(i,j,k,1) + sj(i,j,k-1,1))
                rk2  = volhalf*(sk(i,j,k,2) + sj(i,j,k-1,2))
                rk3  = volhalf*(sk(i,j,k,3) + sj(i,j,k-1,3))
                metterm     = rk1*rk1+rk2*rk2+rk3*rk3
                eps2        = epsval*epsval*cInf2*metterm
                diagPlus(1) =half*(qq_k(i,j,k)+fac*sqrt(qq_k(i,j,k)**2+eps2))
                diagPlus(2) =diagPlus(1)
                diagPlus(3) =diagPlus(1)
                diagPlus(4) =half*(qq_k(i,j,k)+cc_k(i,j,k)   &
                     +fac*sqrt((qq_k(i,j,k)+cc_k(i,j,k))**2+eps2))
                diagPlus(5) =half*(qq_k(i,j,k)-cc_k(i,j,k)   &
                     +fac*sqrt((qq_k(i,j,k)-cc_k(i,j,k))**2+eps2))
                diagMinus(1) =half*(qq_k(i,j,k)-fac*sqrt(qq_k(i,j,k)**2+eps2))
                diagMinus(2) =diagMinus(1)
                diagMinus(3) =diagMinus(1)
                diagMinus(4) =half*(qq_k(i,j,k)+cc_k(i,j,k)   &
                     -fac*sqrt((qq_k(i,j,k)+cc_k(i,j,k))**2+eps2))
                diagMinus(5) =half*(qq_k(i,j,k)-cc_k(i,j,k)   &
                     -fac*sqrt((qq_k(i,j,k)-cc_k(i,j,k))**2+eps2))

                do n=1,5
                   bbk(k+1,n)= -viscTerm1-diagPlus(n)
                   ddk(k-1,n)= -viscTerm3+diagMinus(n)
                   cck(k  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n)
                enddo
             enddo

             do n=1,5
                bbk(ke,n)=zero
                ddk(1 ,n)=zero
                do k=2,kl
                   bbk(k,n)=bbk(k,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   ddk(k,n)=ddk(k,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)
                   cck(k,n)=one+cck(k,n)*dual_dt(i,j,k)*max(real(iblank(i,j,k),realType),zero)+spectral_i(i,j,k)+spectral_j(i,j,k)
                   ffk(k,n)=dw(i,j,k,n)
                enddo
             enddo

             call tridiagsolve(bbk,cck,ddk,ffk,kl)

             do n=1,5
                do k=2,kl
                   dw(i,j,k,n)=ffk(k,n)
                enddo
             enddo

          enddo
       enddo
    endif

    !       Multiply by T_zeta

    do k=2,kl
       do j=2,jl
          do i=2,il
             uvel  = w(i,j,k,ivx)
             vvel  = w(i,j,k,ivy)
             wvel  = w(i,j,k,ivz)

             rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
             rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
             rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))

             rk = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)
             uu = uvel*rk1+vvel*rk2+wvel*rk3

             rk1  = rk1/rk
             rk2  = rk2/rk
             rk3  = rk3/rk

             uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)
             cijkinv = sqrt(w(i,j,k,irho)/gamma(i,j,k)/p(i,j,k))
             alph  = w(i,j,k,irho)*cijkinv*sqrt2inv
             xfact = two/cijkinv

             ge=gamma(i,j,k)*w(i,j,k,irhoE)/w(i,j,k,irho)-  &
                  (gamma(i,j,k)-one)*uvw  

             dw1   = dw(i,j,k,1)
             dw2   = dw(i,j,k,2)
             dw3   = dw(i,j,k,3)
             dw4   = dw(i,j,k,4)*alph
             dw5   = dw(i,j,k,5)*alph

             a1    =  dw1*rk1+dw2*rk2+dw3*rk3+dw4+dw5
             a2    =  half*xfact*(dw4-dw5)
             a3    =  uvw*(rk1*dw1+rk2*dw2+rk3*dw3)

             dw(i,j,k,1)  = a1
             dw(i,j,k,2)  = a1*uvel-w(i,j,k,irho)*(rk3*dw2-rk2*dw3) & 
                  +a2*rk1
             dw(i,j,k,3)  = a1*vvel-w(i,j,k,irho)*(rk1*dw3-rk3*dw1) & 
                  +a2*rk2
             dw(i,j,k,4)  = a1*wvel-w(i,j,k,irho)*(rk2*dw1-rk1*dw2) & 
                  +a2*rk3
             dw(i,j,k,5)  = a3+ w(i,j,k,irho)*       &
                  ((vvel*rk3-wvel*rk2)*dw1            &   
                  +(wvel*rk1-uvel*rk3)*dw2            &         
                  +(uvel*rk2-vvel*rk1)*dw3)           & 
                  +(ge+half*xfact*uu/rk)*dw4 &   
                  +(ge-half*xfact*uu/rk)*dw5       

          enddo
       enddo
    enddo


    !      For consistency with update.

    do k=2,kl
       do j=2,jl
          do i=2,il
             volfact=-one/vol(i,j,k)
             dw(i,j,k,1)=dw(i,j,k,1)*volfact
             dw(i,j,k,2)=dw(i,j,k,2)*volfact
             dw(i,j,k,3)=dw(i,j,k,3)*volfact
             dw(i,j,k,4)=dw(i,j,k,4)*volfact
             dw(i,j,k,5)=dw(i,j,k,5)*volfact
          enddo
       enddo
    enddo


  end subroutine computedwDADI

  subroutine tridiagsolve(bb,cc,dd,ff,nn)
    use precision
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType) :: nn
    real(kind=realType), dimension(nn+1,5) :: bb,cc,dd,ff

    !     local variables

    integer(kind=intType) :: m,n
    real(kind=realType) :: d0,d2

    do n=1,5
       m=2
       d0=1./cc(m,n)
       dd(m,n)=dd(m,n)*d0
       ff(m,n)=ff(m,n)*d0


       do m=3,nn
          d2=bb(m,n)
          d0=1./(cc(m,n)-d2*dd(m-1,n))
          ff(m,n)=(ff(m,n)-d2*ff(m-1,n))*d0
          dd(m,n)=dd(m,n)*d0
       enddo


       do m=nn-1,2,-1
          ff(m,n)=ff(m,n)-dd(m,n)*ff(m+1,n)
       enddo

    enddo


  end subroutine tridiagsolve


  subroutine residualAveraging
    !
    !       Implicit residual smoothing is a simple procedure that         
    !       replaces the residual at each point by a weighted sum of all   
    !       of the residuals in the block (although the residuals that are 
    !       closer to the cell under consideration are weighted more       
    !       heavily). This smoothing can be applied explicitly, but may    
    !       result in zero smoothed residuals for non-zero initial         
    !       residual modes.  For this reason, the smoothing is applied     
    !       implicitly in the following form:                              
    !       -epz R{i+1} + (1 + 2 epz) R{i} -epz R{i-1} = r{i}              
    !       Where r{i} is the original residual at point i, and R{i} is    
    !       the implicitly smoothed residual at point i.  The analysis for 
    !       the 1-D scalar convection-diffusion equation shows that if     
    !       Epz >= (1/4)*((lambda/lambda*)^2 -1),                          
    !       where lambda is the cfl number desired to be used, and         
    !       lambda* is the CFL limit of the unsmoothed scheme, the scheme  
    !       can be made unconditionally stable (arbitrarily large lambda). 
    !       In practice, lambda = 6-8 is common for Euler solutions.  For  
    !       RANS solutions lambda = 3-4 is what we have used in practice   
    !       resulting in a slight improvement in convergence rate, but,    
    !       more importantly, an increase in the robustness of the solver. 
    !       Note that this theoretical result can be shown for infinite    
    !       1-D problems, but also for finite-periodic 1-D problems and    
    !       finite-size 1-D problems (i.e. with block boundaries).  Such   
    !       theoretical results are not available for 3-D cases, where the 
    !       scheme is applied in an ADI fashion:                           
    !       (1 -epzI d_ii)(1 -epzJ d_jj)(1 -epzK d_kk) r = r               
    !       Where d_ii, d_jj, d_kk are second difference operators in each 
    !       of the mesh coordinate directions.                             
    !       For each of the coordinate direction solves, the initial       
    !       matrix problem is of the form:                                 
    !         -                                             - - -   - -    
    !         | (1+2 epz)   -epz                            | |r|   |r|    
    !         |   -epz    (1 +2 epz)    -epz                | |r|   |r|    
    !         |             -epz       (1 + 2 epz)   -epz   | |r| = |r|    
    !         |                 .           .           .   | |.|   |.|    
    !         |                   .            .          . | |.|   |.|    
    !         -                                             - - -   - -    
    !       And after the forward elimination phase a normalization is     
    !       applied so the result looks like:                              
    !         -                   - - -   - -                              
    !         |  1  -d            | |r|   |r|                              
    !         |  0   1  -d        | |r|   |r|                              
    !         |      0   1  -d    | |r| = |r|                              
    !         |       .   .   .   | |.|   |.|                              
    !         |          .  .   . | |.|   |.|                              
    !         -                   - - -   - -                              
    !       Which can then be used with a straightforward backsolve to     
    !       obtain the answer.                                             
    !       It is assumed that the pointers in blockPointers already       
    !       point to the correct block.                                    
    !
    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, dw, p, iBlank
    use inputIteration, only: cfl, cflCoarse, cflLimit, smoop
    use flowVarRefState, only : pInfCorr, nwf
    use iteration, only : currentLevel

    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l

    real(kind=realType), parameter :: b = 2.0_realType

    real(kind=realType) :: currentCfl, rfl0, plim
    real(kind=realType) :: dpi, dpj, dpk, r
    real(kind=realType), dimension(il,max(jl,kl)) :: epz, d, t, rfl


    !      rfl0 is a measure of the ratio lambda/lambda*
    !
    currentCfl = cflCoarse
    if (currentLevel == 1) then
       currentCfl = cfl
    end if

    rfl0  = half*currentCfl/cflLimit

    plim  = 0.001_realType*pInfCorr

    !       Smoothing in the i-direction. Only done when enough cells are  
    !       present in the i-direction.                                    
    !
    if(nx > 1) then

       do k=2,kl

          ! Compute smoothing coeficients

          do j=2,jl
             do i=2,il
                dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                     /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                     /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                     /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                rfl(i,j) = one/(one + b*(dpi  +dpj  +dpk))
             end do
          end do

          do j=2,jl
             do i=2,nx
                r = rfl0*(rfl(i,j) +rfl(i+1,j))
                epz(i,j) = fourth*smoop*dim(r*r,one)*max(real(iblank(i,j,k), realType), zero)
             end do
          end do

          ! Zero out coefficients for boundary condition treatment

          do j=2,jl
             epz(1,j)  = zero
             epz(il,j) = zero
             d(1,j)    = zero
          end do

          ! Compute coefficients for forward elimination process

          do i=2,il
             do j=2,jl
                t(i,j) = one &
                     / (one +epz(i,j) +epz(i-1,j) -epz(i-1,j)*d(i-1,j))
                d(i,j) = t(i,j)*epz(i,j)
             end do
          end do

          ! Apply same transformation to the rhs vector of residuals

          do i=2,il
             do j=2,jl
                do l=1,nwf
                   dw(i,j,k,l) = t(i,j) &
                        * (dw(i,j,k,l) +epz(i-1,j)*dw(i-1,j,k,l))
                end do
             end do
          end do

          ! Backsolve operation.  Smoothed residuals are left
          !  in dw(i,j,k,l)

          do i=nx,2,-1
             do j=2,jl
                do l=1,nwf
                   dw(i,j,k,l) = dw(i,j,k,l) +d(i,j)*dw(i+1,j,k,l)
                end do
             end do
          end do

       enddo
    endif
    !
    !       Smoothing in the j-direction. Only done when enough cells are  
    !       present in the j-direction.                                    
    !
    if(ny > 1) then

       do k=2,kl

          ! Compute smoothing coeficients

          do j=2,jl
             do i=2,il
                dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                     /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                     /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                     /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                rfl(i,j) = one/(one + b*(dpi  +dpj  +dpk))
             end do
          end do

          do j=2,ny
             do i=2,il
                r = rfl0*(rfl(i,j) +rfl(i,j+1))
                epz(i,j) = fourth*smoop*dim(r*r,one)*max(real(iblank(i,j,k), realType), zero)
             end do
          end do

          ! Zero out coefficients for boundary condition treatment

          do i=2,il
             epz(i,1)  = zero
             epz(i,jl) = zero
             d(i,1)    = zero
          end do

          ! Compute coefficients for forward eliMination process

          do j=2,jl
             do i=2,il
                t(i,j) = one &
                     / (one +epz(i,j) +epz(i,j-1) -epz(i,j-1)*d(i,j-1))
                d(i,j) = t(i,j)*epz(i,j)
             end do
          end do

          ! Apply same transformation to the rhs vector of residuals

          do j=2,jl
             do i=2,il
                do l=1,nwf
                   dw(i,j,k,l) = t(i,j) &
                        * (dw(i,j,k,l) +epz(i,j-1)*dw(i,j-1,k,l))
                end do
             end do
          end do

          ! Backsolve operation.  Smoothed residuals are left
          !  in dw(i,j,k,l)

          do j=ny,2,-1
             do i=2,il
                do l=1,nwf
                   dw(i,j,k,l) = dw(i,j,k,l) +d(i,j)*dw(i,j+1,k,l)
                end do
             end do
          end do

       enddo
    endif
    !
    !       Smoothing in the k-direction. Only done when enough cells are  
    !       present in the k-direction.                                    
    !
    if(nz > 1) then

       do j=2,jl

          ! Compute smoothing coeficients

          do k=2,kl
             do i=2,il
                dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                     /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                     /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                     /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                rfl(i,k) = one/(one + b*(dpi  +dpj  +dpk))
             end do
          end do

          do k=2,nz
             do i=2,il
                r = rfl0*(rfl(i,k) +rfl(i,k+1))
                epz(i,k) = fourth*smoop*dim(r*r,one)*max(real(iblank(i,j,k), realType), zero)
             end do
          end do

          ! Zero out coefficients for boundary condition treatment

          do i=2,il
             epz(i,1)  = zero
             epz(i,kl) = zero
             d(i,1)    = zero
          end do

          ! Compute coefficients for forward eliMination process

          do k=2,kl
             do i=2,il
                t(i,k) = one &
                     / (one +epz(i,k) +epz(i,k-1) -epz(i,k-1)*d(i,k-1))
                d(i,k) = t(i,k)*epz(i,k)
             end do
          end do

          ! Apply same transformation to the rhs vector of residuals

          do k=2,kl
             do i=2,il
                do l=1,nwf
                   dw(i,j,k,l) = t(i,k) &
                        * (dw(i,j,k,l) +epz(i,k-1)*dw(i,j,k-1,l))
                end do
             end do
          end do

          ! Backsolve operation.  Smoothed residuals are left
          !  in dw(i,j,k,l)

          do k=nz,2,-1
             do i=2,il
                do l=1,nwf
                   dw(i,j,k,l) = dw(i,j,k,l) +d(i,k)*dw(i,j,k+1,l)
                end do
             end do
          end do

       enddo
    endif

  end subroutine residualAveraging
#endif
end module residuals
