subroutine gridVelocitiesFineLevel(useOldCoor, t, sps)
  !
  ! Shell function to call gridVelocitiesFineLevel on all blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps
  logical,               intent(in) :: useOldCoor
  real(kind=realType), dimension(*), intent(in) :: t  !
  !      Local variables.
  !
  integer(kind=intType) :: nn

  ! Loop over the number of blocks.

  domains: do nn=1,nDom

     ! Set the pointers for this block.

     call setPointers(nn, groundLevel, sps)

     call gridVelocitiesFineLevel_block(useOldCoor, t, sps)

  end do domains

end subroutine gridVelocitiesFineLevel
!
!      ******************************************************************
!      *                                                                *
!      * File:          gridVelocities.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine gridVelocitiesFineLevel_block(useOldCoor, t, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * gridVelocitiesFineLevel computes the grid velocities for       *
  !      * the cell centers and the normal grid velocities for the faces  *
  !      * of moving blocks for the currently finest grid, i.e.           *
  !      * groundLevel. The velocities are computed at time t for         *
  !      * spectral mode sps. If useOldCoor is .true. the velocities      *
  !      * are determined using the unsteady time integrator in           *
  !      * combination with the old coordinates; otherwise the analytic   *
  !      * form is used.                                                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use cgnsGrid
  use flowVarRefState
  use inputMotion
  use inputUnsteady
  use iteration
  use inputPhysics
  use inputTSStabDeriv
  use monitor
  use communication
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps
  logical,               intent(in) :: useOldCoor

  real(kind=realType), dimension(*), intent(in) :: t
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, mm
  integer(kind=intType) :: i, j, k, ii, iie, jje, kke

  real(kind=realType) :: oneOver4dt, oneOver8dt
  real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
  real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

  real(kind=realType), dimension(3) :: sc, xc, xxc
  real(kind=realType), dimension(3) :: rotCenter, rotRate
  real(kind=realType), dimension(3) :: rotRateTemp
  real(kind=realType), dimension(3) :: offsetVector
  real(kind=realType), dimension(3,3) :: rotRateTrans

  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix,&
       derivRotationMatrix

  real(kind=realType) :: tNew, tOld
  real(kind=realType), dimension(:,:), pointer :: sFace

  real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
  real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

  integer(kind=intType) :: liftIndex
  real(kind=realType) :: alpha,beta,intervalMach,alphaTS,alphaIncrement,&
       betaTS,betaIncrement
  real(kind=realType), dimension(3) ::velDir,liftDir,dragDir
  real(kind=realType), dimension(3) :: refDirection
 
  !Function Definitions

  real(kind=realType) :: TSAlpha,TSBeta,TSMach
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Compute the mesh velocity from the given mesh Mach number.

  !  aInf = sqrt(gammaInf*pInf/rhoInf)
  !  velxGrid = aInf*MachGrid(1)
  !  velyGrid = aInf*MachGrid(2)
  !  velzGrid = aInf*MachGrid(3)

  aInf = sqrt(gammaInf*pInf/rhoInf)
  velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
  velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
  velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

  ! Compute the derivative of the rotation matrix and the rotation
  ! point; needed for velocity due to the rigid body rotation of
  ! the entire grid. It is assumed that the rigid body motion of
  ! the grid is only specified if there is only 1 section present.

  call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, t(1))

  !compute the rotation matrix to update the velocities for the time
  !spectral stability derivative case...

  if(TSStability)then
     ! Determine the time values of the old and new time level.
     ! It is assumed that the rigid body rotation of the mesh is only
     ! used when only 1 section is present.

     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)

     if(TSpMode.or. TSqMode .or.TSrMode) then
        ! Compute the rotation matrix of the rigid body rotation as
        ! well as the rotation point; the latter may vary in time due
        ! to rigid body translation.

        call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

        velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
             + rotationMatrix(1,2)*velygrid0 &
             + rotationMatrix(1,3)*velzgrid0
        velygrid0 = rotationMatrix(2,1)*velxgrid0 &
             + rotationMatrix(2,2)*velygrid0 &
             + rotationMatrix(2,3)*velzgrid0
        velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
             + rotationMatrix(3,2)*velygrid0 &
             + rotationMatrix(3,3)*velzgrid0
     elseif(tsAlphaMode)then
        ! get the baseline alpha and determine the liftIndex
        call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)
        !Determine the alpha for this time instance
        alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t(1))

        alphaTS = alpha+alphaIncrement
        !Determine the grid velocity for this alpha
        refDirection(:) = zero
        refDirection(1) = one
        call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

        !do I need to update the lift direction and drag direction as well?
        !set the effictive grid velocity for this time interval
        velxGrid0 = (aInf*machgrid)*(-velDir(1))
        velyGrid0 = (aInf*machgrid)*(-velDir(2))
        velzGrid0 = (aInf*machgrid)*(-velDir(3))

     elseif(tsBetaMode)then
        ! get the baseline alpha and determine the liftIndex
        call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)

        !Determine the alpha for this time instance
        betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
             degreeFourBeta,  omegaFourBeta,     &
             cosCoefFourBeta, sinCoefFourBeta, t(1))

        betaTS = beta+betaIncrement
        !Determine the grid velocity for this alpha
        refDirection(:) = zero
        refDirection(1) = one
        call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

        !do I need to update the lift direction and drag direction as well?
        !set the effictive grid velocity for this time interval
        velxGrid0 = (aInf*machgrid)*(-velDir(1))
        velyGrid0 = (aInf*machgrid)*(-velDir(2))
        velzGrid0 = (aInf*machgrid)*(-velDir(3))
     elseif(TSMachMode)then
        !determine the mach number at this time interval
        IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
             degreeFourMach,  omegaFourMach,     &
             cosCoefFourMach, sinCoefFourMach, t(1))
        !set the effective grid velocity
        velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
        velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
        velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

     elseif(TSAltitudeMode)then
        call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
     else
        call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
     end if
  endif

  testMoving: if( blockIsMoving ) then
     ! Determine the situation we are having here.

     testUseOldCoor: if( useOldCoor ) then
        !
        !            ************************************************************
        !            *                                                          *
        !            * The velocities must be determined via a finite           *
        !            * difference formula using the coordinates of the old      *
        !            * levels.                                                  *
        !            *                                                          *
        !            ************************************************************
        !
        ! Set the coefficients for the time integrator and store
        ! the inverse of the physical nonDimensional time step,
        ! divided by 4 and 8, a bit easier.

        call setCoefTimeIntegrator
        oneOver4dt = fourth*timeRef/deltaT
        oneOver8dt = half*oneOver4dt
        !
        !            ************************************************************
        !            *                                                          *
        !            * Grid velocities of the cell centers, including the       *
        !            * 1st level halo cells.                                    *
        !            *                                                          *
        !            ************************************************************
        !
        ! Loop over the cells, including the 1st level halo's.

        do k=1,ke
           do j=1,je
              do i=1,ie

                 ! The velocity of the cell center is determined
                 ! by a finite difference formula. First store
                 ! the current coordinate, multiplied by 8 and
                 ! coefTime(0) in sc.

                 sc(1) = (x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                      +  x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                      +  x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                      +  x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
                      * coefTime(0)
                 sc(2) = (x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                      +  x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                      +  x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                      +  x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
                      * coefTime(0)
                 sc(3) = (x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                      +  x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                      +  x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                      +  x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
                      * coefTime(0)

                 ! Loop over the older levels to complete the
                 ! finite difference formula.

                 do ii=1,nOldLevels
                    sc(1) = sc(1) + (xOld(ii,i-1,j-1,k-1,1)  &
                         +          xOld(ii,i,  j-1,k-1,1)  &
                         +          xOld(ii,i-1,j,  k-1,1)  &
                         +          xOld(ii,i,  j,  k-1,1)  &
                         +          xOld(ii,i-1,j-1,k,  1)  &
                         +          xOld(ii,i,  j-1,k,  1)  &
                         +          xOld(ii,i-1,j,  k,  1)  &
                         +          xOld(ii,i,  j,  k,  1)) &
                         * coefTime(ii)
                    sc(2) = sc(2) + (xOld(ii,i-1,j-1,k-1,2)  &
                         +          xOld(ii,i,  j-1,k-1,2)  &
                         +          xOld(ii,i-1,j,  k-1,2)  &
                         +          xOld(ii,i,  j,  k-1,2)  &
                         +          xOld(ii,i-1,j-1,k,  2)  &
                         +          xOld(ii,i,  j-1,k,  2)  &
                         +          xOld(ii,i-1,j,  k,  2)  &
                         +          xOld(ii,i,  j,  k,  2)) &
                         * coefTime(ii)
                    sc(3) = sc(3) + (xOld(ii,i-1,j-1,k-1,3)  &
                         +          xOld(ii,i,  j-1,k-1,3)  &
                         +          xOld(ii,i-1,j,  k-1,3)  &
                         +          xOld(ii,i,  j,  k-1,3)  &
                         +          xOld(ii,i-1,j-1,k,  3)  &
                         +          xOld(ii,i,  j-1,k,  3)  &
                         +          xOld(ii,i-1,j,  k,  3)  &
                         +          xOld(ii,i,  j,  k,  3)) &
                         * coefTime(ii)
                 enddo

                 ! Divide by 8 delta t to obtain the correct
                 ! velocities.

                 s(i,j,k,1) = sc(1)*oneOver8dt
                 s(i,j,k,2) = sc(2)*oneOver8dt
                 s(i,j,k,3) = sc(3)*oneOver8dt
              enddo
           enddo
        enddo
        !
        !            ************************************************************
        !            *                                                          *
        !            * Normal grid velocities of the faces.                     *
        !            *                                                          *
        !            ************************************************************
        !
        ! Loop over the three directions.

        loopDir: do mm=1,3

           ! Set the upper boundaries depending on the direction.

           select case (mm)
           case (1_intType)       ! normals in i-direction
              iie = ie; jje = je; kke = ke

           case (2_intType)       ! normals in j-direction
              iie = je; jje = ie; kke = ke

           case (3_intType)       ! normals in k-direction
              iie = ke; jje = ie; kke = je
           end select
           !
           !              **********************************************************
           !              *                                                        *
           !              * Normal grid velocities in generalized i-direction.     *
           !              * Mm == 1: i-direction                                   *
           !              * mm == 2: j-direction                                   *
           !              * mm == 3: k-direction                                   *
           !              *                                                        *
           !              **********************************************************
           !
           do i=0,iie

              ! Set the pointers for the coordinates, normals and
              ! normal velocities for this generalized i-plane.
              ! This depends on the value of mm.

              select case (mm)
              case (1_intType)       ! normals in i-direction
                 xx =>  x(i,:,:,:);  xxOld => xOld(:,i,:,:,:)
                 ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)

              case (2_intType)       ! normals in j-direction
                 xx =>  x(:,i,:,:);  xxOld => xOld(:,:,i,:,:)
                 ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)

              case (3_intType)       ! normals in k-direction
                 xx =>  x(:,:,i,:);  xxOld => xOld(:,:,:,i,:)
                 ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
              end select

              ! Loop over the k and j-direction of this
              ! generalized i-face. Note that due to the usage of
              ! the pointers xx and xxOld an offset of +1 must be
              ! used in the coordinate arrays, because x and xOld
              ! originally start at 0 for the i, j and k indices.

              do k=1,kke
                 do j=1,jje

                    ! The velocity of the face center is determined
                    ! by a finite difference formula. First store
                    ! the current coordinate, multiplied by 4 and
                    ! coefTime(0) in sc.

                    sc(1) = coefTime(0)*(xx(j+1,k+1,1) + xx(j,k+1,1) &
                         +              xx(j+1,k,  1) + xx(j,k,  1))
                    sc(2) = coefTime(0)*(xx(j+1,k+1,2) + xx(j,k+1,2) &
                         +              xx(j+1,k,  2) + xx(j,k,  2))
                    sc(3) = coefTime(0)*(xx(j+1,k+1,3) + xx(j,k+1,3) &
                         +              xx(j+1,k,  3) + xx(j,k,  3))

                    ! Loop over the older levels to complete the
                    ! finite difference.

                    do ii=1,nOldLevels

                       sc(1) = sc(1) + coefTime(ii)         &
                            *         (xxOld(ii,j+1,k+1,1) &
                            +          xxOld(ii,j,  k+1,1) &
                            +          xxOld(ii,j+1,k,  1) &
                            +          xxOld(ii,j,  k,  1))
                       sc(2) = sc(2) + coefTime(ii)         &
                            *         (xxOld(ii,j+1,k+1,2) &
                            +          xxOld(ii,j,  k+1,2) &
                            +          xxOld(ii,j+1,k,  2) &
                            +          xxOld(ii,j,  k,  2))
                       sc(3) = sc(3) + coefTime(ii)         &
                            *         (xxOld(ii,j+1,k+1,3) &
                            +          xxOld(ii,j,  k+1,3) &
                            +          xxOld(ii,j+1,k,  3) &
                            +          xxOld(ii,j,  k,  3))
                    enddo

                    ! Determine the dot product of sc and the normal
                    ! and divide by 4 deltaT to obtain the correct
                    ! value of the normal velocity.

                    sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
                         + sc(3)*ss(j,k,3)
                    sFace(j,k) = sFace(j,k)*oneOver4dt

                 enddo
              enddo
           enddo

        enddo loopDir

     else testUseOldCoor
        !
        !            ************************************************************
        !            *                                                          *
        !            * The velocities must be determined analytically.          *
        !            *                                                          *
        !            ************************************************************
        !
        ! Store the rotation center and determine the
        ! nonDimensional rotation rate of this block. As the
        ! reference length is 1 timeRef == 1/uRef and at the end
        ! the nonDimensional velocity is computed.

        j = nbkGlobal

        rotCenter = cgnsDoms(j)%rotCenter
        offSetVector= (rotCenter-rotPoint)
        rotRate   = timeRef*cgnsDoms(j)%rotRate


        if (useWindAxis)then
           !determine the current angles from the free stream velocity
           call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)

           if (liftIndex == 2) then
              ! different coordinate system for aerosurf
              ! Wing is in z- direction
              rotRateTrans(1,1)=cos(alpha)*cos(beta)
              rotRateTrans(1,2)=-sin(alpha)
              rotRateTrans(1,3)=-cos(alpha)*sin(beta)
              rotRateTrans(2,1)=sin(alpha)*cos(beta)
              rotRateTrans(2,2)=cos(alpha)
              rotRateTrans(2,3)=-sin(alpha)*sin(beta)
              rotRateTrans(3,1)=sin(beta)
              rotRateTrans(3,2)=0.0
              rotRateTrans(3,3)=cos(beta)

           elseif(liftIndex ==3) then
              ! Wing is in y- direction
              !Rotate the rotation rate from the wind axis back to the local body axis
              rotRateTrans(1,1)=cos(alpha)*cos(beta)
              rotRateTrans(1,2)=-cos(alpha)*sin(beta)
              rotRateTrans(1,3)=-sin(alpha)
              rotRateTrans(2,1)=sin(beta)
              rotRateTrans(2,2)=cos(beta)
              rotRateTrans(2,3)=0.0
              rotRateTrans(3,1)=sin(alpha)*cos(beta)
              rotRateTrans(3,2)=-sin(alpha)*sin(beta)
              rotRateTrans(3,3)=cos(alpha)
           else
              call terminate('getDirAngle', 'Invalid Lift Direction')
           endif

           rotRateTemp = rotRate
           rotRate=0.0
           do i=1,3
              do j=1,3
                 rotRate(i)=rotRate(i)+rotRateTemp(j)*rotRateTrans(i,j)
              end do
           end do
        end if


!!$             if (useWindAxis)then
!!$                !determine the current angles from the free stream velocity
!!$                call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)
!!$                !Rotate the rotation rate from the wind axis back to the local body axis
!!$                !checkt he relationship between the differnt degrees of freedom!
!!$                rotRateTrans(1,1)=cos(alpha)*cos(beta)
!!$                rotRateTrans(1,2)=-cos(alpha)*sin(beta)
!!$                rotRateTrans(1,3)=-sin(alpha)
!!$                rotRateTrans(2,1)=sin(beta)
!!$                rotRateTrans(2,2)=cos(beta)
!!$                rotRateTrans(2,3)=0.0
!!$                rotRateTrans(3,1)=sin(alpha)*cos(beta)
!!$                rotRateTrans(3,2)=-sin(alpha)*sin(beta)
!!$                rotRateTrans(3,3)=cos(alpha)
!!$
!!$                rotRateTemp = rotRate
!!$                rotRate=0.0
!!$                do i=1,3
!!$                   do j=1,3
!!$                      rotRate(i)=rotRate(i)+rotRateTemp(j)*rotRateTrans(i,j)
!!$                   end do
!!$                end do
!!$             end if

        !subtract off the rotational velocity of the center of the grid
        ! to account for the added overall velocity.
        !             velxGrid =velxgrid0+ 1*(rotRate(2)*rotCenter(3) - rotRate(3)*rotCenter(2))
        !             velyGrid =velygrid0+ 1*(rotRate(3)*rotCenter(1) - rotRate(1)*rotCenter(3))
        !             velzGrid =velzgrid0+ 1*(rotRate(1)*rotCenter(2) - rotRate(2)*rotCenter(1))

        velxGrid =velxgrid0+ 1*(rotRate(2)*offSetVector(3) &
             - rotRate(3)*offSetVector(2)) &
             + derivRotationMatrix(1,1)*offSetVector(1) &
             + derivRotationMatrix(1,2)*offSetVector(2) &
             + derivRotationMatrix(1,3)*offSetVector(3)
        velyGrid =velygrid0+ 1*(rotRate(3)*offSetVector(1)&
             - rotRate(1)*offSetVector(3))&
             + derivRotationMatrix(2,1)*offSetVector(1) &
             + derivRotationMatrix(2,2)*offSetVector(2) &
             + derivRotationMatrix(2,3)*offSetVector(3)
        velzGrid =velzgrid0+ 1*(rotRate(1)*offSetVector(2) &
             - rotRate(2)*offSetVector(1)) &
             + derivRotationMatrix(3,1)*offSetVector(1) &
             + derivRotationMatrix(3,2)*offSetVector(2) &
             + derivRotationMatrix(3,3)*offSetVector(3)

        !add in rotmatrix*rotpoint....

!
!            ************************************************************
!            *                                                          *
!            * Grid velocities of the cell centers, including the       *
!            * 1st level halo cells.                                    *
!            *                                                          *
!            ************************************************************
!
        ! Loop over the cells, including the 1st level halo's.

        do k=1,ke
           do j=1,je
              do i=1,ie

                 ! Determine the coordinates of the cell center,
                 ! which are stored in xc.

                 xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1) &
                      +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1) &
                      +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1) &
                      +         x(i-1,j,  k,  1) + x(i,j,  k,  1))
                 xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2) &
                      +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2) &
                      +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2) &
                      +         x(i-1,j,  k,  2) + x(i,j,  k,  2))
                 xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3) &
                      +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3) &
                      +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3) &
                      +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

                 ! Determine the coordinates relative to the
                 ! center of rotation.

                 xxc(1) = xc(1) - rotCenter(1)
                 xxc(2) = xc(2) - rotCenter(2)
                 xxc(3) = xc(3) - rotCenter(3)

                 ! Determine the rotation speed of the cell center,
                 ! which is omega*r.

                 sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                 sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                 sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                 ! Determine the coordinates relative to the
                 ! rigid body rotation point.

                 xxc(1) = xc(1) - rotationPoint(1)
                 xxc(2) = xc(2) - rotationPoint(2)
                 xxc(3) = xc(3) - rotationPoint(3)

                 ! Determine the total velocity of the cell center.
                 ! This is a combination of rotation speed of this
                 ! block and the entire rigid body rotation.

                 s(i,j,k,1) = sc(1) + velxGrid           &
                      + derivRotationMatrix(1,1)*xxc(1) &
                      + derivRotationMatrix(1,2)*xxc(2) &
                      + derivRotationMatrix(1,3)*xxc(3)
                 s(i,j,k,2) = sc(2) + velyGrid           &
                      + derivRotationMatrix(2,1)*xxc(1) &
                      + derivRotationMatrix(2,2)*xxc(2) &
                      + derivRotationMatrix(2,3)*xxc(3)
                 s(i,j,k,3) = sc(3) + velzGrid           &
                      + derivRotationMatrix(3,1)*xxc(1) &
                      + derivRotationMatrix(3,2)*xxc(2) &
                      + derivRotationMatrix(3,3)*xxc(3)
              enddo
           enddo
        enddo
!
!            ************************************************************
!            *                                                          *
!            * Normal grid velocities of the faces.                     *
!            *                                                          *
!            ************************************************************
!
        ! Loop over the three directions.

        loopDirection: do mm=1,3

           ! Set the upper boundaries depending on the direction.

           select case (mm)
           case (1_intType)       ! Normals in i-direction
              iie = ie; jje = je; kke = ke

           case (2_intType)       ! Normals in j-direction
              iie = je; jje = ie; kke = ke

           case (3_intType)       ! Normals in k-direction
              iie = ke; jje = ie; kke = je
           end select
!
!              **********************************************************
!              *                                                        *
!              * Normal grid velocities in generalized i-direction.     *
!              * mm == 1: i-direction                                   *
!              * mm == 2: j-direction                                   *
!              * mm == 3: k-direction                                   *
!              *                                                        *
!              **********************************************************
!
           do i=0,iie

              ! Set the pointers for the coordinates, normals and
              ! normal velocities for this generalized i-plane.
              ! This depends on the value of mm.

              select case (mm)
              case (1_intType)       ! normals in i-direction
                 xx =>  x(i,:,:,:)
                 ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)

              case (2_intType)       ! normals in j-direction
                 xx =>  x(:,i,:,:)
                 ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)

              case (3_intType)       ! normals in k-direction
                 xx =>  x(:,:,i,:)
                 ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
              end select

              ! Loop over the k and j-direction of this generalized
              ! i-face. Note that due to the usage of the pointer
              ! xx an offset of +1 must be used in the coordinate
              ! array, because x originally starts at 0 for the
              ! i, j and k indices.

              do k=1,kke
                 do j=1,jje

                    ! Determine the coordinates of the face center,
                    ! which are stored in xc.

                    xc(1) = fourth*(xx(j+1,k+1,1) + xx(j,k+1,1) &
                         +         xx(j+1,k,  1) + xx(j,k,  1))
                    xc(2) = fourth*(xx(j+1,k+1,2) + xx(j,k+1,2) &
                         +         xx(j+1,k,  2) + xx(j,k,  2))
                    xc(3) = fourth*(xx(j+1,k+1,3) + xx(j,k+1,3) &
                         +         xx(j+1,k,  3) + xx(j,k,  3))

                    ! Determine the coordinates relative to the
                    ! center of rotation.

                    xxc(1) = xc(1) - rotCenter(1)
                    xxc(2) = xc(2) - rotCenter(2)
                    xxc(3) = xc(3) - rotCenter(3)

                    ! Determine the rotation speed of the face center,
                    ! which is omega*r.

                    sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                    sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                    sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                    ! Determine the coordinates relative to the
                    ! rigid body rotation point.

                    xxc(1) = xc(1) - rotationPoint(1)
                    xxc(2) = xc(2) - rotationPoint(2)
                    xxc(3) = xc(3) - rotationPoint(3)

                    ! Determine the total velocity of the cell face.
                    ! This is a combination of rotation speed of this
                    ! block and the entire rigid body rotation.

                    sc(1) = sc(1) + velxGrid           &
                         + derivRotationMatrix(1,1)*xxc(1) &
                         + derivRotationMatrix(1,2)*xxc(2) &
                         + derivRotationMatrix(1,3)*xxc(3)
                    sc(2) = sc(2) + velyGrid           &
                         + derivRotationMatrix(2,1)*xxc(1) &
                         + derivRotationMatrix(2,2)*xxc(2) &
                         + derivRotationMatrix(2,3)*xxc(3)
                    sc(3) = sc(3) + velzGrid           &
                         + derivRotationMatrix(3,1)*xxc(1) &
                         + derivRotationMatrix(3,2)*xxc(2) &
                         + derivRotationMatrix(3,3)*xxc(3)

                    ! Store the dot product of grid velocity sc and
                    ! the normal ss in sFace.

                    sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
                         + sc(3)*ss(j,k,3)

                 enddo
              enddo
           enddo

        enddo loopDirection
     endif testUseOldCoor
  endif testMoving

end subroutine gridVelocitiesFineLevel_block
