subroutine gridVelocitiesFineLevelPart1(useOldCoor, t, sps)
  !
  ! Shell function to call gridVelocitiesFineLevel on all blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use iteration
  use utils, only : setPointers
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
     call gridVelocitiesFineLevelPart1_block(useOldCoor, t, sps)

  end do domains

end subroutine gridVelocitiesFineLevelPart1
!
!      ******************************************************************
!      *                                                                *
!      * File:          gridVelocities.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-23-2004                                     *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine gridVelocitiesFineLevelPart1_block(useOldCoor, t, sps)
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
  !      * Now it is split up into two parts.                             *
  !      * First part calculate the grid velocity using FIRST order BDF.  *
  !      * Second part calculate the surface normal and normal velocity.  *
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
  use utils, only : TSAlpha, TSBeta, TSMach, terminate, rotMatrixRigidBody
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

  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix,&
       derivRotationMatrix

  real(kind=realType) :: tNew, tOld
  real(kind=realType), dimension(:,:), pointer :: sFace
  real(kind=realType), dimension(:,:,:),   pointer :: sVelo

  real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
  real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

  integer(kind=intType) :: liftIndex
  real(kind=realType) :: alpha,beta,intervalMach,alphaTS,alphaIncrement,&
       betaTS,betaIncrement
  real(kind=realType), dimension(3) ::velDir
  real(kind=realType), dimension(3) :: refDirection
 
  ! Compute the mesh velocity from the given mesh Mach number.

  ! vel{x,y,z}Grid0 is the ACTUAL velocity you want at the
  ! geometry. 
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
     ! *******************************
     ! REMOVED the rigid body rotation part for simplicity
     ! *******************************
     
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

              ! *******************************
              ! Using FIRST order BDF for all cases
              ! Refer to eq. 11b, found paper by C.Farhat http://dx.doi.org/10.1016/S0021-9991(03)00311-5
              ! Same applies for the velocities of the faces below. Theta(n+1) = 1, Theta(n) = -1 therfore
              ! it becoms a first order scheme.
              ! *******************************

              ! The velocity of the cell center is determined
              ! by a finite difference formula. First store
              ! the current coordinate, multiplied by 8 and
              ! coefTime(0) in sc.

              sc(1) = (x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +  x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +  x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +  x(i-1,j,  k,  1) + x(i,j,  k,  1))
              sc(2) = (x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +  x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +  x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +  x(i-1,j,  k,  2) + x(i,j,  k,  2))
              sc(3) = (x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +  x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +  x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +  x(i-1,j,  k,  3) + x(i,j,  k,  3))

              ! Loop over the older levels to complete the
              ! finite difference formula.

              ii = 1 ! There was a loop over all old levels
              sc(1) = sc(1) + (xOld(ii,i-1,j-1,k-1,1)  &
                   +          xOld(ii,i,  j-1,k-1,1)  &
                   +          xOld(ii,i-1,j,  k-1,1)  &
                   +          xOld(ii,i,  j,  k-1,1)  &
                   +          xOld(ii,i-1,j-1,k,  1)  &
                   +          xOld(ii,i,  j-1,k,  1)  &
                   +          xOld(ii,i-1,j,  k,  1)  &
                   +          xOld(ii,i,  j,  k,  1)) &
                   * (-1.0_realType)
              sc(2) = sc(2) + (xOld(ii,i-1,j-1,k-1,2)  &
                   +          xOld(ii,i,  j-1,k-1,2)  &
                   +          xOld(ii,i-1,j,  k-1,2)  &
                   +          xOld(ii,i,  j,  k-1,2)  &
                   +          xOld(ii,i-1,j-1,k,  2)  &
                   +          xOld(ii,i,  j-1,k,  2)  &
                   +          xOld(ii,i-1,j,  k,  2)  &
                   +          xOld(ii,i,  j,  k,  2)) &
                   * (-1.0_realType)
              sc(3) = sc(3) + (xOld(ii,i-1,j-1,k-1,3)  &
                   +          xOld(ii,i,  j-1,k-1,3)  &
                   +          xOld(ii,i-1,j,  k-1,3)  &
                   +          xOld(ii,i,  j,  k-1,3)  &
                   +          xOld(ii,i-1,j-1,k,  3)  &
                   +          xOld(ii,i,  j-1,k,  3)  &
                   +          xOld(ii,i-1,j,  k,  3)  &
                   +          xOld(ii,i,  j,  k,  3)) &
                   * (-1.0_realType)

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
     !            * Velocities of the faces, vector.                         *
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
        !              * Face velocities in generalized i-direction.            *
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
              xx =>  x(i,:,:,:);  xxOld => xOld(:,i,:,:,:)
              sVelo => sVeloIALE(i,:,:,:)

           case (2_intType)       ! normals in j-direction
              xx =>  x(:,i,:,:);  xxOld => xOld(:,:,i,:,:)
              sVelo => sVeloJALE(:,i,:,:)

           case (3_intType)       ! normals in k-direction
              xx =>  x(:,:,i,:);  xxOld => xOld(:,:,:,i,:)
              sVelo => sVeloKALE(:,:,i,:)
           end select

           ! Loop over the k and j-direction of this
           ! generalized i-face. Note that due to the usage of
           ! the pointers xx and xxOld an offset of +1 must be
           ! used in the coordinate arrays, because x and xOld
           ! originally start at 0 for the i, j and k indices.
! print *, mm
           do k=1,kke
              do j=1,jje

                 ! The velocity of the face center is determined
                 ! by a finite difference formula. First store
                 ! the current coordinate, multiplied by 4 and
                 ! coefTime(0) in sc.

                 sc(1) = (xx(j+1,k+1,1) + xx(j,k+1,1) &
                      +  xx(j+1,k,  1) + xx(j,k,  1))
                 sc(2) = (xx(j+1,k+1,2) + xx(j,k+1,2) &
                      +  xx(j+1,k,  2) + xx(j,k,  2))
                 sc(3) = (xx(j+1,k+1,3) + xx(j,k+1,3) &
                      +  xx(j+1,k,  3) + xx(j,k,  3))

                 ii = 1 ! There was a loop who looped over nOldLevels
                 sc(1) = sc(1) + (xxOld(ii,j+1,k+1,1) &
                      +          xxOld(ii,j,  k+1,1) &
                      +          xxOld(ii,j+1,k,  1) &
                      +          xxOld(ii,j,  k,  1)) &
                      * (-1.0_realType)
                 sc(2) = sc(2) + (xxOld(ii,j+1,k+1,2) &
                      +          xxOld(ii,j,  k+1,2) &
                      +          xxOld(ii,j+1,k,  2) &
                      +          xxOld(ii,j,  k,  2)) &
                      * (-1.0_realType)
                 sc(3) = sc(3) + (xxOld(ii,j+1,k+1,3) &
                      +          xxOld(ii,j,  k+1,3) &
                      +          xxOld(ii,j+1,k,  3) &
                      +          xxOld(ii,j,  k,  3)) &
                      * (-1.0_realType)

                 ! Determine the dot product of sc and the normal
                 ! and divide by 4 deltaT to obtain the correct
                 ! value of the normal velocity.

                 sVelo(j,k,1) = sc(1)*oneOver4dt
                 sVelo(j,k,2) = sc(2)*oneOver4dt
                 sVelo(j,k,3) = sc(3)*oneOver4dt

                 ! if ((i.ge.2) .and. (i.le.3) .and. (j.ge.2) .and. (j.le.3) .and. (k.ge.2) .and. (k.le.3)) then
                 !    print *, i,j,k, sVelo(j,k,:)
                 !    print *, '                                   ', xx(j,k,:)
                 !    print *, '                                   ', xxOld(1,j,k,:)
                 ! end if

              enddo
           enddo
        enddo

     enddo loopDir

  endif testMoving

end subroutine gridVelocitiesFineLevelPart1_block

!
!      ******************************************************************
!      Here begins the second part
!      ******************************************************************
!

subroutine gridVelocitiesFineLevelPart2(useOldCoor, t, sps)
  !
  ! Shell function to call gridVelocitiesFineLevel on all blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use iteration
  use utils, only : setPointers
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
     call gridVelocitiesFineLevelPart2_block(useOldCoor, t, sps)

  end do domains

end subroutine gridVelocitiesFineLevelPart2

subroutine gridVelocitiesFineLevelPart2_block(useOldCoor, t, sps)
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
  real(kind=realType), dimension(3) :: sc, xc, xxc
  real(kind=realType), dimension(:,:), pointer :: sFace
  real(kind=realType), dimension(:,:,:),   pointer :: sVelo
  real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
  real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

  testMoving: if( blockIsMoving ) then
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
              ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)
              sVelo => sVeloIALE(i,:,:,:)

           case (2_intType)       ! normals in j-direction
              ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)
              sVelo => sVeloJALE(:,i,:,:)

           case (3_intType)       ! normals in k-direction
              ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
              sVelo => sVeloKALE(:,:,i,:)
           end select

           ! Loop over the k and j-direction of this
           ! generalized i-face. Note that due to the usage of
           ! the pointers xx and xxOld an offset of +1 must be
           ! used in the coordinate arrays, because x and xOld
           ! originally start at 0 for the i, j and k indices.

           do k=1,kke
              do j=1,jje

                 ! Determine the dot product of sc and the normal
                 ! and divide by 4 deltaT to obtain the correct
                 ! value of the normal velocity.

                 sFace(j,k) = sVelo(j,k,1)*ss(j,k,1) &
                            + sVelo(j,k,2)*ss(j,k,2) &
                            + sVelo(j,k,3)*ss(j,k,3)

              enddo
           enddo
        enddo

     enddo loopDir
  endif testMoving

end subroutine gridVelocitiesFineLevelPart2_block
