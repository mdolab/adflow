subroutine slipVelocitiesFineLevel(useOldCoor, t, sps)
  !
  ! Shell function to call slipVelocitiesFineLevel on all blocks
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

     call slipVelocitiesFineLevel_block(useOldCoor, t, sps)

  end do domains

end subroutine slipVelocitiesFineLevel


!
!      ******************************************************************
!      *                                                                *
!      * File:          slipVelocities.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-12-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine slipVelocitiesFineLevel_block(useOldCoor, t, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * slipVelocitiesFineLevel computes the slip velocities for       *
  !      * viscous subfaces on all viscous boundaries on groundLevel for  *
  !      * the given spectral solution. If useOldCoor is .true. the       *
  !      * velocities are determined using the unsteady time integrator;  *
  !      * otherwise the analytic form is used.                           *
  !      *                                                                *
  !      ******************************************************************
  !
  use BCTypes
  use inputTimeSpectral
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
  integer(kind=intType) :: nn, mm, i, j, level

  real(kind=realType) :: oneOver4dt
  real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
  real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

  real(kind=realType), dimension(3) :: xc, xxc
  real(kind=realType), dimension(3) :: rotCenter, rotRate

  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix,&
       derivRotationMatrix

  real(kind=realType) :: tNew, tOld

  real(kind=realType), dimension(:,:,:),   pointer :: uSlip
  real(kind=realType), dimension(:,:,:),   pointer :: xFace
  real(kind=realType), dimension(:,:,:,:), pointer :: xFaceOld

  integer(kind=intType) :: liftIndex
  real(kind=realType) :: alpha,beta,intervalMach,alphaTS,alphaIncrement,&
       betaTS,betaIncrement
  real(kind=realType), dimension(3) ::velDir
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
  ! Determine the situation we are having here.

  testUseOldCoor: if( useOldCoor ) then

     ! The velocities must be determined via a finite difference
     ! formula using the coordinates of the old levels.

     ! Set the coefficients for the time integrator and store the
     ! inverse of the physical nonDimensional time step, divided
     ! by 4, a bit easier.

     call setCoefTimeIntegrator
     oneOver4dt = fourth*timeRef/deltaT

     ! Loop over the number of viscous subfaces.

     bocoLoop1: do mm=1,nViscBocos

        ! Set the pointer for uSlip to make the code more
        ! readable.

        uSlip => BCData(mm)%uSlip

        ! Determine the grid face on which the subface is located
        ! and set some variables accordingly.

        select case (BCFaceID(mm))

        case (iMin)
           xFace => x(1,:,:,:);  xFaceOld => xOld(:,1,:,:,:)

        case (iMax)
           xFace => x(il,:,:,:); xFaceOld => xOld(:,il,:,:,:)

        case (jMin)
           xFace => x(:,1,:,:);  xFaceOld => xOld(:,:,1,:,:)

        case (jMax)
           xFace => x(:,jl,:,:); xFaceOld => xOld(:,:,jl,:,:)

        case (kMin)
           xFace => x(:,:,1,:);  xFaceOld => xOld(:,:,:,1,:)

        case (kMax)
           xFace => x(:,:,kl,:); xFaceOld => xOld(:,:,:,kl,:)

        end select

        ! Some boundary faces have a different rotation speed than
        ! the corresponding block. This happens e.g. in the tip gap
        ! region of turboMachinary problems where the casing does
        ! not rotate. As the coordinate difference corresponds to
        ! the rotation rate of the block, a correction must be
        ! computed. Therefore compute the difference in rotation
        ! rate and store the rotation center a bit easier. Note that
        ! the rotation center of subface is taken, because if there
        ! is a difference in rotation rate this info for the subface
        ! must always be specified.

        j = nbkGlobal
        i = cgnsSubface(mm)

        rotCenter = cgnsDoms(j)%bocoInfo(i)%rotCenter
        rotRate   = timeRef*(cgnsDoms(j)%bocoInfo(i)%rotRate &
             -          cgnsDoms(j)%rotRate)

        ! Loop over the quadrilateral faces of the viscous subface.
        ! Note that due to the usage of the pointers xFace and
        ! xFaceOld an offset of +1 must be used in the coordinate
        ! arrays, because x and xOld originally start at 0 for the
        ! i, j and k indices.

        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd

              ! Determine the coordinates of the centroid of the
              ! face, multiplied by 4.

              xc(1) = xFace(i+1,j+1,1) + xFace(i+1,j,1) &
                   + xFace(i,  j+1,1) + xFace(i,  j,1)
              xc(2) = xFace(i+1,j+1,2) + xFace(i+1,j,2) &
                   + xFace(i,  j+1,2) + xFace(i,  j,2)
              xc(3) = xFace(i+1,j+1,3) + xFace(i+1,j,3) &
                   + xFace(i,  j+1,3) + xFace(i,  j,3)

              ! Multiply the sum of the 4 vertex coordinates with
              ! coefTime(0) to obtain the contribution for the
              ! current time level. The division by 4*deltaT will
              ! take place later. This is both more efficient and
              ! more accurate for extremely small time steps.

              uSlip(i,j,1) = coefTime(0)*xc(1)
              uSlip(i,j,2) = coefTime(0)*xc(2)
              uSlip(i,j,3) = coefTime(0)*xc(3)

              ! Loop over the older time levels and take their
              ! contribution into account.

              do level=1,nOldLevels

                 uSlip(i,j,1) = uSlip(i,j,1) + coefTime(level)      &
                      *          (xFaceOld(level,i+1,j+1,1) &
                      +           xFaceOld(level,i+1,j,  1) &
                      +           xFaceOld(level,i,  j+1,1) &
                      +           xFaceOld(level,i,  j,  1))

                 uSlip(i,j,2) = uSlip(i,j,2) + coefTime(level)      &
                      *          (xFaceOld(level,i+1,j+1,2) &
                      +           xFaceOld(level,i+1,j,  2) &
                      +           xFaceOld(level,i,  j+1,2) &
                      +           xFaceOld(level,i,  j,  2))

                 uSlip(i,j,3) = uSlip(i,j,3) + coefTime(level)      &
                      *          (xFaceOld(level,i+1,j+1,3) &
                      +           xFaceOld(level,i+1,j,  3) &
                      +           xFaceOld(level,i,  j+1,3) &
                      +           xFaceOld(level,i,  j,  3))
              enddo

              ! Divide by 4 times the time step to obtain the
              ! correct velocity.

              uSlip(i,j,1) = uSlip(i,j,1)*oneOver4dt
              uSlip(i,j,2) = uSlip(i,j,2)*oneOver4dt
              uSlip(i,j,3) = uSlip(i,j,3)*oneOver4dt

              ! Determine the correction due to the difference
              ! in rotation rate between the block and subface.

              ! First determine the coordinates relative to the
              ! rotation center. Remember that 4 times this value
              ! is currently stored in xc.

              xc(1) = fourth*xc(1) - rotCenter(1)
              xc(2) = fourth*xc(2) - rotCenter(2)
              xc(3) = fourth*xc(3) - rotCenter(3)

              ! Compute the velocity, which is the cross product
              ! of rotRate and xc and add it to uSlip.

              uSlip(i,j,1) = uSlip(i,j,1) &
                   + rotRate(2)*xc(3) - rotRate(3)*xc(2)
              uSlip(i,j,2) = uSlip(i,j,2) &
                   + rotRate(3)*xc(1) - rotRate(1)*xc(3)
              uSlip(i,j,3) = uSlip(i,j,3) &
                   + rotRate(1)*xc(2) - rotRate(2)*xc(1)

           enddo
        enddo

     enddo bocoLoop1

  else

     ! The velocities must be determined analytically.

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

     call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, &
          t(1))

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

     ! Loop over the number of viscous subfaces.

     bocoLoop2: do mm=1,nViscBocos

        ! Determine the grid face on which the subface is located
        ! and set some variables accordingly.

        select case (BCFaceID(mm))

        case (iMin)
           xFace => x(1,:,:,:)

        case (iMax)
           xFace => x(il,:,:,:)

        case (jMin)
           xFace => x(:,1,:,:)

        case (jMax)
           xFace => x(:,jl,:,:)

        case (kMin)
           xFace => x(:,:,1,:)

        case (kMax)
           xFace => x(:,:,kl,:)

        end select

        ! Store the rotation center and the rotation rate
        ! for this subface.

        j = nbkGlobal
        i = cgnsSubface(mm)

        rotCenter = cgnsDoms(j)%bocoInfo(i)%rotCenter
        rotRate   = timeRef*cgnsDoms(j)%bocoInfo(i)%rotRate

        ! useWindAxis should go back here!
        velXgrid = velXGrid0
        velYgrid = velYGrid0
        velZgrid = velZGrid0

        ! Loop over the quadrilateral faces of the viscous
        ! subface.

        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd

              ! Compute the coordinates of the centroid of the face.
              ! Normally this would be an average of i-1 and i, but
              ! due to the usage of the pointer xFace and the fact
              ! that x starts at index 0 this is shifted 1 index.

              xc(1) = fourth*(xFace(i+1,j+1,1) + xFace(i+1,j,1) &
                   +         xFace(i,  j+1,1) + xFace(i,  j,1))
              xc(2) = fourth*(xFace(i+1,j+1,2) + xFace(i+1,j,2) &
                   +         xFace(i,  j+1,2) + xFace(i,  j,2))
              xc(3) = fourth*(xFace(i+1,j+1,3) + xFace(i+1,j,3) &
                   +         xFace(i,  j+1,3) + xFace(i,  j,3))

              ! Determine the coordinates relative to the center
              ! of rotation.

              xxc(1) = xc(1) - rotCenter(1)
              xxc(2) = xc(2) - rotCenter(2)
              xxc(3) = xc(3) - rotCenter(3)

              ! Compute the velocity, which is the cross product
              ! of rotRate and xc.

              BCData(mm)%uSlip(i,j,1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
              BCData(mm)%uSlip(i,j,2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
              BCData(mm)%uSlip(i,j,3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

              ! Determine the coordinates relative to the
              ! rigid body rotation point.

              xxc(1) = xc(1) - rotationPoint(1)
              xxc(2) = xc(2) - rotationPoint(2)
              xxc(3) = xc(3) - rotationPoint(3)

              ! Determine the total velocity of the cell center.
              ! This is a combination of rotation speed of this
              ! block and the entire rigid body rotation.

              BCData(mm)% uSlip(i,j,1) = BCData(mm)%uSlip(i,j,1) + velxGrid    &
                   + derivRotationMatrix(1,1)*xxc(1) &
                   + derivRotationMatrix(1,2)*xxc(2) &
                   + derivRotationMatrix(1,3)*xxc(3)
              BCData(mm)%uSlip(i,j,2) = BCData(mm)%uSlip(i,j,2) + velyGrid    &
                   + derivRotationMatrix(2,1)*xxc(1) &
                   + derivRotationMatrix(2,2)*xxc(2) &
                   + derivRotationMatrix(2,3)*xxc(3)
              BCData(mm)%uSlip(i,j,3) = BCData(mm)%uSlip(i,j,3) + velzGrid    &
                   + derivRotationMatrix(3,1)*xxc(1) &
                   + derivRotationMatrix(3,2)*xxc(2) &
                   + derivRotationMatrix(3,3)*xxc(3)
           enddo
        enddo

     enddo bocoLoop2


  endif testUseOldCoor

end subroutine slipVelocitiesFineLevel_block

!      ==================================================================

       subroutine slipVelocitiesCoarseLevels(sps)
!
!      ******************************************************************
!      *                                                                *
!      * slipVelocitiesCoarseLevels determines the slip velocities      *
!      * for the given spectral solution starting from the known        *
!      * velocities on the finer level.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps
!
!      Local variables.
!
       integer(kind=intType) :: i, j, iiMax, jjMax
       integer(kind=intType) :: if1, if2, jf1, jf2
       integer(kind=intType) :: nLevels, level, levm1, nn, mm

       integer(kind=intType), dimension(:,:), pointer :: iFine, jFine

       real(kind=realType), dimension(:,:,:), pointer :: uSlip
       real(kind=realType), dimension(:,:,:), pointer :: uSlipFine
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of multigrid levels.

       nLevels = ubound(flowDoms,2)

       ! Loop over coarser grid levels, where ground level is considered
       ! as the finest grid.

       levelLoop: do level=(groundLevel+1),nLevels

         ! Set levm1 for the finer level.

         levm1 = level - 1

         ! Loop over the number of local blocks.

         domains: do nn=1,nDom

           ! Set the pointers to the coarse block.

           call setPointers(nn, level, sps)

           ! Loop over the number of viscous subfaces.

           bocoLoop: do mm=1,nViscBocos

             ! Set the pointers for uSlip and uSlipFine to make the
             ! code more readable.

             uSlip     => BCData(mm)%uSlip
             uSlipFine => flowDoms(nn,levm1,sps)%BCData(mm)%uSlip

             ! Determine the grid face on which the subface is located
             ! and set some variables accordingly.

             select case (BCFaceID(mm))

               case (iMin,iMax)
                 iiMax = jl; jjMax = kl
                 iFine => mgJFine; jFine => mgKFine

               case (jMin,jMax)
                 iiMax = il; jjMax = kl
                 iFine => mgIFine; jFine => mgKFine

               case (kMin,kMax)
                 iiMax = il; jjMax = jl
                 iFine => mgIFine; jFine => mgJFine

             end select

             ! Loop over the number of faces of the viscous subface.
             ! First in the generalized j-direction.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd

               ! Determine the two children in this direction.
               ! Take care of the halo's, as this info is only
               ! available for owned cells.

               if(j < 2) then
                 jf1 = 1; jf2 = 1
               else if(j > jjMax) then
                 jf1 = jFine(jjMax,2) +1; jf2 = jf1
               else
                 jf1 = jFine(j,1); jf2 = jFine(j,2)
               endif

               ! Loop in the generalized i-direction.

               do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ! Determine the two children in this direction.
                 ! Same story as in j-direction.

                 if(i < 2) then
                   if1 = 1; if2 = 1
                 else if(i > iiMax) then
                   if1 = iFine(iiMax,2) +1; if2 = if1
                 else
                   if1 = iFine(i,1); if2 = iFine(i,2)
                 endif

                 ! Average the fine grid velocities to the
                 ! coarse grid velocities.

                 uSlip(i,j,1) = fourth*(uSlipFine(if1,jf1,1) &
                              +         uSlipFine(if2,jf1,1) &
                              +         uSlipFine(if1,jf2,1) &
                              +         uSlipFine(if2,jf2,1))

                 uSlip(i,j,2) = fourth*(uSlipFine(if1,jf1,2) &
                              +         uSlipFine(if2,jf1,2) &
                              +         uSlipFine(if1,jf2,2) &
                              +         uSlipFine(if2,jf2,2))

                 uSlip(i,j,3) = fourth*(uSlipFine(if1,jf1,3) &
                              +         uSlipFine(if2,jf1,3) &
                              +         uSlipFine(if1,jf2,3) &
                              +         uSlipFine(if2,jf2,3))
               enddo
             enddo

           enddo bocoLoop
         enddo domains
       enddo levelLoop

       end subroutine slipVelocitiesCoarseLevels
