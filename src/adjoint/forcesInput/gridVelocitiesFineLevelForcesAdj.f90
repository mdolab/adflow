!
!      ******************************************************************
!      *                                                                *
!      * File:          gridVelocitiesFineLevelForcesAdj.f90            *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 10-22-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine gridVelocitiesFineLevelForcesAdj(useOldCoor, t, sps,xAdj,sAdj,&
            iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End,mm,&
            sFaceIAdj,sFaceJAdj,sFaceKAdj,&
            machGridAdj,velDirFreestreamAdj,&
            liftDirectionAdj,alphaAdj,betaAdj,liftindex,&
            rotPointAdj,pointRefAdj,&
            rotCenterAdj, rotRateAdj,siAdj,sjAdj,skAdj)
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
       use BCTypes     !imin,imax,etc.
       use inputTSStabDeriv
       use inputTimeSpectral !nTimeIntervalsSpectral
       use monitor           !Timeunsteady, timeunsteadyrestart
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!

       integer(kind=intType) :: iiBeg,iiEnd,jjBeg,jjEnd
       integer(kind=intType) :: i2Beg,i2End,j2Beg,j2End

       integer(kind=intType), intent(in) :: sps
       logical,               intent(in) :: useOldCoor

       real(kind=realType), dimension(*), intent(in) :: t
       real(kind=realType), dimension(3):: rotCenterAdj, rotRateAdj
       real(kind=realType), dimension(3):: pointRefAdj

       real(kind=realType), dimension(0:ie,0:je,0:ke,3),intent(in) :: xAdj
       real(kind=realType), dimension(0:ie,0:je,0:ke,3),intent(out) :: sAdj

       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3),intent(in) :: siAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3),intent(in) :: sjAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3),intent(in) :: skAdj

       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd) :: sFaceiAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd) :: sFacejAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2) :: sFacekAdj
       
       real(kind=realType),intent(in):: machGridAdj
       real(kind=realType),dimension(3),intent(in):: velDirFreestreamAdj,liftDirectionAdj
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: i, j, k, ii, iie, jje, kke

       real(kind=realType) :: oneOver4dt, oneOver8dt
       real(kind=realType) :: velxGrid, velyGrid, velzGrid,aInf
       real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

       real(kind=realType), dimension(3) :: sc, xc, xxc

       real(kind=realType), dimension(3)   :: rotationPointAdj,rotPointAdj
       real(kind=realType), dimension(3,3) :: derivrotationMatrixAdj,rotationMatrixAdj
       real(kind=realType), dimension(3) :: rotRateTemp
       real(kind=realType), dimension(3) :: offsetVector
       real(kind=realType), dimension(3,3) :: rotRateTrans
       real(kind=realType) :: alpha,beta
       
       real(kind=realType), dimension(iiBeg-1:iiEnd,jjBeg-1:jjEnd,3) :: xxAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: ssAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd) :: sFaceAdj
      
       real(kind=realType) ::tNew,tOld,intervalMach
       real(kind=realType), dimension(3)::liftDir,velDir,dragDir
       real(kind=realType)::alphaAdj,betaAdj
       real(kind=realType)::alphaIncrement,betaIncrement
       real(kind=realType)::alphaTS,betaTS
       integer(kind=intType) :: liftIndex

       !function definitions
       real(kind=realType) ::TSAlpha,TSBeta,TSMach
       !real(kind=realType), dimension(:,:), pointer :: sFace

!       real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
!       real(kind=realType), dimension(:,:,:,:), pointer :: xxOld
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the mesh velocity from the given mesh Mach number.
       !print *,'in gridvelocities',mm

    !  aInf = sqrt(gammaInf*pInf/rhoInf)
    !  velxGrid = aInf*MachGrid(1)
    !  velyGrid = aInf*MachGrid(2)
    !  velzGrid = aInf*MachGrid(3)

      ! velxGrid = zero
      ! velyGrid = zero
      ! velzGrid = zero

       aInf = sqrt(gammaInf*pInf/rhoInf)
       velxGrid0 = (aInf*machgridAdj)*(-velDirFreestreamAdj(1))
       velyGrid0 = (aInf*machgridAdj)*(-velDirFreestreamAdj(2))
       velzGrid0 = (aInf*machgridAdj)*(-velDirFreestreamAdj(3))


       ! Compute the derivative of the rotation matrix and the rotation
       ! point; needed for velocity due to the rigid body rotation of
       ! the entire grid. It is assumed that the rigid body motion of
       ! the grid is only specified if there is only 1 section present.

      

       call derivativeRotMatrixRigidForcesAdj(derivrotationMatrixAdj, rotationPointAdj,rotPointAdj, t(1))
       !compute the rotation matrix to update the velocities for the time
       !spectral stability derivative case...
       if(TSStability)then
          ! Determine the time values of the old and new time level.
          ! It is assumed that the rigid body rotation of the mesh is only
          ! used when only 1 section is present.
          
          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)
          !print *,'Time',t(1)
          if(TSpMode.or.TSqMode .or.TSrMode) then
             ! Compute the rotation matrix of the rigid body rotation as
             ! well as the rotation point; the latter may vary in time due
             ! to rigid body translation.
             
             call rotMatrixRigidBodyForcesAdj(tNew, tOld, rotationMatrixAdj, rotationPointAdj,rotPointAdj)
             velxgrid0 = rotationMatrixAdj(1,1)*velxgrid0 &
                  + rotationMatrixAdj(1,2)*velygrid0 &
                  + rotationMatrixAdj(1,3)*velzgrid0
             velygrid0 = rotationMatrixAdj(2,1)*velxgrid0 &
                  + rotationMatrixAdj(2,2)*velygrid0 &
                  + rotationMatrixAdj(2,3)*velzgrid0
             velzgrid0 = rotationMatrixAdj(3,1)*velxgrid0 &
                  + rotationMatrixAdj(3,2)*velygrid0 &
                  + rotationMatrixAdj(3,3)*velzgrid0
          elseif(tsAlphaMode)then
             !print*,'liftindex alpha',liftindex,velDirFreestreamAdj,liftDir,liftIndex
             ! get the baseline alpha and determine the liftIndex
!             call getDirAngleForces(velDirFreestreamAdj,liftDir,liftIndex,alpha,beta)
             !print *,'liftindex alpha',liftindex
             !Determine the alpha for this time instance
             alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
                             degreeFourAlpha,  omegaFourAlpha,     &
                             cosCoefFourAlpha, sinCoefFourAlpha, t(1))
             alphaTS = alphaAdj+alphaIncrement
             !Determine the grid velocity for this alpha
             call adjustInflowAngleForcesAdj(alphaTS,betaAdj,velDir,liftDir,dragDir,&
                  liftIndex)
             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgridAdj)*(-velDir(1))
             velyGrid0 = (aInf*machgridAdj)*(-velDir(2))
             velzGrid0 = (aInf*machgridAdj)*(-velDir(3))
             !print *,'base velocity',machgridadj, velxGrid0 , velyGrid0 , velzGrid0 

          elseif(tsBetaMode)then
             ! get the baseline alpha and determine the liftIndex
             !call getDirAngleForces(velDirFreestreamAdj,liftDir,liftIndex,alpha,beta)
             
             !Determine the alpha for this time instance
             betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
                             degreeFourBeta,  omegaFourBeta,     &
                             cosCoefFourBeta, sinCoefFourBeta, t(1))
             betaTS = betaAdj+betaIncrement
             !Determine the grid velocity for this alpha
             call adjustInflowAngleForcesAdj(alphaAdj,betaTS,velDir,liftDir,dragDir,&
                  liftIndex)
             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgridAdj)*(-velDir(1))
             velyGrid0 = (aInf*machgridAdj)*(-velDir(2))
             velzGrid0 = (aInf*machgridAdj)*(-velDir(3))
          elseif(TSMachMode)then
             !determine the mach number at this time interval
             IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
                             degreeFourMach,  omegaFourMach,     &
                             cosCoefFourMach, sinCoefFourMach, t(1))
             !set the effective grid velocity
             velxGrid0 = (aInf*(IntervalMach+machgridAdj))*(-velDirFreestreamAdj(1))
             velyGrid0 = (aInf*(IntervalMach+machgridAdj))*(-velDirFreestreamAdj(2))
             velzGrid0 = (aInf*(IntervalMach+machgridAdj))*(-velDirFreestreamAdj(3))
             
          elseif(TSAltitudeMode)then
             call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
          else
             call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
          end if
       endif
!moved outside
!!$       ! Loop over the number of local blocks.
!!$
!!$       domains: do nn=1,nDom
!!$
!!$         ! Set the pointers for this block.
!!$
!!$         call setPointers(nn, groundLevel, sps)

         ! Check for a moving block.

         testMoving: if( blockIsMoving ) then

           ! Determine the situation we are having here.

           testUseOldCoor: if( useOldCoor ) then
              
              print *,'arbitrary mesh movement not yet supported in the ADjoint'
              call terminate()
!!$!
!!$!            ************************************************************
!!$!            *                                                          *
!!$!            * The velocities must be determined via a finite           *
!!$!            * difference formula using the coordinates of the old      *
!!$!            * levels.                                                  *
!!$!            *                                                          *
!!$!            ************************************************************
!!$!
!!$             ! Set the coefficients for the time integrator and store
!!$             ! the inverse of the physical nonDimensional time step,
!!$             ! divided by 4 and 8, a bit easier.
!!$
!!$             call setCoefTimeIntegrator
!!$             oneOver4dt = fourth*timeRef/deltaT
!!$             oneOver8dt = half*oneOver4dt
!!$!
!!$!            ************************************************************
!!$!            *                                                          *
!!$!            * Grid velocities of the cell centers, including the       *
!!$!            * 1st level halo cells.                                    *
!!$!            *                                                          *
!!$!            ************************************************************
!!$!
!!$             ! Loop over the cells, including the 1st level halo's.
!!$
!!$             do k=1,ke
!!$               do j=1,je
!!$                 do i=1,ie
!!$
!!$                   ! The velocity of the cell center is determined
!!$                   ! by a finite difference formula. First store
!!$                   ! the current coordinate, multiplied by 8 and
!!$                   ! coefTime(0) in sc.
!!$
!!$                   sc(1) = (x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
!!$                         +  x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
!!$                         +  x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
!!$                         +  x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
!!$                         * coefTime(0)
!!$                   sc(2) = (x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
!!$                         +  x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
!!$                         +  x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
!!$                         +  x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
!!$                         * coefTime(0)
!!$                   sc(3) = (x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
!!$                         +  x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
!!$                         +  x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
!!$                         +  x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
!!$                         * coefTime(0)
!!$
!!$                   ! Loop over the older levels to complete the
!!$                   ! finite difference formula.
!!$
!!$                   do ii=1,nOldLevels
!!$                     sc(1) = sc(1) + (xOld(ii,i-1,j-1,k-1,1)  &
!!$                           +          xOld(ii,i,  j-1,k-1,1)  &
!!$                           +          xOld(ii,i-1,j,  k-1,1)  &
!!$                           +          xOld(ii,i,  j,  k-1,1)  &
!!$                           +          xOld(ii,i-1,j-1,k,  1)  &
!!$                           +          xOld(ii,i,  j-1,k,  1)  &
!!$                           +          xOld(ii,i-1,j,  k,  1)  &
!!$                           +          xOld(ii,i,  j,  k,  1)) &
!!$                           * coefTime(ii)
!!$                     sc(2) = sc(2) + (xOld(ii,i-1,j-1,k-1,2)  &
!!$                           +          xOld(ii,i,  j-1,k-1,2)  &
!!$                           +          xOld(ii,i-1,j,  k-1,2)  &
!!$                           +          xOld(ii,i,  j,  k-1,2)  &
!!$                           +          xOld(ii,i-1,j-1,k,  2)  &
!!$                           +          xOld(ii,i,  j-1,k,  2)  &
!!$                           +          xOld(ii,i-1,j,  k,  2)  &
!!$                           +          xOld(ii,i,  j,  k,  2)) &
!!$                           * coefTime(ii)
!!$                     sc(3) = sc(3) + (xOld(ii,i-1,j-1,k-1,3)  &
!!$                           +          xOld(ii,i,  j-1,k-1,3)  &
!!$                           +          xOld(ii,i-1,j,  k-1,3)  &
!!$                           +          xOld(ii,i,  j,  k-1,3)  &
!!$                           +          xOld(ii,i-1,j-1,k,  3)  &
!!$                           +          xOld(ii,i,  j-1,k,  3)  &
!!$                           +          xOld(ii,i-1,j,  k,  3)  &
!!$                           +          xOld(ii,i,  j,  k,  3)) &
!!$                           * coefTime(ii)
!!$                   enddo
!!$
!!$                   ! Divide by 8 delta t to obtain the correct
!!$                   ! velocities.
!!$
!!$                   s(i,j,k,1) = sc(1)*oneOver8dt
!!$                   s(i,j,k,2) = sc(2)*oneOver8dt
!!$                   s(i,j,k,3) = sc(3)*oneOver8dt
!!$
!!$                 enddo
!!$               enddo
!!$             enddo
!!$!
!!$!            ************************************************************
!!$!            *                                                          *
!!$!            * Normal grid velocities of the faces.                     *
!!$!            *                                                          *
!!$!            ************************************************************
!!$!
!!$             ! Loop over the three directions.
!!$
!!$             loopDir: do mm=1,3
!!$
!!$               ! Set the upper boundaries depending on the direction.
!!$
!!$               select case (mm)
!!$                 case (1_intType)       ! normals in i-direction
!!$                   iie = ie; jje = je; kke = ke
!!$
!!$                 case (2_intType)       ! normals in j-direction
!!$                   iie = je; jje = ie; kke = ke
!!$
!!$                 case (3_intType)       ! normals in k-direction
!!$                   iie = ke; jje = ie; kke = je
!!$               end select
!!$!
!!$!              **********************************************************
!!$!              *                                                        *
!!$!              * Normal grid velocities in generalized i-direction.     *
!!$!              * Mm == 1: i-direction                                   *
!!$!              * mm == 2: j-direction                                   *
!!$!              * mm == 3: k-direction                                   *
!!$!              *                                                        *
!!$!              **********************************************************
!!$!
!!$               do i=0,iie
!!$
!!$                 ! Set the pointers for the coordinates, normals and
!!$                 ! normal velocities for this generalized i-plane.
!!$                 ! This depends on the value of mm.
!!$
!!$                 select case (mm)
!!$                   case (1_intType)       ! normals in i-direction
!!$                     xx =>  x(i,:,:,:);  xxOld => xOld(:,i,:,:,:)
!!$                     ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)
!!$
!!$                   case (2_intType)       ! normals in j-direction
!!$                     xx =>  x(:,i,:,:);  xxOld => xOld(:,:,i,:,:)
!!$                     ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)
!!$
!!$                   case (3_intType)       ! normals in k-direction
!!$                     xx =>  x(:,:,i,:);  xxOld => xOld(:,:,:,i,:)
!!$                     ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
!!$                 end select
!!$
!!$                 ! Loop over the k and j-direction of this
!!$                 ! generalized i-face. Note that due to the usage of
!!$                 ! the pointers xx and xxOld an offset of +1 must be
!!$                 ! used in the coordinate arrays, because x and xOld
!!$                 ! originally start at 0 for the i, j and k indices.
!!$
!!$                 do k=1,kke
!!$                   do j=1,jje
!!$
!!$                     ! The velocity of the face center is determined
!!$                     ! by a finite difference formula. First store
!!$                     ! the current coordinate, multiplied by 4 and
!!$                     ! coefTime(0) in sc.
!!$
!!$                     sc(1) = coefTime(0)*(xx(j+1,k+1,1) + xx(j,k+1,1) &
!!$                           +              xx(j+1,k,  1) + xx(j,k,  1))
!!$                     sc(2) = coefTime(0)*(xx(j+1,k+1,2) + xx(j,k+1,2) &
!!$                           +              xx(j+1,k,  2) + xx(j,k,  2))
!!$                     sc(3) = coefTime(0)*(xx(j+1,k+1,3) + xx(j,k+1,3) &
!!$                           +              xx(j+1,k,  3) + xx(j,k,  3))
!!$
!!$                     ! Loop over the older levels to complete the
!!$                     ! finite difference.
!!$
!!$                     do ii=1,nOldLevels
!!$
!!$                       sc(1) = sc(1) + coefTime(ii)         &
!!$                             *         (xxOld(ii,j+1,k+1,1) &
!!$                             +          xxOld(ii,j,  k+1,1) &
!!$                             +          xxOld(ii,j+1,k,  1) &
!!$                             +          xxOld(ii,j,  k,  1))
!!$                       sc(2) = sc(2) + coefTime(ii)         &
!!$                             *         (xxOld(ii,j+1,k+1,2) &
!!$                             +          xxOld(ii,j,  k+1,2) &
!!$                             +          xxOld(ii,j+1,k,  2) &
!!$                             +          xxOld(ii,j,  k,  2))
!!$                       sc(3) = sc(3) + coefTime(ii)         &
!!$                             *         (xxOld(ii,j+1,k+1,3) &
!!$                             +          xxOld(ii,j,  k+1,3) &
!!$                             +          xxOld(ii,j+1,k,  3) &
!!$                             +          xxOld(ii,j,  k,  3))
!!$                     enddo
!!$
!!$                     ! Determine the dot product of sc and the normal
!!$                     ! and divide by 4 deltaT to obtain the correct
!!$                     ! value of the normal velocity.
!!$
!!$                     sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
!!$                                + sc(3)*ss(j,k,3)
!!$                     sFace(j,k) = sFace(j,k)*oneOver4dt
!!$
!!$                   enddo
!!$                 enddo
!!$               enddo
!!$
!!$             enddo loopDir

           else testUseOldCoor
!
!            ************************************************************
!            *                                                          *
!            * The velocities must be determined analytically.          *
!            *                                                          *
!            ************************************************************
!

             !!! Pass these in, set them in copyADjointStencil.f90
!!$
!!$             ! Store the rotation center and determine the
!!$             ! nonDimensional rotation rate of this block. As the
!!$             ! reference length is 1 timeRef == 1/uRef and at the end
!!$             ! the nonDimensional velocity is computed.
!!$
!!$             j = nbkGlobal
!!$
!!$             rotCenter = cgnsDoms(j)%rotCenter
!!$             rotRate   = timeRef*cgnsDoms(j)%rotRate
              if (useWindAxis)then
                 alpha =alphaAdj
                 beta= betaAdj
                 !Rotate the rotation rate from the wind axis back to the local body axis 
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
                    !Rotate the rotation rate from the winsd axis back to the local body axis
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
!!$                 rotRateTrans(1,1)=cos(alpha)*cos(beta)
!!$                 rotRateTrans(1,2)=-cos(alpha)*sin(beta)
!!$                 rotRateTrans(1,3)=-sin(alpha)
!!$                 rotRateTrans(2,1)=sin(beta)
!!$                 rotRateTrans(2,2)=cos(beta)
!!$                 rotRateTrans(2,3)=0.0
!!$                 rotRateTrans(3,1)=sin(alpha)*cos(beta)
!!$                 rotRateTrans(3,2)=-sin(alpha)*sin(beta)
!!$                 rotRateTrans(3,3)=cos(alpha)
                 
                 rotRateTemp = rotRateAdj
                 rotRateAdj=0.0
                 do i=1,3
                    do j=1,3
                       rotRateAdj(i)=rotRateAdj(i)+rotRateTemp(j)*rotRateTrans(i,j)
                    end do
                 end do
              end if

             offSetVector= (rotCenterAdj-pointRefAdj)
             !subtract off the rotational velocity of the center of the grid
             ! to account for the added overall velocity.
             velxGrid =velxgrid0+ 1*(rotRateAdj(2)*offSetVector(3)&
                                - rotRateAdj(3)*offSetVector(2)) &
                                + derivRotationMatrixAdj(1,1)*rotPointAdj(1) &
                                + derivRotationMatrixAdj(1,2)*rotPointAdj(2) &
                                + derivRotationMatrixAdj(1,3)*rotPointAdj(3)
             velyGrid =velygrid0+ 1*(rotRateAdj(3)*offSetVector(1)&
                                - rotRateAdj(1)*offSetVector(3))&
                              + derivRotationMatrixAdj(2,1)*rotPointAdj(1) &
                              + derivRotationMatrixAdj(2,2)*rotPointAdj(2) &
                              + derivRotationMatrixAdj(2,3)*rotPointAdj(3)
             velzGrid =velzgrid0+ 1*(rotRateAdj(1)*offSetVector(2)&
                                - rotRateAdj(2)*offSetVector(1)) &
                              + derivRotationMatrixAdj(3,1)*rotPointAdj(1) &
                              + derivRotationMatrixAdj(3,2)*rotPointAdj(2) &
                              + derivRotationMatrixAdj(3,3)*rotPointAdj(3)


!
!            ************************************************************
!            *                                                          *
!            * Grid velocities of the cell centers, including the       *
!            * 1st level halo cells.                                    *
!            *                                                          *
!            ************************************************************
!
             ! Loop over the cells, including the 1st level halo's.
             ! print *,'calculating grid velocities',ie,je,ke,shape(xadj)
             do k=1,ke
               do j=1,je
                 do i=1,ie

                   ! Determine the coordinates of the cell center,
                   ! which are stored in xc.

                   xc(1) = eighth*(xAdj(i-1,j-1,k-1,1) + xAdj(i,j-1,k-1,1) &
                         +         xAdj(i-1,j,  k-1,1) + xAdj(i,j,  k-1,1) &
                         +         xAdj(i-1,j-1,k,  1) + xAdj(i,j-1,k,  1) &
                         +         xAdj(i-1,j,  k,  1) + xAdj(i,j,  k,  1))
                   xc(2) = eighth*(xAdj(i-1,j-1,k-1,2) + xAdj(i,j-1,k-1,2) &
                         +         xAdj(i-1,j,  k-1,2) + xAdj(i,j,  k-1,2) &
                         +         xAdj(i-1,j-1,k,  2) + xAdj(i,j-1,k,  2) &
                         +         xAdj(i-1,j,  k,  2) + xAdj(i,j,  k,  2))
                   xc(3) = eighth*(xAdj(i-1,j-1,k-1,3) + xAdj(i,j-1,k-1,3) &
                         +         xAdj(i-1,j,  k-1,3) + xAdj(i,j,  k-1,3) &
                         +         xAdj(i-1,j-1,k,  3) + xAdj(i,j-1,k,  3) &
                         +         xAdj(i-1,j,  k,  3) + xAdj(i,j,  k,  3))

                   ! Determine the coordinates relative to the
                   ! center of rotation.

                   xxc(1) = xc(1) - rotCenterAdj(1)
                   xxc(2) = xc(2) - rotCenterAdj(2)
                   xxc(3) = xc(3) - rotCenterAdj(3)

                   ! Determine the rotation speed of the cell center,
                   ! which is omega*r.

                   sc(1) = rotRateAdj(2)*xxc(3) - rotRateAdj(3)*xxc(2)
                   sc(2) = rotRateAdj(3)*xxc(1) - rotRateAdj(1)*xxc(3)
                   sc(3) = rotRateAdj(1)*xxc(2) - rotRateAdj(2)*xxc(1)

                   ! Determine the coordinates relative to the
                   ! rigid body rotation point.

                   xxc(1) = xc(1) - rotationPointAdj(1)
                   xxc(2) = xc(2) - rotationPointAdj(2)
                   xxc(3) = xc(3) - rotationPointAdj(3)

                   ! Determine the total velocity of the cell center.
                   ! This is a combination of rotation speed of this
                   ! block and the entire rigid body rotation.

                   sAdj(i,j,k,1) = sc(1) + velxGrid           &
                              + derivrotationMatrixAdj(1,1)*xxc(1) &
                              + derivrotationMatrixAdj(1,2)*xxc(2) &
                              + derivrotationMatrixAdj(1,3)*xxc(3)
                   sAdj(i,j,k,2) = sc(2) + velyGrid           &
                              + derivrotationMatrixAdj(2,1)*xxc(1) &
                              + derivrotationMatrixAdj(2,2)*xxc(2) &
                              + derivrotationMatrixAdj(2,3)*xxc(3)
                   sAdj(i,j,k,3) = sc(3) + velzGrid           &
                              + derivrotationMatrixAdj(3,1)*xxc(1) &
                              + derivrotationMatrixAdj(3,2)*xxc(2) &
                              + derivrotationMatrixAdj(3,3)*xxc(3)
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
! this loop is equivalent to the BC loop. Thus mm is determined by the global mm.
!!$             ! Loop over the three directions.
!!$
!!$             loopDirection: do mm=1,3
             
             

               ! Set the upper boundaries depending on the direction.
             !print *,'mm1',mm
               select case (BCFaceID(mm))
               !case (1_intType)
               case (imin,imax)   ! Normals in i-direction
                  iie = ie; jje = je; kke = ke
                  
                  !case (2_intType)
               case (jmin,jmax)       ! Normals in j-direction
                  iie = je; jje = ie; kke = ke
                  
               !case (3_intType)
               case (kmin,kmax)      ! Normals in k-direction
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
               ii=1
               !do i=0,iie
               do i=1,iie,iie-1 ! 1 to iie in steps of iie-1
                 ! print *,'iie',i,iie,mm,jje,kke
                 ! Set the pointers for the coordinates, normals and
                 ! normal velocities for this generalized i-plane.
                 ! This depends on the value of mm.
                 !print *,'xindices',iibeg,iiend,0,ie
                 select case (BCFaceID(mm))
                 case (imin,imax)       ! normals in i-direction
                      !print *,'case1',i,shape(xxadj),shape(xadj)
                     !xxAdj =  xAdj(i,:,:,:)
                     xxAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,:)=xAdj(i,iiBeg-1:iiEnd,jjBeg-1:jjEnd,:)
                     ssAdj = siAdj(ii,:,:,:)!;  sFaceAdj = sFaceIAdj(i,:,:)

                  case (jmin,jmax)      ! normals in j-direction
                      !print *,'i2',i,shape(xxadj),shape(xadj)
                     xxAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,:) =  xAdj(iiBeg-1:iiEnd,i,jjBeg-1:jjEnd,:)
                     !print *,'indices',iiBeg,iiEnd,jjBeg,jjEnd
                     ssAdj = sjAdj(:,ii,:,:)!;  sFaceAdj = sFaceJAdj(:,i,:)

                  case (kmin,kmax)       ! normals in k-direction
                     !print *,'i3',i,shape(xxadj),shape(xadj)
                     xxAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,:)=xAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,i,:)
                     !xxAdj =  xAdj(:,:,i,:)
                     ssAdj = skAdj(:,:,ii,:)!;  sFaceAdj = sFaceKAdj(:,:,i)
                 end select
               
                 ! Loop over the k and j-direction of this generalized
                 ! i-face. Note that due to the usage of the pointer
                 ! xx an offset of +1 must be used in the coordinate
                 ! array, because x originally starts at 0 for the
                 ! i, j and k indices.

                 do k=jjBeg,jjEnd!1,kke
                   do j=iiBeg,iiEnd!1,jje
                      !print *,'j',j,jje,'k',k,kke
                     ! Determine the coordinates of the face center,
                     ! which are stored in xc.

!!$                     xc(1) = fourth*(xxAdj(j+1,k+1,1) + xxAdj(j,k+1,1) &
!!$                           +         xxAdj(j+1,k,  1) + xxAdj(j,k,  1))
!!$                     xc(2) = fourth*(xxAdj(j+1,k+1,2) + xxAdj(j,k+1,2) &
!!$                           +         xxAdj(j+1,k,  2) + xxAdj(j,k,  2))
!!$                     xc(3) = fourth*(xxAdj(j+1,k+1,3) + xxAdj(j,k+1,3) &
!!$                           +         xxAdj(j+1,k,  3) + xxAdj(j,k,  3))
                     xc(1) = fourth*(xxAdj(j,k,1) + xxAdj(j-1,k,1) &
                           +         xxAdj(j,k-1,  1) + xxAdj(j-1,k-1,  1))
                     xc(2) = fourth*(xxAdj(j,k,2) + xxAdj(j-1,k,2) &
                           +         xxAdj(j,k-1,  2) + xxAdj(j-1,k-1,  2))
                     xc(3) = fourth*(xxAdj(j,k,3) + xxAdj(j-1,k,3) &
                           +         xxAdj(j,k-1,  3) + xxAdj(j-1,k-1,  3))

                     ! Determine the coordinates relative to the
                     ! center of rotation.

                     xxc(1) = xc(1) - rotCenterAdj(1)
                     xxc(2) = xc(2) - rotCenterAdj(2)
                     xxc(3) = xc(3) - rotCenterAdj(3)

                     ! Determine the rotation speed of the face center,
                     ! which is omega*r.

                     sc(1) = rotRateAdj(2)*xxc(3) - rotRateAdj(3)*xxc(2)
                     sc(2) = rotRateAdj(3)*xxc(1) - rotRateAdj(1)*xxc(3)
                     sc(3) = rotRateAdj(1)*xxc(2) - rotRateAdj(2)*xxc(1)

                     ! Determine the coordinates relative to the
                     ! rigid body rotation point.

                     xxc(1) = xc(1) - rotationPointAdj(1)
                     xxc(2) = xc(2) - rotationPointAdj(2)
                     xxc(3) = xc(3) - rotationPointAdj(3)

                     ! Determine the total velocity of the cell face.
                     ! This is a combination of rotation speed of this
                     ! block and the entire rigid body rotation.

                     sc(1) = sc(1) + velxGrid           &
                           + derivrotationMatrixAdj(1,1)*xxc(1) &
                           + derivrotationMatrixAdj(1,2)*xxc(2) &
                           + derivrotationMatrixAdj(1,3)*xxc(3)
                     sc(2) = sc(2) + velyGrid           &
                           + derivrotationMatrixAdj(2,1)*xxc(1) &
                           + derivrotationMatrixAdj(2,2)*xxc(2) &
                           + derivrotationMatrixAdj(2,3)*xxc(3)
                     sc(3) = sc(3) + velzGrid           &
                           + derivrotationMatrixAdj(3,1)*xxc(1) &
                           + derivrotationMatrixAdj(3,2)*xxc(2) &
                           + derivrotationMatrixAdj(3,3)*xxc(3)

                     ! Store the dot product of grid velocity sc and
                     ! the normal ss in sFace.

                     sFaceAdj(j,k) = sc(1)*ssAdj(j,k,1) + sc(2)*ssAdj(j,k,2) &
                                + sc(3)*ssAdj(j,k,3)

                   enddo
                 enddo

                 select case (BCFaceID(mm))
                 case (imin,imax)       ! normals in i-direction
                    sFaceIAdj(ii,:,:) =  sFaceAdj
                 case (jmin,jmax)       ! normals in j-direction
                    sFaceJAdj(:,ii,:) = sFaceAdj
                 case (kmin,kmax)       ! normals in k-direction
                    sFaceKAdj(:,:,ii )= sFaceAdj 
                 end select

                 !increment counter for s(i,j,k)adj, sFace(i,j,k)Adj
                 ii=ii+1
                 
               enddo

  !           enddo loopDirection
           endif testUseOldCoor
         endif testMoving
 !      enddo domains
        ! print *,'finished grid velocities'
       end subroutine gridVelocitiesFineLevelForcesAdj
