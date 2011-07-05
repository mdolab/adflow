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
       subroutine gridVelocitiesFineLevel(useOldCoor, t, sps)
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
       real(kind=realType) :: velxGrid, velyGrid, velzGrid

       real(kind=realType), dimension(3) :: sc, xc, xxc
       real(kind=realType), dimension(3) :: rotCenter, rotRate

       real(kind=realType), dimension(3)   :: rotationPoint
       real(kind=realType), dimension(3,3) :: rotationMatrix

       real(kind=realType), dimension(:,:), pointer :: sFace

       real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
       real(kind=realType), dimension(:,:,:,:), pointer :: xxOld
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

       velxGrid = zero
       velyGrid = zero
       velzGrid = zero

       ! Compute the derivative of the rotation matrix and the rotation
       ! point; needed for velocity due to the rigid body rotation of
       ! the entire grid. It is assumed that the rigid body motion of
       ! the grid is only specified if there is only 1 section present.

       call derivativeRotMatrixRigid(rotationMatrix, rotationPoint, t(1))

       ! Loop over the number of local blocks.

       domains: do nn=1,nDom

         ! Set the pointers for this block.

         call setPointers(nn, groundLevel, sps)

         ! Check for a moving block.

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
             rotRate   = timeRef*cgnsDoms(j)%rotRate
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
                              + rotationMatrix(1,1)*xxc(1) &
                              + rotationMatrix(1,2)*xxc(2) &
                              + rotationMatrix(1,3)*xxc(3)
                   s(i,j,k,2) = sc(2) + velyGrid           &
                              + rotationMatrix(2,1)*xxc(1) &
                              + rotationMatrix(2,2)*xxc(2) &
                              + rotationMatrix(2,3)*xxc(3)
                   s(i,j,k,3) = sc(3) + velzGrid           &
                              + rotationMatrix(3,1)*xxc(1) &
                              + rotationMatrix(3,2)*xxc(2) &
                              + rotationMatrix(3,3)*xxc(3)
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
                           + rotationMatrix(1,1)*xxc(1) &
                           + rotationMatrix(1,2)*xxc(2) &
                           + rotationMatrix(1,3)*xxc(3)
                     sc(2) = sc(2) + velyGrid           &
                           + rotationMatrix(2,1)*xxc(1) &
                           + rotationMatrix(2,2)*xxc(2) &
                           + rotationMatrix(2,3)*xxc(3)
                     sc(3) = sc(3) + velzGrid           &
                           + rotationMatrix(3,1)*xxc(1) &
                           + rotationMatrix(3,2)*xxc(2) &
                           + rotationMatrix(3,3)*xxc(3)

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
       enddo domains

       end subroutine gridVelocitiesFineLevel
