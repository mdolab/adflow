!
!      ******************************************************************
!      *                                                                *
!      * File:          gridVelocities.f90                              *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 10-22-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine gridVelocitiesFineLevelAdj(useOldCoor, t, sps,xAdj,&
            siAdj, sjAdj, skAdj,rotCenterAdj, rotRateAdj,sAdj,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj, iCell, jCell, kCell)
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
       real(kind=realType), dimension(3), intent(in) :: rotCenterAdj, rotRateAdj
       
       !real(kind=realType), dimension(:,:), intent(out) :: sFace
       real(kind=realType), dimension(-2:2,-2:2,-2:2,3),intent(out) :: sAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2), intent(out) ::sFaceIAdj,sFaceJAdj,sFaceKAdj

       !new ADjoint variables
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3), intent(in) :: xAdj
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3), intent(in) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(-3:2,-3:2,3)::xxAdj,ssAdj
       !real(kind=realType) :: volAdj
       !real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(out) :: normAdj
       integer(kind=intType), intent(in) :: iCell, jCell, kCell
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: i, j, k, ii, iie, jje, kke
       integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

       real(kind=realType) :: oneOver4dt, oneOver8dt
       real(kind=realType) :: velxGrid, velyGrid, velzGrid

       real(kind=realType), dimension(3) :: sc, xc, xxc
       

       real(kind=realType), dimension(3)   :: rotationPointAdj,rotPointAdj
       real(kind=realType), dimension(3,3) :: rotationMatrixAdj

       !real(kind=realType), dimension(:,:), pointer :: sFace
       real(kind=realType), dimension(-2:2,-2:2)::sFaceAdj 

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

       call derivativeRotMatrixRigidAdj(rotationMatrixAdj, rotationPointAdj,rotPointAdj, t(1))

       ! Loop over the number of local blocks.

       domains: do nn=1,nDom

         ! Set the pointers for this block.

         call setPointersAdj(nn, groundLevel, sps)

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

!
!            ************************************************************
!            *                                                          *
!            * Grid velocities of the cell centers, including the       *
!            * 1st level halo cells.                                    *
!            *                                                          *
!            ************************************************************
!
             ! Loop over the cells, including the 1st level halo's.
              kStart=-2; kEnd=2
              jStart=-2; jEnd=2
              iStart=-2; iEnd=2
              do k=kStart,kEnd !-2,2
                 do j=jStart,jEnd !-2,2
                    do i=iStart,iEnd !-2,2
!!$             do k=1,ke
!!$               do j=1,je
!!$                 do i=1,ie

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
                              + rotationMatrixAdj(1,1)*xxc(1) &
                              + rotationMatrixAdj(1,2)*xxc(2) &
                              + rotationMatrixAdj(1,3)*xxc(3)
                   sAdj(i,j,k,2) = sc(2) + velyGrid           &
                              + rotationMatrixAdj(2,1)*xxc(1) &
                              + rotationMatrixAdj(2,2)*xxc(2) &
                              + rotationMatrixAdj(2,3)*xxc(3)
                   sAdj(i,j,k,3) = sc(3) + velzGrid           &
                              + rotationMatrixAdj(3,1)*xxc(1) &
                              + rotationMatrixAdj(3,2)*xxc(2) &
                              + rotationMatrixAdj(3,3)*xxc(3)
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

                kStart=-2; kEnd=2
                jStart=-2; jEnd=2
                iStart=-2; iEnd=2
                
               ! Set the upper boundaries depending on the direction.
                
               select case (mm)
                 case (1_intType)       ! Normals in i-direction
                   iie = ie; jje = je; kke = ke

!!$                   
!!$                   if(iCell==il) iEnd=1
!!$                   if(iCell==2) iStart=-2
!!$                   
!!$                   if(jCell==2)  jStart=-1
!!$                   if(jCell==jl) jEnd=1 
!!$                   
!!$                   if(kCell==2) kStart=-1
!!$                   if(kCell==kl) kEnd=1

                 case (2_intType)       ! Normals in j-direction
                   iie = je; jje = ie; kke = ke
                   
!!$                   if(iCell==2)  iStart=-1
!!$                   if(iCell==il) iEnd=1
!!$                   
!!$                   if(jCell==jl) jEnd=1 
!!$                   if(jCell==2) jStart=-2
!!$                   
!!$                   if(kCell==2) kStart=-1
!!$                   if(kCell==kl) kEnd=1

                 case (3_intType)       ! Normals in k-direction
                   iie = ke; jje = ie; kke = je
                   
                   !if(iCell==2)  iStart=-1
                   !if(iCell==il) iEnd=1
                   
                   !if(jCell==2)  jStart=-1
                   !if(jCell==jl) jEnd=1 
                   
                   !if(kCell==kl) kEnd=1
                   !if(kCell==2)  kStart=-2

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
               !do i=0,iie
               do i=istart,iend

                 ! Set the pointers for the coordinates, normals and
                 ! normal velocities for this generalized i-plane.
                 ! This depends on the value of mm.

                 select case (mm)
                   case (1_intType)       ! normals in i-direction
                     xxAdj =  xAdj(i,:,:,:)
                     ssAdj = siAdj(i,:,:,:)

                   case (2_intType)       ! normals in j-direction
                     xxAdj =  xAdj(:,i,:,:)
                     ssAdj = sjAdj(:,i,:,:)

                   case (3_intType)       ! normals in k-direction
                     xxAdj =  xAdj(:,:,i,:)
                     ssAdj = skAdj(:,:,i,:)
                 end select

                 ! Loop over the k and j-direction of this generalized
                 ! i-face. Note that due to the usage of the pointer
                 ! xx an offset of +1 must be used in the coordinate
                 ! array, because x originally starts at 0 for the
                 ! i, j and k indices.

                 !do k=1,kke
                  ! do j=1,jje
                 do k=kstart,kend
                    do j=jstart,jend

                     ! Determine the coordinates of the face center,
                     ! which are stored in xc.

                     !xc(1) = fourth*(xxAdj(j+1,k+1,1) + xxAdj(j,k+1,1) &
                     !      +         xxAdj(j+1,k,  1) + xxAdj(j,k,  1))
                     !xc(2) = fourth*(xxAdj(j+1,k+1,2) + xxAdj(j,k+1,2) &
                     !      +         xxAdj(j+1,k,  2) + xxAdj(j,k,  2))
                     !xc(3) = fourth*(xxAdj(j+1,k+1,3) + xxAdj(j,k+1,3) &
                     !      +         xxAdj(j+1,k,  3) + xxAdj(j,k,  3))

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
                           + rotationMatrixAdj(1,1)*xxc(1) &
                           + rotationMatrixAdj(1,2)*xxc(2) &
                           + rotationMatrixAdj(1,3)*xxc(3)
                     sc(2) = sc(2) + velyGrid           &
                           + rotationMatrixAdj(2,1)*xxc(1) &
                           + rotationMatrixAdj(2,2)*xxc(2) &
                           + rotationMatrixAdj(2,3)*xxc(3)
                     sc(3) = sc(3) + velzGrid           &
                           + rotationMatrixAdj(3,1)*xxc(1) &
                           + rotationMatrixAdj(3,2)*xxc(2) &
                           + rotationMatrixAdj(3,3)*xxc(3)

                     ! Store the dot product of grid velocity sc and
                     ! the normal ss in sFace.

                     sFaceAdj(j,k) = sc(1)*ssAdj(j,k,1) + sc(2)*ssAdj(j,k,2) &
                                + sc(3)*ssAdj(j,k,3)

                   enddo
                 enddo
                 select case (mm)
                 case (1_intType)       ! normals in i-direction
                    sFaceIAdj(i,:,:) = sFaceAdj
                    
                 case (2_intType)       ! normals in j-direction
                    sFaceJAdj(:,i,:) = sFaceAdj
                    
                 case (3_intType)       ! normals in k-direction
                    sFaceKAdj(:,:,i) = sFaceAdj
                 end select
                 
               enddo

             enddo loopDirection
           endif testUseOldCoor
         endif testMoving
       enddo domains

     end subroutine gridVelocitiesFineLevelAdj
