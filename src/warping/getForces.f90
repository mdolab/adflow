subroutine getForceSize(size)

  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use BCTypes
  use blockPointers
  
  implicit none

  integer(kind=intType),intent(out) :: size
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd
  size = 0_intType
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
        end if
     end do bocos
  end do domains
end subroutine getForceSize

! subroutine getForces()

!   use BCTypes
!   use blockPointers
!   use cgnsGrid
!   use communication
!   use flowVarRefState
!   use warpingPetsc
!   use block
!   implicit none
!   !
!   !      Local variables.
!   !
!   integer :: ierr

!   integer(kind=intType) :: mm, nn, i, j
!   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

!   real(kind=realType) :: scaleDim, fact, pp, f(3)

!   real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
!   real(kind=realType), dimension(:,:,:), pointer :: ss
!   integer(kind=intType) :: iset(3)

!   !      ******************************************************************
!   !      *                                                                *
!   !      * Begin execution                                                *
!   !      *                                                                *
!   !      ******************************************************************

!   ! Compute the scaling factor to create the correct dimensional
!   ! force in newton. As the coordinates are already in meters,
!   ! this scaling factor is pRef.

!   scaleDim = pRef

!   ! Compute the local forces. Take the scaling factor into
!   ! account to obtain the forces in SI-units, i.e. Newton.

!   call vecZeroEntries(sumbGridVec,ierr)
!   call vecZeroEntries(cgnsGridVec,ierr)

!   domains: do nn=1,nDom
!      call setPointers(nn,1_intType,1_intType)

!      ! Loop over the number of boundary subfaces of this block.
!      bocos: do mm=1,nBocos

!         if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
!              BCType(mm) == NSWallIsothermal) then
!            select case (BCFaceID(mm))

!            case (iMin)
!               pp2 => p( 2,1:,1:); pp1 => p( 1,1:,1:); ss => si0( 1,:,:,:)
!               fact = -one

!            case (iMax)
!               pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:); ss => si0(il,:,:,:)
!               fact = one

!            case (jMin)
!               pp2 => p(1:, 2,1:); pp1 => p(1:, 1,1:); ss => sj0(:, 1,:,:)
!               fact = -one

!            case (jMax)
!               pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:); ss => sj0(:,jl,:,:)
!               fact = one

!            case (kMin)
!               pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1); ss => sk0(:,:, 1,:)
!               fact = -one

!            case (kMax)
!               pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke); ss => sk0(:,:,kl,:)
!               fact = one

!            end select

!            ! Store the cell range of the subfaces a bit easier.
!            ! As only owned faces must be considered the nodal range
!            ! in BCData must be used to obtain this data.

!            jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
!            iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

!            ! Compute the inviscid force on each of the faces and
!            ! scatter it to the 4 nodes, whose forces are updated.

!            do j=jBeg, jEnd ! This is a face loop
!               do i=iBeg, iEnd ! This is a face loop 

!                  ! Compute the pressure in the center of the boundary
!                  ! face, which is an average between pp2 and pp1. The
!                  ! value of pp is multiplied by 1/4 (the factor to
!                  ! scatter to its 4 nodes, by scaleDim (to obtain the
!                  ! correct dimensional value) and by fact (which takes
!                  ! the possibility of inward or outward pointing normal
!                  ! into account).

!                  pp = half*(pp2(i,j) + pp1(i,j))-Pinf
!                  pp = fourth*fact*scaleDim*pp

!                  ! Compute the corresponding force.

!                  f(1) = pp*ss(i,j,1)
!                  f(2) = pp*ss(i,j,2)
!                  f(3) = pp*ss(i,j,3)

!                  ! Now its just an indexing problem...we need to
!                  ! set each of the 4 entries in the sumbGridVec
!                  ! with AddValues

!                  select case(BCFaceID(mm))
!                  case(imin)
!                     call getIset(nn,1  ,i-1,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,1  ,i  ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,1  ,i-1,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,1  ,i  ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  case(imax)
!                     call getIset(nn,il ,i-1,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,il ,i  ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,il ,i-1,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,il ,i  ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  case(jmin)
!                     call getIset(nn,i-1,1  ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,1  ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i-1,1  ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,1  ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  case(jmax)
!                     call getIset(nn,i-1,jl ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,jl ,j-1,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i-1,jl ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,jl ,j  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  case(kmin)
!                     call getIset(nn,i-1,j-1,1  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,j-1,1  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i-1,j  ,1  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,j  ,1  ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  case(kmax)
!                     call getIset(nn,i-1,j-1,kl ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i  ,j-1,kl ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                     call getIset(nn,i-1,j  ,kl ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)          
!                     call getIset(nn,i  ,j  ,kl ,il,jl,kl,iset)
!                     call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
!                  end select
!               end do
!            end do
!         end if
!      end do bocos
!   end do domains

!   call VecAssemblyBegin(sumbGridVec,ierr)
!   call VecAssemblyEnd(sumbGridVec,ierr)

!   call VecScatterBegin(sumbTOcgnsForce,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
!   call VecScatterEnd  (sumbTOcgnsForce,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)

! end subroutine getForces


subroutine getForces2(forces,pts,npts)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use warpingPetsc
  use communication
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts
  real(kind=realType), intent(in)  :: pts(3,npts)
  real(kind=realType), intent(out) :: forces(3,npts)

  integer :: ierr

  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  real(kind=realType) :: scaleDim, fact, fact2, pp, fx,fy,fz

  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType) :: ss(3),v2(3),v1(3)
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  ! Compute the scaling factor to create the correct dimensional
  ! force in newton. As the coordinates are already in meters,
  ! this scaling factor is pRef.

  scaleDim = pRef

  ! Compute the local forces. Take the scaling factor into
  ! account to obtain the forces in SI-units, i.e. Newton.
  forces = 0.0
  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     if (flowDoms(nn,1_intType,1_intType)%rightHanded) then
        fact2 = half
     else
        fact2 = -half
     end if

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           select case (BCFaceID(mm))

           case (iMin)
              pp2 => p( 2,1:,1:); pp1 => p( 1,1:,1:)
              fact = -one

           case (iMax)
              pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:)
              fact = one

           case (jMin)
              pp2 => p(1:, 2,1:); pp1 => p(1:, 1,1:)
              fact = -one

           case (jMax)
              pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:)
              fact = one

           case (kMin)
              pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1)
              fact = -one

           case (kMax)
              pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke)
              fact = one

           end select

           ! Store the cell range of the subfaces a bit easier.
           ! As only owned faces must be considered the nodal range
           ! in BCData must be used to obtain this data.

           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

           ! Compute the inviscid force on each of the faces and
           ! scatter it to the 4 nodes, whose forces are updated.

           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 

                 ! Compute the pressure in the center of the boundary
                 ! face, which is an average between pp2 and pp1. The
                 ! value of pp is multiplied by 1/4 (the factor to
                 ! scatter to its 4 nodes, by scaleDim (to obtain the
                 ! correct dimensional value) and by fact (which takes
                 ! the possibility of inward or outward pointing normal
                 ! into account).

                 pp = half*(pp2(i,j) + pp1(i,j))-Pinf
                 pp = fourth*fact*scaleDim*pp

                 ! Compute Normal
                                  
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 v1(:) = pts(:,upper_right)-pts(:,lower_left)
                 v2(:) = pts(:,upper_left)-pts(:,lower_right)
                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact2 is
                 ! either -0.5 or 0.5.

                 ss(1) = fact2*(v1(2)*v2(3) - v1(3)*v2(2))
                 ss(2) = fact2*(v1(3)*v2(1) - v1(1)*v2(3))
                 ss(3) = fact2*(v1(1)*v2(2) - v1(2)*v2(1))

                 fx = pp*ss(1)
                 fy = pp*ss(2)
                 fz = pp*ss(3)

                 forces(1,lower_left)  = forces(1,lower_left)  + fx
                 forces(2,lower_left)  = forces(2,lower_left)  + fy
                 forces(3,lower_left)  = forces(3,lower_left)  + fz

                 forces(1,lower_right) = forces(1,lower_right) + fx
                 forces(2,lower_right) = forces(2,lower_right) + fy
                 forces(3,lower_right) = forces(3,lower_right) + fz

                 forces(1,upper_left)  = forces(1,upper_left)  + fx
                 forces(2,upper_left)  = forces(2,upper_left)  + fy
                 forces(3,upper_left)  = forces(3,upper_left)  + fz

                 forces(1,upper_right) = forces(1,upper_right) + fx
                 forces(2,upper_right) = forces(2,upper_right) + fy
                 forces(3,upper_right) = forces(3,upper_right) + fz

             !     if (myid == 0 .and. ii ==0 .and. i == 2 .and. j == 2) then
!                     print *,'getForces,ss,pp:',i,j
!                     print *,ss
!                     print *,pp
!                     print *,lower_left,lower_right,upper_left,upper_right
!                     print *,fx,fy,fz
!                  end if

              end do
           end do
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)

        end if
     end do bocos
  end do domains
end subroutine getForces2

subroutine getForcePoints(points,npts)

  use BCTypes
  use blockPointers

  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts
  real(kind=realType), intent(out) :: points(3,npts)
  integer :: ierr

  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
         
           ! NODE Based
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           ! Compute the inviscid force on each of the faces and
           ! scatter it to the 4 nodes, whose forces are updated.

           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 ii = ii +1
                 select case(BCFaceID(mm))

                 case(imin)
                    points(:,ii) = x(1,i,j,:)
                 case(imax)
                    points(:,ii) = x(il,i,j,:)
                 case(jmin) 
                    points(:,ii) = x(i,1,j,:)
                 case(jmax) 
                    points(:,ii) = x(i,jl,j,:)
                 case(kmin) 
                    points(:,ii) = x(i,j,1,:)
                 case(kmax) 
                    points(:,ii) = x(i,j,kl,:)
                 end select
              end do
           end do
        end if
     end do bocos
  end do domains
end subroutine getForcePoints


! subroutine getCGNSData(ndofcgns,data)

!   use BCTypes
!   use blockPointers
!   use cgnsGrid
!   use communication
!   use flowVarRefState
!   use warpingPetsc
!   use block
!   implicit none
!   !
!   !      Subroutine arguments.
!   !
!   integer(kind=intType),intent(in)  :: ndofcgns
!   real(kind=realType) ,intent(out)  :: data(ndofcgns)

!   !
!   !      Local variables.
!   !
!   integer(kind=intType) :: lowInd,highInd,ierr,size,i
!   integer(kind=intType) ,allocatable,dimension(:) :: indices

!   call VecGetOwnershipRange(cgnsGridVec,lowInd,HighInd,ierr)
  
!   size = highInd-lowInd
!   allocate(indices(size))
!   do i=1,size
!      indices(i) = lowInd + i - 1
!   end do
  
!   call VecGetValues(cgnsGridVec,size,indices,data,ierr)
!   deallocate(indices)
   
!  end subroutine getCGNSData

subroutine getIset(nn,i,j,k,il,jl,kl,iset)
  use precision 
  use communication
  use block
  use warpingPETSc
  implicit none

  integer(kind=intType) :: i,j,k,il,jl,kl,iset(3),nn

  iset(1) = cumdofproc(myID) + cumdofblock(nn) + &
       (k-1)*jl*il*3 + &
       (j-1)*il*3    + &
       (i-1)*3

  iset(2) = iset(1) + 1
  iset(3) = iset(1) + 2

end subroutine getIset


