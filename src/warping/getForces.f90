subroutine getForces()

  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use flowVarRefState
  use mdDataLocal
  use warpingPetsc
  use block
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: ii, jj, mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  real(kind=realType) :: scaleDim, fact, pp, f(3)
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz

  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType), dimension(:,:,:), pointer :: ss
  logical :: storeSubface
  integer(kind=intType) :: iset(3),size,highInd,lowInd,ndof

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

  ii = 0
  call vecZeroEntries(sumbGridVec,ierr)
  call vecZeroEntries(cgnsGridVec,ierr)

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           select case (BCFaceID(mm))

           case (iMin)
              pp2 => p( 2,1:,1:); pp1 => p( 1,1:,1:); ss => si0( 1,:,:,:)
              fact = -one

           case (iMax)
              pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:); ss => si0(il,:,:,:)
              fact = one

           case (jMin)
              pp2 => p(1:, 2,1:); pp1 => p(1:, 1,1:); ss => sj0(:, 1,:,:)
              fact = -one

           case (jMax)
              pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:); ss => sj0(:,jl,:,:)
              fact = one

           case (kMin)
              pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1); ss => sk0(:,:, 1,:)
              fact = -one

           case (kMax)
              pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke); ss => sk0(:,:,kl,:)
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

                 ! Compute the corresponding force.

                 f(1) = pp*ss(i,j,1)
                 f(2) = pp*ss(i,j,2)
                 f(3) = pp*ss(i,j,3)

                 ! Now its just an indexing problem...we need to
                 ! set each of the 4 entries in the sumbGridVec
                 ! with AddValues

                 select case(BCFaceID(mm))
                 case(imin)
                    call getIset(nn,1  ,i-1,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,1  ,i  ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,1  ,i-1,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,1  ,i  ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 case(imax)
                    call getIset(nn,il ,i-1,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,il ,i  ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,il ,i-1,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,il ,i  ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 case(jmin)
                    call getIset(nn,i-1,1  ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,1  ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i-1,1  ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,1  ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 case(jmax)
                    call getIset(nn,i-1,jl ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,jl ,j-1,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i-1,jl ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,jl ,j  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 case(kmin)
                    call getIset(nn,i-1,j-1,1  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,j-1,1  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i-1,j  ,1  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,j  ,1  ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 case(kmax)
                    call getIset(nn,i-1,j-1,kl ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i  ,j-1,kl ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                    call getIset(nn,i-1,j  ,kl ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)          
                    call getIset(nn,i  ,j  ,kl ,il,jl,kl,iset)
                    call VecSetValues(sumbGridVec,3,iset,f,ADD_VALUES,ierr)
                 end select
              end do
           end do
        end if
     end do bocos
  end do domains

  call VecAssemblyBegin(sumbGridVec,ierr)
  call VecAssemblyEnd(sumbGridVec,ierr)

  call VecScatterBegin(sumbTOcgnsForce,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd  (sumbTOcgnsForce,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)

end subroutine getForces

subroutine getCGNSData(ndofcgns,data)

  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use flowVarRefState
  use mdDataLocal
  use warpingPetsc
  use block
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType),intent(in)  :: ndofcgns
  real(kind=realType) ,intent(out)  :: data(ndofcgns)

  !
  !      Local variables.
  !
  integer(kind=intType) :: lowInd,highInd,ierr,size,i
  integer(kind=intType) ,allocatable,dimension(:) :: indices

  call VecGetOwnershipRange(cgnsGridVec,lowInd,HighInd,ierr)
  
  size = highInd-lowInd
   allocate(indices(size))
   do i=1,size
      indices(i) = lowInd + i - 1
   end do
  
   call VecGetValues(cgnsGridVec,size,indices,data,ierr)
   deallocate(indices)
   
 end subroutine getCGNSData

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

  iset(2) = cumdofproc(myID) + cumdofblock(nn) + &
       (k-1)*jl*il*3 + &
       (j-1)*il*3    + &
       (i-1)*3 + 1

  iset(3) = cumdofproc(myID) + cumdofblock(nn) + &
       (k-1)*jl*il*3 + &
       (j-1)*il*3    + &
       (i-1)*3 + 2

end subroutine getIset

