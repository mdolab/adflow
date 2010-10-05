subroutine getNumberLocalNodes(ndof)
  use block
  use blockPointers
  use cgnsGrid
  implicit none
  integer(kind=intType) :: nn,ndof
  ndof = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     ndof = ndof + il*jl*kl*3
  end do
end subroutine getNumberLocalNodes

subroutine getNumberLocalForceNodes(ndof)
  use block
  use blockPointers
  use bctypes
  use communication
  implicit none
  integer(kind=intType) :: nn,ndof,mm,nni,nnj,nnk
  ndof = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if((BCType(mm) == EulerWall       .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal)) then
           nni = inEnd(mm) - inBeg(mm) + 1
           nnj = jnEnd(mm) - jnBeg(mm) + 1
           nnk = knEnd(mm) - knBeg(mm) + 1

           ndof = ndof + nni*nnj*nnk*3
        end if
     end do bocos
  end do domains
end subroutine getNumberLocalForceNodes

subroutine getCGNSMeshIndices(ndof,indices)
  use block
  use blockPointers
  use cgnsGrid
  implicit none

  ! subroutine arguments
  integer(kind=intType), intent(in) :: ndof
  integer(kind=intType), intent(out):: indices(ndof)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,ii,indx,indy,indz,il_cg,jl_cg,kl_cg
  integer(kind=intType) ,allocatable,dimension(:) :: dof_offset

  allocate(dof_offset(cgnsNDom))
  dof_offset(1) = 0
  do nn=2,cgnsNDom
     dof_offset(nn) = dof_offset(nn-1) + cgnsDoms(nn-1)%il*cgnsDoms(nn-1)%jl*cgnsDoms(nn-1)%kl*3
  end do

  ii = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     do k=1,kl
        do j=1,jl
           do i=1,il
              il_cg = cgnsDoms(nbkGlobal)%il
              jl_cg = cgnsDoms(nbkGlobal)%jl
              kl_cg = cgnsDoms(nbkGlobal)%kl

              ii = ii + 1

              indx = iBegOr + i - 1
              indy = jBegOr + j - 1
              indz = kBegOr + k - 1

              indices(ii*3-2) = dof_offset(nbkGlobal)  + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3

              indices(ii*3-1) = dof_offset(nbkGlobal)  + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3 + 1

              indices(ii*3 ) = dof_offset(nbkGlobal)   + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3 + 2

           end do ! k loop
        end do ! j loop
     end do ! i loop
  end do ! domain loop
end subroutine getCGNSMeshIndices

subroutine getCGNSForceIndices(ndof,indices)
  use block
  use blockPointers
  use cgnsGrid
  use bctypes
  implicit none

  ! subroutine arguments
  integer(kind=intType), intent(in) :: ndof
  integer(kind=intType), intent(out):: indices(ndof)

  ! Local Variables
  integer(kind=intType) :: nn,mm,i,j,k,ii,indx,indy,indz,il_cg,jl_cg,kl_cg
  integer(kind=intType) ,allocatable,dimension(:) :: dof_offset

  allocate(dof_offset(cgnsNDom))
  dof_offset(1) = 0
  do nn=2,cgnsNDom
     dof_offset(nn) = dof_offset(nn-1) + cgnsDoms(nn-1)%il*cgnsDoms(nn-1)%jl*cgnsDoms(nn-1)%kl*3
  end do

  ii = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     il_cg = cgnsDoms(nbkGlobal)%il
     jl_cg = cgnsDoms(nbkGlobal)%jl
     kl_cg = cgnsDoms(nbkGlobal)%kl

     bocos: do mm=1,nBocos
        if((BCType(mm) == EulerWall       .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal)) then

           do k=knBeg(mm),knEnd(mm)
              do j=jnBeg(mm),jnEnd(mm)
                 do i=inBeg(mm),inEnd(mm)
                    ii = ii + 1

                    indx = iBegOr + i - 1
                    indy = jBegOr + j - 1
                    indz = kBegOr + k - 1

                    indices(ii*3-2) = dof_offset(nbkGlobal)  + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3

                    indices(ii*3-1) = dof_offset(nbkGlobal)  + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3 + 1

                    indices(ii*3 ) = dof_offset(nbkGlobal)   + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3 + 2
                 end do ! i loop
              end do ! k loop
           end do ! k loop
        end if ! BC Face
     end do bocos
  end do ! Domain Loop
end subroutine getCGNSForceIndices
