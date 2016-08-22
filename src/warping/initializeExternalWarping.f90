subroutine getCGNSMeshIndices(ndof,indices)

  use constants
  use blockPointers, only : nDom, nBKGlobal, il, jl, kl, iBegOr, jBegOr, kBegOr
  use cgnsGrid, only : cgnsDoms, cgnsnDom
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

