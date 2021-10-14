module warping

  ! This module cotains the required inferface functions for using an
  ! external mesh warping utility with ADflow

contains

  subroutine getCGNSMeshIndices(ndof,indices)

    use constants
    use blockPointers, only : nDom, nBKGlobal, il, jl, kl, iBegOr, jBegOr, kBegOr
    use cgnsGrid, only : cgnsDoms, cgnsnDom
    use utils, only : setPointers
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
       dof_offset(nn) = dof_offset(nn-1) + &
            cgnsDoms(nn-1)%il*cgnsDoms(nn-1)%jl*cgnsDoms(nn-1)%kl*3
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

             end do ! i loop
          end do ! j loop
       end do ! k loop
    end do ! domain loop
    deallocate(dof_offset)
  end subroutine getCGNSMeshIndices

  subroutine setGrid(grid,ndof)

    ! The purpose of this routine is to set the grid dof as returned by
    ! the external warping. This function takes the "Base" grid at the
    ! first time instance and does rotation/translation operations to
    ! get the grid at subsequent time instances
    use constants
    use blockPointers, only : nDom, il, jl, kl, x
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use section, only : sections, nSections
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use monitor, only : timeUnsteadyRestart, timeUnsteady
    use inputPhysics, only : equationMode
    use utils, only : setPointers, rotMatrixRigidBody
    use preprocessingAPI, only : xhalo
    implicit none

    integer(kind=intType),intent(in) :: ndof
    real(kind=realType) ,intent(in) :: grid(ndof)

    ! Local Variables

    integer(kind=intType) :: nn,i,j,k,counter,sps
    real(kind=realType) :: t(nSections),dt(nSections)
    real(kind=realType) :: displ(3)
    real(kind=realType) :: tOld,tNew

    real(kind=realType), dimension(3)   :: rotationPoint,r
    real(kind=realType), dimension(3,3) :: rotationMatrix


    if (equationMode == steady .or. equationMode == TimeSpectral) then
       timeUnsteady = zero

       ! This is very straight forward...loop over all domains and set all elements
       do nn=1,nSections
          dt(nn) = sections(nn)%timePeriod &
               / real(nTimeIntervalsSpectral,realType)
       enddo

       do sps = 1,nTimeIntervalsSpectral
          do nn=1,nSections
             t(nn) = (sps-1)*dt(nn)
          enddo

          ! Compute the displacements due to the rigid motion of the mesh.

          displ(:) = zero

          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
          counter = 0
          do nn=1,nDom
             call setPointers(nn,1_intType,sps)
             do k=1,kl
                do j=1,jl
                   do i=1,il
                      ! r is distance from grid point to rotationPoint
                      r = grid(3*counter+1:3*counter+3) - rotationPoint

                      X(i,j,k,:) = rotationPoint + matmul(rotationMatrix,r) + displ
                      counter = counter + 1

                   end do
                end do
             end do
          end do
          call xhalo(1_intType)
       end do
    else
       counter = 0
       sps = 1
       do nn=1,nDom
          call setPointers(nn,1_intType,sps)
          do k=1,kl
             do j=1,jl
                do i=1,il
                   X(i,j,k,:) = grid(3*counter+1:3*counter+3)
                   counter = counter + 1
                end do
             end do
          end do
       end do
       call xhalo(1_intType)
    end if

  end subroutine setGrid

  subroutine setGridForOneInstance(grid,sps)

    ! The purpose of this routine is to set the grid dof as returned by
    ! the external warping. This routine will take in the deformed mesh
    ! and set it to "sps"th time instance
    use constants
    use blockPointers, only : nDom, il, jl, kl, x
    use section, only : sections, nSections
    use inputPhysics, only : equationMode
    use utils, only : setPointers
    use preprocessingAPI, only : xhalo
    implicit none

    real(kind=realType) ,dimension(:), intent(in) :: grid
    integer, intent(in) :: sps

    ! Local Variables

    integer(kind=intType) :: nn,i,j,k,counter

    counter = 0
    do nn=1,nDom
       call setPointers(nn,1_intType,sps)
       do k=1,kl
          do j=1,jl
             do i=1,il
                X(i,j,k,:) = grid(3*counter+1:3*counter+3)
                counter = counter + 1
             end do
          end do
       end do
    end do
    call xhalo(1_intType)

  end subroutine setGridForOneInstance


  subroutine getGrid(grid,ndof)

    ! Opposite of setGrid. This is ONLY a debugging function. NOT used
    ! in regular usage. Really only useful for direct mesh manipulation
    ! on single block and a single processor. s

    use constants
    use blockPointers, only : nDom, il, jl, kl, x
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    implicit none
    integer(kind=intType),intent(in) :: ndof
    real(kind=realType) ,intent(out) :: grid(ndof)

    ! Local Variables
    integer(kind=intType) :: nn,i,j,k,l,counter,sps

    ! This is very straight forward...loop over all domains and copy out
    counter = 1
    do sps = 1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn,1_intType,sps)
          do k=1,kl
             do j=1,jl
                do i=1,il
                   do l=1,3
                      grid(counter) = X(i,j,k,l)
                      counter = counter + 1
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine getGrid

  subroutine getStatePerturbation(randVec, nRand, randState, nRandState)

    use constants
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use blockPointers, only : nDom, il, jl, kl, nx, ny, nz, x, nbkglobal, iBegOr, jBegOr, kBegOr
    use communication, only : adflow_comm_world, myid
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use adjointVars, only : nCellsLocal
    use flowVarRefState, only : nw
    use utils, only : setPointers, EChk
    implicit none

    ! Input Parameters
    real(kind=realType), intent(in), dimension(nRand) :: randVec
    integer(kind=intType), intent(in) :: nRand, nRandState

    ! Ouput Parameters
    real(kind=realType), intent(out), dimension(nRandState) :: randState

    ! Working parameters
    integer(kind=intType) :: i, j, k, ierr, l, nx_cg, ny_cg, nz_cg
    integer(kind=intType) :: sps, ii, indx, indy, indz, nn, cgnsInd
    integer(kind=intType) :: dofCGNSPerInstance
    integer(kind=intType) ,allocatable, dimension(:) :: dof_offset

    allocate(dof_offset(cgnsNDom))
    dof_offset(1) = 0
    do nn=2,cgnsNDom
       dof_offset(nn) = dof_offset(nn-1) + &
            cgnsDoms(nn-1)%nx*cgnsDoms(nn-1)%ny*cgnsDoms(nn-1)%nz*nw
    end do

    dofCGNSPerInstance = nRand/nTimeIntervalsSpectral

    ii = 0
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, 1, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il

                   nx_cg = cgnsDoms(nbkGlobal)%nx
                   ny_cg = cgnsDoms(nbkGlobal)%ny
                   nz_cg = cgnsDoms(nbkGlobal)%nz

                   indx = iBegOr + i - 2
                   indy = jBegOr + j - 2
                   indz = kBegOr + k - 2

                   do l=1, nw
                      cgnsInd =  (sps-1)*dofCGNSPerInstance + &
                           dof_offset(nbkGlobal)  + &
                           (indz-1)*ny_cg*nx_cg*nw + &
                           (indy-1)*nx_cg*nw       + &
                           (indx-1)*nw + l
                      randState(nw*ii + l) = randVec(cgnsInd)
                   end do

                   ii = ii + 1

                end do
             end do
          end do
       end do
    end do
  end subroutine getStatePerturbation

  subroutine getSurfacePerturbation(xRand, nRand, randSurface, nRandSurface, famList, nFamList, sps)

    use constants
    use blockPointers, only : nDom, BCData, nBocos, BCFaceID, il, jl ,kl
    use communication, only : adflow_comm_world, myid
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers, EChk
    use sorting, only : famInList
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
    use surfaceFamilies, only : BCFamGroups, familyExchange, BCFamExchange
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Parameters
    real(kind=realType), intent(in), dimension(nRand) :: xRand
    integer(kind=intType), intent(in) :: nRand, nRandSurface
    integer(kind=intType), intent(in) :: famList(nFamList), nFamList, sps
    ! Ouput Parameters
    real(kind=realType), intent(inout), dimension(3*nRandSurface) :: randSurface

    ! Working parameters
    integer(kind=intType) :: i, j, k, ierr, iDim, iBeg, iEnd, jBeg, jEnd, nn, mm
    integer(kind=intType) :: ii, jj, indI, indJ, indK, jjInd, iBCGroup
    type(zipperMesh), pointer :: zipper
    type(familyexchange), pointer :: exch
    logical :: BCGroupNeeded
    real(kind=realType), dimension(:), pointer :: localPtr

    ii = 0
    jj = 0
    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, sps)

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,nBocos

          ! NODE Based
          jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
          iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

          famInclude: if (famInList(BCdata(mm)%famID, famList)) then

             do j=jBeg, jEnd ! This is a node loop
                do i=iBeg, iEnd ! This is a node loop
                   select case(BCFaceID(mm))
                   case(imin)
                      indI = 1
                      indJ = i
                      indK = j
                   case(imax)
                      indI = il
                      indJ = i
                      indK = j
                   case(jmin)
                      indI = i
                      indJ = 1
                      indK = j
                   case(jmax)
                      indI = i
                      indJ = jl
                      indK = j
                   case(kmin)
                      indI = i
                      indJ = j
                      indK = 1
                   case(kmax)
                      indI = i
                      indJ = j
                      indK = kl
                   end select

                   do iDim=1,3
                      jjInd = jj + (indK-1)*il*jl + (indJ-1)*il + (indI-1) + iDim
                      randSurface(3*ii+iDim) = xRand(jjInd)
                   end do
                   ii = ii +1
                end do
             end do
          end if famInclude
       end do bocos

       ! jj is the counter through xRand. Increment it by the full
       ! block.
       jj = jj + il*jl*kl*3
    end do domains

    ! No overset or not zipper, return
    if (.not. oversetPresent) then ! .or. .not. includeZipper) then
       return
    end if


    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    do iBCGroup=1, nFamExchange

       zipper => zipperMeshes(iBCGroup)

       if (.not. zipper%allocated) then
          cycle
       end if

       exch => BCFamExchange(iBCGroup, sps)
       BCGroupNeeded = .False.
       BCGroupFamLoop: do i=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList)) then
             BCGroupNeeded = .True.
             exit BCGroupFamLoop
          end if
       end do BCGroupFamLoop

       if (.not. BCGroupNeeded) then
          cycle
       end if

       ! Now we know we *actually* need something from this BCGroup.

       ! Loop over each dimension individually since we have a scalar
       ! scatter.
       dimLoop: do iDim=1,3

          call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! The local Pointer is just the localRandSurface we've set
          ! above.
          do j=1, size(localPtr)
             localPtr(i) = randSurface(3*(j-1) + iDim)
          end do

          ! Restore the pointer
          call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Now scatter this to the zipper
          call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! The values we need are precisely what is in zipper%localVal
          call vecGetArrayF90(zipper%localVal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Just copy the received seeds into the random aray
          do j=1, size(localPtr)
             ! Careful here becuase we have to interlate the dim
             randSurface(3*ii + 3*(j-1) + iDim) = localPtr(j)
          end do

       end do dimLoop

       ! Increcment the running ii counter.
       ii = ii + size(localPtr)
    end do
  end subroutine getSurfacePerturbation

end module warping
