
subroutine getForces(forces, npts, sps)
  use constants
  use communication, only : myid
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use inputPhysics, only : forcesAsTractions
  use utils, only : setPointers, terminate, EChk
  use surfaceIntegrations, only : integrateSurfaces
  use surfaceFamilies, only : fullfamList
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(inout) :: forces(3,npts)

  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: sss(3),v2(3),v1(3), qa, sepSensor, Cavitation
  real(kind=realType) :: sepSensorAvg(3)
  real(kind=realType) :: Fp(3), Fv(3), Mp(3), Mv(3), yplusmax, qf(3)
  real(kind=realType) :: localValues(nLocalValues)
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch
  real(kind=realType), dimension(:), pointer :: localPtr
  real(kind=realType), dimension(nCostFunction) :: funcValues
  ! Make sure *all* forces are computed. Sectioning will be done
  ! else-where.

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     localValues = zero
     call integrateSurfaces(localValues, fullFamList)
  end do domains

  if (forcesAsTractions) then
     ! Compute tractions if necessary
     call computeNodalTractions(sps)
  else
     call computeNodalForces(sps)
  end if

  ii = 0
  domains2: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    Forces(:, ii) = bcData(mm)%Tp(i, j, :) + bcData(mm)%Tv(i, j, :)
                 else
                    Forces(:, ii) = bcData(mm)%F(i, j, :)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  ! We know must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  if (.not. oversetPresent .or. .not. zipper%allocated) then
     return
  end if

  if (.not. forcesAsTractions) then
     ! We have a zipper and regular forces are requested. This is not yet supported.
     call terminate('getForces', 'getForces() is not implmented for zipper meshes and '&
          &'forcesAsTractions=False')
  end if

  ! Loop over each dimension individually since we have a scalar
  ! scatter.
  dimLoop: do iDim=1,3

     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Copy in the values we already have to the exchange.
     ii = size(LocalPtr)
     localPtr = forces(iDim, 1:ii)

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

     ! Just copy the received data into the forces array. Just on root proc:
     if (myid == 0) then
        forces(iDim, ii+1:ii+size(localPtr)) = localPtr
     end if

     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do dimLoop

end subroutine getForces

subroutine getForces_d(forces, forcesd, npts, sps)

  ! This routine performs the forward mode linearization getForces. It
  ! takes in perturbations defined on bcData(mm)%Fp, bcData(mm)%Fv and
  ! bcData(mm)%area and computes either the nodal forces or nodal
  ! tractions.
  use constants
  use communication, only : myid
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use utils, only : setPointers, setPointers_d, EChk, terminate
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out), dimension(3, npts) :: forces, forcesd
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qa, qad, qf, qfd
  real(kind=realType), dimension(:), pointer :: localPtr, localPtrd
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch

  if (forcesAsTractions) then
     call computeNodalTractions_d(sps)
  else
     call computeNodalForces_d(sps)
  end if

  ! Extract the values out into the output derivative array
  ii = 0
  domains2: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    Forcesd(:, ii) = bcDatad(mm)%Tp(i, j, :) + bcDatad(mm)%Tv(i, j, :)
                 else
                    Forcesd(:, ii) = bcDatad(mm)%F(i, j, :)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  ! We know must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  if (.not. oversetPresent .or. .not. zipper%allocated) then
     return
  end if

  if (.not. forcesAsTractions) then
     ! We have a zipper and regular forces are requested. This is not yet supported.
     call terminate('getForces', 'getForces() is not implmented for zipper meshes and '&
          &'forcesAsTractions=False')
  end if

  ! Loop over each dimension individually since we have a scalar
  ! scatter.
  dimLoop: do iDim=1,3

     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Copy in the values we already have to the exchange.
     ii = size(LocalPtr)
     localPtr = forcesd(iDim, 1:ii)

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

     ! Just copy the received data into the forces array. Just on root proc:
     if (myid == 0) then
        forcesd(iDim, ii+1:ii+size(localPtr)) = localPtr
     end if

     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do dimLoop

end subroutine getForces_d

subroutine getForces_b(forcesd, npts, sps)

  ! This routine performs the reverse of getForces. It takes in
  ! forces_b and perfroms the reverse of the nodal averaging procedure
  ! in getForces to compute bcDatad(mm)%Fp, bcDatad(mm)%Fv and
  ! bcDatad(mm)%area.
  use constants
  use communication, only : myid
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use utils, only : EChk, setPointers, setPointers_d
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(inout) :: forcesd(3, npts)
  integer(kind=intType) :: mm, nn, i, j, ii, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch
  real(kind=realType), dimension(:), pointer :: localPtr
  real(kind=realType), dimension(3, npts) :: forces

  ! Run nonlinear code to make sure that all intermediate values are
  ! updated.
  call getForces(forces, npts, sps)

  ! We know must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  zipperReverse: if (oversetPresent .and. zipper%allocated) then

     ! Loop over each dimension individually since we have a scalar
     ! scatter.
     dimLoop: do iDim=1,3

        call vecGetArrayF90(zipper%localVal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ii = exch%nNodes
        ! Just copy the received data into the forces array. Just on root proc:
        if (myid == 0) then
           do i=1, size(localPtr)
              localPtr(i) = forcesd(iDim, ii+i)
              forcesd(iDim, ii+i) = zero
           end do
        end if

        call vecGetArrayF90(zipper%localVal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Zero the vector we are scatting into:
        call VecSet(exch%nodeValLocal, zero, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Scatter values from the root using the zipper scatter.
        call VecScatterBegin(zipper%scatter, zipper%localVal, &
             exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call VecScatterEnd(zipper%scatter, zipper%localVal, &
             exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Accumulate the scatted values onto forcesd
        ii = size(localPtr)
        forcesd(iDim, 1:ii) = forcesd(iDim, 1:ii) + localPtr

        ! Restore the pointer
        call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

     end do dimLoop
  end if zipperReverse

  ! Set the incoming derivative values
  ii = 0
  domains2: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    bcDatad(mm)%Tp(i, j, :) = forcesd(:, ii)
                    bcDatad(mm)%Tv(i, j, :) = forcesd(:, ii)
                 else
                    bcDatad(mm)%F(i, j, :) = forcesd(:, ii)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  if (.not. forcesAsTractions) then
     ! For forces, we can accumulate the nodal seeds on the Fp and Fv
     ! values. The area seed is zeroed.
     call computeNodalForces_b(sps)
  else
     call computeNodalTractions_b(sps)
  end if

end subroutine getForces_b


subroutine surfaceCellCenterToNode(exch)

  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : familyExchange
  use utils, only : setPointers, EChk
  use sorting, only : famInList
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

  type(familyExchange) :: exch
  integer(kind=intType) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qv
  real(kind=realType), dimension(:), pointer :: localPtr

  ! We assume that normalization factor is already computed
  sps = exch%sps
  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  localPtr = zero

  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0,nj-2
              do i=0,ni-2
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 qv = fourth * BCData(mm)%cellVal(i+1, j+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtr(ind(jj)) = localPtr(ind(jj)) + qv
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if famInclude
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the current face based value
  call vecSet(exch%nodeValGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now divide by the weighting. We can do this with a vecpointwisemult
  call vecPointwiseMult(exch%nodeValGlobal, exch%nodeValGlobal, &
       exch%sumGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Push back to the local values
  call VecScatterBegin(exch%scatter, exch%nodeValGlobal, &
       exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValGlobal, &
       exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=1,nj
              do i=1,ni
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 ii = ii + 1
                 BCData(mm)%nodeVal(i, j) = localPtr(ii)
              end do
           end do
        end if famInclude2
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine surfaceCellCenterToNode

subroutine computeWeighting(exch)

  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : familyExchange
  use utils, only : setPointers, EChk
  use sorting, only : famInList
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  type(familyExchange) :: exch
  integer(kind=intType) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf, qa
  real(kind=realType), dimension(:), pointer :: localPtr

  sps = exch%sps

  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = zero
  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0,nj-2
              do i=0,ni-2

                 ! Scatter a quarter of the face value to each node:
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 qa = fourth*BCData(mm)%cellVal(i+1, j+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtr(ind(jj)) = localPtr(ind(jj)) + qa
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if famInclude
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the face value
  call vecSet(exch%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing. Note that we check dividing by zero and just
  ! set those to zero.

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  do i=1, size(localPtr)
     if (localPtr(i) == zero) then
        localPtr(i) = zero
     else
        localPtr(i) = one/localPtr(i)
     end if
  end do

  call vecRestoreArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine computeWeighting

subroutine computeNodalTractions(sps)
  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : BCFamExchange, familyExchange
  use utils, only : setPointers, EChk, isWallType
  implicit none

  integer(kind=intType), intent(in) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf, qa
  real(kind=realType), dimension(:), pointer :: localPtr
  type(familyExchange), pointer :: exch

  ! Set the pointer to the wall exchange:
  exch => BCfamExchange(iBCGroupWalls, sps)

  ! Set the weighting factors. In this case, area
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

        bocoType1: if(isWallType(BCType(mm))) then
           BCData(mm)%cellVal => BCData(mm)%area(:, :)
        end if bocoType1
     end do
  end do
  call computeWeighting(exch)

  FpFvLoop: do iDim=1, 6
     ! ii is the running counter through the pointer array.
     ii = 0
     do nn=1, nDom
        call setPointers(nn, 1_intType, sps)
        do mm=1, nBocos
           bocoType2: if(isWallType(BCType(mm))) then
              if (iDim <= 3) then
                 BCData(mm)%cellVal => BCData(mm)%Fp(:, :, iDim)
                 BCData(mm)%nodeVal => BCData(mm)%Tp(:, :, iDim)
              else
                 BCData(mm)%cellVal => BCData(mm)%Fv(:, :, iDim-3)
                 BCData(mm)%nodeVal => BCData(mm)%Tv(:, :, iDim-3)
              end if
           end if bocoType2
        end do
     end do

     call surfaceCellCenterToNode(exch)

  end do FpFVLoop

end subroutine computeNodalTractions

subroutine computeNodalTractions_d(sps)

  ! Forward mode lineariation of nodal tractions

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use utils, only : setPointers, setPointers_d, EChk
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  integer(kind=intType), intent(in) :: sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qa, qad, qf, qfd
  real(kind=realType), dimension(:), pointer :: localPtr, localPtrd
  type(familyExchange), pointer :: exch
  Vec nodeValLocald, nodeValGlobald, sumGlobald, tmp

  exch => BCFamExchange(iBCGroupWalls, sps)

  call VecDuplicate(exch%nodeValLocal, nodeValLocald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%nodeValGlobal, nodeValGlobald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%sumGlobal, sumGlobald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%sumGlobal, tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(nodeValLocald, localPtrd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtrd = zero
  localPtr = zero
  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers_d(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
        ni = iEnd - iBeg + 1
        nj = jEnd - jBeg + 1

        if(BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           do j=0,nj-2
              do i=0,ni-2

                 ! Scatter a quarter of the area to each node:
                 qa = fourth*BCData(mm)%area(i+iBeg+1, j+jBeg+1)
                 qad = fourth*BCDatad(mm)%area(i+iBeg+1, j+jBeg+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtrd(ind(jj)) = localPtrd(ind(jj)) + qad
                    localPtr(ind(jj)) = localPtr(ind(jj)) + qa
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if
     end do
  end do
  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecRestoreArrayF90(nodeValLocald, localPtrd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the area
  call vecSet(exch%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the area derivative
  call vecSet(sumGlobald, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, nodeValLocald, &
       sumGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, nodeValLocald, &
       sumGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing. Here we need the original value too:

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(sumGlobald, localPtrd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtrd = -(localPtrd/localPtr**2)
  localPtr = one/localPtr

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecRestoreArrayF90(sumGlobald, localPtrd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now do each of the three dimensions for the pressure and viscous forces
  dimLoop: do iDim=1, 6

     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecGetArrayF90(nodeValLocald, localPtrd, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     localPtr  = zero
     localPtrd = zero

     ! ii is the running counter through the pointer array.
     ii = 0
     do nn=1, nDom
        call setPointers_d(nn, 1_intType, sps)
        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           if(BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              do j=0,nj-2
                 do i=0,ni-2
                    if (iDim <= 3) then
                       qf  = fourth*BCData (mm)%Fp(i+iBeg+1, j+jBeg+1, iDim)
                       qfd = fourth*BCDatad(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim)
                    else
                       qf  = fourth*BCData (mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3)
                       qfd = fourth*BCDatad(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3)
                    end if

                    ind(1) = ii + (j  )*ni + i + 1
                    ind(2) = ii + (j  )*ni + i + 2
                    ind(3) = ii + (j+1)*ni + i + 2
                    ind(4) = ii + (j+1)*ni + i + 1
                    do jj=1,4
                       localPtr (ind(jj)) = localPtr (ind(jj)) + qf
                       localPtrd(ind(jj)) = localPtrd(ind(jj)) + qfd
                    end do
                 end do
              end do
              ii = ii + ni*nj
           end if
        end do
     end do

     call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecRestoreArrayF90(nodeValLocald, localPtrd, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Globalize the current force
     call vecSet(exch%nodeValGlobal, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Globalize the current force derivative
     call vecSet(nodeValGlobald, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterBegin(exch%scatter, nodeValLocald, &
          nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, nodeValLocald, &
          nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! The product rule here: (since we are multiplying)
     ! nodeValGlobal = nodeValGlobal * invArea
     ! nodeValGlobald = nodeValGlobald*invArea + nodeValGlobal*invAread

     ! First term:  nodeValGlobald = nodeValGlobald*invArea
     call vecPointwiseMult(nodeValGlobald, nodeValGlobald, &
          exch%sumGlobal, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Second term:, tmp = nodeValGlobal*invAread
     call vecPointwiseMult(tmp, exch%nodeValGlobal, sumGlobald, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Sum the second term into the first
     call VecAXPY(nodeValGlobald, one, tmp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Push back to the local values
     call VecScatterBegin(exch%scatter, nodeValGlobald, &
          nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, nodeValGlobald, &
          nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecGetArrayF90(nodeValLocald, localPtrd, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ii = 0
     do nn=1, nDom
        call setPointers_d(nn, 1_intType, sps)
        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

           if(BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    ii = ii + 1
                    if (iDim <= 3) then
                       bcDatad(mm)%Tp(i, j, iDim) = localPtrd(ii)
                    else
                       bcDatad(mm)%Tv(i, j, iDim-3) = localPtrd(ii)
                    end if
                 end do
              end do
           end if
        end do
     end do

     call vecRestoreArrayF90(nodeValLocald, localPtrd, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do dimLoop

  call VecDestroy(nodeValLocald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(nodeValGlobald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(sumGlobald, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine computeNodalTractions_d

subroutine computeNodalTractions_b(sps)

  ! This routine performs the reverse of computeNodalTractions. Tt
  ! takes in bcDatad%Tv and bcDatad%Tp and perfroms the reverse of the
  ! nodal averaging procedure in getForces to compute bcDatad(mm)%Fp,
  ! bcDatad(mm)%Fv and bcDatad(mm)%area.

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use communication
  use utils, only : EChk, setPointers, setPointers_d
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

  integer(kind=intType), intent(in) :: sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf_b, qf, qa, qa_b
  real(kind=realType), dimension(:), pointer :: localPtr, localPtr_b
  type(familyExchange), pointer :: exch
  Vec nodeValLocal_b, nodeValGlobal_b, sumGlobal_b, tmp, tmp_b, T_b

  ! For better readibility
  exch => BCFamExchange(iBCGroupWalls, sps)

  call VecDuplicate(exch%nodeValLocal, nodeValLocal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%nodeValGlobal, nodeValGlobal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%sumGlobal, sumGlobal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(exch%sumGlobal, tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(tmp, T_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! For tractions it's (a lot) more difficult becuase we have to do
  ! the scatter/gather operation.

  ! ==================================
  !  Recompute the dual area
  ! ==================================

  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = zero
  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
        ni = iEnd - iBeg + 1
        nj = jEnd - jBeg + 1

        if(BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           do j=0,nj-2
              do i=0,ni-2

                 ! Scatter a quarter of the area to each node:
                 qa = fourth*BCData(mm)%area(i+iBeg+1, j+jBeg+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtr(ind(jj)) = localPtr(ind(jj)) + qa
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the area
  call vecSet(exch%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing.

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = one/localPtr

  call vecRestoreArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! ==================================
  ! Now trace through the computeNodalTractions() routine
  ! backwards. All the scatters flip direction and INSERT_VALUES
  ! becomes ADD_VALUES and vice-versa
  ! ==================================
  dimLoop: do iDim=1, 6

     ! ====================
     ! Do the forward pass:
     ! ====================
     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     localPtr = zero
     ! ii is the running counter through the pointer array.
     ii = 0
     do nn=1, nDom
        call setPointers(nn, 1_intType, sps)
        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           if(BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              do j=0,nj-2
                 do i=0,ni-2
                    if (iDim <= 3) then
                       qf = fourth*BCData(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim)
                    else
                       qf = fourth*BCData(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3)
                    end if

                    ind(1) = ii + (j  )*ni + i + 1
                    ind(2) = ii + (j  )*ni + i + 2
                    ind(3) = ii + (j+1)*ni + i + 2
                    ind(4) = ii + (j+1)*ni + i + 1
                    do jj=1,4
                       localPtr(ind(jj)) = localPtr(ind(jj)) + qf
                    end do
                 end do
              end do
              ii = ii + ni*nj
           end if
        end do
     end do

     call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)


     ! Globalize the current force
     call vecSet(exch%nodeValGlobal, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! ====================
     ! Do the reverse pass:
     ! ====================

     ! Copy the reverse seed into the local values
     call vecGetArrayF90(nodeValLocal_b, localPtr_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ii = 0
     domains: do nn=1,nDom
        call setPointers_d(nn, 1_intType, sps)

        ! Loop over the number of boundary subfaces of this block.
        bocos: do mm=1, nBocos
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then

              ! This is easy, just copy out F or T in continuous ordering.
              do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
                 do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                    ii = ii + 1
                    if (iDim <= 3) then
                       localPtr_b(ii) = BCDatad(mm)%Tp(i, j, iDim)
                    else
                       localPtr_b(ii) = BCDatad(mm)%Tv(i, j, iDim-3)
                    end if
                 end do
              end do
           end if
        end do bocos
     end do domains

     call vecRestoreArrayF90(nodeValLocal_b, localPtr_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecSet(T_b, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     ! Push up to the global values
     call VecScatterBegin(exch%scatter, nodeValLocal_b, &
          T_b, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, nodeValLocal_b, &
          T_b, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! this is particularly nasty. This is why you don't do
     ! derivatives by hand, kids.
     ! exch%nodeValGlobal = F
     ! nodeValGlobal_b = F_b
     ! T_b = reverse seed for tractions
     ! sumGlobal_b =  inverseDualarea_b
     ! exch%sumGlobal = invDualArea

     ! Basically what we have to compute here is:
     ! Fb = invDualArea * T_b
     ! invDualAreab = invDualAreab + F*T_b

     call vecPointwiseMult(nodeValGlobal_b, exch%sumGlobal, T_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecPointwiseMult(tmp, exch%nodeValGlobal, T_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Accumulate seed on adflowGlobal_b
     call vecAXPY(sumGlobal_b, one, tmp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now communicate F_b back to the local patches

     call VecScatterBegin(exch%scatter, nodeValGlobal_b, &
          nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, nodeValGlobal_b, &
          nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! ============================
     ! Copy the values into patches
     ! ============================

     call vecGetArrayF90(nodeValLocal_b, localPtr_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! ii is the running counter through the pointer array.
     ii = 0
     do nn=1, nDom
        call setPointers_d(nn, 1_intType, sps)
        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           if(BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then

              do j=0,nj-2
                 do i=0,ni-2

                    ind(1) = ii + (j  )*ni + i + 1
                    ind(2) = ii + (j  )*ni + i + 2
                    ind(3) = ii + (j+1)*ni + i + 2
                    ind(4) = ii + (j+1)*ni + i + 1
                    qf_b = zero
                    do jj=1,4
                       qf_b = qf_b + localPtr_b(ind(jj))
                    end do
                    qf_b = qf_b*fourth

                    if (iDim <= 3) then
                       BCDatad(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim) = &
                            BCDatad(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim) + qf_b
                    else
                       BCDatad(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3) = &
                            BCDatad(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3) + qf_b
                    end if

                 end do
              end do
              ii = ii + ni*nj
           end if
        end do
     end do

     call vecRestoreArrayF90(nodeValLocal_b, localPtr_b, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do dimLoop

  ! ============================
  ! Finish the dual area sensitivity.
  ! ============================

  ! On the forward pass we computed:
  ! sumGlobal = one/sumGlobal
  ! So on the reverse pass we need:
  ! sumGlobalb = -(sumGlobalb/sumGlobal**2)

  ! We will do this by getting pointers

  call vecGetArrayF90(sumGlobal_b, localPtr_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Keep in mind localPtr points to sumGlobal which already has
  ! been inversed so we just multiply.
  localPtr_b = -localPtr_b*localPtr**2

  call vecRestoreArrayF90(sumGlobal_b, localPtr_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecRestoreArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Push back to the local patches
  call VecScatterBegin(exch%scatter, sumGlobal_b, &
       nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, sumGlobal_b, &
       nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(nodeValLocal_b, localPtr_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers_d(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
        ni = iEnd - iBeg + 1
        nj = jEnd - jBeg + 1

        if(BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           do j=0,nj-2
              do i=0,ni-2

                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 qa_b = zero
                 do jj=1,4
                    qa_b = qa_b + localPtr_b(ind(jj))
                 end do
                 qa_b = fourth*qa_b
                 BCDatad(mm)%area(i+iBeg+1, j+jBeg+1) = &
                      BCDatad(mm)%area(i+iBeg+1, j+jBeg+1) + qa_b
              end do
           end do
           ii = ii + ni*nj
        end if
     end do
  end do

  call vecRestoreArrayF90(nodeValLocal_b, localPtr_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Remove temporary petsc vecs
  call VecDestroy(nodeValLocal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(nodeValGlobal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(sumGlobal_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(T_b, ierr)
  call EChk(ierr,__FILE__,__LINE__)





end subroutine computeNodalTractions_b

subroutine computeNodalForces(sps)

  ! This subroutine averages the cell based forces and tractions to
  ! node based values. There is no need for communication since we are
  ! simplying summing a quarter of each value to each corner.

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use utils, only : setPointers
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qf(3)

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           BCData(mm)%F = zero
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 qf = fourth*(BCData(mm)%Fp(i,j,:) + BCData(mm)%Fv(i,j,:))
                 BCData(mm)%F(i  , j,   :) = BCData(mm)%F(i  , j,   :) + qf
                 BCData(mm)%F(i-1, j,   :) = BCData(mm)%F(i-1, j  , :) + qf
                 BCData(mm)%F(i  , j-1, :) = BCData(mm)%F(i  , j-1, :) + qf
                 BCData(mm)%F(i-1, j-1, :) = BCData(mm)%F(i-1, j-1, :) + qf
              end do
           end do
        end if
     end do
  end do
end subroutine computeNodalForces

subroutine computeNodalForces_d(sps)

  ! Forward mode linearization of nodalForces

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
  use utils, only : setPointers
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qfd(3)

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           BCDatad(mm)%F = zero
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 qfd = fourth*(BCDatad(mm)%Fp(i,j,:) + BCDatad(mm)%Fv(i,j,:))
                 BCDatad(mm)%F(i  , j,   :) = BCDatad(mm)%F(i  , j,   :) + qfd
                 BCDatad(mm)%F(i-1, j,   :) = BCDatad(mm)%F(i-1, j  , :) + qfd
                 BCDatad(mm)%F(i  , j-1, :) = BCDatad(mm)%F(i  , j-1, :) + qfd
                 BCDatad(mm)%F(i-1, j-1, :) = BCDatad(mm)%F(i-1, j-1, :) + qfd
              end do
           end do
        end if
     end do
  end do
end subroutine computeNodalForces_d

subroutine computeNodalForces_b(sps)

  ! Reverse mode linearization of nodalForces

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
  use utils, only : setPointers_d
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qf_b(3)

  domains: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 qf_b = fourth*(BCDatad(mm)%F(i, j, :) + BCdatad(mm)%F(i-1, j, :) + &
                      BCDatad(mm)%F(i, j-1, :) + BCDatad(mm)%F(i-1, j-1, :))

                 ! Fp and Fv are face-based values
                 BCDatad(mm)%Fp(i, j, :) = BCDatad(mm)%Fp(i, j, :) + qf_b
                 BCDatad(mm)%Fv(i, j, :) = BCDatad(mm)%Fv(i, j, :) + qf_b
              end do
           end do
           ! this needs to be after the update to be the reverse of the forward mode.
           BCDatad(mm)%F = zero
        end if
     end do
  end do domains
end subroutine computeNodalForces_b

subroutine getHeatFlux(hflux, npts, sps)
  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use surfaceFamilies, only : BCFamExchange, familyExchange, &
       zeroCellVal, zeroNodeVal
  use utils, only : setPointers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out) :: hflux(npts)

  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  type(familyExchange), pointer :: exch

  exch => BCFamExchange(iBCGroupWalls, sps)
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     call heatFluxes()

     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

        bocoType1: if (BCType(mm) == NSWallIsoThermal) then
           BCData(mm)%cellVal => BCData(mm)%area(:, :)
        else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
           BCData(mm)%cellVal => zeroCellVal
           BCData(mm)%nodeVal => zeroNodeVal
        end if bocoType1
     end do
  end do

  call computeWeighting(exch)

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        bocoType2: if (BCType(mm) == NSWallIsoThermal) then
           BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
           BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)
        end if bocoType2
     end do
  end do

  call surfaceCellCenterToNode(exch)

  ! Now extract into the flat array:
  ii = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,sps)

     ! Loop over the number of viscous boundary subfaces of this block.
     ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
     ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
     do mm=1,nBocos
        bocoType3: if (BCType(mm) == NSWallIsoThermal) then
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 hflux(ii) = BCData(mm)%nodeHeatFlux(i, j)
              end do
           end do
           ! Simply put in zeros for the other wall BCs
        else if (BCType(mm) == NSWallAdiabatic .or. BCType(mm) == EulerWall) then
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 hflux(ii) = zero
              end do
           end do
        end if bocoType3
     end do
  end do
end subroutine getHeatFlux

subroutine heatFluxes
  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType, BCFaceID, viscSubFace
  use BCPointers, only : ssi
  use flowVarRefState, only : pRef, rhoRef
  use utils, only : setPointers, setBCPointers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, ii, mm
  real(kind=realType) :: fact, scaleDim
  real(kind=realType) :: qw, qA
  logical :: heatedSubface

  ! Set the actual scaling factor such that ACTUAL heat flux is computed
  ! The factor is determined from stanton number
  scaleDim = pRef*sqrt(pRef/rhoRef)

  ! Loop over the boundary subfaces of this block.
  bocos: do mm=1, nBocos

     ! Only do this on isoThermalWalls
     if (BCType(mm) == NSWallIsoThermal) then

        ! Set a bunch of pointers depending on the face id to make
        ! a generic treatment possible. The routine setBcPointers
        ! is not used, because quite a few other ones are needed.
        call setBCPointers(mm, .True.)

        select case (BCFaceID(mm))
        case (iMin, jMin, kMin)
           fact = -one
        case (iMax, jMax, kMax)
           fact = one
        end select

        ! Loop over the quadrilateral faces of the subface. Note that
        ! the nodal range of BCData must be used and not the cell
        ! range, because the latter may include the halo's in i and
        ! j-direction. The offset +1 is there, because inBeg and jnBeg
        ! refer to nodal ranges and not to cell ranges.
        !
        do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
           do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

              ! Compute the normal heat flux on the face. Inward positive.
              BCData(mm)%cellHeatFlux(i,j) = -fact*scaleDim* &
                   sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2) * &
                   ( viscSubface(mm)%q(i,j,1)*BCData(mm)%norm(i,j,1) &
                   + viscSubface(mm)%q(i,j,2)*BCData(mm)%norm(i,j,2) &
                   + viscSubface(mm)%q(i,j,3)*BCData(mm)%norm(i,j,3))
           enddo
        end do
     end if
  enddo bocos
end subroutine heatFluxes

subroutine setTNSWall(tnsw, npts, sps)

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType
  use flowVarRefState, only : TRef
  use utils, only : setPointers
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(in) :: tnsw(npts)

  ! Local Variables
  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  ii = 0
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     ! Loop over the number of viscous boundary subfaces of this block.
     bocos: do mm=1,nBocos
        isoWall: if (BCType(mm) == NSWallIsoThermal) then
           jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd
           do j=jBeg,jEnd
              do i=iBeg, iEnd
                 ii = ii + 1
                 BCData(mm)%TNS_Wall(i,j) = tnsw(ii)/TRef
              end do
           end do
        end if isoWall
     end do bocos
  end do domains
  ! TODO: The temperature must be interpolated to the coarse meshes.
  !
  ! The following lines are extracted from BCData/setBCDataCoarseGrid
  ! By the design of the subroutine, TNSWall shall be interpolated during this process.
  ! Yet, the subroutine requires an internal subroutine interpolateBcData.
  ! It remains a question whether interpolateBcData shall become a normal subroutine.
  !
  ! use blockPointers, only : flowDoms
  ! use inputTimeSpectral, only : nTimeIntervalsSpectral
  ! use iteration, only : groundLevel
  ! implicit none
  ! !
  ! !      Local variables.
  ! !
  ! integer(kind=intType) :: nLevels, level, levm1

  ! ! Determine the number of grid levels.

  ! nLevels = ubound(flowDoms,2)

  ! ! Loop over the coarser grid levels. It is assumed that the
  ! ! bc data of groundLevel is set correctly.

  ! coarseLevelLoop: do level=(groundLevel+1),nLevels

  !    ! Store the fine grid level a bit easier.

  !    levm1 = level - 1

  !    ! Loop over the number of spectral solutions and local blocks.

  !    spectralLoop: do sps=1,nTimeIntervalsSpectral
  !       domainsLoop: do i=1,nDom

  !          ! Set the pointers to the coarse block.

  !          call setPointers(i, level, sps)

  !          ! Loop over the boundary subfaces and interpolate the
  !          ! prescribed boundary data for this grid level.

  !          bocoLoop: do j=1,nBocos

  !             ! Interpolate the data for the possible prescribed boundary
  !             ! data.

  !             call interpolateBcData(BCData(j)%TNS_Wall, &
  !                  flowDoms(i,levm1,sps)%BCData(j)%TNS_Wall)

  !          enddo bocoLoop
  !       enddo domainsLoop
  !    enddo spectralLoop
  ! enddo coarseLevelLoop

end subroutine setTNSWall

subroutine getTNSWall(tnsw, npts, sps)

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType
  use flowVarRefState, only : TRef
  use utils, only : setPointers
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out) :: tnsw(npts)

  ! Local Variables
  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  ii = 0
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     ! Loop over the number of viscous boundary subfaces of this block.
     bocos: do mm=1,nBocos
        isoWall: if (BCType(mm) == NSWallIsoThermal) then
           do j=jBeg,jEnd
              do i=iBeg, iEnd
                 ii = ii + 1
                 tnsw(ii) = BCData(mm)%TNS_Wall(i,j)*Tref
              end do
           end do
        end if isoWall
     end do bocos
  end do domains
end subroutine getTNSWall
