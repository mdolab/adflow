subroutine getSurfaceSize(size, sizeCell, famList, n)
  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use constants
  use blockPointers
  use inputTimeSpectral
  use surfaceFamilies
  use utils, only : setPointers
  use sorting, only : bsearchIntegers
  implicit none

  integer(kind=intType),intent(out) :: size, sizeCell
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd
  integer(kind=intType), intent(in) :: famList(n), n

  size = 0_intType
  sizeCell = 0_intType

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        ! Check if this surface should be included or not:
        famInclude: if (bsearchIntegers(BCdata(mm)%famID, famList, n) > 0) then 

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           sizeCell = sizeCell + (iEnd - iBeg)*(jEnd - jBeg)
        end if famInclude
     end do bocos
  end do domains
end subroutine getSurfaceSize

subroutine getSurfaceConnectivity(conn, ncell)
  ! Return the connectivity list for the each of the patches
  use constants
  use blockPointers
  use inputPhysics
  use surfaceFamilies
  use utils, only : setPointers
  use sorting, only : bsearchIntegers
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ncell
  integer(kind=intType), intent(inout) :: conn(4*ncell)

  ! Working
  integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  logical regularOrdering
  cellCount = 0
  nodeCount = 0

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)
     bocos: do mm=1,nBocos
        famInclude: if (bsearchIntegers(BCdata(mm)%famID, &
             famGroups, size(famGroups)) > 0) then 

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1

           ! We want to ensure that all the normals of the faces are
           ! consistent. To ensure this, we enforce that all normals
           ! are "into" the domain. Therefore we must treat difference
           ! faces of a block differently. For example for an iLow
           ! face, when looping over j-k in the regular way, results
           ! in in a domain inward pointing normal for iLow but
           ! outward pointing normal for iHigh. The same is true for
           ! kMin and kMax. However, it is reverse for the J-faces:
           ! This is becuase the way the pointers are extracted i then
           ! k is the reverse of what "should" be for consistency. The
           ! other two, the pointers are cyclic consistent: i,j->k,
           ! j,k (wrap) ->i, but for the j-direction is is i,k->j when
           ! to be consistent with the others it should be
           ! k,i->j. Hope that made sense. 

           select case(BCFaceID(mm))
           case(iMin, jMax, kMin)
              regularOrdering = .True.
           case default
              regularOrdering = .False.
           end select

           ! Now this can be reversed *again* if we have a block that
           ! is left handed. 
           if (.not. rightHanded) then 
              regularOrdering = .not. (regularOrdering)
           end if

           if (regularOrdering) then 
              ! Do regular ordering.

              ! Loop over generic face size...Note we are doing zero
              ! based ordering!

              ! This cartoon of a generic cell might help:
              !
              ! i, j+1 +-----+ i+1, j+1
              !   n4   |     | n3
              !        +-----+
              !       i,j    i+1, j
              !       n1     n2
              !

              do j=0,nj-2
                 do i=0,ni-2
                    conn(4*cellCount+1) = nodeCount + (j  )*ni + i       + 1! n1
                    conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1   + 1! n2
                    conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1   + 1! n3
                    conn(4*cellCount+4) = nodeCount + (j+1)*ni + i       + 1! n4
                    cellCount = cellCount + 1
                 end do
              end do
           else
              ! Do reverse ordering:
              do j=0,nj-2
                 do i=0,ni-2
                    conn(4*cellCount+1) = nodeCount + (j  )*ni + i        + 1! n1
                    conn(4*cellCount+2) = nodeCount + (j+1)*ni + i        + 1! n4
                    conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1    + 1! n3
                    conn(4*cellCount+4) = nodeCount + (j  )*ni + i + 1    + 1! n2
                    cellCount = cellCount + 1
                 end do
              end do
           end if
           nodeCount = nodeCount + ni*nj
        end if famInclude
     end do bocos
  end do domains
end subroutine getSurfaceConnectivity


subroutine getSurfaceFamily(elemFam, ncell)
  ! Return the connectivity list for the each of the patches
  use constants
  use blockPointers
  use inputPhysics
  use surfaceFamilies
  use utils, only : setPointers
  use sorting, only : bsearchIntegers
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ncell
  integer(kind=intType), intent(inout) :: elemFam(nCell)

  ! Working
  integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  logical regularOrdering
  cellCount = 0
  nodeCount = 0

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)
     bocos: do mm=1,nBocos
        famInclude: if (bsearchIntegers(BCdata(mm)%famID, &
             famGroups, size(famGroups)) > 0) then 

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0,nj-2
              do i=0,ni-2
                 cellCount = cellCount + 1
                 elemFam(cellCount) = BCdata(mm)%famID
              end do
           end do
           nodeCount = nodeCount + ni*nj
        end if famInclude
     end do bocos
  end do domains
end subroutine getSurfaceFamily


subroutine getSurfacePoints(points, npts, sps_in)
  use constants
  use blockPointers
  use inputTimeSpectral
  use surfaceFamilies
  use utils, only : setPointers
  use sorting, only : bsearchIntegers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts,sps_in
  real(kind=realType), intent(inout) :: points(3,npts)

  integer(kind=intType) :: mm, nn, i, j, ii,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd



  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos

        famInclude: if (bsearchIntegers(BCdata(mm)%famID, &
             famGroups, size(famGroups)) > 0) then 

           ! NODE Based
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           do j=jBeg, jEnd ! This is a node loop
              do i=iBeg, iEnd ! This is a node loop
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
        end if famInclude
     end do bocos
  end do domains

end subroutine getSurfacePoints

subroutine getForces(forces, npts, sps)
  use constants
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  use surfaceFamilies
  use utils, only : setPointers
  use sorting, only : bsearchIntegers
  implicit none

  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(inout) :: forces(3*npts)

  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: sss(3),v2(3),v1(3), qa, sepSensor, Cavitation
  real(kind=realType) :: sepSensorAvg(3)
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yplusmax, qf(3)

  ! Make sure *all* forces are computed. Sectioning will be done
  ! else-where.
  call setFullFamilyList()
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, &
          sepSensor, sepSensorAvg, Cavitation)
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
                    Forces(3*ii-2:3*ii) = bcData(mm)%T(i, j, :)
                 else
                    Forces(3*ii-2:3*ii) = bcData(mm)%F(i, j, :) 
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2
end subroutine getForces

subroutine computeNodalTractions(sps)
  use constants
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  use surfaceFamilies, only: wallExchange, wallFamilies, totalWallFamilies
  use utils, only : setPointers, EChk
  use sorting, only : bsearchIntegers
  implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  integer(kind=intType), intent(in) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf, qa
  real(kind=realType), dimension(:), pointer :: localPtr

  call vecGetArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)
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
 
  call vecRestoreArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the area
  call vecSet(wallExchange(sps)%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecScatterBegin(wallExchange(sps)%scatter, wallExchange(sps)%nodeValLocal, &
       wallExchange(sps)%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(wallExchange(sps)%scatter, wallExchange(sps)%nodeValLocal, &
       wallExchange(sps)%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing.

  call vecGetArrayF90(wallExchange(sps)%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = one/localPtr
  call vecRestoreArrayF90(wallExchange(sps)%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)


  ! Now do each of the three dimensions for the pressure and viscous forces
  dimLoop: do iDim=1, 6

     call vecGetArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)

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

     call vecRestoreArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Globalize the current force
     call vecSet(wallExchange(sps)%nodeValGlobal, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterBegin(wallExchange(sps)%scatter, wallExchange(sps)%nodeValLocal, &
          wallExchange(sps)%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(wallExchange(sps)%scatter, wallExchange(sps)%nodeValLocal, &
          wallExchange(sps)%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now divide by the dual area. We can do this with a vecpointwisemult
     call vecPointwiseMult(wallExchange(sps)%nodeValGlobal, wallExchange(sps)%nodeValGlobal, &
          wallExchange(sps)%sumGlobal, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Push back to the local values
     call VecScatterBegin(wallExchange(sps)%scatter, wallExchange(sps)%nodeValGlobal, &
          wallExchange(sps)%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(wallExchange(sps)%scatter, wallExchange(sps)%nodeValGlobal, &
          wallExchange(sps)%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecGetArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ii = 0
     do nn=1, nDom
        call setPointers(nn, 1_intType, sps)
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
                       bcData(mm)%Tp(i, j, iDim) = localPtr(ii) 
                    else
                       bcData(mm)%Tv(i, j, iDim-3) = localPtr(ii) 
                    end if
                 end do
              end do
           end if
        end do
     end do

     call vecRestoreArrayF90(wallExchange(sps)%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do dimLoop

  ! Finally sum the Tp and Tv together
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           do iDim=1,3
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    bcData(mm)%T(i, j, iDim) = bcData(mm)%Tp(i, j, iDim) + &
                         bcData(mm)%Tv(i, j, iDim)
                 end do
              end do
           end do
        end if
     end do
  end do

end subroutine computeNodalTractions

subroutine computeNodalForces(sps)

  ! This subroutine averages the cell based forces and tractions to
  ! node based values. There is no need for communication since we are
  ! simplying summing a quarter of each value to each corner. 
  use constants
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  use surfaceFamilies
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
           BCData(mM)%F = zero
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

subroutine getForces_b(forces_b, npts, sps)

  ! This routine performs the reverse of getForces. It takes in
  ! forces_b and perfroms the reverse of the nodal averaging procedure
  ! in getForces to compute bcDatad(mm)%Fp, bcDatad(mm)%Fv and
  ! bcDatad(mm)%area.
  use constants
  use blockPointers
  use inputPhysics
  use surfaceFamilies, only: wallExchange, familyExchange
  use communication
  use utils, only : EChk, setPointers

  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(in) :: forces_b(3,npts)
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf_b, qf, qa, qa_b
  real(kind=realType), dimension(:), pointer :: localPtr, localPtr_b
  type(familyExchange), pointer :: exch
  Vec nodeValLocal_b, nodeValGlobal_b, sumGlobal_b, tmp, tmp_b, T_b

  ! For better readibility
  exch => wallExchange(sps)

  if (.not. forcesAsTractions) then 
     ! For forces, we can accumulate the nodal seeds on the Fp and Fv
     ! values. The area seed is zeroed. 
 
     ii = 0
     domains: do nn=1,nDom
        call setPointers_d(nn, 1_intType, sps)

        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              BCDatad(mm)%Fp= zero
              BCDatad(mm)%Fv = zero
              do j=0,nj-2
                 do i=0,ni-2
                    do iDim=1,3
                       
                       ind(1) = ii + (j  )*ni + i + 1
                       ind(2) = ii + (j  )*ni + i + 2 
                       ind(3) = ii + (j+1)*ni + i + 2 
                       ind(4) = ii + (j+1)*ni + i + 1
                       qf_b = zero
                       do jj=1,4
                          qf_b = qf_b + forces_b(iDim, ind(jj))
                       end do
                       qf_b = qf_b*fourth
                       
                       BCDatad(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim) = & 
                            BCDatad(mm)%Fp(i+iBeg+1, j+jBeg+1, iDim) + qf_b
                       BCDatad(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim-3) = & 
                            BCDatad(mm)%Fv(i+iBeg+1, j+jBeg+1, iDim) + qf_b
                    end do
                 end do
              end do
           end if
        end do
     end do domains
  else

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

        do i=1, nPts
           if (iDim <= 3) then 
              localPtr_b(i) = forces_b(iDim, i)
           else
              localPtr_b(i) = forces_b(iDim-3, i)
           end if
        end do

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
        
        ! Accumulate seed on sumbGlobal_b
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

                 ! Zero the accumulation:
                 if (iDim <= 3) then 
                    BCDatad(mm)%Fp(:, :, iDim) = zero
                 else
                    BCDatad(mm)%Fv(:, :, iDim-3) = zero
                 end if

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
              BCDatad(mm)%area = zero
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

  end if

end subroutine getForces_b




