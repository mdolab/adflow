

subroutine testRev(dwbar, wbar, m)
#ifndef USE_COMPLEX
  use BCTypes
  use blockPointers
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use communication
  use diffSizes
  use ADjointVars
  use cgnsGrid
  use block
  use inputiteration
  use adjointpetsc, only : psi_like3
  use bcroutines_b
  use samodule_fast_b
  !use samodule_b
  use paramturb

  implicit none

  ! Input Variables
  real(kind=realType), dimension(m), intent(in) :: dwbar
  real(kind=realType), dimension(m), intent(out) :: wbar
  integer(kind=intTYpe) :: m

#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  real(kind=realType),dimension(:),allocatable :: vec1, vec2 

  ! Local variables.
  integer(kind=intType) :: i, j, k, l, nn, ii, ierr, jj, irow, sps
  integer(kind=intType) :: nState, level
  real(kind=realType) :: timea, timeb, fwdTime, revTime, ovol
  logical :: resetToRANS

  ! call VecPlaceArray(psi_like3, wbar, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif
  ! Assembling matrix on coarser levels is not entirely implemented yet. 
  level = 1
  currentLevel = level
  groundLevel = level

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Allocate the memory we need for this block to do the forward
  ! mode derivatives and copy reference values
  if (.not. derivVarsAllocated) then 
     call alloc_derivative_values(level)
  end if
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        call zeroADSeeds(nn,level, sps)
     end do
  end do

  ii = 0
  sps = 1
  timeA = mpi_wtime()
  do nn=1,nDom
     ! Set pointers and derivative pointers
     call setPointers_d(nn, level, 1)

     do k=2,kl
        do j=2,jl
           do i=2,il
              ovol = one/vol(i,j,k)
              do l=1,nwf
                 dwd(i,j,k,l) = dwbar(ii+ l)*ovol
                 fwd(i,j,k,l) = dwd(i,j,k,l)
              end do
              do l=nt1,nt2
                 dwd(i,j,k,l) = dwbar(ii+ l)*ovol*turbresscale
                 fwd(i,j,k,l) = dwd(i,j,k,l)
              end do
              ii = ii + nw
           end do
        end do
     end do

     cv13 = rsacv1**3
     kar2inv = one/rsak**2
     cw36 = rsacw3**6
     cb3inv = one/rsacb3
     !call viscousFlux_fast_b

     !call allnodalgradients_fast_b
     !call computespeedofsoundsquared_fast_b
     call inviscidDissFluxScalar_fast_b
     call inviscidcentralflux_fast_b

     ! ! Turblent sa stuff
     ! call saresscale_fast_b()
     ! call saviscous_fast_b()
     ! call turbadvection_fast_b(1_inttype, 1_inttype, itu1-1, qq)
     ! call sasource_fast_b()

     ! select case  (turbprod) 
     ! case (strain) 
     !    call prodsmag2_fast_b()
     ! case (vorticity) 
     !    call prodwmag2_fast_b()
     ! case (katolaunder) 
     !    call prodkatolaunder_fast_b()
     ! end select

     call pushreal8array(radk, size(radk, 1)*size(radk, 2)*size(radk, 3))
     call pushreal8array(radj, size(radj, 1)*size(radj, 2)*size(radj, 3))
     call pushreal8array(radi, size(radi, 1)*size(radi, 2)*size(radi, 3))

     if (nn==1) then 
        !print *,'radi before:', sum(radi), sum(w)
     end if
     !call timestep_block(.False.)
     if (nn==1) then 
        !print *,'radi middle:', sum(radi), sum(w)
     end if
     call timestep_block_fast_b(.False.)
     if (nn==1) then 
        !print *,'radi after:', sum(radi), sum(w)
     end if
     !call timestep_block(.False.)
     !call applyallturbbcthisblock_b(.true.)
     !call bcturbtreatment_b()
     call applyallbc_block_b(.True.)
  call popreal8array(radi, size(radi, 1)*size(radi, 2)*size(radi, 3))
  call popreal8array(radj, size(radj, 1)*size(radj, 2)*size(radj, 3))
  call popreal8array(radk, size(radk, 1)*size(radk, 2)*size(radk, 3))

  end do

  call whalo2_b(1, 1, nw, .True., .True., .True.)
  ii = 0
  do nn=1,nDom
     call setPointers_d(nn, level, 1)
     !call saeddyviscosity_b
     !call computelamviscosity_fast_b
     call computepressuresimple_fast_b

     ! We can put stuff directly into wbar with no assembly; the
     ! whalo_b already takes care of it. 

     do k=2, kl
        do j=2,jl
           do i=2,il
              do l=1,nw
                 ii =ii + 1
                 wbar(ii) = flowdomsd(nn, level, sps)%w(i,j,k,l)
              end do
           end do
        end do
     end do
  end do

  timeB = mpi_wtime()
  !print *,'Fortran Time:', myid, timeB-timeA

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if

#endif
end subroutine testRev



subroutine testRevSpatial(dwbar, xbar, stateSize, spatialSize)
#ifndef USE_COMPLEX
  use BCTypes
  use blockPointers
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use communication
  use diffSizes
  use ADjointVars
  use cgnsGrid
  use block
  use inputiteration
  use adjointpetsc, only : psi_like3, x_like
  use bcroutines_b
  use samodule_b!spatial_b
  !use samodule_b
  use paramturb
  use walldistancedata, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input Variables
  integer(kind=intTYpe) :: stateSize, SpatialSize
  real(kind=realType), dimension(stateSize), intent(in) :: dwbar
  real(kind=realType), dimension(spatialSize), intent(out) :: xbar


  ! Local variables.
  integer(kind=intType) :: i, j, k, l, nn, ii, ierr, jj, irow, sps
  integer(kind=intType) :: nState, level
  real(kind=realType) :: timea, timeb, fwdTime, revTime, ovol, ovol2, ovold
  logical :: resetToRANS
  real(kind=realType), dimension(:), allocatable :: xSurfbSum

  call VecPlaceArray(x_like, xbar, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecSet(x_like, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif
  ! Assembling matrix on coarser levels is not entirely implemented yet. 
  level = 1
  currentLevel = level
  groundLevel = level

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if 
  call alloc_derivative_values(level)


  ! Now extract the vector of the surface data we need
  call VecGetArrayF90(xSurfVec(level), xSurf, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And it's derivative
  call VecGetArrayF90(xSurfVecd, xSurfd, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  allocate(xSurfbSum(size(xSurfd)))
  xSurfbSum = zero

  ! Allocate the memory we need for this block to do the forward
  ! mode derivatives and copy reference values
  
  call mpi_barrier(sumb_comm_world, ierr)
  ii = 0
  sps = 1
  timeA = mpi_wtime()
  call computeResidualNK
  do nn=1,nDom
     ! Set pointers and derivative pointers
     call setPointers_d(nn, level, 1)
     do k=2,kl
        do j=2,jl
           do i=2,il
              ovol = one/vol(i,j,k)
              ovold = zero
              do l=1,nwf
                 dwd(i,j,k,l) = dwbar(ii+ l)*ovol
                 fwd(i,j,k,l) = dwd(i,j,k,l)
                 ovold = ovold + (dw(i,j,k,l))*dwbar(ii+l)
              end do
              do l=nt1,nt2
                 dwd(i,j,k,l) = dwbar(ii+ l)*ovol*turbresscale
                 fwd(i,j,k,l) = dwd(i,j,k,l)
              end do
              ii = ii + nw

              vold(i,j,k) = vold(i,j,k) - ovold*(ovol*ovol)
           end do
        end do
     end do
     
     cv13 = rsacv1**3
     kar2inv = one/rsak**2
     cw36 = rsacw3**6
     cb3inv = one/rsacb3
     
     call viscousFlux_b!spatial_b
     call allnodalgradients_b!spatial_b
     call inviscidDissFluxScalar_b!spatial_b
     call inviscidcentralflux_b!spatial_b
     ! ! Turblent sa stuff
     ! call saresscale_b()!spatial_b()
     ! call saviscous_b!spatial_b()
     ! call turbadvection_b(1_inttype, 1_inttype, itu1-1, qq)
     ! call sasource_b()

     ! select case  (turbprod) 
     ! case (strain) 
     !    call prodsmag2_b!spatial_b()
     ! case (vorticity) 
     !    call prodwmag2_b!spatial_b()
     ! case (katolaunder) 
     !    call prodkatolaunder_b!spatial_b()
     ! end select

     call timestep_block_b(.False.)!spatial_b(.False.)
     ! call applyallturbbcthisblock_b(.true.)
     ! call bcturbtreatment_b()
     
     call applyallbc_block_b(.True.)
     !call computeeddyviscosity_b()
     !call computelamviscosity_b()
     !call computepressuresimple_b()
     !call updatewalldistancesquickly_b(nn, 1, 1)
     call boundarynormals_b()
     call metric_block_b()
     call volume_block_b()
     call xhalo_block_b()

     ! We can assemble stuff into x_like 

     do k=0,ke
        do j=0,je
           do i=0,ie
              do l=1,3
                 irow = flowDoms(nn, 1, 1)%globalNode(i,j,k)*3 + l -1
                 if (irow >= 0) then 
                    call VecSetValues(x_like, 1, (/irow/), &
                         (/flowdomsd(nn, level, sps)%x(i, j, k, l)/), ADD_VALUES, ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                 end if
              end do
           end do
        end do
     end do

     ! We need to acculumate the contribution from this block into xSurfbSum
     !xSurfbSum = xSurfbSum + xSurfd
  end do
  ! And perform assembly on the x vectors
  call VecAssemblyBegin(x_like, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecAssemblyEnd(x_like, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! call VecScatterBegin(wallScatter(1), xSurfVecd, x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)
  
  ! call VecScatterEnd(wallScatter(1), xSurfVecd, x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)
  
  call VecResetArray(x_like, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  timeB = mpi_wtime()
  !print *,'Fortran Time:', myid, timeB-timeA
  call dealloc_derivative_values(level)

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if

#endif
end subroutine testRevSpatial

! subroutine testRev
!   use BCTypes
!   use blockPointers
!   use inputDiscretization 
!   use inputTimeSpectral 
!   use inputPhysics
!   use iteration         
!   use flowVarRefState     
!   use inputAdjoint       
!   use communication
!   use diffSizes
!   use ADjointVars
!   use inputDiscretization
!   use cgnsGrid
!   use block
!   implicit none

!   ! Input Variables
! #define PETSC_AVOID_MPIF_H
! #include "include/finclude/petsc.h"
!   real(kind=realType),dimension(:),allocatable :: vec1, vec2 

!   ! Local variables.
!   integer(kind=intType) :: i, j, k, l, nn, ii, ierr, jj
!   integer(kind=intType) :: nState, level
!   real(kind=realType) :: timea, timeb, fwdTime, revTime
!   logical :: resetToRANS

! #ifndef USE_COMPLEX

!   ! Setup number of state variable based on turbulence assumption
!   if ( frozenTurbulence ) then
!      nState = nwf
!   else
!      nState = nw
!   endif
!   ! Assembling matrix on coarser levels is not entirely implemented yet. 
!   level = 1
!   currentLevel = level
!   groundLevel = level

!   ! If we are computing the jacobian for the RANS equations, we need
!   ! to make block_res think that we are evauluating the residual in a
!   ! fully coupled sense.  This is reset after this routine is
!   ! finished.
!   if (equations == RANSEquations) then
!      nMGVar = nw
!      nt1MG = nt1
!      nt2MG = nt2

!      turbSegregated = .False.
!      turbCoupled = .True.
!   end if

!   ! Determine if we want to use frozenTurbulent Adjoint
!   resetToRANS = .False. 
!   if (frozenTurbulence .and. equations == RANSEquations) then
!      equations = NSEquations 
!      resetToRANS = .True.
!   end if


!   ! Allocate the memory we need for this block to do the forward
!   ! mode derivatives and copy reference values
!   !call alloc_derivative_values( level)
!   call alloc_derivative_values(level)
!   call mpi_barrier(sumb_comm_world, ierr)
!   timeA = mpi_wtime()
!   do ii=1,10
!      call whalo2(1_intType, 1_intType, nw, .true., &
!           .true., .true.)

!      do nn=1,nDom
!         ! Run the regular forward call:

!         ! Set pointers to the first timeInstance...just to getSizes
!         call setPointers(nn, level, 1)

!         call computeLamViscosity
!         ! call computeEddyViscosity
!         ! call applyAllBC_block(.True.)
!         ! if (equations == RANSequations) &
!         !    call applyAllTurbBCThisBLock(.True.)

!         call timeStep_block(.False.)

!         ! if (equations == RANSEquations)  then
!         !    call sa_block(.true.)
!         !   end if

!         call inviscidCentralFlux
!         call inviscidDissFluxScalar
!         call computespeedofsoundsquared
!         call allnodalgradients
!         call viscousFlux
!      end do
!   end do

!   call mpi_barrier(sumb_comm_world, ierr)
!   timeB = mpi_wtime()
!   fwdTime = timeB-timeA
!   if (myid == 0) then   
!      print *,'Forward Time',fwdTime
!   end if
!   call mpi_barrier(sumb_comm_world, ierr)
!   timeA = mpi_wtime()
!   do ii=1,10
!      call whalo2_b(1_intType, 1_intType, nw, .true., &
!           .true., .true.)

!      do nn=1,nDom
!         ! Set pointers and derivative pointers
!         call setPointers_d(nn, level, 1)
!         fwd = one
!         dwd = one
!         call viscousFlux_fast_b
!         call allnodalgradients_fast_b
!         call computespeedofsoundsquared_fast_b
!         call inviscidDissFluxScalar_fast_b
!         call inviscidcentralflux_fast_b
!         call timestep_block_fast_b(.False.)
!         call computelamviscosity_fast_b
!      end do
!   end do
!   print *,sum(flowDomsd(1,1,1)%w)
!   call mpi_barrier(sumb_comm_world, ierr)
!   timeB = mpi_wtime()
!   revTime = timeB-timeA
!   if (myid == 0) then   
!      print *,'Reverse Time:', revTime
!      print *,'Ratio:', revTime/fwdTime
!   end if
!   call dealloc_derivative_values(level)

!   ! Reset the correct equation parameters if we were useing the frozen
!   ! Turbulent 
!   if (resetToRANS) then
!      equations = RANSEquations
!   end if

!   ! Reset the paraters to use segrated turbulence solve. 
!   if (equations == RANSEquations) then
!      nMGVar = nwf
!      nt1MG = nwf + 1
!      nt2MG = nwf

!      turbSegregated = .True.
!      turbCoupled = .False.
!      restrictEddyVis = .false.
!      if( eddyModel ) restrictEddyVis = .true.
!   end if

! #endif
! end subroutine testRev
