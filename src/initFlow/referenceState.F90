!
!      ******************************************************************
!      *                                                                *
!      * File:          referenceState.f90                              *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 05-29-2003                                      *
!      * Last modified: 04-22-2006                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine referenceState
  !
  !      ******************************************************************
  !      *                                                                *
  !      * referenceState computes the reference state values in case     *
  !      * these have not been specified. A distinction is made between   *
  !      * internal and external flows. In case nothing has been          *
  !      * specified for the former a dimensional computation will be     *
  !      * made. For the latter the reference state is set to an          *
  !      * arbitrary state for an inviscid computation and computed for a *
  !      * viscous computation. Furthermore for internal flows an average *
  !      * velocity direction is computed from the boundary conditions,   *
  !      * which is used for initialization.                              *
  !      *                                                                *
  !      * The original version has been nuked since the computations are *
  !      * no longer necessary when calling from python                   *
  !      ******************************************************************
  !
  use BCTypes
  use block
  use communication
  use constants
#ifndef USE_TAPENADE
  use couplerParam
#endif
  use flowVarRefState
  use inputMotion
  use inputPhysics
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: sps, nn, mm

  real(kind=realType) :: gm1, ratio, tmp
  real(kind=realType) :: mx, my, mz, Re, v, TInfDim

  ! The following values MUST be set:
  ! pInfDim, reynolds, tempFreeStream

  ! pRef, rhoRef and Tref may be optionally set of if less than 0
  ! will take free stream values.

  TInfDim = tempFreestream
  rhoInfDim = pInfDim/(RGasDim*TInfDim)
  muDim     = muSuthDim                                    &
       * ((TSuthDim + SSuthDim)/(TInfDim + SSuthDim)) &
       * ((TInfDim/TSuthDim)**1.5_realType)

  ! External flow. Compute the value of gammaInf.

  call computeGamma(tempFreestream, gammaInf, 1)

  ! In case of a viscous problem, compute the
  ! dimensional free stream density and pressure.

  if(equations == NSEquations .or. &
       equations == RANSEquations) then

     ! Compute the x, y, and z-components of the Mach number
     ! relative to the body; i.e. the mesh velocity must be
     ! taken into account here.

     mx = MachCoef*velDirFreestream(1)
     my = MachCoef*velDirFreestream(2)
     mz = MachCoef*velDirFreestream(3)

     ! Reynolds number per meter, the viscosity using sutherland's
     ! law and the free stream velocity relative to the body.

     Re = Reynolds/ReynoldsLength
     muDim = muSuthDim                    &
          * ((TSuthDim + SSuthDim)       &
          /  (tempFreestream + SSuthDim)) &
          * ((tempFreestream/TSuthDim)**1.5)
     v  = sqrt((mx*mx + my*my + mz*mz) &
          *      gammaInf*RGasDim*tempFreestream)
     
     ! Compute the free stream density and pressure.
     ! Set TInfDim to tempFreestream.

     rhoInfDim = Re*muDim/v
     pInfDim   = rhoInfDim*RGasDim*tempFreestream

  endif

  ! In case the reference pressure, density and temperature were
  ! not specified, set them to the infinity values.

  if(pRef   <= zero) then 
     pRef   = pInfDim
  end if
  if(rhoRef <= zero) then 
     rhoRef = rhoInfDim
  end if
  if(TRef   <= zero) then 
     TRef   = TInfDim
  end if



  ! Compute the value of muRef, such that the nonDimensional
  ! equations are identical to the dimensional ones.
  ! Note that in the non-dimensionalization of muRef there is
  ! a reference length. However this reference length is 1.0
  ! in this code, because the coordinates are converted to
  ! meters.

  muRef = sqrt(pRef*rhoRef)

  ! Compute timeRef for a correct nonDimensionalization of the
  ! unsteady equations. Some story as for the reference viscosity
  ! concerning the reference length.

  timeRef = sqrt(rhoRef/pRef)

  ! Compute the nonDimensional pressure, density, velocity,
  ! viscosity and gas constant.

  pInf   = pInfDim/pRef
  rhoInf = rhoInfDim/rhoRef
  uInf   = Mach*sqrt(gammaInf*pInf/rhoInf)
  RGas   = RGasDim*rhoRef*TRef/pRef
  muInf  = muDim/muRef

end subroutine referenceState

!=================================================================

subroutine velMagnAndDirectionSubface(vmag, dir, BCData, mm)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * VelMagnAndDirectionSubface determines the maximum value    *
  !      * of the magnitude of the velocity as well as the sum of the     *
  !      * flow directions for the currently active subface.              *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: mm

  real(kind=realType), intent(out) :: vmag
  real(kind=realType), dimension(3), intent(inout) :: dir

  type(BCDataType), dimension(:), pointer :: BCData
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j
  real(kind=realType)   :: vel
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Initialize vmag to -1.0.

  vmag = -one

  ! Check if the velocity is prescribed.

  if( associated(BCData(mm)%velx) .and. &
       associated(BCData(mm)%vely) .and. &
       associated(BCData(mm)%velz) ) then

     ! Loop over the owned faces of the subface. As the cell range
     ! may contain halo values, the nodal range is used.

     do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
        do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

           ! Compute the magnitude of the velocity and compare it
           ! with the current maximum. Store the maximum of the two.

           vel  = sqrt(BCData(mm)%velx(i,j)**2 &
                +      BCData(mm)%vely(i,j)**2 &
                +      BCData(mm)%velz(i,j)**2)
           vmag = max(vmag, vel)

           ! Compute the unit vector of the velocity and add it to dir.

           vel    = one/max(eps,vel)
           dir(1) = dir(1) + vel*BCData(mm)%velx(i,j)
           dir(2) = dir(2) + vel*BCData(mm)%vely(i,j)
           dir(3) = dir(3) + vel*BCData(mm)%velz(i,j)

        enddo
     enddo
  endif

  ! Check if the velocity direction is prescribed.

  if( associated(BCData(mm)%flowXdirInlet) .and. &
       associated(BCData(mm)%flowYdirInlet) .and. &
       associated(BCData(mm)%flowZdirInlet) ) then

     ! Add the unit vectors to dir by looping over the owned
     ! faces of the subfaces. Again the nodal range must be
     ! used for this.

     do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
        do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

           dir(1) = dir(1) + BCData(mm)%flowXdirInlet(i,j)
           dir(2) = dir(2) + BCData(mm)%flowYdirInlet(i,j)
           dir(3) = dir(3) + BCData(mm)%flowZdirInlet(i,j)

        enddo
     enddo

  endif

end subroutine velMagnAndDirectionSubface
