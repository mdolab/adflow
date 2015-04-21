subroutine slipVelocitiesFineLevel_ALE(useOldCoor, t, sps)
  !
  ! Shell function to call slipVelocitiesFineLevel on all blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps
  logical,               intent(in) :: useOldCoor
  real(kind=realType), dimension(*), intent(in) :: t  !
  !      Local variables.
  !
  integer(kind=intType) :: nn

  ! Loop over the number of blocks.

  domains: do nn=1,nDom

     ! Set the pointers for this block.

     call setPointers(nn, groundLevel, sps)

     call slipVelocitiesFineLevelALE_block(useOldCoor, t, sps)

  end do domains

end subroutine slipVelocitiesFineLevel_ALE


!
!      ******************************************************************
!      *                                                                *
!      * File:          slipVelocities.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-12-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine slipVelocitiesFineLevelALE_block(useOldCoor, t, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * slipVelocitiesFineLevel computes the slip velocities for       *
  !      * viscous subfaces on all viscous boundaries on groundLevel for  *
  !      * the given spectral solution. If useOldCoor is .true. the       *
  !      * velocities are determined using the unsteady time integrator;  *
  !      * otherwise the analytic form is used.                           *
  !      *                                                                *
  !      * Calculates the surface normal and normal velocity on BC using  *
  !      * FIRST order BDF.                                               *
  !      *                                                                *
  !      ******************************************************************
  !
  use BCTypes
  use inputTimeSpectral
  use blockPointers
  use cgnsGrid
  use flowVarRefState
  use inputMotion
  use inputUnsteady
  use iteration
  use inputPhysics
  use inputTSStabDeriv
  use monitor
  use communication
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps
  logical,               intent(in) :: useOldCoor

  real(kind=realType), dimension(*), intent(in) :: t
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, mm, i, j, level

  real(kind=realType) :: oneOver4dt
  real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
  real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

  real(kind=realType), dimension(3) :: xc, xxc
  real(kind=realType), dimension(3) :: rotCenter, rotRate

  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix,&
       derivRotationMatrix

  real(kind=realType) :: tNew, tOld

  real(kind=realType), dimension(:,:,:),   pointer :: uSlip
  real(kind=realType), dimension(:,:,:),   pointer :: xFace
  real(kind=realType), dimension(:,:,:,:), pointer :: xFaceOld

  integer(kind=intType) :: liftIndex
  real(kind=realType) :: alpha,beta,intervalMach,alphaTS,alphaIncrement,&
       betaTS,betaIncrement
  real(kind=realType), dimension(3) ::velDir
  real(kind=realType), dimension(3) :: refDirection

  !Function Definitions

  real(kind=realType) :: TSAlpha,TSBeta,TSMach
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the situation we are having here.

     ! *******************************
     ! REMOVED the rigid body rotation part for simplicity
     ! *******************************
     
     ! The velocities must be determined via a finite difference
     ! formula using the coordinates of the old levels.

     ! Set the coefficients for the time integrator and store the
     ! inverse of the physical nonDimensional time step, divided
     ! by 4, a bit easier.

     call setCoefTimeIntegrator_ALE
     oneOver4dt = fourth*timeRef/deltaT

     ! Loop over the number of viscous subfaces.

     bocoLoop1: do mm=1,nViscBocos

        ! Set the pointer for uSlip to make the code more
        ! readable.

        uSlip => BCData(mm)%uSlip

        ! Determine the grid face on which the subface is located
        ! and set some variables accordingly.

        select case (BCFaceID(mm))

        case (iMin)
           xFace => x(1,:,:,:);  xFaceOld => xOld(:,1,:,:,:)

        case (iMax)
           xFace => x(il,:,:,:); xFaceOld => xOld(:,il,:,:,:)

        case (jMin)
           xFace => x(:,1,:,:);  xFaceOld => xOld(:,:,1,:,:)

        case (jMax)
           xFace => x(:,jl,:,:); xFaceOld => xOld(:,:,jl,:,:)

        case (kMin)
           xFace => x(:,:,1,:);  xFaceOld => xOld(:,:,:,1,:)

        case (kMax)
           xFace => x(:,:,kl,:); xFaceOld => xOld(:,:,:,kl,:)

        end select

        ! Some boundary faces have a different rotation speed than
        ! the corresponding block. This happens e.g. in the tip gap
        ! region of turboMachinary problems where the casing does
        ! not rotate. As the coordinate difference corresponds to
        ! the rotation rate of the block, a correction must be
        ! computed. Therefore compute the difference in rotation
        ! rate and store the rotation center a bit easier. Note that
        ! the rotation center of subface is taken, because if there
        ! is a difference in rotation rate this info for the subface
        ! must always be specified.

        j = nbkGlobal
        i = cgnsSubface(mm)

        ! Loop over the quadrilateral faces of the viscous subface.
        ! Note that due to the usage of the pointers xFace and
        ! xFaceOld an offset of +1 must be used in the coordinate
        ! arrays, because x and xOld originally start at 0 for the
        ! i, j and k indices.

        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd

              ! Determine the coordinates of the centroid of the
              ! face, multiplied by 4.

              uSlip(i,j,1) = (xFace(i+1,j+1,1) + xFace(i+1,j,1) &
                            + xFace(i,  j+1,1) + xFace(i,  j,1))
              uSlip(i,j,2) = (xFace(i+1,j+1,2) + xFace(i+1,j,2) &
                            + xFace(i,  j+1,2) + xFace(i,  j,2))
              uSlip(i,j,3) = (xFace(i+1,j+1,3) + xFace(i+1,j,3) &
                            + xFace(i,  j+1,3) + xFace(i,  j,3))

              ! Loop over the older time levels and take their
              ! contribution into account.

              level = 1 ! There was a loop over all old levels
              uSlip(i,j,1) = uSlip(i,j,1) &
                   +          (xFaceOld(level,i+1,j+1,1) &
                   +           xFaceOld(level,i+1,j,  1) &
                   +           xFaceOld(level,i,  j+1,1) &
                   +           xFaceOld(level,i,  j,  1)) &
                   * (-1.0_realType)
              uSlip(i,j,2) = uSlip(i,j,2) &
                   +          (xFaceOld(level,i+1,j+1,2) &
                   +           xFaceOld(level,i+1,j,  2) &
                   +           xFaceOld(level,i,  j+1,2) &
                   +           xFaceOld(level,i,  j,  2)) &
                   * (-1.0_realType)
              uSlip(i,j,3) = uSlip(i,j,3) &
                   +          (xFaceOld(level,i+1,j+1,3) &
                   +           xFaceOld(level,i+1,j,  3) &
                   +           xFaceOld(level,i,  j+1,3) &
                   +           xFaceOld(level,i,  j,  3)) &
                   * (-1.0_realType)

              ! Divide by 4 times the time step to obtain the
              ! correct velocity.

              uSlip(i,j,1) = uSlip(i,j,1)*oneOver4dt
              uSlip(i,j,2) = uSlip(i,j,2)*oneOver4dt
              uSlip(i,j,3) = uSlip(i,j,3)*oneOver4dt
           enddo
        enddo

     enddo bocoLoop1

end subroutine slipVelocitiesFineLevelALE_block




!
!      ******************************************************************
!      Assuming coarse level grid has nothing to do with ALE
!      So the original one will be used
!      ******************************************************************
