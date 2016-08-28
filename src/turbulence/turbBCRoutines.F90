module turbBCRoutines

contains
#ifndef USE_TAPENADE
  subroutine applyAllTurbBC(secondHalo)
    !
    !       applyAllTurbBC applies all boundary conditions to the          
    !       turbulent transport equations for the all blocks on the grid   
    !       level currentLevel.                                            
    !
    use constants
    use blockPointers
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: secondHalo
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, sps

    ! Loop over the number of spectral modes and local blocks.

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom

          ! Set the pointers to this block. The min function is present
          ! because this routine can be called from movfin.

          call setPointers(nn, min(currentLevel,groundLevel), sps)

          ! Set the arrays for the boundary condition treatment
          ! and set the turbulent halo values.

          call bcTurbTreatment
          call applyAllTurbBCThisBlock(secondHalo)

       enddo
    enddo

  end subroutine applyAllTurbBC
#endif
  !      ==================================================================

  subroutine applyAllTurbBCThisBlock(secondHalo)
    !
    !       applyAllTurbBCThisBlock sets the halo values of the            
    !       turbulent variables and eddy viscosity for the block the       
    !       variables in blockPointers currently point to.                 
    !
    use constants
    use blockPointers
    use flowVarRefState
    use inputPhysics

    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: secondHalo
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, i, j, l, m

    real(kind=realType), dimension(:,:,:,:), pointer :: bmt
    real(kind=realType), dimension(:,:,:),   pointer :: bvt, ww1, ww2

    ! Loop over the boundary condition subfaces of this block.

    bocos: do nn=1,nBocos

       ! Loop over the faces and set the state in
       ! the turbulent halo cells.

       if( wallFunctions ) then
#ifndef USE_TAPENADE        
          ! Determine the block face on which this subface is located
          ! and set some pointers accordingly.

          select case (BCFaceID(nn))
          case (iMin)
             bmt => bmti1; bvt => bvti1
             ww1 => w(1 ,1:,1:,:); ww2 => w(2 ,1:,1:,:)

          case (iMax)
             bmt => bmti2; bvt => bvti2
             ww1 => w(ie,1:,1:,:); ww2 => w(il,1:,1:,:)

          case (jMin)
             bmt => bmtj1; bvt => bvtj1
             ww1 => w(1:,1 ,1:,:); ww2 => w(1:,2 ,1:,:)

          case (jMax)
             bmt => bmtj2; bvt => bvtj2
             ww1 => w(1:,je,1:,:); ww2 => w(1:,jl,1:,:)

          case (kMin)
             bmt => bmtk1; bvt => bvtk1
             ww1 => w(1:,1:,1 ,:); ww2 => w(1:,1:,2 ,:)

          case (kMax)
             bmt => bmtk2; bvt => bvtk2
             ww1 => w(1:,1:,ke,:); ww2 => w(1:,1:,kl,:)
          end select

          ! Write an approximate value into the halo cell for
          ! postprocessing (it is not used in computation).

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                do l=nt1,nt2
                   ww1(i,j,l) = bvt(i,j,l) - bmt(i,j,l,l)*ww2(i,j,l)
                   do m=nt1,nt2
                      if(m /= l .and. bmt(i,j,l,m) /= zero) &
                           ww1(i,j,l) = ww2(i,j,l)
                   enddo
                enddo
             enddo
          enddo
#endif
       else

          select case (BCFaceID(nn))
          case (iMin)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(1,i,j,l) = bvti1(i,j,l)
                      do m=nt1,nt2
                         w(1,i,j,l) = w(1,i,j,l) - bmti1(i,j,l,m)*w(2,i,j,m)
                      enddo
                   enddo
                enddo
             enddo

          case (iMax)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(ie,i,j,l) = bvti2(i,j,l)
                      do m=nt1,nt2
                         w(ie,i,j,l) = w(ie,i,j,l) - bmti2(i,j,l,m)*w(il,i,j,m)
                      enddo
                   enddo
                enddo
             enddo

          case (jMin)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(i,1,j,l) = bvtj1(i,j,l)
                      do m=nt1,nt2
                         w(i,1,j,l) = w(i,1,j,l) - bmtj1(i,j,l,m)*w(i,2,j,m)
                      enddo
                   enddo
                enddo
             enddo

          case (jMax)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(i,je,j,l) = bvtj2(i,j,l)
                      do m=nt1,nt2
                         w(i,je,j,l) = w(i,je,j,l) - bmtj2(i,j,l,m)*w(i,jl,j,m)
                      enddo
                   enddo
                enddo
             enddo

          case (kMin)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(i,j,1,l) = bvtk1(i,j,l)
                      do m=nt1,nt2
                         w(i,j,1,l) = w(i,j,1,l) - bmtk1(i,j,l,m)*w(i,j,2,m)
                      enddo
                   enddo
                enddo
             enddo

          case (kMax)
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   do l=nt1,nt2
                      w(i,j,ke,l) = bvtk2(i,j,l)
                      do m=nt1,nt2
                         w(i,j,ke,l) = w(i,j,ke,l) - bmtk2(i,j,l,m)*w(i,j,kl,m)
                      enddo
                   enddo
                enddo
             enddo
          end select

       endif

       ! Set the value of the eddy viscosity, depending on the type of
       ! boundary condition. Only if the turbulence model is an eddy
       ! viscosity model of course.

       if( eddyModel ) then

          if(BCType(nn) == NSWallAdiabatic .or. &
               BCType(nn) == NSWallIsothermal) then

             ! Viscous wall boundary condition. Eddy viscosity is
             ! zero at the wall.

             call bcEddyWall(nn)

          else

             ! Any boundary condition but viscous wall. A homogeneous
             ! Neumann condition is applied to the eddy viscosity.

             call bcEddyNoWall(nn)

          endif

       endif

       ! Extrapolate the turbulent variables in case a second halo
       ! is needed.

       if( secondHalo ) call turb2ndHalo(nn)

    enddo bocos

  end subroutine applyAllTurbBCThisBlock

  subroutine bcEddyNoWall(nn)
    !
    !       bcEddyNoWall sets the eddy viscosity in the halo cells of      
    !       subface nn of the block given in blockPointers. The boundary   
    !       condition on the subface can be anything but a viscous wall.   
    !       A homogeneous neumann condition is applied, which means that   
    !       the eddy viscosity is simply copied from the interior cell.    
    !
    use constants
    use blockPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j


    ! Determine the face id on which the subface and copy

    select case (BCFaceid(nn))
    case (iMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(1,i,j) = rev(2,i,j)
          enddo
       enddo

    case (iMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(ie,i,j) = rev(il,i,j)
          enddo
       enddo

    case (jMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,1,j) = rev(i,2,j)
          enddo
       enddo

    case (jMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,je,j) = rev(i,jl,j)
          enddo
       enddo

    case (kMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,j,1) = rev(i,j,2)
          enddo
       enddo

    case (kMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,j,ke) = rev(i,j,kl)
          enddo
       enddo
    end select

  end subroutine bcEddyNoWall
  subroutine bcEddyWall(nn)
    !
    !       bcEddyWall sets the eddy viscosity in the halo cells of        
    !       viscous subface nn of the block given in blockPointers.        
    !       As the eddy viscosity is zero at the wall, the value in the    
    !       halo is simply the negative value of the first interior cell.  
    !
    use constants
    use blockPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j


    ! Determine the face id on which the subface is located and
    ! loop over the faces of the subface and set the eddy viscosity
    ! in the halo cells.

    select case (BCFaceid(nn))
    case (iMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(1,i,j) = -rev(2,i,j)
          enddo
       enddo

    case (iMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(ie,i,j) = -rev(il,i,j)
          enddo
       enddo

    case (jMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,1,j) = -rev(i,2,j)
          enddo
       enddo

    case (jMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,je,j) = -rev(i,jl,j)
          enddo
       enddo

    case (kMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,j,1) = -rev(i,j,2)
          enddo
       enddo

    case (kMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             rev(i,j,ke) = -rev(i,j,kl)
          enddo
       enddo
    end select

  end subroutine bcEddyWall
  subroutine bcTurbFarfield(nn)
    !
    !       bcTurbFarfield applies the implicit treatment of the           
    !       farfield boundary condition to subface nn. As the farfield     
    !       boundary condition is independent of the turbulence model,     
    !       this routine is valid for all models. It is assumed that the   
    !       pointers in blockPointers are already set to the correct       
    !       block on the correct grid level.                               
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l

    real(kind=realType) :: nnx, nny, nnz, dot

    ! Loop over the faces of the subfaces and set the values of
    ! bmt and bvt for an implicit treatment.

    do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
       do i=BCData(nn)%icBeg, BCData(nn)%icEnd

          ! Determine the dot product between the outward pointing
          ! normal and the free stream velocity direction and add the
          ! possible grid velocity.

          dot = BCData(nn)%norm(i,j,1)*wInf(ivx) + &
               BCData(nn)%norm(i,j,2)*wInf(ivy) + &
               BCData(nn)%norm(i,j,3)*wInf(ivz) - BCData(nn)%rface(i,j)

          ! Determine whether we are dealing with an inflow or
          ! outflow boundary here.

          if(dot > zero) then

             ! Outflow. Simply extrapolation or zero Neumann BC
             ! of the turbulent variables.

             do l=nt1,nt2
                select case (BCFaceID(nn))
                case (iMin)
                   bmti1(i,j,l,l) = -one
                case (iMax)
                   bmti2(i,j,l,l) = -one
                case (jMin)
                   bmtj1(i,j,l,l) = -one
                case (jMax)
                   bmtj2(i,j,l,l) = -one
                case (kMin)
                   bmtk1(i,j,l,l) = -one
                case (kMax)
                   bmtk2(i,j,l,l) = -one
                end select
             end do

          else

             ! Inflow. Turbulent variables are prescribed.

             do l=nt1,nt2
                select case(BCFaceID(nn))
                case (iMin)
                   bvti1(i,j,l) = wInf(l)
                case (iMax)
                   bvti2(i,j,l) = wInf(l)
                case (jMin)
                   bvtj1(i,j,l) = wInf(l)
                case (jMax)
                   bvtj2(i,j,l) = wInf(l)
                case (kMin)
                   bvtk1(i,j,l) = wInf(l)
                case (kMax)
                   bvtk2(i,j,l) = wInf(l)
                end select
             enddo
          endif
       enddo
    enddo
  end subroutine bcTurbFarfield
  subroutine bcTurbInflow(nn)
    !
    !       bcTurbInflow applies the implicit treatment of the inflow      
    !       boundary conditions to subface nn. As the inflow boundary      
    !       condition is independent of the turbulence model, this routine 
    !       is valid for all models. It is assumed that the pointers in    
    !       blockPointers are already set to the correct block on the      
    !       correct grid level.                                            
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l


    ! Loop over the faces of the subfaces and set the values of
    ! bvt and bmt such that the inflow state is linearly extrapolated
    ! with a fixed state at the face.

    do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
       do i=BCData(nn)%icBeg, BCData(nn)%icEnd

          ! Loop over the number of turbulent variables.

          do l=nt1,nt2
             select case (BCFaceID(nn))
             case (iMin)
                bvti1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmti1(i,j,l,l) = one
             case (iMax)
                bvti2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmti2(i,j,l,l) = one
             case (jMin)
                bvtj1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmtj1(i,j,l,l) = one
             case (jMax)
                bvtj2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmtj2(i,j,l,l) = one
             case (kMin)
                bvtk1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmtk1(i,j,l,l) = one
             case (kMax)
                bvtk2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
                bmtk2(i,j,l,l) = one
             end select
          end do
       enddo
    enddo
  end subroutine bcTurbInflow
  subroutine bcTurbInterface(nn)
    !
    !       bcTurbInterface applies the halo treatment for interface halo  
    !       cells, sliding mesh interface and domain interface. As these   
    !       are not really boundary conditions, the variable bvt is simply 
    !       set to keep the current value.                                 
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l


    ! Loop over the faces of the subfaces and set the values of
    ! bvt to keep the current value.

    do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
       do i=BCData(nn)%icBeg, BCData(nn)%icEnd
          do l=nt1,nt2
             select case (BCFaceID(nn))
             case (iMin)
                bvti1(i,j,l) = w(1,i,j,l)
             case (iMax)
                bvti2(i,j,l) = w(ie,i,j,l)
             case (jMin)
                bvtj1(i,j,l) = w(i,1,j,l)
             case (jMax)
                bvtj2(i,j,l) = w(i,je,j,l)
             case (kMin)
                bvtk1(i,j,l) = w(i,j,1,l)
             case (kMax)
                bvtk2(i,j,l) = w(i,j,ke,l)
             end select
          enddo
       enddo
    enddo

    ! Note that the original code had an error in the pointers...they
    ! were pointing to {il,jl,kl} and not {ie, je, ke}.

  end subroutine bcTurbInterface
  subroutine bcTurbOutflow(nn)
    !
    !       bcTurbOutflow applies the implicit treatment of the outflow    
    !       boundary conditions to subface nn. As the outflow boundary     
    !       condition is independent of the turbulence model, either       
    !       extrapolation or zero Neumann, this routine is valid for all   
    !       models. It is assumed that the pointers in blockPointers are   
    !       already set to the correct block on the correct grid level.    
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l

    ! Loop over the faces of the subfaces and set the values of bmt
    ! for an implicit treatment. For an outflow the turbulent variable
    ! variable is either extrapolated or zero Neumann. As constant
    ! extrapolation is used this leads to an identical treatment, i.e.
    ! the halo value is identical to the value of the internal cell.

    do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
       do i=BCData(nn)%icBeg, BCData(nn)%icEnd
          do l=nt1,nt2
             select case (BCFaceID(nn))
             case (iMin)
                bmti1(i,j,l,l) = -one
             case (iMax)
                bmti2(i,j,l,l) = -one
             case (jMin)
                bmtj1(i,j,l,l) = -one
             case (jMax)
                bmtj2(i,j,l,l) = -one
             case (kMin)
                bmtk1(i,j,l,l) = -one
             case (kMax)
                bmtk2(i,j,l,l) = -one
             end select
          enddo
       enddo
    enddo

  end subroutine bcTurbOutflow
  subroutine bcTurbSymm(nn)
    !
    !       bcTurbSymm applies the implicit treatment of the symmetry      
    !       boundary condition (or inviscid wall) to subface nn. As the    
    !       symmetry boundary condition is independent of the turbulence   
    !       model, this routine is valid for all models. It is assumed     
    !       that the pointers in blockPointers are already set to the      
    !       correct block on the correct grid level.                       
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l

    ! Loop over the faces of the subfaces and set the values of bmt
    ! for an implicit treatment. For a symmetry face this means
    ! that the halo value is set to the internal value.

    do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
       do i=BCData(nn)%icBeg, BCData(nn)%icEnd
          do l=nt1,nt2
             select case (BCFaceID(nn))
             case (iMin)
                bmti1(i,j,l,l) = -one
             case (iMax)
                bmti2(i,j,l,l) = -one
             case (jMin)
                bmtj1(i,j,l,l) = -one
             case (jMax)
                bmtj2(i,j,l,l) = -one
             case (kMin)
                bmtk1(i,j,l,l) = -one
             case (kMax)
                bmtk2(i,j,l,l) = -one
             end select
          enddo
       enddo
    enddo
  end subroutine bcTurbSymm

  subroutine bcTurbTreatment
    !
    !       bcTurbTreatment sets the arrays bmti1, bvti1, etc, such that   
    !       the physical boundary conditions are treated correctly.        
    !       It is assumed that the variables in blockPointers already      
    !       point to the correct block.                                    
    !       The turbulent variable in the halo is computed as follows:     
    !       wHalo = -bmt*wInternal + bvt for every block facer. As it is   
    !       possible to have a coupling in the boundary conditions bmt     
    !       actually are matrices. If there is no coupling between the     
    !       boundary conditions of the turbulence equations bmt is a       
    !       diagonal matrix.                                               
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Local variable.
    !
    integer(kind=intType) :: nn, i, j, k, l, m

    ! Initialize the arrays for the boundary condition treatment
    ! to zero, such that internal block boundaries are solved
    ! correctly (i.e. explicitly).

    do k=1,ke
       do j=1,je
          do l=nt1,nt2
             do m=nt1,nt2
                bmti1(j,k,l,m) = zero
                bmti2(j,k,l,m) = zero
             enddo
             bvti1(j,k,l) = zero
             bvti2(j,k,l) = zero
          enddo
       enddo
    enddo

    do k=1,ke
       do i=1,ie
          do l=nt1,nt2
             do m=nt1,nt2
                bmtj1(i,k,l,m) = zero
                bmtj2(i,k,l,m) = zero
             enddo
             bvtj1(i,k,l) = zero
             bvtj2(i,k,l) = zero
          enddo
       enddo
    enddo

    do j=1,je
       do i=1,ie
          do l=nt1,nt2
             do m=nt1,nt2
                bmtk1(i,j,l,m) = zero
                bmtk2(i,j,l,m) = zero
             enddo
             bvtk1(i,j,l) = zero
             bvtk2(i,j,l) = zero
          enddo
       enddo
    enddo

    ! Loop over the boundary condition subfaces of this block.

    bocos: do nn=1,nBocos

       ! Determine the kind of boundary condition for this subface.

       typeBC: select case (BCType(nn))

       case (NSWallAdiabatic, NSWallIsothermal)

          ! Viscous wall. There is no difference between an adiabatic
          ! and an isothermal wall for the turbulent equations.
          ! Set the implicit treatment of the wall boundary conditions.

          call bcTurbWall(nn)

          !=============================================================
#ifndef USE_TAPENADE
       case (SubsonicInflow, SupersonicInflow, MassBleedInflow)

          ! Inflow. Subsonic, supersonic or mass bleed inflow is
          ! identical for the turbulent transport equations.

          call bcTurbInflow(nn)

          !=============================================================

       case (SubsonicOutflow,  SupersonicOutflow, &
            MassBleedOutflow, Extrap)

          ! Outflow. Subsonic, supersonic or mass bleed outflow is
          ! identical for the turbulent transport equations. The
          ! extrapolation boundary is also an outflow.

          call bcTurbOutflow(nn)
#endif
          !=============================================================

       case (Symm, SymmPolar, EulerWall)

          ! Symmetry, polar symmetry or inviscid wall. Treatment of
          ! the turbulent equations is identical.

          call bcTurbSymm(nn)

          !=============================================================

       case (FarField)

          ! Farfield. The kind of boundary condition to be applied,
          ! inflow or outflow, depends on the local conditions.

          call bcTurbFarfield(nn)

          !=============================================================

       case (SlidingInterface,   OversetOuterBound,     &
            DomainInterfaceAll, DomainInterfaceRhoUVW, &
            DomainInterfaceP,   DomainInterfaceRho,    &
            DomainInterfaceTotal)

          ! Sliding mesh interface, overset outer boudaries, and 
          ! domain interface with another code are not really boundary
          ! condition and therefore the values are kept.

          call bcTurbInterface(nn)

       end select typeBC

    enddo bocos

  end subroutine bcTurbTreatment
  subroutine bcTurbWall(nn)
    !
    !       bcTurbWall applies the implicit treatment of the viscous       
    !       wall boundary condition for the turbulence model used to the   
    !       given subface nn.                                              
    !       It is assumed that the pointers in blockPointers are           
    !       already set to the correct block.                              
    !
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use constants
    use paramTurb
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, ii, jj, iiMax, jjMax

    real(kind=realType) :: tmpd, tmpe, tmpf, nu

    real(kind=realType), dimension(:,:,:,:), pointer :: bmt
    real(kind=realType), dimension(:,:,:),   pointer :: bvt, ww2
    real(kind=realType), dimension(:,:),     pointer :: rlv2, dd2Wall



    ! Determine the turbulence model used and loop over the faces
    ! of the subface and set the values of bmt and bvt for an
    ! implicit treatment.

    select case (turbModel)

    case (spalartAllmaras, spalartAllmarasEdwards)

       ! Spalart-allmaras type of model. Value at the wall is zero,
       ! so simply negate the internal value.
       select case (BCFaceID(nn))
       case (iMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmti1(i,j,itu1,itu1) = one
             enddo
          enddo
       case (iMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmti2(i,j,itu1,itu1) = one
             enddo
          enddo
       case (jMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtj1(i,j,itu1,itu1) = one
             enddo
          enddo
       case (jMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtj2(i,j,itu1,itu1) = one
             enddo
          enddo

       case (kMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtk1(i,j,itu1,itu1) = one
             enddo
          enddo

       case (kMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtk2(i,j,itu1,itu1) = one
             enddo
          enddo
       end select

       !        ================================================================
    case (komegaWilcox, komegaModified, menterSST)

       ! K-omega type of models. K is zero on the wall and thus the
       ! halo value is the negative of the first internal cell.
       ! For omega the situation is a bit more complicated.
       ! Theoretically omega is infinity, but it is set to a large
       ! value, see menter's paper. The halo value is constructed
       ! such that the wall value is correct. Make sure that i and j
       ! are limited to physical dimensions of the face for the wall
       ! distance. Due to the usage of the dd2Wall pointer and the
       ! fact that the original d2Wall array starts at 2, there is
       ! an offset of -1 present in dd2Wall.

       select case (BCFaceID(nn))
       case (iMin)
          iiMax = jl; jjMax = kl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(2,i,j)/w(2,i,j,irho)
                tmpd = one/(rkwBeta1*(d2Wall(2,ii,jj)**2))

                bmti1(i,j,itu1,itu1) = one
                bmti1(i,j,itu2,itu2) = one

                bvti1(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo

       case (iMax)
          iiMax = jl; jjMax = kl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(jl,i,j)/w(il,i,j,irho)
                tmpd = one/(rkwBeta1*(d2Wall(il,ii,jj)**2))

                bmti2(i,j,itu1,itu1) = one
                bmti2(i,j,itu2,itu2) = one

                bvti2(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo

       case (jMin)
          iiMax = il; jjMax = kl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(i,2,j)/w(i,2,j,irho)
                tmpd = one/(rkwBeta1*(d2Wall(ii,2,jj)**2))

                bmtj1(i,j,itu1,itu1) = one
                bmtj1(i,j,itu2,itu2) = one

                bvtj1(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo

       case (jMax)
          iiMax = il; jjMax = kl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(i,jl,j)/w(i,jl,j,irho)
                tmpd = one/(rkwBeta1*(d2Wall(ii,jl,jj)**2))

                bmtj2(i,j,itu1,itu1) = one
                bmtj2(i,j,itu2,itu2) = one

                bvtj2(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo

       case (kMin)
          iiMax = il; jjMax = jl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(i,j,2)/w(i,j,2,irho)
                tmpd = one/(rkwBeta1*(d2Wall(ii,jj,2)**2))

                bmtk1(i,j,itu1,itu1) = one
                bmtk1(i,j,itu2,itu2) = one

                bvtk1(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo

       case (kMax)
          iiMax = il; jjMax = jl

          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                ii = max(2,min(i,iiMax))

                nu   = rlv(i,j,kl)/w(i,j,kl,irho)
                tmpd = one/(rkwBeta1*(d2Wall(ii,jj,kl)**2))

                bmtk2(i,j,itu1,itu1) = one
                bmtk2(i,j,itu2,itu2) = one

                bvtk2(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
          enddo
       end select

       !        ================================================================

    case (ktau)

       ! K-tau model. Both k and tau are zero at the wall, so the
       ! negative value of the internal cell is taken for the halo.
       select case (BCFaceID(nn))
       case (iMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmti1(i,j,itu1,itu1) = one
                bmti1(i,j,itu2,itu2) = one
             enddo
          enddo
       case (iMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmti2(i,j,itu1,itu1) = one
                bmti2(i,j,itu2,itu2) = one
             enddo
          enddo
       case (jMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtj1(i,j,itu1,itu1) = one
                bmtj1(i,j,itu2,itu2) = one
             enddo
          enddo
       case (jMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtj2(i,j,itu1,itu1) = one
                bmtj2(i,j,itu2,itu2) = one
             enddo
          enddo

       case (kMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtk1(i,j,itu1,itu1) = one
                bmtk1(i,j,itu2,itu2) = one
             enddo
          enddo

       case (kMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                bmtk2(i,j,itu1,itu1) = one
                bmtk2(i,j,itu2,itu2) = one
             enddo
          enddo
       end select

       !        ================================================================
#ifndef USE_TAPENADE
    case (v2f)

       ! Set some variables depending on the block face on which the
       ! subface is located. Needed for a general treatment.

       select case (BCFaceID(nn))
       case (iMin)
          iiMax = jl; jjMax = kl
          bmt => bmti1; bvt => bvti1; ww2 => w(2 ,1:,1:,:)
          rlv2 => rlv(2, 1:,1:); dd2Wall => d2Wall(2, :,:)

       case (iMax)
          iiMax = jl; jjMax = kl
          bmt => bmti2; bvt => bvti2; ww2 => w(il,1:,1:,:)
          rlv2 => rlv(il,1:,1:); dd2Wall => d2Wall(il,:,:)

       case (jMin)
          iiMax = il; jjMax = kl
          bmt => bmtj1; bvt => bvtj1; ww2 => w(1:,2 ,1:,:)
          rlv2 => rlv(1:,2 ,1:); dd2Wall => d2Wall(:,2 ,:)

       case (jMax)
          iiMax = il; jjMax = kl
          bmt => bmtj2; bvt => bvtj2; ww2 => w(1:,jl,1:,:)
          rlv2 => rlv(1:,jl,1:); dd2Wall => d2Wall(:,jl,:)

       case (kMin)
          iiMax = il; jjMax = jl
          bmt => bmtk1; bvt => bvtk1; ww2 => w(1:,1:,2 ,:)
          rlv2 => rlv(1:,1:,2 ); dd2Wall => d2Wall(:,:,2 )

       case (kMax)
          iiMax = il; jjMax = jl
          bmt => bmtk2; bvt => bvtk2; ww2 => w(1:,1:,kl,:)
          rlv2 => rlv(1:,1:,kl); dd2Wall => d2Wall(:,:,kl)
       end select

       ! V2f turbulence model. Same story for the wall distance as
       ! for k-omega. For this model there is a coupling between the
       ! equations via the boundary conditions.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          jj = max(2,min(j,jjMax))

          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             ii = max(2,min(i,iiMax))

             nu   = rlv2(i,j)/ww2(i,j,irho)
             tmpd = one/(dd2Wall(ii-1,jj-1)**2)
             tmpe = two*nu*tmpd
             tmpf =-20.0_realType*(nu*tmpd)**2 &
                  / abs(tmpe*ww2(i,j,itu1))
             if(rvfN == 6) tmpf = zero

             bmt(i,j,itu1,itu1) = one
             bmt(i,j,itu2,itu2) = one
             bmt(i,j,itu3,itu3) = one
             bmt(i,j,itu4,itu4) = one

             bmt(i,j,itu2,itu1) = -two*tmpe
             bmt(i,j,itu4,itu3) = -two*tmpf
          enddo
       enddo
#endif
    end select
  end subroutine bcTurbWall

  subroutine turb2ndHalo(nn)
    !
    !       turb2ndHalo sets the turbulent variables in the second halo    
    !       cell for the given subface. Simple constant extrapolation is   
    !       used to avoid problems.                                        
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l

    ! Determine the face on which this subface is located and set
    ! some pointers accordingly.

    ! Loop over the turbulent variables and set the second halo
    ! value. If this is an eddy model, also set the eddy viscosity.

    select case (BCFaceID(nn))
    case (iMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(0,i,j,l) = w(1,i,j,l)
             enddo
             if( eddyModel ) rev(0,i,j) = rev(1,i,j)
          enddo
       enddo

       !===============================================================

    case (iMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(ib,i,j,l) = w(ie,i,j,l)
             enddo
             if( eddyModel ) rev(ib,i,j) = rev(ie,i,j)
          enddo
       enddo

       !===============================================================

    case (jMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(i,0,j,l) = w(i,1,j,l)
             enddo
             if( eddyModel ) rev(i,0,j) = rev(i,1,j)
          enddo
       enddo

       !===============================================================

    case (jMax)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(i,jb,j,l) = w(i,je,j,l)
             enddo
             if( eddyModel ) rev(i,jb,j) = rev(i,je,j)
          enddo
       enddo

       !===============================================================

    case (kMin)
       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(i,j,0,l) = w(i,j,1,l)
             enddo
             if( eddyModel ) rev(i,j,0) = rev(i,j,1)
          enddo
       enddo

       !===============================================================

    case (kMax)

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
          do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
                w(i,j,kb,l) = w(i,j,ke,l)
             enddo
             if( eddyModel ) rev(i,j,kb) = rev(i,j,ke)
          enddo
       enddo

    end select

  end subroutine turb2ndHalo

  subroutine turbBCNSWall(secondHalo)
    !
    !       turbBCNSWall applies the viscous wall boundary conditions      
    !       of the turbulent transport equations to a block. It is assumed 
    !       that the pointers in blockPointers are already set to the      
    !       correct block on the correct grid level.                       
    !
    use constants
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Subroutine argument.
    !
    logical, intent(in) :: secondHalo
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, i, j, l, m


    ! Loop over the viscous subfaces of this block.

    bocos: do nn=1,nViscBocos

       ! Set the corresponding arrays.

       call BCTurbWall(nn)

       ! Loop over the faces and set the state in
       ! the turbulent halo cells.

       select case (BCFaceID(nn))
       case (iMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(1,i,j,l) = bvti1(i,j,l)
                   do m=nt1,nt2
                      w(1,i,j,l) = w(1,i,j,l) - bmti1(i,j,l,m)*w(2,i,j,m)
                   enddo
                   if (secondHalo) w(0,i,j,l) = w(1,i,j,l)
                end do

                if (eddyModel) then 
                   rev(1,i,j) = -rev(2,i,j)
                   if (secondHalo) then 
                      rev(0,i,j) = rev(1,i,j)
                   end if
                end if
             end do
          end do
       case (iMax) 
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(ie,i,j,l) = bvti2(i,j,l)
                   do m=nt1,nt2
                      w(ie,i,j,l) = w(ie,i,j,l) - bmti2(i,j,l,m)*w(il,i,j,m)
                   enddo
                   if (secondHalo) w(ib,i,j,l) = w(ie,i,j,l)
                end do

                if (eddyModel) then 
                   rev(ie,i,j) = -rev(il,i,j)
                   if (secondHalo) then 
                      rev(ib,i,j) = rev(ie,i,j)
                   end if
                end if
             end do
          end do
       case (jMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(i,1,j,l) = bvtj1(i,j,l)
                   do m=nt1,nt2
                      w(i,1,j,l) = w(i,1,j,l) - bmtj1(i,j,l,m)*w(i,2,j,m)
                   enddo
                   if (secondHalo) w(i,0,j,l) = w(i,1,j,l)
                end do

                if (eddyModel) then 
                   rev(i,1,j) = -rev(i,2,j)
                   if (secondHalo) then 
                      rev(i,0,j) = rev(i,1,j)
                   end if
                end if
             end do
          end do
       case (jMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(i,je,j,l) = bvtj2(i,j,l)
                   do m=nt1,nt2
                      w(i,je,j,l) = w(i,je,j,l) - bmtj2(i,j,l,m)*w(i,jl,j,m)
                   enddo
                   if (secondHalo) w(i,jb,j,l) = w(i,je,j,l)
                end do

                if (eddyModel) then 
                   rev(i,je,j) = -rev(i,jl,j)
                   if (secondHalo) then 
                      rev(i,jb,j) = rev(i,je,j)
                   end if
                end if
             end do
          end do
       case (kMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(i,j,1,l) = bvtk1(i,j,l)
                   do m=nt1,nt2
                      w(i,j,1,l) = w(i,j,1,l) - bmtk1(i,j,l,m)*w(i,j,2,m)
                   enddo
                   if (secondHalo) w(i,j,0,l) = w(i,j,1,l)
                end do

                if (eddyModel) then 
                   rev(i,j,1) = -rev(i,j,2)
                   if (secondHalo) then 
                      rev(i,j,0) = rev(i,j,1)
                   end if
                end if
             end do
          end do
       case (kMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                do l=nt1,nt2
                   w(i,j,ke,l) = bvtk2(i,j,l)
                   do m=nt1,nt2
                      w(i,j,ke,l) = w(i,j,ke,l) - bmtk2(i,j,l,m)*w(i,j,kl,m)
                   enddo
                   if (secondHalo) w(i,j,kb,l) = w(i,j,ke,l)
                end do

                if (eddyModel) then 
                   rev(i,j,ke) = -rev(i,j,kl)
                   if (secondHalo) then 
                      rev(i,j,kb) = rev(i,j,ke)
                   end if
                end if
             end do
          end do
       end select
    enddo bocos
  end subroutine turbBCNSWall

end module turbBCRoutines
