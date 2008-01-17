!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllTurbBC.f90                              *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      * Starting date: 05-01-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine applyAllTurbBC(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllTurbBC applies all boundary conditions to the          *
!      * turbulent transport equations for the all blocks on the grid   *
!      * level currentLevel.                                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo
!
!      Local variables.
!
       integer(kind=intType) :: nn, sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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

!      ==================================================================

       subroutine applyAllTurbBCThisBlock(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllTurbBCThisBlock sets the halo values of the            *
!      * turbulent variables and eddy viscosity for the block the       *
!      * variables in blockPointers currently point to.                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

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

         ! Loop over the faces and set the state in
         ! the turbulent halo cells.

         if( wallFunctions ) then

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

         else

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               do l=nt1,nt2
                 ww1(i,j,l) = bvt(i,j,l)
                 do m=nt1,nt2
                   ww1(i,j,l) = ww1(i,j,l) - bmt(i,j,l,m)*ww2(i,j,m)
                 enddo
               enddo
             enddo
           enddo

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
