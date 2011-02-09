!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbFarfield.f90                              *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-15-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbFarfield(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbFarfield applies the implicit treatment of the           *
!      * farfield boundary condition to subface nn. As the farfield     *
!      * boundary condition is independent of the turbulence model,     *
!      * this routine is valid for all models. It is assumed that the   *
!      * pointers in blockPointers are already set to the correct       *
!      * block on the correct grid level.                               *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
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

       real(kind=realType), dimension(:,:,:,:), pointer :: bmt
       real(kind=realType), dimension(:,:,:),   pointer :: bvt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointers for bmt and bvt, depending on the block face
       ! on which the subface is located.

       select case (BCFaceID(nn))
         case (iMin)
           bmt => bmti1; bvt => bvti1
         case (iMax)
           bmt => bmti2; bvt => bvti2
         case (jMin)
           bmt => bmtj1; bvt => bvtj1
         case (jMax)
           bmt => bmtj2; bvt => bvtj2
         case (kMin)
           bmt => bmtk1; bvt => bvtk1
         case (kMax)
           bmt => bmtk2; bvt => bvtk2
       end select

       ! Loop over the faces of the subfaces and set the values of
       ! bmt and bvt for an implicit treatment.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd

           ! Store the three components of the unit normal a bit easier.

           nnx = BCData(nn)%norm(i,j,1)
           nny = BCData(nn)%norm(i,j,2)
           nnz = BCData(nn)%norm(i,j,3)

           ! Determine the dot product between the outward pointing
           ! normal and the free stream velocity direction and add the
           ! possible grid velocity.

           dot = nnx*wInf(ivx) + nny*wInf(ivy) + nnz*wInf(ivz) &
               - BCData(nn)%rface(i,j)

           ! Determine whether we are dealing with an inflow or
           ! outflow boundary here.

           if(dot > zero) then

             ! Outflow. Simply extrapolation or zero Neumann BC
             ! of the turbulent variables.

             do l=nt1,nt2
               bmt(i,j,l,l) = -one
             enddo

           else

             ! Inflow. Turbulent variables are prescribed.

             do l=nt1,nt2
               bvt(i,j,l) = wInf(l)
             enddo

           endif

         enddo
       enddo

       end subroutine bcTurbFarfield
