!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbInflow.f90                                *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbInflow(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbInflow applies the implicit treatment of the inflow      *
!      * boundary conditions to subface nn. As the inflow boundary      *
!      * condition is independent of the turbulence model, this routine *
!      * is valid for all models. It is assumed that the pointers in    *
!      * blockPointers are already set to the correct block on the      *
!      * correct grid level.                                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
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

       real(kind=realType), dimension(:,:,:),   pointer :: bvt
       real(kind=realType), dimension(:,:,:,:), pointer :: bmt
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
       ! bvt and bmt such that the inflow state is linearly extrapolated
       ! with a fixed state at the face.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd

           ! Loop over the number of turbulent variables.

           do l=nt1,nt2
             bvt(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
             bmt(i,j,l,l) = one
           enddo

         enddo
       enddo

       end subroutine bcTurbInflow
