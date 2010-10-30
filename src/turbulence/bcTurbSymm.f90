!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbSymm.F90                                  *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbSymm(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbSymm applies the implicit treatment of the symmetry      *
!      * boundary condition (or inviscid wall) to subface nn. As the    *
!      * symmetry boundary condition is independent of the turbulence   *
!      * model, this routine is valid for all models. It is assumed     *
!      * that the pointers in blockPointers are already set to the      *
!      * correct block on the correct grid level.                       *
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

       real(kind=realType), dimension(:,:,:,:), pointer :: bmt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointer for bmt, depending on the block face on which
       ! the subface is located.

       select case (BCFaceID(nn))
         case (iMin)
           bmt => bmti1
         case (iMax)
           bmt => bmti2
         case (jMin)
           bmt => bmtj1
         case (jMax)
           bmt => bmtj2
         case (kMin)
           bmt => bmtk1
         case (kMax)
           bmt => bmtk2
       end select

       ! Loop over the faces of the subfaces and set the values of bmt
       ! for an implicit treatment. For a symmetry face this means
       ! that the halo value is set to the internal value.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           do l=nt1,nt2
             bmt(i,j,l,l) = -one
           enddo
         enddo
       enddo

       end subroutine bcTurbSymm
