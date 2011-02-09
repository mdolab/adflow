!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingPlaneInterpol.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-28-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mixingPlaneInterpol(level, color, donorSlideID, &
                                      haloSlideID, nSlices, commPattern)
!
!      ******************************************************************
!      *                                                                *
!      * mixingPlaneInterpol determines the data for the                *
!      * communication/interpolation pattern of the given side of a     *
!      * mixing plane.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use commMixing
       use mixingData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, color
       integer(kind=intType), intent(in) :: donorSlideID
       integer(kind=intType), intent(in) :: haloSlideID
       integer(kind=intType), intent(in) :: nSlices

       type(commMixingType), intent(out) :: commPattern
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine distribution used as interpolation intervals.
       ! This is done using the grid information of the donor.

       call mixingIntervals(donorSlideID, level, color, nSlices)

       ! Set the number of interpolation intervals. As these intervals
       ! correspond to cell centered data, this is the number of
       ! interpolation points - 1.

       commPattern%nInter = nMixingPoints - 1

       ! Determine the donor interpolation data.

       call mixingDonorInterpol(level, donorSlideID, color, commPattern)

       ! Determine the halo interpolation data.

       call mixingHaloInterpol(level, haloSlideID, color, commPattern)

       ! Release the memory of mixingPoints and mixingCells, which are
       ! stored in the module mixingData.

       deallocate(mixingPoints, mixingCells, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("mixingPlaneInterpol", &
                        "Deallocation failure for mixingPoints and &
                        &mixingCells.")

       end subroutine mixingPlaneInterpol
