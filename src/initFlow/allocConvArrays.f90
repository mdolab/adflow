!
!      ******************************************************************
!      *                                                                *
!      * File:          allocConvArrays.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-23-2003                                      *
!      * Last modified: 07-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine allocConvArrays(nIterTot)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * allocConvArrays allocates the memory for the convergence       *
  !      * arrays. The number of iterations allocated, nIterTot, is       *
  !      * enough to store the maximum number of iterations specified     *
  !      * plus possible earlier iterations read from the restart file.   *
  !      * This routine MAY be called with data already inside of         *
  !      * convArray and this will be saved.                              *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use inputIO
  use inputTimeSpectral
  use monitor
  implicit none
  !
  !      Subroutine argument.
  !
  integer(kind=intType) :: nIterTot,ll,mm,nn,i,j,k
  !
  !      Local variables.
  !
  integer :: ierr
  real(kind=realType), dimension(:,:,:),allocatable ::  tempArray
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Return immediately if the convergence history (of the inner
  ! iterations) does not need to be stored. This logical can
  ! only be .false. for an unsteady computation.
  if(.not. storeConvInnerIter) return 

  if (allocated(convArray)) then
     ! Its already allocated, so copy the data out first:

     ll = ubound(convArray,1)
     mm = ubound(convArray,2)
     nn = ubound(convArray,3)
     allocate(tempArray(0:ll,1:mm,1:nn))

     do k=1,nn
        do j=1,mm
           do i=0,ll
              tempArray(i,j,k) = convArray(i,j,k)
           end do
        end do
     end do

     deallocate(convArray)

     allocate(convArray(0:nIterTot,nTimeIntervalsSpectral,nMon), stat=ierr)
     if(ierr /= 0) then
        call terminate("allocConvArrays", &
             "Memory allocation failure for convArray")
     end if
     
     ! Copy Data back in:
     
     do k=1,nn
        do j=1,mm
           do i=1,ll
              convArray(i,j,k) = tempArray(i,j,k)
           end do
        end do
     end do

     deallocate(tempArray)

  else ! Just allocate:

     allocate(convArray(0:nIterTot,nTimeIntervalsSpectral,nMon), stat=ierr)
     if(ierr /= 0) then
        call terminate("allocConvArrays", &
             "Memory allocation failure for convArray")
     end if

     ! Zero Array:
     convArray = zero
  end if

end subroutine allocConvArrays

