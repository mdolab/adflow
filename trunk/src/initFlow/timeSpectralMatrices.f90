!
!      ******************************************************************
!      *                                                                *
!      * File:          timeSpectralMatrices.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-06-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine timeSpectralMatrices
!
!      ******************************************************************
!      *                                                                *
!      * timeSpectralMatrices computes the matrices for the time        *
!      * derivative of the time spectral method for all sections. For   *
!      * scalar quantities these matrices only differ if sections have  *
!      * different periodic times. For vector quantities, such as       *
!      * momentum, these matrices can be different depending on whether *
!      * the section is rotating or not and the number of slices        *
!      * present.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use section
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, ll, kk, ii
       integer(kind=intType) :: i, j

       real(kind=realType), dimension(3,3) :: tmpMat

       real(kind=realType), dimension(:,:), allocatable :: coefSpectral
       real(kind=realType), dimension(:,:,:,:), allocatable :: &
                                               matrixCoefSpectral
       real(kind=realType), dimension(:,:,:), allocatable :: &
                                               diagMatCoefSpectral
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! This routine is only used for the spectral solutions. Return
       ! immediately if a different mode is solved.

       if(equationMode /= timeSpectral) return

       ! Allocate the memory for the matrices as well as the help
       ! variables needed to construct these matrices.

       nn = nTimeIntervalsSpectral
       mm = 3*nn
       kk = nn - 1
       allocate(dscalar(nSections,nn,nn),             &
                dvector(nSections,mm,mm),             &
                coefSpectral(nSections,kk),           &
                matrixCoefSpectral(nSections,kk,3,3), &
                diagMatCoefSpectral(nSections,3,3), stat=ierr)
       if(ierr /= 0)                              &
         call terminate("timeSpectralMatrices", &
                        "Memory allocation failure for the matrices of &
                        &the spectral time derivatives.")

       ! Determine the help variables needed to construct the
       ! actual matrices.

       call timeSpectralCoef(coefSpectral, matrixCoefSpectral, &
                             diagMatCoefSpectral)
!
!      ******************************************************************
!      *                                                                *
!      * Determine the time derivative matrices for the sections.       *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of sections.

       sectionLoop: do ii=1,nSections
!
!        ****************************************************************
!        *                                                              *
!        * Matrix for scalar quantities.                                *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the number of rows.

         do nn=1,nTimeIntervalsSpectral

           ! Set the diagonal element to zero, i.e. there is no
           ! contribution to the own time derivative.

           dscalar(ii,nn,nn) = zero

           ! Loop over the rest of the columns.

           do mm=1,(nTimeIntervalsSpectral - 1)

             ! Determine the corresponding column index.

             ll = nn + mm
             if(ll > nTimeIntervalsSpectral) &
               ll = ll - nTimeIntervalsSpectral

             ! Store the corresponding coefficient in dscalar.

             dscalar(ii,nn,ll) = coefSpectral(ii,mm)

           enddo
         enddo
!
!        ****************************************************************
!        *                                                              *
!        * Matrices for vector quantities.                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the number of time intervals; the number of rows
         ! is 3 times this number.

         rowLoop: do nn=1,nTimeIntervalsSpectral

           ! Initialize the diagonal block to diagMatCoefSpectral,
           ! the additional diagonal entry needed for the rotational
           ! periodicity.

           kk = 3*(nn-1)
           do j=1,3
             do i=1,3
               dvector(ii,kk+i,kk+j) = diagMatCoefSpectral(ii,i,j)
             enddo
           enddo

           ! Loop over the other time intervals, which contribute to
           ! the time derivative.

           columnLoop: do mm=1,(nTimeIntervalsSpectral - 1)

             ! Determine the corresponding column index and check the
             ! situation we are having here.

             ll = nn + mm
             if(ll > nTimeIntervalsSpectral) then

               ! Index is outside the range and a shift must be applied.

               ll = ll - nTimeIntervalsSpectral

               ! The vector must be rotated. This effect is incorporated
               ! directly in the matrix of time derivatives.

               do j=1,3
                 do i=1,3
                   tmpMat(i,j) = matrixCoefSpectral(ii,mm,i,1) &
                               * rotMatrixSpectral(ii,1,j)     &
                               + matrixCoefSpectral(ii,mm,i,2) &
                               * rotMatrixSpectral(ii,2,j)     &
                               + matrixCoefSpectral(ii,mm,i,3) &
                               * rotMatrixSpectral(ii,3,j)
                 enddo
               enddo

             else

               ! Index is in the range. Copy the matrix coefficient
               ! into tmpMat.

               do j=1,3
                 do i=1,3
                   tmpMat(i,j) = matrixCoefSpectral(ii,mm,i,j)
                 enddo
               enddo

             endif

             ! Determine the offset for the column index and store
             ! this submatrix in the correct place of dvector.

             ll = 3*(ll-1)
             do j=1,3
               do i=1,3
                 dvector(ii,kk+i,ll+j) = tmpMat(i,j)
               enddo
             enddo

           enddo columnLoop
         enddo rowLoop
       enddo sectionLoop

       ! Release the memory of the help variables needed to construct
       ! the matrices of the time derivatives.

       deallocate(coefSpectral, matrixCoefSpectral, &
                  diagMatCoefSpectral, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("timeSpectralMatrices", &
                        "Deallocation failure for the help variables.")

       end subroutine timeSpectralMatrices
