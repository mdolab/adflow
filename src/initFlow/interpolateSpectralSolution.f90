!
!      ******************************************************************
!      *                                                                *
!      * File:          interpolateSpectralSolution.f90                 *
!      * Author:        Edwin van der Weide, Arathi K. Gopinath.        *
!      * Starting date: 08-08-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine interpolateSpectralSolution
!
!      ******************************************************************
!      *                                                                *
!      * interpolateSpectralSolution uses a spectral interpolation to   *
!      * determine the initialization of the flow solution.             *
!      * The solution is interpolated from the solution read, which     *
!      * contains a different number of time instances and is stored in *
!      * IOVar()%w. This variable can be found in IOModule.             *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputTimeSpectral
       use IOModule
       use section
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: jj, nn, ll, sps, i, j, k, l

       real(kind=realType) :: t

       real(kind=realType), dimension(nSolsRead) :: alpScal
       real(kind=realType), dimension(nSections,nSolsRead,3,3) :: alpMat
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions to be interpolated.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Determine the ratio of the time of this solution and the
         ! periodic time.

         nn = sps - 1
         t  = real(nn,realType)/real(nTimeIntervalsSpectral,realType)

         ! Determine the interpolation coefficients for both the scalar
         ! and the vector quantities.

         call spectralInterpolCoef(nSolsRead, t, alpScal, alpMat)

         ! Loop over the local number of blocks.

         domains: do nn=1,nDom

           ! Set the pointers for this block to the finest grid level.

           call setPointers(nn, 1_intType, sps)

           ! Loop over the number of variables to be interpolated.

           varLoop: do l=1,nw

             ! Check if this is a velocity variable.

             velTest: if(l == ivx .or. l == ivy .or. l == ivz) then

               ! Velocity variable. Set ll, which is the row in the
               ! matrix coefficients.

               select case (l)
                 case (ivx)
                   ll = 1
                 case (ivy)
                   ll = 2
                 case (ivz)
                   ll = 3
               end select

               ! Loop over the owned cells to interpolate the variable.

               do k=2,kl
                 do j=2,jl
                   do i=2,il

                     ! Initialization to zero and loop over the number of
                     ! spectral solutions used in the interpolation and
                     ! update the variable accordingly. Note that for the
                     ! vector variables the matrix coefficients must be
                     ! used; these matrices can be different for the
                     ! sections present in the grid.

                     w(i,j,k,l) = zero
                     do jj=1,nSolsRead
                       w(i,j,k,l) = w(i,j,k,l)                &
                                  + alpMat(sectionID,jj,ll,1) &
                                  * IOVar(nn,jj)%w(i,j,k,ivx) &
                                  + alpMat(sectionID,jj,ll,2) &
                                  * IOVar(nn,jj)%w(i,j,k,ivy) &
                                  + alpMat(sectionID,jj,ll,3) &
                                  * IOVar(nn,jj)%w(i,j,k,ivz)
                     enddo

                   enddo
                 enddo
               enddo

             else velTest

               ! Scalar variable.
               ! Loop over the owned cells to interpolate the variable.

               do k=2,kl
                 do j=2,jl
                   do i=2,il

                     ! Initialization to zero and loop over the number of
                     ! spectral solutions used in the interpolation and
                     ! update the variable accordingly.

                     w(i,j,k,l) = zero
                     do jj=1,nSolsRead
                       w(i,j,k,l) = w(i,j,k,l)  &
                                  + alpScal(jj)*IOVar(nn,jj)%w(i,j,k,l)
                     enddo

                   enddo
                 enddo
               enddo

             endif velTest

           enddo varLoop

         enddo domains

       enddo spectralLoop

       ! Release the memory of w of IOVar.

       do sps=1,nSolsRead
         do nn=1,nDom
           deallocate(IOVar(nn,sps)%w, stat=ierr)
           if(ierr /= 0) &
             call terminate("interpolateSpectralSolution", &
                            "Deallocation failure for w.")
         enddo
       enddo

       end subroutine interpolateSpectralSolution
