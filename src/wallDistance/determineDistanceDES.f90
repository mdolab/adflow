!
!      ******************************************************************
!      *                                                                *
!      * File:          determineDistance.F90                           *
!      * Author:        Edwin van der Weide, Eran Arad                  *
!      * Starting date: 03-03-2003                                      *
!      * Last modified: 08-10-2005                                      *
!      *                                                                *
!      * This routine is supposed to be called after determineDistance  *
!      ******************************************************************
!
       subroutine determineDistanceDES(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * DetermineDistance determines the distance from the center      *
!      * of the cell to the nearest viscous wall for owned cells.       *
!      * In addition, for des model ,calculate dTilda                   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputPhysics
       use section
       use viscSurface
       use inputDES

       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer(kind=intType) ::  nn, i, j, k, km1, jm1, im1

       real(kind=realType) ::  Lengthscale
!
!
       real(kind=realType), dimension(3) ::  tmp, distd

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       
       ! Loop over the domains.


       domains: do nn=1,nDom

         ! Set the pointers for this block and store the section id
         ! to which this block a bit easier in mm.

         call setPointers(nn, level, sps)

         ! Loop over the cell centers

         kLoop: do k=2,kl
            km1 = k-1

            jLoop: do j=2,jl
               jm1 = j-1

               iLoop: do i=2,il

                  im1 = i-1

! determine dx(k)

                  tmp(1) = x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                        + x(i-1,j,  k,  1) + x(i,j,  k,  1)  &
                        - x(i-1,j-1,km1,1) - x(i,j-1,km1,1)  &
                        - x(i-1,j,  km1,1) - x(i,j,  km1,1)

                  tmp(2) = x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                        + x(i-1,j,  k,  2) + x(i,j,  k,  2)  &
                        - x(i-1,j-1,km1,2) - x(i,j-1,km1,2)  &
                        - x(i-1,j,  km1,2) - x(i,j,  km1,2)

                  tmp(3) = x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                        + x(i-1,j,  k,  3) + x(i,j,  k,  3)  &
                        - x(i-1,j-1,km1,3) - x(i,j-1,km1,3)  &
                        - x(i-1,j,  km1,3) - x(i,j,  km1,3)

                 ! Compute the sum of the distance to the k-neighbor

                  distd(3)  = fourth*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
 
!
! determine dx(j)
!

                  tmp(1) = x(i-1,j  ,k-1,1) + x(i,j  ,k-1,1)  &
                        + x(i-1,j,  k,  1) + x(i,j,  k,  1)  &
                        - x(i-1,jm1,k-1,1) - x(i,jm1,k-1,1)  &
                        - x(i-1,jm1,k,  1) - x(i,jm1,k,  1)

                  tmp(2) = x(i-1,j  ,k-1,2) + x(i,j  ,k-1,2)  &
                        + x(i-1,j,  k,  2) + x(i,j,  k,  2)  &
                        - x(i-1,jm1,k-1,2) - x(i,jm1,k-1,2)  &
                        - x(i-1,jm1,k,  2) - x(i,jm1,k,  2)

                  tmp(3) = x(i-1,j  ,k-1,3) + x(i,j  ,k-1,3)  &
                        + x(i-1,j,  k,  3) + x(i,j,  k,  3)  &
                        - x(i-1,jm1,k-1,3) - x(i,jm1,k-1,3)  &
                        - x(i-1,jm1,k,  3) - x(i,jm1,k,  3)

                 ! Compute the sum of the distance to the j-neighbor

                  distd(2)  = fourth*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
!
! determine dx(i)
!

                  tmp(1) = x(i  ,j-1,k-1,1) + x(i  ,j,k-1,1)  &
                        + x(i  ,j-1,k,  1) + x(i  ,j,k,  1)  &
                        - x(im1,j-1,k-1,1) - x(im1,j,k-1,1)  &
                        - x(im1,j-1,k,  1) - x(im1,j,k,  1)

                  tmp(2) = x(i  ,j-1,k-1,2) + x(i  ,j,k-1,2)  &
                        + x(i  ,j-1,k,  2) + x(i  ,j,k,  2)  &
                        - x(im1,j-1,k-1,2) - x(im1,j,k-1,2)  &
                        - x(im1,j-1,k,  2) - x(im1,j,k,  2)

                  tmp(3) = x(i  ,j-1,k-1,3) + x(i  ,j,k-1,3)  &
                        + x(i  ,j-1,k,  3) + x(i  ,j,k,  3)  &
                        - x(im1,j-1,k-1,3) - x(im1,j,k-1,3)  &
                        - x(im1,j-1,k,  3) - x(im1,j,k,  3)

                 ! Compute the sum of the distance to the i-neighbor

                  distd(1)  = fourth*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
              
!
! find the max(delta xI)
!
              !    filterDES(i,j,k)  = cDes*max(distd(1),distd(2),distd(3))
                  filterDES(i,j,k)  = cDes*(distd(1)*distd(2)*distd(3))**0.3333


               enddo iLoop

            enddo jLoop

         enddo kLoop

         if(turbModel == spalartAllmaras .or. turbModel ==  spalartAllmarasEdwards)then

            kLoopsa: do k=2,kl

               jLoopsa: do j=2,jl

                  iLoopsa: do i=2,il
!
!--- avoid cell size calculations for areas out of DES zone
!

                     if (x(i,j,k,1) <=  xDESmin .or.  x(i,j,k,1) >=  xDESmax  )then
                        filterDES(i,j,k) = d2wall(i,j,k)
                     end if

!
! des: Calculate dTilde
!
                     if(d2wall(i,j,k) <= distDESmin  .or. d2wall(i,j,k) >= distDESmax )then

                        filterDES(i,j,k) = d2wall(i,j,k)

                     else
                      
                        if (d2Wall(i,j,k) <= distRANSmax)then
                           if(.not.applyDDES)&
                                filterDES(i,j,k) = min(filterDES(i,j,k),d2wall(i,j,k))
! The selection is done here for DES97; for DDES the selection is done
! in the solution loop (saSolve, and SSTSolve)
                        end if
                      
                     end if

                  enddo iLoopsa

               enddo jLoopsa

            enddo kLoopsa

         end if ! turbModel

      enddo domains

       write(*,*)'determineDistanceDes was called'

       end subroutine determineDistanceDes 
