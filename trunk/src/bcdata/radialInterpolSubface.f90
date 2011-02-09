!
!      ******************************************************************
!      *                                                                *
!      * File:          radialInterpolSubface.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-02-2004                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, nbcVar, &
                                        cgnsBoco, blockFaceId, indCoor, &
                                        ind, bcVarPresent, bcVarArray,  &
                                        axAssumed)
!
!      ******************************************************************
!      *                                                                *
!      * radialInterpolSubface interpolates the prescribed variables    *
!      * in the data set of the given cgns subface onto the subface     *
!      * indicated by iBeg, iEnd, jBeg, jEnd and blockFaceId.           *
!      * This routine performs a 1d interpolation in radial direction   *
!      * assuming that there is no variation in the other directions.   *
!      * The variables in blockPointers are already set.                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use section
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType), intent(in) :: nbcVar, cgnsBoco
       integer(kind=intType), intent(in) :: blockFaceId

       integer(kind=intType), dimension(2), intent(in) :: indCoor
       integer(kind=intType), dimension(2,nbcVar), intent(in) :: ind

       real(kind=realType), dimension(iBeg:iEnd,jBeg:jEnd,nbcVar), &
                                            intent(out) :: bcVarArray

       logical, intent(inout) :: axAssumed
       logical, dimension(nbcVar), intent(in) :: bcVarPresent
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l, n1, n2, nn, ii
       integer(kind=intType) :: nDim, nPoints, start

       real(kind=realType) :: t, length, fact, rad, ww1, ww2

       real(kind=realType), dimension(3) :: rotAxis, xc, xaxis

       real(kind=realType), dimension(:),     pointer :: rr
       real(kind=realType), dimension(:,:,:), pointer :: xf

       character(len=maxStringLen) :: errorMessage

       type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the rotation axis of the section a bit easier.

       rotAxis = sections(sectionId)%rotAxis

       ! Check if a rotation axis could be constructed for the section
       ! to which this block belongs. If not it is assumed that the
       ! x-axis is the axial direction.

       length = rotAxis(1)**2 + rotAxis(2)**2 + rotAxis(3)**2

       if(length < half) then

         ! No axis could be constructed from the rotation info.
         ! Assume this is the x-axis and set axAssumed to .True.

         rotAxis(1) = one; rotAxis(2) = zero; rotAxis(3) = zero
         axAssumed = .true.

       endif

       ! Set the pointer for dataSet to make the code more readable.

       dataSet => cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%dataSet

       ! Check the number of dimensions of the specified data set.
       ! This should be 1, because a 1d interpolation in radial
       ! direction is performed.

       k    = indCoor(1)
       l    = indCoor(2)
       nDim = dataSet(k)%dirichletArrays(l)%nDimensions

       if(nDim > 1) then
         write(errorMessage,101) &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
         call terminate("radialInterpolSubface", errorMessage)
 101     format("Zone",1X,A,", subface",1X,A,": Multidimensional &
                &radially varying data specified. Only 1d data possible")
       endif

       ! Set the pointer for the radial coordinate and abbreviate the
       ! number of interpolation points a bit easier.

       rr     => dataSet(k)%dirichletArrays(l)%dataArr
       nPoints = dataSet(k)%dirichletArrays(l)%dataDim(1)

       ! Check if the data is specified for increasing radius.

       do i=2,nPoints
         if(rr(i) < rr(i-1)) exit
       enddo

       if(i <= nPoints) then
         write(errorMessage,102) &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
         call terminate("radialInterpolSubface", errorMessage)
 102     format("Zone",1X,A,", subface",1X,A,": Data should be &
                &specified for increasing radius.")
       endif

       ! Set the pointer for the coordinates of the block face on which
       ! the boundary subface is located.

       select case (blockFaceId)
         case (iMin)
           xf => x(1,:,:,:)
         case (iMax)
           xf => x(il,:,:,:)
         case (jMin)
           xf => x(:,1,:,:)
         case (jMax)
           xf => x(:,jl,:,:)
         case (kMin)
           xf => x(:,:,1,:)
         case (kMax)
           xf => x(:,:,kl,:)
       end select

       ! Compute the factor needed to compute the coordinates in the
       ! original units. The fourth comes from the averaging of the 4
       ! nodal coordinates.

       fact = fourth/cgnsDoms(nbkGlobal)%LRef

       ! Loop over the range of the subface.

       jloop: do j=jBeg,jEnd
         iloop: do i=iBeg,iEnd

           ! Determine the coordinates of the face center. Normally this
           ! is an average of i-1, i, j-1, j, but due to the usage of
           ! the pointer xf and the fact that x originally starts at 0,
           ! an offset of 1 is introduced and thus the average should
           ! be taken of i, i+1, j and j+1.

           xc(1) = fact*(xf(i,j,  1) + xf(i+1,j,  1) &
                 +       xf(i,j+1,1) + xf(i+1,j+1,1))
           xc(2) = fact*(xf(i,j,  2) + xf(i+1,j,  2) &
                 +       xf(i,j+1,2) + xf(i+1,j+1,2))
           xc(3) = fact*(xf(i,j,  3) + xf(i+1,j,  3) &
                 +       xf(i,j+1,3) + xf(i+1,j+1,3))

           ! Determine the parameter, which defines the closest point
           ! on the rotation axis. Note that rotAxis is a unit-vector.

           t = rotAxis(1)*(xc(1) - sections(sectionId)%rotCenter(1)) &
             + rotAxis(2)*(xc(2) - sections(sectionId)%rotCenter(2)) &
             + rotAxis(3)*(xc(3) - sections(sectionId)%rotCenter(3))

           ! Determine the coordinates of this point.

           xaxis(1) = sections(sectionId)%rotCenter(1) + t*rotAxis(1)
           xaxis(2) = sections(sectionId)%rotCenter(2) + t*rotAxis(2)
           xaxis(3) = sections(sectionId)%rotCenter(3) + t*rotAxis(3)

           ! Determine the radius of this point.

           rad = sqrt((xc(1)-xaxis(1))**2 + (xc(2)-xaxis(2))**2 &
               +      (xc(3)-xaxis(3))**2)

           ! Determine the interpolation interval and set the
           ! interpolation weights. Take care of the exceptions.

           checkInterpol: if(rad <= rr(1)) then

             ! Radius is less than the minimum value specified.
             ! Use constant extrapolation.

             n1  = 1;   n2  = 1
             ww1 = one; ww2 = zero

           else if(rad >= rr(nPoints)) then checkInterpol

             ! Radius is larger than the maximum value specified.
             ! Use constant extrapolation.

             n1  = nPoints; n2  = nPoints
             ww1 = one;     ww2 = zero

           else checkInterpol

             ! Radius is in the range. Determine the correct interval
             ! using a binary search algorithm.

             ii    = nPoints - 1
             start = 1
             interval: do

               ! Next guess for the interval and determine the new
               ! situation.

               nn = start + ii/2
               if(rad > rr(nn+1)) then

                 ! Rad is larger than the upper boundary of the
                 ! current interval. Update the lower boundary.

                 start = nn + 1; ii = ii - 1

               else if(rad >= rr(nn)) then

                 ! This is the correct range. Exit the loop.

                 exit

               endif

               ! Modify ii for the next branch to search.

               ii = ii/2

             enddo interval

             ! Rad is in the interval nn, nn+1. Store this and
             ! determine the interpolation weight.

             n1 = nn
             n2 = nn + 1
             ww1 = (rr(nn+1) - rad)/(rr(nn+1) - rr(nn))
             ww2 = one - ww1

           endif checkInterpol
 
           ! Interpolate the values the values present for this face.

           do nn=1,nbcVar
             if( bcVarPresent(nn) ) then

               ! Easier storage of the indices in the data set.

               k = ind(1,nn)
               l = ind(2,nn)

               ! Interpolate this variable.

               bcVarArray(i,j,nn) =                                 &
                     ww1*dataSet(k)%dirichletArrays(l)%dataArr(n1) &
                   + ww2*dataSet(k)%dirichletArrays(l)%dataArr(n2)
             endif
           enddo

         enddo iloop
       enddo jloop

       end subroutine radialInterpolSubface
