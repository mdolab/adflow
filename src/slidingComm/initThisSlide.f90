!
!      ******************************************************************
!      *                                                                *
!      * File:          initThisSlide.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-13-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initThisSlide(color)
!
!      ******************************************************************
!      *                                                                *
!      * initThisSlide determines the communicator, the rotation axis   *
!      * and the radial scale for this sliding mesh interface.          *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: color
!
!      Local parameters. The value of tolDotmin corresponds
!      to 0.1 degrees.
!
       real(kind=realType), parameter :: notPresent = -1.e+6_realType
       real(kind=realType), parameter :: tolDotmin  = 0.9999985_realType
!
!      Local variables.
!
       integer :: ierr, groupColor
       integer :: nProcSlide, myIDSlide, commSlide

       integer(kind=intType) :: nn, mm, ii, i, j, k

       real(kind=realType) :: dot, length, xx, yy, zz
       real(kind=realType) :: rMin, rMax, axMin, axMax

       real(kind=realType), dimension(3) :: tmpCenter, center
       real(kind=realType), dimension(3) :: tmp, tmpAxis
       real(kind=realType), dimension(3) :: axis, radVec1, radVec2

       logical :: firstTime

       character(len=7)            :: integerString
       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the value of groupColor to the current interface ID if
       ! this processor contributes and to mpi_undefined otherwise.

       if( myInterfaces(color)%procContributes ) then
         groupColor = myInterfaces(color)%globalSlideID
       else
         groupColor = mpi_undefined
       endif

       ! Create the communicator. Also store it in myInterfaces.

       call mpi_comm_split(SUmb_comm_world, groupColor, myID, &
                          commSlide, ierr)
       myInterfaces(color)%commSlide = commSlide

       ! Check if this processor should participate in a search.
       ! If not, return.

       if(.not. myInterfaces(color)%procContributes) return

       ! Determine the number of processors and my processor number
       ! in this new group. Store it in myInterfaces too.

       call mpi_comm_rank(commSlide, myIDSlide,  ierr)
       call mpi_comm_size(commSlide, nProcSlide, ierr)

       myInterfaces(color)%myIDSlide  = myIDSlide
       myInterfaces(color)%nProcSlide = nProcSlide

       ! Determine for the locally stored blocks the rotation axis
       ! and the rotation center.

       firstTime = .true.
       domains1: do nn=1,nDom

         ! Set the pointers for this block; note that only the 1st
         ! time spectral interval needs to be considered.

         call setPointers(nn,1_intType,1_intType)

         ! Check if this block is part of the current sliding mesh
         ! interface.

         bocos1: do mm=1,nBocos
           if((BCType(mm) == slidingInterface) .and. &
              (abs(groupNum(mm)) ==                  &
               myInterfaces(color)%globalSlideID)) then

             ! Block is part of the sliding mesh. Check if a rotating
             ! frame was specified for this block.

             ii = nbkGlobal
             if( cgnsDoms(ii)%rotatingFrameSpecified ) then

               ! A rotating frame was specified for this block.
               ! Determine the situation we are having here.

               if( firstTime ) then

                 ! First time that a rotating block is encountered
                 ! on this processor. Copy the data specified.

                 tmpAxis   = cgnsDoms(ii)%rotRate
                 tmpCenter = cgnsDoms(ii)%rotCenter

                 ! Create a unit vector out of tmpAxis.

                 length = sqrt(tmpAxis(1)**2 + tmpAxis(2)**2 &
                        +      tmpAxis(3)**2)
                 dot = one/max(eps,length)
                 tmpAxis(1) = tmpAxis(1)*dot
                 tmpAxis(2) = tmpAxis(2)*dot
                 tmpAxis(3) = tmpAxis(3)*dot

                 ! Set to firstTime to .false. to indicate that a
                 ! rotating frame has been specified.

                 firstTime = .false.

               else

                 ! A rotating frame has been specified before. Check
                 ! if the data is correct. If not exit the code.

                 length = sqrt(cgnsDoms(ii)%rotRate(1)**2 &
                        +      cgnsDoms(ii)%rotRate(2)**2 &
                        +      cgnsDoms(ii)%rotRate(3)**2)
                 dot = (tmpAxis(1)*cgnsDoms(ii)%rotRate(1)  &
                     +  tmpAxis(2)*cgnsDoms(ii)%rotRate(2)  &
                     +  tmpAxis(3)*cgnsDoms(ii)%rotRate(3)) &
                     /  max(eps,length)

                 if(abs(dot) < tolDotmin) then

                   ! The rotation axis are not aligned. Print
                   ! an error message.

                   write(integerString,"(i6)") &
                         myInterfaces(color)%globalSlideID
                   integerString = adjustl(integerString)
                   write(errorMessage,100) trim(integerString)

                   call terminate("initThisSlide", errorMessage)

                 endif

                 ! Also check the rotation center.

                 if(tmpCenter(1) /= cgnsDoms(ii)%rotCenter(1) .or. &
                    tmpCenter(2) /= cgnsDoms(ii)%rotCenter(2) .or. &
                    tmpCenter(3) /= cgnsDoms(ii)%rotCenter(3)) then

                    write(integerString,"(i6)") &
                          myInterfaces(color)%globalSlideID
                    integerString = adjustl(integerString)
                    write(errorMessage,101) trim(integerString)

                    call terminate("initThisSlide", errorMessage)

                 endif
               endif
             endif
           endif
         enddo bocos1

       enddo domains1

       ! Set the value of axis and center, depending on the situation.

       if( firstTime) then
         axis   = notPresent
         center = notPresent
       else
         axis   = tmpAxis
         center = tmpCenter

         ! Make sure that the largest component of vector is positive.
         ! In this way the axial direction is uniquely defined.
         ! Use dot and length as temporary storage.

         dot = axis(1); length = abs(dot)
         if(abs(axis(2)) > length) then
           dot = axis(2); length = abs(dot)
         endif
         if(abs(axis(3)) > length) then
           dot = axis(3); length = abs(dot)
         endif

         if(dot < zero) then
           axis(1) = -axis(1)
           axis(2) = -axis(2)
           axis(3) = -axis(3)
         endif

       endif

       ! Determine the maximum value of axis of all processors and
       ! check if everything was done correctly.

       call mpi_allreduce(axis, tmp, 3, sumb_real, mpi_max, &
                          commSlide, ierr)

       if(tmp(1) == notPresent) then

         ! Rotation rate not specified for any of the blocks of this
         ! interface. Local processor 0 prints an error message,
         ! while the rest wait to be killed.

         if(myIDSlide == 0) then

           write(integerString,"(i6)") myInterfaces(color)%globalSlideID
           integerString = adjustl(integerString)
           write(errorMessage,102) trim(integerString)

           call terminate("initThisSlide", errorMessage)

         endif

         call mpi_barrier(commSlide, ierr)

       endif

       ! Check if the value is the same for all processors.
       ! This must only be checked if this processor had a locally
       ! stored rotation axis.

       if(.not. firstTime) then
         if((tmp(1) /= axis(1)) .or. (tmp(2) /= axis(2)) .or. &
            (tmp(3) /= axis(3))) then

           ! Local value of the rotation rate differs from the maximum for
           ! this sliding interface. Print an error message and exit.

           write(integerString,"(i6)") myInterfaces(color)%globalSlideID
           integerString = adjustl(integerString)
           write(errorMessage,100) trim(integerString)

           call terminate("initThisSlide", errorMessage)

         endif
       endif

       ! Store the rotation axis in axis.

       axis = tmp

       ! Idem for the rotation center.

       call mpi_allreduce(center, tmp, 3, sumb_real, mpi_max, &
                          commSlide, ierr)

       if(tmp(1) == notPresent) then

         ! Rotation center not specified for any of the blocks of this
         ! interface. Local processor 0 prints an error message,
         ! while the rest wait to be killed.

         if(myIDSlide == 0) then

           write(integerString,"(i6)") myInterfaces(color)%globalSlideID
           integerString = adjustl(integerString)
           write(errorMessage,103) trim(integerString)

           call terminate("initThisSlide", errorMessage)

         endif

         call mpi_barrier(commSlide, ierr)

       endif

       ! Check if the value is the same for all processors.
       ! This must only be checked if this processor had a locally
       ! stored rotation axis.

       if(.not. firstTime) then
         if((tmp(1) /= center(1)) .or. (tmp(2) /= center(2)) .or. &
            (tmp(3) /= center(3))) then

           ! Local value of the rotation center differs from the maximum
           ! for this sliding interface. Print an error message and exit.

           write(integerString,"(i6)") myInterfaces(color)%globalSlideID
           integerString = adjustl(integerString)
           write(errorMessage,101) trim(integerString)

           call terminate("initThisSlide", errorMessage)

         endif
       endif

       ! Store the rotation center in center and store both the center
       ! and the rotation axis for this sliding mesh interface.

       center = tmp
       myInterfaces(color)%rotAxis   = axis
       myInterfaces(color)%rotCenter = center

       ! For the transformation to cylindrical coordinates two unit
       ! vectors are needed which span the radial plane. These vectors
       ! must be normal to the axial vector as well as normal to each
       ! other. Note that this defines the vectors up to an arbitrary
       ! angle. First try the y-axis. If this is not good enough use
       ! the z-axis. One of them must be approximately right.

       if(abs(axis(2)) < 0.707107_realType) then
         radVec1(1) = zero
         radVec1(2) = one
         radVec1(3) = zero
       else
         radVec1(1) = zero
         radVec1(2) = zero
         radVec1(3) = one
       endif

       ! Make sure that radVec1 is normal to axis. Create a unit
       ! vector again.

       dot = radVec1(1)*axis(1) + radVec1(2)*axis(2) &
           + radVec1(3)*axis(3)
       radVec1(1) = radVec1(1) - dot*axis(1)
       radVec1(2) = radVec1(2) - dot*axis(2)
       radVec1(3) = radVec1(3) - dot*axis(3)

       dot = one/sqrt(radVec1(1)**2 + radVec1(2)**2 + radVec1(3)**2)
       radVec1(1) = radVec1(1)*dot
       radVec1(2) = radVec1(2)*dot
       radVec1(3) = radVec1(3)*dot

       ! Create the second vector which spans the radial plane. This must
       ! be normal to both axis and radVec1, i.e. the cross-product.

       radVec2(1) = axis(2)*radVec1(3) - axis(3)*radVec1(2)
       radVec2(2) = axis(3)*radVec1(1) - axis(1)*radVec1(3)
       radVec2(3) = axis(1)*radVec1(2) - axis(2)*radVec1(1)

       ! Store both vectors in the data structure for this sliding
       ! mesh interface.

       myInterfaces(color)%radVec1 = radVec1
       myInterfaces(color)%radVec2 = radVec2

       ! Determine the local values of the minimum and maximum axial
       ! and radial coordinates. First the initialization.

       rMin  =  large
       rMax  = -large
       axMin =  large
       axMax = -large

       ! Loop over the locally stored blocks and its possible subfaces
       ! on this sliding mesh interface.

       domains2: do nn=1,nDom

         ! Set the pointers for this block; note that only the 1st
         ! time spectral interval needs to be considered.

         call setPointers(nn,1_intType,1_intType)

         bocos2: do mm=1,nBocos

           ! Check if this boundary subface is part of the current
           ! sliding mesh interface.

           if((BCType(mm) == slidingInterface) .and. &
              (abs(groupNum(mm)) ==                  &
               myInterfaces(color)%globalSlideID)) then

             ! Loop over the nodes of this subface.

             do k=knBeg(mm), knEnd(mm)
               do j=jnBeg(mm), jnEnd(mm)
                 do i=inBeg(mm), inEnd(mm)

                   ! Compute the coordinates relative to the rotation
                   ! center.

                   xx = x(i,j,k,1) - center(1)
                   yy = x(i,j,k,2) - center(2)
                   zz = x(i,j,k,3) - center(3)

                   ! Compute the axial coordinate and check if this is
                   ! either a minimum or a maximum.

                   dot   = xx*axis(1) + yy*axis(2) + zz*axis(3)
                   axMin = min(axMin,dot)
                   axMax = max(axMax,dot)

                   ! Substract the axial component such that a vector
                   ! in the radial plane is kept.

                   xx = xx - dot*axis(1)
                   yy = yy - dot*axis(2)
                   zz = zz - dot*axis(3)

                   ! Compute the radius squared (stored in dot) and
                   ! check if this is either a minimum or maximum.

                   dot  = xx*xx + yy*yy + zz*zz
                   rMin = min(rMin,dot)
                   rMax = max(rMax,dot)

                 enddo
               enddo
             enddo

           endif

         enddo bocos2

       enddo domains2

       ! Take the square root of rMin and rMax to obtain the correct
       ! values.

       rMin = sqrt(rMin)
       rMax = sqrt(rMax)

       ! Determine the global minima and maxima of the axial and radial
       ! component for the entire sliding mesh.

       call mpi_allreduce(axMin, myInterfaces(color)%axMin, 1, &
                          sumb_real, mpi_min, commSlide, ierr)
       call mpi_allreduce(axMax, myInterfaces(color)%axMax, 1, &
                          sumb_real, mpi_max, commSlide, ierr)

       call mpi_allreduce(rMin, myInterfaces(color)%rMin, 1, &
                          sumb_real, mpi_min, commSlide, ierr)
       call mpi_allreduce(rMax, myInterfaces(color)%rMax, 1, &
                          sumb_real, mpi_max, commSlide, ierr)

       ! Determine the scale factor for the radial and axial
       ! coordinates for this sliding mesh interface.

       axMax = myInterfaces(color)%axMax - myInterfaces(color)%axMin
       rMax  = myInterfaces(color)%rMax  - myInterfaces(color)%rMin

       myInterfaces(color)%scale = one/max(axMax,rMax)

       ! Format statements.

 100   format("Sliding mesh", 1x,a,": Rotation axis for certain blocks &
              &are not aligned.")
 101   format("Sliding mesh", 1x,a,": Conflicting values &
              &for the rotation center.")
 102   format("Sliding mesh", 1x,a,": Rotation rate not specified.")
 103   format("Sliding mesh", 1x,a,": Rotation center not specified.")

       end subroutine initThisSlide
