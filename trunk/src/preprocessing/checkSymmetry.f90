!
!      ******************************************************************
!      *                                                                *
!      * File:          checkSymmetry.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-21-2003                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkSymmetry(level)
!
!      ******************************************************************
!      *                                                                *
!      * checkSymmetry checks whether or not the symmetry planes are    *
!      * really planar (within a certain tolerance). If this is not the *
!      * case for the finest level, a warning is printed. In all cases  *
!      * the unit normals are replaced by the face averaged unit        *
!      * normal.                                                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use cgnsGrid
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local parameter, tolerance for planar, 0.1 degrees.
!
       real(kind=realType), parameter :: tolDotmin = 0.9999985_realType
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, sps, i, j

       real(kind=realType) :: fact, dotMin, dot, mult

       real(kind=realType), dimension(3) :: faceNorm
       real(kind=realType), dimension(:,:,:), pointer :: ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and local domains.

       spectral: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, level, sps)

           ! Loop over the number of boundary subfaces for this block.

           bocos: do mm=1,nBocos

             ! Check for symmetry boundary condition.

             symmetry: if(BCType(mm) == symm) then

               ! Determine the block face on which this subface is
               ! located and set some variables accordingly.

               select case (BCFaceID(mm))

                 case (iMin)
                   mult = -one; ss => si(1,:,:,:)
                 case (iMax)
                   mult = one; ss => si(il,:,:,:)
                 case (jMin)
                   mult = -one; ss => sj(:,1,:,:)
                 case (jMax)
                   mult = one; ss => sj(:,jl,:,:)
                 case (kMin)
                   mult = -one; ss => sk(:,:,1,:)
                 case (kMax)
                   mult = one; ss => sk(:,:,kl,:)

               end select

               ! Loop over the range of the subface compute the face
               ! normal. The halo cells should not be taken into account,
               ! which explains why the nodal range of BCData is used.
               ! As the starting index of the cell range is shifted 1,
               ! (inBeg+1) and (jnBeg+1) are the starting indices for the
               ! owned cell range.

               faceNorm = zero

               do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
                 do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd
                   faceNorm(1) = faceNorm(1) + ss(i,j,1)
                   faceNorm(2) = faceNorm(2) + ss(i,j,2)
                   faceNorm(3) = faceNorm(3) + ss(i,j,3)
                 enddo
               enddo

               ! Create the unit normal for faceNorm. Make sure it
               ! is outward pointing by multiplying it by mult;
               ! mult is either 1.0 or -1.0.

               fact = sqrt(faceNorm(1)*faceNorm(1) &
                    +      faceNorm(2)*faceNorm(2) &
                    +      faceNorm(3)*faceNorm(3))
               if(fact > zero) fact = mult/fact

               faceNorm(1) = faceNorm(1)*fact
               faceNorm(2) = faceNorm(2)*fact
               faceNorm(3) = faceNorm(3)*fact

               ! Check if the symmetry plane is really planar. This is
               ! only done on the finest mesh and for the 1st spectral
               ! solution, because it is only to inform the user.
               ! Afterwards the normals will be reset to the unit
               ! normal of the face anyway.

               fineLevelTest: if(level == 1 .and. sps == 1) then

                 ! Initialize dotMin such that it will always
                 ! be overwritten.

                 dotMin = one

                 ! Loop over the physical faces of the symmetry plane,
                 ! i.e. no halo's.

                 do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
                   do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

                     ! Compute the dot product between the normal of
                     ! this face and the averaged normal of the plane.

                     dot = BCData(mm)%norm(i,j,1)*faceNorm(1) &
                         + BCData(mm)%norm(i,j,2)*faceNorm(2) &
                         + BCData(mm)%norm(i,j,3)*faceNorm(3)

                     ! And determine the minimum of dot and dotMin

                     dotMin = min(dot,dotMin)
                   enddo
                 enddo

                 ! Test if the minimum dot product is smaller than the
                 ! tolerance. If so, the plane is considered as not
                 ! planar.

                 if(dotMin < tolDotmin) then

                   ! Determine the corresponding angle in degrees of
                   ! dotmin.

                   fact = acos(dotMin)*180.0_realType/pi

                   ! Store the corresponding cgns block id and the
                   ! subface in this block a bit easier.

                   i = nbkGlobal
                   j = cgnsSubface(mm)

                   ! Print a warning.

                   print "(a)", "#"
                   print "(a)", "#                      Warning"
                   print 100,                              &
                     trim(cgnsDoms(i)%bocoInfo(j)%bocoName), &
                     trim(cgnsDoms(i)%zonename)
                   print 110, fact
                   print "(a)", "#"
 100               format("# Symmetry boundary face",1X,A,1X,"of zone", &
                          1x,a,1x, "is not planar.")
 110               format("# Maximum deviation from the mean normal: ", &
                          e12.5, " degrees")

                 endif

               endif fineLevelTest

               ! Set the unit normals to the unit normal of the entire
               ! plane. All the cells, also possible halo's, are treated.

               do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   BCData(mm)%norm(i,j,1) = faceNorm(1)
                   BCData(mm)%norm(i,j,2) = faceNorm(2)
                   BCData(mm)%norm(i,j,3) = faceNorm(3)

                 enddo
               enddo

             endif symmetry
           enddo bocos
         enddo domains
       enddo spectral

       end subroutine checkSymmetry
