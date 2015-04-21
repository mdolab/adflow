       module checkVolBlock
!
!      ******************************************************************
!      *                                                                *
!      * Local module, which contains the definition of the derived     *
!      * datatype used to test for negative volumes in the grid.        *
!      *                                                                *
!      ******************************************************************
!
       implicit none
       save

       type checkVolBlockType

         ! blockHasNegVol:              Whether or not the block
         !                              contains negative volumes.
         ! volumeIsNeg(2:il,2:jl,2:kl): Whether or not the owned volumes
         !                              are negative.

         logical :: blockHasNegVol
         logical, dimension(:,:,:), pointer :: volumeIsNeg

       end type checkVolBlockType

       end module checkVolBlock

       subroutine metric_ALE(level, idxALE)
!
!      ******************************************************************
!      *                                                                *
!      * metric computes the face normals and the volume for the given  *
!      * grid level for all spectral solutions. First the volumes are   *
!      * computed assuming that the block is right handed. Then the     *
!      * number of positive and negative volumes are determined. If all *
!      * volumes are positive the block is indeed right handed; if all  *
!      * volumes are negative the block is left handed and both the     *
!      * volumes and the normals must be negated (for the normals this  *
!      * is done by the introduction of fact, which is either -0.5 or   *
!      * 0.5); if there are both positive and negative volumes the mesh *
!      * is not valid.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use inputTimeSpectral
       use checkVolBlock
       use inputIteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
       integer(kind=intType), intent(in) :: idxALE
!
!      Local parameter.
!
       real(kind=realType), parameter :: thresVolume = 1.e-2_realType
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, n, m, l
       integer(kind=intType) :: nn, mm, sps
       integer(kind=intType) :: nVolNeg,   nVolPos
       integer(kind=intType) :: nVolBad,   nVolBadGlobal
       integer(kind=intType) :: nBlockBad, nBlockBadGlobal

       real(kind=realType) :: fact, mult
       real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6

       real(kind=realType), dimension(3) :: v1, v2

       real(kind=realType), dimension(:,:,:), pointer :: ss

       character(len=10) :: integerString

       logical :: checkK, checkJ, checkI, checkAll
       logical :: badVolume

       logical, dimension(:,:,:), pointer :: volumeIsNeg

       type(checkVolBlockType), &
                  dimension(nDom,nTimeIntervalsSpectral) :: checkVolDoms
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the number of bad volumes and bad blocks to 0.

       nVolBad   = 0
       nBlockBad = 0

       ! Loop over the number of spectral solutions and local blocks.

       spectral: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers to this block and allocate the memory for
           ! volumeIsNeg. Set a pointer to this entry afterwards to make
           ! the code more readable.

           call setPointers(nn, level, sps)

           allocate(checkVolDoms(nn,sps)%volumeIsNeg(2:il,2:jl,2:kl), &
                    stat=ierr)
           if(ierr /= 0)              &
             call terminate("metric", &
                            "Memory allocation failure for volumeIsNeg")
           volumeIsNeg => checkVolDoms(nn,sps)%volumeIsNeg
!
!          **************************************************************
!          *                                                            *
!          * Volume and block orientation computation.                  *
!          *                                                            *
!          **************************************************************
!
           ! Initialize the number of positive and negative volumes for
           ! this block to 0.

           nVolNeg = 0
           nVolPos = 0

           ! Compute the volumes. The hexahedron is split into 6 pyramids
           ! whose volumes are computed. The volume is positive for a
           ! right handed block.
           ! Initialize the volumes to zero. The reasons is that the second
           ! level halo's must be initialized to zero and for convenience
           ! all the volumes are set to zero.

           vol = zero

           do k=1,ke
             n = k -1

             checkK = .true.
             if(k == 1 .or. k == ke) checkK = .false.

             do j=1,je
               m = j -1

               checkJ = .true.
               if(j == 1 .or. j == je) checkJ = .false.

               do i=1,ie
                 l = i -1

                 checkI = .true.
                 if(i == 1 .or. i == ie) checkI = .false.

                 ! Determine whether or not the voluem must be checked for
                 ! quality. Only owned volumes are checked, not halo's.

                 checkAll = .false.
                 if(checkK .and. checkJ .and. checkI) checkAll = .true.

                 ! Compute the coordinates of the center of gravity.

                 xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
                    +         x(i,m,n,1) + x(i,j,n,1) &
                    +         x(l,j,k,1) + x(l,m,k,1) &
                    +         x(l,m,n,1) + x(l,j,n,1))
                 yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
                    +         x(i,m,n,2) + x(i,j,n,2) &
                    +         x(l,j,k,2) + x(l,m,k,2) &
                    +         x(l,m,n,2) + x(l,j,n,2))
                 zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
                    +         x(i,m,n,3) + x(i,j,n,3) &
                    +         x(l,j,k,3) + x(l,m,k,3) &
                    +         x(l,m,n,3) + x(l,j,n,3))

                 ! Compute the volumes of the 6 sub pyramids. The
                 ! arguments of volpym must be such that for a (regular)
                 ! right handed hexahedron all volumes are positive.

                 vp1 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                              x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                              x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                              x(i,m,k,1), x(i,m,k,2), x(i,m,k,3))

                 vp2 = volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                              x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                              x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                              x(l,j,n,1), x(l,j,n,2), x(l,j,n,3))

                 vp3 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                              x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                              x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                              x(i,j,n,1), x(i,j,n,2), x(i,j,n,3))

                 vp4 = volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                              x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                              x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                              x(l,m,k,1), x(l,m,k,2), x(l,m,k,3))

                 vp5 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                              x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                              x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                              x(l,j,k,1), x(l,j,k,2), x(l,j,k,3))

                 vp6 = volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                              x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                              x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                              x(i,m,n,1), x(i,m,n,2), x(i,m,n,3))

                 ! Set the volume to 1/6 of the sum of the volumes of the
                 ! pyramid. Remember that volpym computes 6 times the
                 ! volume.

                 vol(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)

                 ! Check the volume and update the number of positive
                 ! and negative volumes if needed.

                 if( checkAll ) then

                   ! Update either the number of negative or positive
                   ! volumes. Negative volumes should only occur for left
                   ! handed blocks. This is checked later.
                   ! Set the logical volumeIsNeg accordingly.

                   if(vol(i,j,k) < zero) then
                     nVolNeg            = nVolNeg + 1
                     volumeIsNeg(i,j,k) = .true.
                   else
                     nVolPos            = nVolPos + 1
                     volumeIsNeg(i,j,k) = .false.
                   endif

                   ! Set the threshold for the volume quality.

                   fact = thresVolume*abs(vol(i,j,k))

                   ! Check the quality of the volume.

                   badVolume = .false.
                   if(vp1*vol(i,j,k) < zero .and. &
                      abs(vp1)       > fact) badVolume = .true.
                   if(vp2*vol(i,j,k) < zero .and. &
                      abs(vp2)       > fact) badVolume = .true.
                   if(vp3*vol(i,j,k) < zero .and. &
                      abs(vp3)       > fact) badVolume = .true.
                   if(vp4*vol(i,j,k) < zero .and. &
                      abs(vp4)       > fact) badVolume = .true.
                   if(vp5*vol(i,j,k) < zero .and. &
                      abs(vp5)       > fact) badVolume = .true.
                   if(vp6*vol(i,j,k) < zero .and. &
                      abs(vp6)       > fact) badVolume = .true.

                   ! Update nVolBad if this is a bad volume.

                   if( badVolume ) nVolBad = nVolBad + 1

                 endif

                 ! Set the volume to the absolute value.

                 vol(i,j,k) = abs(vol(i,j,k))

               enddo
             enddo
           enddo

           ! Some additional safety stuff for halo volumes.

           do k=2,kl
             do j=2,jl
               if(vol(1, j,k) <= eps) vol(1, j,k) = vol(2, j,k)
               if(vol(ie,j,k) <= eps) vol(ie,j,k) = vol(il,j,k)
             enddo
           enddo

           do k=2,kl
             do i=1,ie
               if(vol(i,1, k) <= eps) vol(i,1, k) = vol(i,2, k)
               if(vol(i,je,k) <= eps) vol(i,je,k) = vol(i,jl,k)
             enddo
           enddo

           do j=1,je
             do i=1,ie
               if(vol(i,j,1)  <= eps) vol(i,j,1)  = vol(i,j,2)
               if(vol(i,j,ke) <= eps) vol(i,j,ke) = vol(i,j,kl)
             enddo
           enddo

           ! Determine the orientation of the block. For the fine level
           ! this is based on the number of positive and negative
           ! volumes; on the coarse levels the corresponding fine level
           ! value is taken. If both positive and negative volumes are
           ! present it is assumed that the block was intended to be
           ! right handed. The code will terminate later on anyway.

           if(level == 1) then
             if(nVolPos == 0) then       ! Left handed block.
               flowDoms(nn,level,sps)%rightHanded = .false.
             else                        ! Right handed (or bad) block.
               flowDoms(nn,level,sps)%rightHanded = .true.
             endif
           else
             flowDoms(nn,level,sps)%rightHanded = &
                 flowDoms(nn,1,sps)%rightHanded
           endif

           ! Set the factor in the surface normals computation. For a
           ! left handed block this factor is negative, such that the
           ! normals still point in the direction of increasing index.
           ! The formulae used later on assume a right handed block
           ! and fact is used to correct this for a left handed block,
           ! as well as the scaling factor of 0.5

           if( flowDoms(nn,level,sps)%rightHanded ) then
             fact =  half
           else
             fact = -half
           endif

           ! Check if both positive and negative volumes occur. If so,
           ! the block is bad and the counter nBlockBad is updated.

           if(nVolNeg > 0 .and. nVolPos > 0) then
             checkVolDoms(nn,sps)%blockHasNegVol = .true.
             nBlockBad = nBlockBad + 1
           else
             checkVolDoms(nn,sps)%blockHasNegVol = .false.
           endif
!
!          **************************************************************
!          *                                                            *
!          * Computation of the face normals in i-, j- and k-direction. *
!          * Formula's are valid for a right handed block; for a left   *
!          * handed block the correct orientation is obtained via fact. *
!          * The normals point in the direction of increasing index.    *
!          * The absolute value of fact is 0.5, because the cross       *
!          * product of the two diagonals is twice the normal vector.   *
!          *                                                            *
!          * Note that also the normals of the first level halo cells   *
!          * are computed. These are needed for the viscous fluxes.     *
!          *                                                            *
!          **************************************************************
!
           ! Projected areas of cell faces in the i direction.

           do k=1,ke
             n = k -1
             do j=1,je
               m = j -1
               do i=0,ie

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = x(i,j,n,1) - x(i,m,k,1)
                 v1(2) = x(i,j,n,2) - x(i,m,k,2)
                 v1(3) = x(i,j,n,3) - x(i,m,k,3)

                 v2(1) = x(i,j,k,1) - x(i,m,n,1)
                 v2(2) = x(i,j,k,2) - x(i,m,n,2)
                 v2(3) = x(i,j,k,3) - x(i,m,n,3)

                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                 

                 sIALE(idxALE,i,j,k,1) = si(i,j,k,1)
                 sIALE(idxALE,i,j,k,2) = si(i,j,k,2)
                 sIALE(idxALE,i,j,k,3) = si(i,j,k,3)

               enddo
             enddo
           enddo

           ! Projected areas of cell faces in the j direction.

           do k=1,ke
             n = k -1
             do j=0,je
               do i=1,ie
                 l = i -1

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = x(i,j,n,1) - x(l,j,k,1)
                 v1(2) = x(i,j,n,2) - x(l,j,k,2)
                 v1(3) = x(i,j,n,3) - x(l,j,k,3)

                 v2(1) = x(l,j,n,1) - x(i,j,k,1)
                 v2(2) = x(l,j,n,2) - x(i,j,k,2)
                 v2(3) = x(l,j,n,3) - x(i,j,k,3)

                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

 
                 sJALE(idxALE,i,j,k,1) = sj(i,j,k,1)
                 sJALE(idxALE,i,j,k,2) = sj(i,j,k,2)
                 sJALE(idxALE,i,j,k,3) = sj(i,j,k,3)

               enddo
             enddo
           enddo

           ! Projected areas of cell faces in the k direction.

           do k=0,ke
             do j=1,je
               m = j -1
               do i=1,ie
                 l = i -1

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = x(i,j,k,1) - x(l,m,k,1)
                 v1(2) = x(i,j,k,2) - x(l,m,k,2)
                 v1(3) = x(i,j,k,3) - x(l,m,k,3)

                 v2(1) = x(l,j,k,1) - x(i,m,k,1)
                 v2(2) = x(l,j,k,2) - x(i,m,k,2)
                 v2(3) = x(l,j,k,3) - x(i,m,k,3)
                 
                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

                 sKALE(idxALE,i,j,k,1) = sk(i,j,k,1)
                 sKALE(idxALE,i,j,k,2) = sk(i,j,k,2)
                 sKALE(idxALE,i,j,k,3) = sk(i,j,k,3)

               enddo
             enddo
           enddo
!
!          **************************************************************
!          *                                                            *
!          * The unit normals on the boundary faces. These always point *
!          * out of the domain, so a multiplication by -1 is needed for *
!          * the iMin, jMin and kMin boundaries.                        *
!          *                                                            *
!          **************************************************************
!
           ! Loop over the boundary subfaces of this block.

           bocoLoop: do mm=1,nBocos

             ! Determine the block face on which this subface is located
             ! and set ss and mult accordingly.

             select case (BCFaceID(mm))

               case (iMin)
                 mult = -one; ss => si(1,:,:,:)

               case (iMax)
                 mult = one;  ss => si(il,:,:,:)

               case (jMin)
                 mult = -one; ss => sj(:,1,:,:)

               case (jMax)
                 mult = one;  ss => sj(:,jl,:,:)

               case (kMin)
                 mult = -one; ss => sk(:,:,1,:)

               case (kMax)
                 mult = one;  ss => sk(:,:,kl,:)

             end select

             ! Loop over the boundary faces of the subface.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
               do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ! Compute the inverse of the length of the normal vector
                 ! and possibly correct for inward pointing.

                 xp = ss(i,j,1);  yp = ss(i,j,2);  zp = ss(i,j,3)
                 fact = sqrt(xp*xp + yp*yp + zp*zp)
                 if(fact > zero) fact = mult/fact

                 ! Compute the unit normal.

                 BCData(mm)%norm(i,j,1) = fact*xp
                 BCData(mm)%norm(i,j,2) = fact*yp
                 BCData(mm)%norm(i,j,3) = fact*zp

               enddo
             enddo

           enddo bocoLoop
!
!          **************************************************************
!          *                                                            *
!          * Check in debug mode the sum of the normals of the cells.   *
!          * If everything is correct this should sum up to zero.       *
!          *                                                            *
!          **************************************************************
!
           debugging: if( debug ) then

             ! Loop over the cells including the 1st level halo's.

             do k=2,kl
               n = k -1
               do j=2,jl
                 m = j -1
                 do i=2,il
                   l = i -1

                   ! Store the sum of the outward pointing surrounding
                   ! normals in v1. Due to the outward convention the
                   ! normals with the lowest index get a negative sign;
                   ! normals point in the direction of the higher index.

                   v1(1) = si(i,j,k,1) + sj(i,j,k,1) + sk(i,j,k,1) &
                         - si(l,j,k,1) - sj(i,m,k,1) - sk(i,j,n,1)
                   v1(2) = si(i,j,k,2) + sj(i,j,k,2) + sk(i,j,k,2) &
                         - si(l,j,k,2) - sj(i,m,k,2) - sk(i,j,n,2)
                   v1(3) = si(i,j,k,3) + sj(i,j,k,3) + sk(i,j,k,3) &
                         - si(l,j,k,3) - sj(i,m,k,3) - sk(i,j,n,3)

                   ! Store the inverse of the sum of the areas of the
                   ! six faces in fact.

                   fact = one/(sqrt(si(i,j,k,1)*si(i,j,k,1)  &
                        +           si(i,j,k,2)*si(i,j,k,2)  &
                        +           si(i,j,k,3)*si(i,j,k,3)) &
                        +      sqrt(si(l,j,k,1)*si(l,j,k,1)  &
                        +           si(l,j,k,2)*si(l,j,k,2)  &
                        +           si(l,j,k,3)*si(l,j,k,3)) &
                        +      sqrt(sj(i,j,k,1)*sj(i,j,k,1)  &
                        +           sj(i,j,k,2)*sj(i,j,k,2)  &
                        +           sj(i,j,k,3)*sj(i,j,k,3)) &
                        +      sqrt(sj(i,m,k,1)*sj(i,m,k,1)  &
                        +           sj(i,m,k,2)*sj(i,m,k,2)  &
                        +           sj(i,m,k,3)*sj(i,m,k,3)) &
                        +      sqrt(sk(i,j,k,1)*sk(i,j,k,1)  &
                        +           sk(i,j,k,2)*sk(i,j,k,2)  &
                        +           sk(i,j,k,3)*sk(i,j,k,3)) &
                        +      sqrt(sk(i,j,n,1)*sk(i,j,n,1)  &
                        +           sk(i,j,n,2)*sk(i,j,n,2)  &
                        +           sk(i,j,n,3)*sk(i,j,n,3)))

                   ! Multiply v1 by fact to obtain a nonDimensional
                   ! quantity and take tha absolute value of it.

                   v1(1) = abs(v1(1)*fact)
                   v1(2) = abs(v1(2)*fact)
                   v1(3) = abs(v1(3)*fact)

                   ! Check if the control volume is closed.

                   if(v1(1) > thresholdReal .or. &
                      v1(2) > thresholdReal .or. &
                      v1(3) > thresholdReal)     &
                     call terminate("metric", &
                                    "Normals do not sum up to 0")

                 enddo
               enddo
             enddo

           endif debugging

         enddo domains
       enddo spectral

       ! Determine the global number of bad blocks. The result must be
       ! known on all processors and thus an allreduce is needed.

       call mpi_allreduce(nBlockBad, nBlockBadGlobal, 1, sumb_integer, &
                          mpi_sum, SUmb_comm_world, ierr)

       ! Test if bad blocks are present in the grid. If so, the action
       ! taken depends on the grid level.

       if(nBlockBadGlobal > 0) then
         if(level == 1) then

           ! Negative volumes present on the fine grid level. Print a
           ! list of the bad volumes and terminate executation.

           call writeNegVolumes(checkVolDoms)

           if(myID == 0) &
             call terminate("metric", "Negative volumes present in grid")
           call mpi_barrier(SUmb_comm_world, ierr)

         else

           ! Coarser grid level. The fine grid is okay, but due to the
           ! coarsening negative volumes are introduced. Print a warning.

           if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print 100, level
 100         format("#* Negative volumes present on coarse grid &
                    &level",1x,i1,".")
             print "(a)", "#* Computation continues, &
                          &but be aware of this"
             print "(a)", "#"
           endif

         endif
       endif

       ! Determine the global number of bad volumes. The result will
       ! only be known on processor 0. The quality volume check will
       ! only be done for the finest grid level.

       if(level == 1) then
         call mpi_reduce(nVolBad, nVolBadGlobal, 1, sumb_integer, &
                         mpi_sum, 0, SUmb_comm_world, ierr)

         ! Print a warning in case bad volumes were found. Only processor
         ! 0 prints this warning.

         if(myID == 0 .and. nVolBadGlobal > 0 .and. printWarnings) then
           write(integerString,"(i10)") nVolBadGlobal
           integerString = adjustl(integerString)
           integerString = trim(integerString)
           print "(a)", "#"
           print "(a)", "#                      Warning"
           print 101, trim(integerString)
           print 102
           print "(a)", "#"
 101       format("# ",a," bad quality volumes found.")
 102       format("# Computation will continue, but be aware of this")
         endif
       endif

       ! Release the memory of volumeIsNeg of all local blocks again.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom
           deallocate(checkVolDoms(nn,sps)%volumeIsNeg, stat=ierr)
           if(ierr /= 0)              &
             call terminate("metric", &
                            "Deallocation failure for volumeIsNeg")
         enddo
       enddo

       contains

!        ================================================================

         function volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd)
!
!        ****************************************************************
!        *                                                              *
!        * volpym computes 6 times the volume of a pyramid. Node p,     *
!        * whose coordinates are set in the subroutine metric itself,   *
!        * is the top node and a-b-c-d is the quadrilateral surface.    *
!        * It is assumed that the cross product vCa * vDb points in     *
!        * the direction of the top node. Here vCa is the diagonal      *
!        * running from node c to node a and vDb the diagonal from      *
!        * node d to node b.                                            *
!        *                                                              *
!        ****************************************************************
!
         use precision
         implicit none
!
!        Function type.
!
         real(kind=realType) :: volpym
!
!        Function arguments.
!
         real(kind=realType), intent(in) :: xa, ya, za, xb, yb, zb
         real(kind=realType), intent(in) :: xc, yc, zc, xd, yd, zd
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         volpym = (xp - fourth*(xa + xb  + xc + xd))              &
                * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
                  (yp - fourth*(ya + yb  + yc + yd))              &
                * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
                  (zp - fourth*(za + zb  + zc + zd))              &
                * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

         end function volpym

       end subroutine metric_ALE
