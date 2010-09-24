!
!      ******************************************************************
!      *                                                                *
!      * File:          volumeInterpol.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-13-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine volumeInterpol(intInfo, realInfo, donorInfo, &
                                 nHalo, level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * volumeInterpol computes the actual interpolation weights for   *
!      * the halo's and stores these weights and the connectivity in    *
!      * donorInfo. The arguments sps is passed for consistency         *
!      * reasons, but it is not really needed. The info used from       *
!      * flowDoms is identical for all spectral solutions.              *
!      *                                                                *
!      ******************************************************************
!
       use updateComm
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nHalo, level, sps

       integer(kind=intType), dimension(*), intent(in) :: intInfo
       real(kind=realType),   dimension(*), intent(in) :: realInfo

       type(updateCommType), intent(out) :: donorInfo
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory of donorInfo, for which the amount is
       ! guessed to be 4 times the number of values to be interpolated.

       nn = 4*nHalo
       allocate(donorInfo%indBuf(nn),    donorInfo%block(nn),  &
                donorInfo%indices(nn,3), donorInfo%weight(nn), &
                stat=ierr)
       if(ierr /= 0)                      &
         call terminate("volumeInterpol", &
                        "Memory allocation failure for the member &
                        &variables of donorInfo.")

       ! First determine the interpolation weights and cell indices for
       ! the 1st level halo's.

       nn = 0
       call interpolHalos(1_intType)

       donorInfo%nCopy1st = nn

       ! Determine the interpolation weights and cell indices for the
       ! 2nd level halos's, which are added after the first level halo's.

       call interpolHalos(2_intType)

       donorInfo%nCopy2nd = nn

       !=================================================================

       contains

         !===============================================================

         subroutine interpolHalos(haloLevel)
!
!        ****************************************************************
!        *                                                              *
!        * interpolHalos determines the donor cells and the             *
!        * corresponding interpolation weights for the halo cells       *
!        * stored in intInfo and realInfo.                              *
!        *                                                              *
!        ****************************************************************
!
         use BCTypes
         use block
         use localSubfacesMod
         use thisSlide
         implicit none
!
!        Subroutine argument.
!
         integer(kind=intType), intent(in) :: haloLevel
!
!        Local variables.
!
         integer(kind=intType) :: ii, jj, mm, nDonor, nindDonor, nip
         integer(kind=intType) :: subface, quadID, donorQuad
         integer(kind=intType) :: blockID, faceID

         integer(kind=intType), dimension(4) :: ip

         integer(kind=intType), dimension(:,:), pointer :: indHalo
         integer(kind=intType), dimension(:,:), pointer :: conn

         real(kind=realType) :: u, v, uv
         real(kind=realType), dimension(4) :: weight
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the number of halo cells stored in intInfo
         ! and realInfo.

         loopDir: do mm=1,nHalo

           ! Determine the offset in intInfo for this halo cell.

           ii = 2*(mm-1)

           ! Check if the halo level of the cell to be interpolated is
           ! the correct level.

           testHaloLevel: if(intInfo(ii+2) == haloLevel) then

             ! Cell has must be stored. Abbreviate the other information
             ! stored in intInfo and realInfo a bit easier.

             donorQuad = intInfo(ii+1)

             u = realInfo(ii+1)
             v = realInfo(ii+2)

             ! Determine the block, face ID and quadrilateral offset and
             ! set the pointers for the connectivity and the indices.
             ! The convention is such that if donorQuad is negative the
             ! subface belongs to part 1 of the sliding interface and
             ! if it is positive it belongs to part 2 of that interface.
             ! Note that indHalo always point to indHalo1 even if the
             ! halo level is not 1. The reason is that indHalo is used
             ! to find the donor indices and the indices in the face
             ! are identical for both the halo and the donor.

             if(donorQuad < 0) then
               subface = subface1(-donorQuad)
               quadID  = quadID1(-donorQuad)

               blockID  = mySubfaces1(subface)%blockID
               faceID   = mySubfaces1(subface)%faceID
               conn    => mySubfaces1(subface)%connDual
               indHalo => mySubfaces1(subface)%indHalo1
             else
               subface = subface2(donorQuad)
               quadID  = quadID2(donorQuad)

               blockID  = mySubfaces2(subface)%blockID
               faceID   = mySubfaces2(subface)%faceID
               conn    => mySubfaces2(subface)%connDual
               indHalo => mySubfaces2(subface)%indHalo1
             endif

             ! Determine the donor index normal to the face and
             ! its corresponding value. This depends on the block face
             ! on which the subface is located.

             select case (faceID)

               case (iMin)
                 nDonor    = haloLevel + 1
                 nindDonor = 1
               case (iMax)
                 nDonor    = flowDoms(blockID,level,sps)%ie - haloLevel
                 nindDonor = 1
               case (jMin)
                 nDonor    = haloLevel + 1
                 nindDonor = 2
               case (jMax)
                 nDonor    = flowDoms(blockID,level,sps)%je - haloLevel
                 nindDonor = 2
               case (kMin)
                 nDonor    = haloLevel + 1
                 nindDonor = 3
               case (kMax)
                 nDonor    = flowDoms(blockID,level,sps)%ke - haloLevel
                 nindDonor = 3

             end select

             ! Store the cell center ID's and weights a bit easier.

             nip = 4; uv = u*v

             weight(1) = one - u - v + uv; ip(1) = conn(quadID,1)
             weight(2) =       u     - uv; ip(2) = conn(quadID,2)
             weight(3) =               uv; ip(3) = conn(quadID,3)
             weight(4) =           v - uv; ip(4) = conn(quadID,4)

             ! Update the donor info.

             do jj=1,nip
               if(weight(jj) > zero) then

                 nn = nn + 1
                 donorInfo%indBuf(nn) = mm
                 donorInfo%block(nn)  = blockID
                 donorInfo%weight(nn) = weight(jj)

                 donorInfo%indices(nn,1) = indHalo(ip(jj),1)
                 donorInfo%indices(nn,2) = indHalo(ip(jj),2)
                 donorInfo%indices(nn,3) = indHalo(ip(jj),3)

                 donorInfo%indices(nn,nindDonor) = nDonor

               endif
             enddo

           endif testHaloLevel
         enddo loopDir

         end subroutine interpolHalos

       end subroutine volumeInterpol
