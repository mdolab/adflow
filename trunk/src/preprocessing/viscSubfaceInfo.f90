!
!      ******************************************************************
!      *                                                                *
!      * File:          viscSubfaceInfo.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine viscSubfaceInfo(level)
!
!      ******************************************************************
!      *                                                                *
!      * viscSubfaceInfo allocates the memory for the storage of the    *
!      * stress tensor and heat flux vector of viscous subfaces for the *
!      * given multigrid level and all spectral solutions. Furthermore  *
!      * the pointers viscIminPointer, etc. Are allocated and set.      *
!      * These pointers contain info to which viscous subface the faces *
!      * of the block faces possibly belong. If not part of a viscous   *
!      * subface these values are set to 0. Note that these pointers    *
!      * are only allocated and determined for the 1st spectral         *
!      * solution, because the info is the same for all of them.        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputTimeSpectral
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, sps, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

       integer(kind=intType), dimension(:,:), pointer :: viscPointer
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number blocks stored on this processor.

       domains: do nn=1,nDom

         ! Set the pointers to the block of the 1st spectral solution.

         call setPointers(nn, level, 1_intType)

         ! Allocate the memory for viscSubface and the pointers
         ! viscIminPointer, etc. ViscSubface must be allocated for
         ! all spectral solutions, the pointers only for the 1st.

         do sps=1,nTimeIntervalsSpectral
           allocate(flowDoms(nn,level,sps)%viscSubface(nViscBocos), &
                    stat=ierr)
           if(ierr /= 0)                         &
             call terminate("viscSubfaceInfo", &
                            "Memory allocation failure for viscSubface")
         enddo
 
         allocate(flowDoms(nn,level,1)%viscIminPointer(2:jl,2:kl), &
                  flowDoms(nn,level,1)%viscImaxPointer(2:jl,2:kl), &
                  flowDoms(nn,level,1)%viscJminPointer(2:il,2:kl), &
                  flowDoms(nn,level,1)%viscJmaxPointer(2:il,2:kl), &
                  flowDoms(nn,level,1)%viscKminPointer(2:il,2:jl), &
                  flowDoms(nn,level,1)%viscKmaxPointer(2:il,2:jl), &
                  stat=ierr)
         if(ierr /= 0)                         &
           call terminate("viscSubfaceInfo", &
                          "Memory allocation failure for subface info")

         ! Reset the pointers viscIminPointer, etc. to make it more
         ! readable and initialize them to 0. This indicates that
         ! the faces are not part of a viscous wall subfaces.

         viscIminPointer => flowDoms(nn,level,1)%viscIminPointer
         viscImaxPointer => flowDoms(nn,level,1)%viscImaxPointer
         viscJminPointer => flowDoms(nn,level,1)%viscJminPointer
         viscJmaxPointer => flowDoms(nn,level,1)%viscJmaxPointer
         viscKminPointer => flowDoms(nn,level,1)%viscKminPointer
         viscKmaxPointer => flowDoms(nn,level,1)%viscKmaxPointer

         viscIminPointer = 0
         viscImaxPointer = 0
         viscJminPointer = 0
         viscJmaxPointer = 0
         viscKminPointer = 0
         viscKmaxPointer = 0

         ! Loop over the viscous subfaces to allocate the memory for the
         ! stress tensor and the heat flux vector and to set the range
         ! in viscIminPointer, etc.

         viscSubfaces: do mm=1,nViscBocos

           ! Store the cell range in iBeg, iEnd, etc. As the viscous data
           ! do not allow for an overlap, the nodal range of the
           ! subface must be used.

           iBeg = BCData(mm)%inBeg + 1
           iEnd = BCData(mm)%inEnd

           jBeg = BCData(mm)%jnBeg + 1
           jEnd = BCData(mm)%jnEnd

           ! Loop over the spectral solutions and allocate the memory
           ! for the stress tensor, heat flux and friction velocity.

           do sps=1,nTimeIntervalsSpectral

             ! Set the pointer for viscSubface to make the code
             ! more readable and allocate the memory.

             viscSubface => flowDoms(nn,level,sps)%viscSubface

             allocate(viscSubface(mm)%tau( iBeg:iEnd,jBeg:jEnd,6), &
                      viscSubface(mm)%q(   iBeg:iEnd,jBeg:jEnd,3), &
                      viscSubface(mm)%utau(iBeg:iEnd,jBeg:jEnd),   &
                      stat=ierr)
             if(ierr /= 0)                       &
               call terminate("viscSubfaceInfo", &
                              "Memory allocation failure for tau, q &
                              &and utau.")
           enddo

           ! Set the pointer viscPointer, depending on the block face
           ! on which the subface is located.

           select case (BCFaceID(mm))
             case (iMin)
               viscPointer => viscIminPointer

             case (iMax)
               viscPointer => viscImaxPointer

             case (jMin)
               viscPointer => viscJminPointer

             case (jMax)
               viscPointer => viscJmaxPointer

             case (kMin)
               viscPointer => viscKminPointer

             case (kMax)
               viscPointer => viscKmaxPointer
           end select

           ! Set this range in viscPointer to viscous subface mm.

           do j=jBeg,jEnd
             do i=iBeg,iEnd
               viscPointer(i,j) = mm
             enddo
           enddo

         enddo viscSubfaces

       enddo domains

       end subroutine viscSubfaceInfo
