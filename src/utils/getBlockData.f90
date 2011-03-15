!
!      ******************************************************************
!      *                                                                *
!      * FILE:          getBlockData.f90                                *
!      * AUTHOR:        C.A.(Sandy) Mader                               *
!      * STARTING DATE: Oct. 11,2007                                    *
!      * LAST MODIFIED: Oct. 11,2007                                    *
!      *                                                                *
!      ******************************************************************
!
       SUBROUTINE getBlockCoords(BLOCKNUM,imax,jmax,kmax,XYZ)
!
!      ******************************************************************
!      *                                                                *
!      * This subroutine returns the block information from the block   *
!      * module in SUmbVertex to python for use in the mesh warping     *
!      * routine.                                                       *
!      * Each of the components is return with its original shape and   *
!      * size.                                                          *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       IMPLICIT NONE
!
!      SUBROUTINE ARGUMENTS
!
       INTEGER(KIND=intType), INTENT(IN) :: BLOCKNUM
       integer(kind=intType), intent(in):: imax,jmax,kmax
       REAL(KIND=REALTYPE), DIMENSION(3,imax,jmax,kmax), INTENT(OUT) :: XYZ
!
!      LOCAL VARIABLES
!
       INTEGER(KIND=INTTYPE) :: I, J, K, COUNT,N
!
       !print *,'coords',blocknum,imax,jmax,kmax
       call setPointers(blocknum,1,1) !(Blocknumber,MG level,Sectral Solution)
       
       DO K=1,kmax
          DO J=1,jmax
             DO I=1,imax
                do N = 1,3
                   XYZ(n,i,j,k) = X(I,J,K,n)
                   !print *,'xyz',XYZ(i,j,k,n),i,j,k,n
                enddo
             ENDDO
          ENDDO
       ENDDO
       !print *,'end function',shape(xyz)
       RETURN
       END SUBROUTINE GETBLOCKCOORDS

!!$       SUBROUTINE getBlockNodes(BLOCKNUM,imax,jmax,kmax,XYZ)
!!$!
!!$!      ******************************************************************
!!$!      *                                                                *
!!$!      * This subroutine returns the block information from the block   *
!!$!      * module in SUmbVertex to python for use in the mesh warping     *
!!$!      * routine.                                                       *
!!$!      * Each of the components is return with its original shape and   *
!!$!      * size.                                                          *
!!$!      *                                                                *
!!$!      ******************************************************************
!!$!
!!$       use blockPointers
!!$       IMPLICIT NONE
!!$!
!!$!      SUBROUTINE ARGUMENTS
!!$!
!!$       INTEGER(KIND=intType), INTENT(IN) :: BLOCKNUM
!!$       integer(kind=intType), intent(in):: imax,jmax,kmax
!!$       REAL(KIND=REALTYPE), DIMENSION(3,imax,jmax,kmax), INTENT(OUT) :: XYZ
!!$!
!!$!      LOCAL VARIABLES
!!$!
!!$       INTEGER(KIND=INTTYPE) :: I, J, K, COUNT,N
!!$!
!!$       !print *,'coords',blocknum,imax,jmax,kmax
!!$       call setPointers(blocknum,1,1) !(Blocknumber,MG level,Sectral Solution)
!!$       
!!$       DO K=1,kmax
!!$          DO J=1,jmax
!!$             DO I=1,imax
!!$                do N = 1,3
!!$                   !XYZ(n,i,j,k) = X(I,J,K,n)
!!$                   XYZ(n,i,j,k) = globalNode(I,J,K)*3+n-1
!!$                   !print *,'xyz',XYZ(i,j,k,n),i,j,k,n
!!$                enddo
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$       !print *,'end function',shape(xyz)
!!$       RETURN
!!$       END SUBROUTINE GETBLOCKNodes


       subroutine setBlockCoords(BLOCKNUM,imax,jmax,kmax,XYZ)
!
!      ******************************************************************
!      *                                                                *
!      * This subroutine returns the block information from the block   *
!      * module in SUmbVertex to python for use in the mesh warping     *
!      * routine.                                                       *
!      * Each of the components is return with its original shape and   *
!      * size.                                                          *
!      *                                                                *
!      ******************************************************************
!
         use blockPointers
         implicit none
!
!      SUBROUTINE ARGUMENTS
!
         integer(kind=intType), intent(in) :: blocknum
         integer(kind=intType), intent(in):: imax,jmax,kmax
         real(kind=realtype), dimension(3,imax,jmax,kmax), intent(in) :: xyz
!
!      LOCAL VARIABLES
!
         integer(kind=inttype) :: i, j, k, count,n
!
         !print *,'coords',blocknum,imax,jmax,kmax
         call setPointers(blocknum,1,1) !(Blocknumber,MG level,Sectral Solution)
         !print *,'fortran compare',xyz(1,1,1,1)
         do K=1,kmax
            do J=1,jmax
               do I=1,imax
                  do N = 1,3
                     !print *,'ijk',i,j,k,n
                     !print *,'xyz',XYZ(i,j,k,n)
                     !print *,'x',X(i,j,k,n)
                     X(i,j,k,n) = xyz(n,I,J,K)
                     !print *,'xyz',XYZ(i,j,k,n)!,globalnode(i,j,k),i,j,k,n
                  enddo
               enddo
            enddo
         enddo
         !print *,'end function',shape(xyz)
         return
       end subroutine setblockcoords



       SUBROUTINE getBlockDims(BLOCKNUM,IMAX,JMAX,KMAX)
!
!      ******************************************************************
!      *                                                                *
!      * THIS SUBROUTINE ALLOWS US TO GET THE DIMENSIONS OF X FOR A     *
!      * PARTICULAR block in the SUmbVertex mesh.                       *
!      *                                                                *
!      ******************************************************************
!
       use blockpointers
       IMPLICIT NONE
!
!      SUBROUTINE ARGUMENTS.
!
       INTEGER(KIND=INTTYPE), INTENT(IN) :: BLOCKNUM
       INTEGER(KIND=INTTYPE), INTENT(OUT) :: IMAX, JMAX, KMAX

       !Begin Execution
       
       call setPointers(blocknum,1,1) !(Blocknumber,MG level,Spectral Solution)
       
       !print *,'blocknum',blocknum,ndom
!
       IF (BLOCKNUM > nDom .OR. BLOCKNUM < 1) THEN
           IMAX = -1
           JMAX = -1
           KMAX = -1
       ELSE
           IMAX = IL
           JMAX = JL
           KMAX = KL
       END IF
 
       RETURN
       END SUBROUTINE GETBLOCKDIMS

       SUBROUTINE getBlockCGNSID(BLOCKNUM,CGNS_BLOCK_ID)
!
!      ******************************************************************
!      *                                                                *
!      * THIS SUBROUTINE ALLOWS US TO GET THE CGNS_BLOCK_ID FOR A       *
!      * PARTICULAR block in SUmbVertex                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockpointers
       IMPLICIT NONE
!
!      SUBROUTINE ARGUMENTS
!
       INTEGER(KIND=INTTYPE), INTENT(IN) :: BLOCKNUM
       INTEGER(KIND=INTTYPE), INTENT(OUT) :: CGNS_BLOCK_ID
!
       CGNS_BLOCK_ID = flowdoms(blocknum,1,1)%cgnsblockid

       RETURN
       END SUBROUTINE GETBLOCKCGNSID

       
       subroutine getnBlocksLocal(nBlocks)

!
!      ******************************************************************
!      *                                                                *
!      * THIS SUBROUTINE ALLOWS US TO GET THE number of blocks on this  *
!      * processor after splitting                                      *
!      *                                                                *
!      ******************************************************************
!
       use blockpointers
       IMPLICIT NONE
!
!      SUBROUTINE ARGUMENTS
!
       INTEGER(KIND=INTTYPE), INTENT(out) :: nBlocks

       !Begin execution

       nBlocks = nDom

       return
  
       end subroutine getnBlocksLocal
       
       subroutine getnSubfacesBlock(blocknum,nsubfaces, n1to1internals, n1to1s, nnonmatchs,numBocos)
         use blockpointers
         implicit none
         !Subroutine Variables
         integer(kind=intType), intent(in)::blocknum
         integer(kind=intType), intent(out)::nSubfaces, n1to1Internals, n1to1s, nNonMatchs, numBocos
         !Begin Execution
         call setPointers(blocknum,1,1)

         nsubfaces=flowDoms(blocknum,1,1)%nSubface
         n1to1Internals = 0!flowDoms(blocknum,1,1)% n1to1Internal
         n1to1s = flowDoms(blocknum,1,1)% n1to1
         nNonMatchs = 0!flowDoms(blocknum,1,1)% nNonMatch
         numBocos = flowDoms(blocknum,1,1)%nBocos

         return
       end subroutine getnSubfacesBlock

       subroutine getBlockCommunicationInfo(blocknum,nSub,nOne,nNon,i1,i2,j1,j2,k1,k2,&
            di1,di2,dj1,dj2,dk1,dk2,neighbourBlock,neighbourProc,trans1,&
            trans2,trans3,bcinfo,faceinfo)
         
         use blockpointers
         implicit none
         
         !Subroutine Variables
         integer(kind=intType),intent(in)::blocknum,nSub,nOne,nNon

         integer(kind=intType),dimension(nSub),intent(out)::i1,i2,j1,j2,k1,k2
         integer(kind=intType),dimension(nSub),intent(out)::di1,di2,dj1,dj2,dk1,dk2
         integer(kind=intType),dimension(nSub),intent(out):: neighbourBlock,&
              &neighbourProc
         integer(kind=intType),dimension(nSub),intent(out):: trans1,trans2,trans3
         integer(kind=inttype),dimension(nSub),intent(out)::bcinfo,faceinfo
         !Local Variables
         integer(kind=intType):: i
         
         !begin execution

         call setPointers(blocknum,1,1)

         do i = 1,nSub
            !print *,'bctype',BCType(i)
            i1(i) = flowDoms(blocknum,1,1)%inBeg(i)
            !print *,'i',shape(ibeg)
            i2(i) = flowDoms(blocknum,1,1)%inEnd(i)
            j1(i) = flowDoms(blocknum,1,1)%jnBeg(i)
            j2(i) = flowDoms(blocknum,1,1)%jnEnd(i)
            k1(i) = flowDoms(blocknum,1,1)%knBeg(i)
            k2(i) = flowDoms(blocknum,1,1)%knEnd(i)
            bcinfo(i)=BCType(i)
            faceinfo(i) = BCFaceID(i)

            if (BCType(i)== -16) then
               di1(i) = dinbeg(i)
               !print *,'di',di1(i),dibeg(i),flowDoms(blocknum,1,1)%diBeg(i)
               di2(i) = dinEnd(i)
               dj1(i) = djnBeg(i)
               dj2(i) = djnEnd(i)
               dk1(i) = dknBeg(i)
               dk2(i) = dknEnd(i)

               neighbourBlock(i) =  neighBlock(i)
               !print *,'nebl',neighbourBlock(i),i
               neighbourProc(i)=  neighProc(i)
               !print *,'nproc',neighProc(i),i
               trans1(i)=  l1(i)
               trans2(i)=  l2(i)
               trans3(i)=  l3(i)
            else
               di1(i) = -1
               !print *,'di',di1(i),dibeg(i),flowDoms(blocknum,1,1)%diBeg(i)
               di2(i) = -1
               dj1(i) = -1
               dj2(i) = -1
               dk1(i) = -1
               dk2(i) = -1

               neighbourBlock(i) =  -1
               !print *,'nebl',neighbourBlock(i),flowDoms(blocknum,1,1)% neighBlock(i),flowDoms(blocknum,1,1)% neighProc(i)
               neighbourProc(i)=  -1
               !print *,'nproc2',shape(neighProc)
               trans1(i)=  -1
               trans2(i)=  -1
               trans3(i)=  -1
            endif


         enddo

!!$         do i=1,nOne
!!$            
!!$            di1(i) = flowDoms(blocknum,1,1)%diBeg(i)
!!$            di2(i) = flowDoms(blocknum,1,1)% diEnd(i)
!!$            dj1(i) = flowDoms(blocknum,1,1)% djBeg(i)
!!$            dj2(i) = flowDoms(blocknum,1,1)% djEnd(i)
!!$            dk1(i) = flowDoms(blocknum,1,1)% dkBeg(i)
!!$            dk2(i) = flowDoms(blocknum,1,1)% dkEnd(i)
!!$            di1(i) = dibeg(i+4)
!!$            print *,'di',di1(i),dibeg(i),flowDoms(blocknum,1,1)%diBeg(i)
!!$            di2(i) = diEnd(i)
!!$            dj1(i) = djBeg(i)
!!$            dj2(i) = djEnd(i)
!!$            dk1(i) = dkBeg(i)
!!$            dk2(i) = dkEnd(i)
!!$            !print *,'di',shape(dibeg),dibeg,di1,flowDoms(blocknum,1,1)%diBeg(i),i
!!$            !print *,'nebl1',neighbourBlock(i),flowDoms(blocknum,1,1)% neighBlock(i),blocknum
!!$            neighbourBlock(i) = flowDoms(blocknum,1,1)% neighBlock(i)
!!$            !print *,'nebl',neighbourBlock(i),flowDoms(blocknum,1,1)% neighBlock(i),flowDoms(blocknum,1,1)% neighProc(i)
!!$            neighbourProc(i)= flowDoms(blocknum,1,1)% neighProc(i)
!!$            !print *,'nproc2',shape(neighProc)
!!$            trans1(i)= flowDoms(blocknum,1,1)% l1(i)
!!$            trans2(i)= flowDoms(blocknum,1,1)% l2(i)
!!$            trans3(i)= flowDoms(blocknum,1,1)% l3(i)
!!$            !print *,'l2',shape(l2)
!!$         end do
         !print *,'neighbour blocks', shape(neighbourBlock(:)),shape(flowDoms(blocknum,1,1)% neighBlock(:))
        
       end subroutine getBlockCommunicationInfo
         
!!$       subroutine getnDonorNonMatch(blocknum,nNon,nDonor)
!!$         use blockpointers
!!$         implicit none
!!$         
!!$         !Subroutine Variables
!!$         integer(kind=inttype), intent(in)::blocknum,nNon
!!$         integer(kind=intType),dimension(nNon),intent(out):: nDonor
!!$         
!!$         !Local Variables
!!$         integer(kind=intType):: i
!!$         !Begin Execution
!!$         
!!$         call setPointers(blocknum,1,1)
!!$         
!!$         do i=1,nNon
!!$            nDonor(i)=flowDoms(blocknum,1,1)%nonMatchNeigh(i)%nDonorBlocks
!!$         enddo 
!!$         
!!$       end subroutine getnDonorNonMatch

!!$       subroutine getNonMatchingSubfaceInfo(blocknum,nNon,nDonor,blocks,proc,faceid)
!!$         !get donor block info for nonmatching subfaces
!!$         use blockpointers
!!$         implicit none
!!$
!!$         !Subroutine Variables
!!$         integer(kind=intType),intent(in)::blocknum,nNon,nDonor
!!$         integer(kind=intType),dimension(nNon,nDonor),intent(out)::blocks,proc,faceid
!!$         !Local Variables
!!$         
!!$         integer(kind=intType)::i,j
!!$
!!$         call setPointers(blocknum,1,1)
!!$
!!$         !begin execution
!!$         do i = 1,nNon
!!$            do j=1,nDonor
!!$               blocks(i,j)=flowDoms(blocknum,1,1)%nonMatchNeigh(i)%donorBlocks(j)
!!$               proc(i,j)=flowDoms(blocknum,1,1)%nonMatchNeigh(i)%donorProcs(j)
!!$               faceid(i,j)=flowDoms(blocknum,1,1)%nonMatchNeigh(i)%donorFaceIDs(j)
!!$            enddo
!!$         enddo
!!$       
!!$         return
!!$       end subroutine getNonMatchingSubfaceInfo


!!$       subroutine getGlobalNodeNumbering(blocknum,i,j,k,nodeNumbers)
!!$         !return the global indices for the current block
!!$         use blockpointers
!!$         implicit none
!!$
!!$         INTEGER(KIND=INTTYPE), INTENT(IN) :: BLOCKNUM,i,j,k
!!$         REAL(KIND=REALTYPE), DIMENSION(i,j,k), INTENT(OUT) :: nodeNumbers
!!$         
!!$         print *,'blocknum',blocknum
!!$         print *,'ijk',i,j,k
!!$         call setPointers(blocknum,1,1)
!!$         print *,'pointers set'
!!$         print *,'shapes',shape(nodenumbers),shape(flowDoms(blocknum,1,1)%globalNode)
!!$         nodenumbers = flowDoms(blocknum,1,1)%globalNode
!!$
!!$         return
!!$       end subroutine getGlobalNodeNumbering


       subroutine setSingleState(blocknum,i,j,k,l,state)
         
         use blockpointers
         implicit none
         
         !Subroutine Variables
         integer(kind=intType),intent(in)::blocknum,i,j,k,l

         real(kind=realType),intent(in)::state
         
         !Local Variables
         

         !begin Execution

         call setPointersAdj(blocknum,1,1)

         w(i,j,k,l)=state

       end subroutine setSingleState


       subroutine getSingleState(blocknum,i,j,k,l,state)
         
         use blockpointers
         implicit none
         
         !Subroutine Variables
         integer(kind=intType),intent(in)::blocknum,i,j,k,l

         real(kind=realType),intent(out)::state
         
         !Local Variables
         

         !begin Execution

         call setPointersAdj(blocknum,1,1)

         state= w(i,j,k,l)
         
       end subroutine getSingleState
