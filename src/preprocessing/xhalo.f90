!
!      ******************************************************************
!      *                                                                *
!      * File:          xhalo.f90                                       *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader            *
!      * Starting date: 02-23-2003                                      *
!      * Last modified: 08-12-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine xhalo(level)
!
!      ******************************************************************
!      *                                                                *
!      * xhalo determines the coordinates of the nodal halo's.          *
!      * First it sets all halo coordinates by simple extrapolation,    *
!      * then the symmetry planes are treated (also the unit normal of  *
!      * symmetry planes are determined) and finally an exchange is     *
!      * made for the internal halo's.                                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType) :: level
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, sps, i, j, k,ii,jj
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iiMax, jjMax
       !integer(kind=intType) :: iBeg2, iEnd2, jBeg2, jEnd2

       real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2

       real(kind=realType) :: length, dot

       real(kind=realType), dimension(3) :: v1, v2, norm

       !File IO
       integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart,jstart,kstart,iend2,jend2,kend2,n
       integer ::unitx = 12
       !logical ::isopen 

       LOGICAL :: opened=.false., named
       CHARACTER(LEN=80) :: fname
       
       INQUIRE (UNIT=unitx, NAMED=named,  OPENED=opened, NAME=fname)
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and the local
       ! number of blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers to this block.

           call setPointers(nn, level, sps)
           
!!$           if(opened)then
!!$              write(unitx,*) 'block Num',nn
!!$           endif
!
!          **************************************************************
!          *                                                            *
!          * Extrapolation of the coordinates. First extrapolation in   *
!          * i-direction, without halo's, followed by extrapolation in  *
!          * j-direction, with i-halo's and finally extrapolation in    *
!          * k-direction, with both i- and j-halo's. In this way also   *
!          * the indirect halo's get a value, albeit a bit arbitrary.   *
!          *                                                            *
!          **************************************************************
!
           ! Extrapolation in i-direction.

           do k=1,kl
             do j=1,jl
               x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
               x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
               x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)

               x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
               x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
               x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
             enddo
           enddo

           ! Extrapolation in j-direction.

           do k=1,kl
             do i=0,ie
               x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
               x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
               x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)

               x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
               x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
               x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
             enddo
           enddo

           ! Extrapolation in k-direction.

           do j=0,je
             do i=0,ie
               x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
               x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
               x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)

               x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
               x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
               x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
             enddo
           enddo
!
!          **************************************************************
!          *                                                            *
!          * Mirror the halo coordinates adjacent to the symmetry       *
!          * planes                                                     *
!          *                                                            *
!          **************************************************************
!
           ! Loop over boundary subfaces.

           loopBocos: do mm=1,nBocos
              !if(opened) write(unitx,*)'loopBocos',mm,nbocos
             ! The actual correction of the coordinates only takes
             ! place for symmetry planes.

             testSymmetry: if(BCType(mm) == Symm) then
               ! if(opened) write(unitx,*)'testSymmetry',bcfaceID(mm)
               ! Set some variables, depending on the block face on
               ! which the subface is located.
 
               select case (BCFaceID(mm))
                 case (iMin)
                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(0,:,:,:); x1 => x(1,:,:,:); x2 => x(2,:,:,:)

                 case (iMax)
                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(ie,:,:,:); x1 => x(il,:,:,:); x2 => x(nx,:,:,:)

                 case (jMin)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(:,0,:,:); x1 => x(:,1,:,:); x2 => x(:,2,:,:)

                 case (jMax)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(:,je,:,:); x1 => x(:,jl,:,:); x2 => x(:,ny,:,:)

                 case (kMin)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   x0 => x(:,:,0,:); x1 => x(:,:,1,:); x2 => x(:,:,2,:)

                 case (kMax)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   x0 => x(:,:,ke,:); x1 => x(:,:,kl,:); x2 => x(:,:,nz,:)
               end select
 
!!$               if(opened)then
!!$ !                 do nnnn=1,ndom
!!$ !                    call setPointersAdj(nnnn,1,sps)
!!$                     do iii = 2,il
!!$                        do jjj = 2,jl
!!$                           do kkk = 2,kl
!!$                              if(opened)then
!!$                                 write(unitx,*)  'kl',kkk,kl,BCFaceID(mm)
!!$                              endif
!!$                              istart = -3
!!$                              jstart = -3
!!$                              !kstart = -3
!!$                              iend2 = 2
!!$                              jend2 = 2
!!$                              !kend2 = 2 
!!$                              select case (BCFaceID(mm))
!!$                              case (iMin)
!!$                                 if(jjj==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(jjj==jl) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (iMax)
!!$                                 if(jjj==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(jjj==jl) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (jMin)
!!$                                 if(iii==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (jMax)
!!$                                 if(iii==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (kMin)
!!$                                 if(iii==2) istart=-2
!!$                                 if(jjj==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(jjj==jl) jend2=1
!!$                                  
!!$                              case (kMax)
!!$                                 if(iii==2) istart=-2
!!$                                 if(jjj==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(jjj==jl) jend2=1
!!$                                  
!!$                              end select
!!$                             
!!$                              
!!$                              do jjjj = jstart,jend2
!!$                                 do iiii = istart,iend2
!!$                                    !do kkkk = kstart,kend2
!!$                                    do n = 1,3
!!$                                       i = iii+iiii
!!$                                       j = jjj+jjjj
!!$                                       !k = kkk+kkkk
!                                          write(unitx,10) i,j,k,n,nn,v1(1),v1(2),v1(3)
!10                                        format(1x,'res',5I8,3f20.14)
!!$                                       write(unitx,10) i,j,n,iii,jjj,x0(i+1,j+1,n), x1(i+1,j+1,n), x2(i+1,j+1,n)!, dot*norm(n)
!!$10                                     format(1x,'res',5I8,3f20.14)
!!$                                    enddo
!!$                                    !enddo
!!$                                 enddo
!!$                              enddo
!!$                           enddo
!!$                        enddo
!!$                     enddo
!                     call setPointersAdj(nn,1,sps)
!                  enddo
!!$                  endif
!CAM commented out on Aug. 11, 2009 by C.A.Mader. Computation modified to fit 
!CAM into a single residual stencil...
!!$
!!$               ! Determine the sum of all face normals and store this
!!$               ! in norm. The sum of all faces is taken instead of the
!!$               ! cross product of the diagonals, because for some c-type
!!$               ! grids the diagonals of the subface are aligned.
!!$
!!$               norm = zero
!!$
!!$               do j=(jBeg+1),jEnd
!!$                 do i=(iBeg+1),iEnd
!!$
!!$                   ! Determine the vector from the lower left corner to
!!$                   ! the upper right corner. Due to the usage of pointers
!!$                   ! an offset of +1 must be used, because the original
!!$                   ! array x start at 0.
!!$
!!$                   v1(1) = x1(i+1,j+1,1) - x1(i,j,1)
!!$                   v1(2) = x1(i+1,j+1,2) - x1(i,j,2)
!!$                   v1(3) = x1(i+1,j+1,3) - x1(i,j,3)
!!$                   
!!$                   ! And the vector from the upper left corner to the
!!$                   ! lower right corner.
!!$
!!$                   v2(1) = x1(i+1,j,1) - x1(i,j+1,1)
!!$                   v2(2) = x1(i+1,j,2) - x1(i,j+1,2)
!!$                   v2(3) = x1(i+1,j,3) - x1(i,j+1,3)
!!$                   
!!$                   ! Determine the normal of the face by taking the cross
!!$                   ! product of v1 and v2 and add it to norm.
!!$
!!$                   norm(1) = norm(1) + v1(2)*v2(3) - v1(3)*v2(2)
!!$                   norm(2) = norm(2) + v1(3)*v2(1) - v1(1)*v2(3)
!!$                   norm(3) = norm(3) + v1(1)*v2(2) - v1(2)*v2(1)
!!$                   print *,'norm', x1(i+1,j+1,1) ,x1(i,j,1) ,norm(1),norm(2),norm(3),i,j,BCFaceID(mm)
!!$                 enddo
!!$               enddo
!!$
!!$               ! Compute the length of the normal and test if this is
!!$               ! larger than eps. If this is the case this means that
!!$               ! it is a nonsingular subface and the coordinates are
!!$               ! corrected.
!!$
!!$               length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
!!$
!!$               testSingular: if(length > eps) then
!!$
!!$                 ! Compute the unit normal of the subface.
!!$
!!$                 norm(1) = norm(1)/length
!!$                 norm(2) = norm(2)/length
!!$                 norm(3) = norm(3)/length
!!$
!!$                 ! Add an overlap to the symmetry subface if the
!!$                 ! boundaries coincide with the block boundaries.
!!$                 ! This way the indirect halo's are treated properly.
!!$
!!$                 if(iBeg == 1)     iBeg = 0
!!$                 if(iEnd == iiMax) iEnd = iiMax + 1
!!$
!!$                 if(jBeg == 1)     jBeg = 0
!!$                 if(jEnd == jjMax) jEnd = jjMax + 1
!!$
!!$                 ! Loop over the nodes of the subface and set the
!!$                 ! corresponding halo coordinates.
!!$
!!$                 do j=jBeg,jEnd
!!$                   do i=iBeg,iEnd
!!$
!!$                     ! Determine the vector from the internal node to the
!!$                     ! node on the face. Again an offset of +1 must be
!!$                     ! used, due to the usage of pointers.
!!$
!!$                     v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
!!$                     v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
!!$                     v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)
!!$
!!$                     ! Determine two times the normal component of this
!!$                     ! vector; this vector must be added to the
!!$                     ! coordinates of the internal node to obtain the
!!$                     ! halo coordinates. Again the offset of +1.
!!$
!!$                     dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
!!$                         +      v1(3)*norm(3))
!!$
!!$                     x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
!!$                     x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
!!$                     x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)
!!$                     print *,'xhalo', x0(i+1,j+1,1) ,x2(i+1,j+1,1) , dot,norm(1),i,j,BCFaceID(mm)
!!$
!!$                   enddo
!!$                 enddo
!!$
!!$                endif testSingular
               !compute a norm from the 4 corners of the subface. This
               ! will serve to check the subface for singularity and
               !reduce the number of points to be dealt with in the derivatives

               ! Determine the vector from the lower left corner to
               ! the upper right corner. Due to the usage of pointers
               ! an offset of +1 must be used, because the original
               ! array x start at 0.
               
               v1(1) = x1(iimax+1,jjmax+1,1) - x1(1+1,1+1,1)
               v1(2) = x1(iimax+1,jjmax+1,2) - x1(1+1,1+1,2)
               v1(3) = x1(iimax+1,jjmax+1,3) - x1(1+1,1+1,3)
               !print *,'v1',v1,x1(iimax+1,jjmax+1,1), x1(1+1,1+1,1),iimax,jjmax
               ! And the vector from the upper left corner to the
               ! lower right corner.
               
               v2(1) = x1(iimax+1,1+1,1) - x1(1+1,jjmax+1,1)
               v2(2) = x1(iimax+1,1+1,2) - x1(1+1,jjmax+1,2)
               v2(3) = x1(iimax+1,1+1,3) - x1(1+1,jjmax+1,3)
               !print *,'v2',v2
               ! Determine the normal of the face by taking the cross
               ! product of v1 and v2 and add it to norm.
                           
               norm(1) = v1(2)*v2(3) - v1(3)*v2(2)
               norm(2) = v1(3)*v2(1) - v1(1)*v2(3)
               norm(3) = v1(1)*v2(2) - v1(2)*v2(1)
               !print *,'norm',norm

               ! Compute the length of the normal and test if this is
               ! larger than eps. If this is the case this means that
               ! it is a nonsingular subface and the coordinates are
               ! corrected.
               
               length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
               
               testSingular: if(length > eps) then
                  
                  ! Compute the unit normal of the subface.
                  
                  norm(1) = norm(1)/length
                  norm(2) = norm(2)/length
                  norm(3) = norm(3)/length

                  ! Add an overlap to the symmetry subface if the
                  ! boundaries coincide with the block boundaries.
                  ! This way the indirect halo's are treated properly.
                  
                  if(iBeg == 1)     iBeg = 0
                  if(iEnd == iiMax) iEnd = iiMax + 1
                  
                  if(jBeg == 1)     jBeg = 0
                  if(jEnd == jjMax) jEnd = jjMax + 1
                  
                  ! Loop over the nodes of the subface and set the
                  ! corresponding halo coordinates.
                  
                  do j=jBeg,jEnd
                     do i=iBeg,iEnd
                        
                        ! Determine the vector from the internal node to the
                        ! node on the face. Again an offset of +1 must be
                        ! used, due to the usage of pointers.
                        
                        v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
                        v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
                        v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)
                        
                        ! Determine two times the normal component of this
                        ! vector; this vector must be added to the
                        ! coordinates of the internal node to obtain the
                        ! halo coordinates. Again the offset of +1.
                        
                        dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                             +      v1(3)*norm(3))
                        
                        x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
                        x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
                        x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)
                        !print *,'xhalo', x0(i+1,j+1,1) ,x2(i+1,j+1,1) , dot,norm(1),i,j,BCFaceID(mm)

                     enddo
                  enddo
!!$                  print *,'opened',opened,nn
!!$                  if(opened)then
!                  
!                  
!                  do nnnn=1,ndom
!                     call setPointersAdj(nnnn,1,sps)
!!$                     do iii = 2,il
!!$                        do jjj = 2,jl
!!$                           do kkk = 2,kl
!!$                              if(opened)then
!!$                                 write(unitx,*)  'kl',kkk,kl
!!$                              endif
!!$                              istart = -3
!!$                              jstart = -3
!!$                              !kstart = -3
!!$                              iend2 = 2
!!$                              jend2 = 2
!!$                              !kend2 = 2 
!!$                              select case (BCFaceID(mm))
!!$                              case (iMin)
!!$                                 if(jjj==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(jjj==jl) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (iMax)
!!$                                 if(jjj==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(jjj==jl) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (jMin)
!!$                                 if(iii==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (jMax)
!!$                                 if(iii==2) istart=-2
!!$                                 if(kkk==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(kkk==kl) jend2=1 
!!$                              case (kMin)
!!$                                 if(iii==2) istart=-2
!!$                                 if(jjj==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(jjj==jl) jend2=1
!!$                                  
!!$                              case (kMax)
!!$                                 if(iii==2) istart=-2
!!$                                 if(jjj==2) jstart=-2
!!$                                 if(iii==il) iend2=1
!!$                                 if(jjj==jl) jend2=1
!!$                                  
!!$                              end select
!!$                             
!!$                              
!!$                              do jjjj = jstart,jend2
!!$                                 do iiii = istart,iend2
!!$                                    !do kkkk = kstart,kend2
!!$                                    do n = 1,3
!!$                                       i = iii+iiii
!!$                                       j = jjj+jjjj
!!$                                       !k = kkk+kkkk
!                                          write(unitx,10) i,j,k,n,nn,v1(1),v1(2),v1(3)
!10                                        format(1x,'res',5I8,3f20.14)
!!$                                       write(unitx,10) i,j,n,iii,jjj,x0(i+1,j+1,n), x2(i+1,j+1,n)!, dot*norm(n)
!!$10                                     format(1x,'res',5I8,2f20.14)
!!$                                    enddo
!!$                                    !enddo
!!$                                 enddo
!!$                              enddo
!!$                           enddo
!!$                        enddo
!!$                     enddo
!                     call setPointersAdj(nn,1,sps)
!                  enddo
!!$                  endif
               endif testSingular
            endif testSymmetry
         enddo loopBocos
      enddo domains
      
!!$print *,'opened',opened
!!$   
!!$   if(opened)then                  
!!$      do nnnn=1,ndom
!!$         call setPointersAdj(nnnn,1,sps)
!!$         do iii = 2,il
!!$            do jjj = 2,jl
!!$               do kkk = 2,kl
!!$                  istart = -3
!!$                  jstart = -3
!!$                  kstart = -3
!!$                  iend2 = 2
!!$                  jend2 = 2
!!$                  kend2 = 2
!!$                  if(iii==2) istart=-2
!!$                  if(jjj==2) jstart=-2
!!$                  if(kkk==2) kstart=-2
!!$                  if(iii==il) iend2=1
!!$                  if(jjj==jl) jend2=1
!!$                  if(kkk==kl) kend2=1
!!$                  do iiii = istart,iend2
!!$                     do jjjj = jstart,jend2
!!$                        do kkkk = kstart,kend2
!!$                           do n = 1,3
!!$                              i = iii+iiii
!!$                              j = jjj+jjjj
!!$                              k = kkk+kkkk
!!$                              write(unitx,10) i,j,k,n,x(i,j,k,n)
!!$10                            format(1x,'res',4I8,f20.14)
!!$                           enddo
!!$                        enddo
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$            enddo
!!$         enddo
!!$         call setPointersAdj(nn,1,sps)
!!$      enddo
!!$   endif
   enddo spectralLoop 
  
!
!      ******************************************************************
!      *                                                                *
!      * Exchange the coordinates for the internal halo's.              *
!      *                                                                *
!      ******************************************************************
!
   call exchangeCoor(level)
!!$sps=1
!!$print *,'opened',opened
!!$   
!!$   if(opened)then                  
!!$      do nnnn=1,ndom
!!$         call setPointersAdj(nnnn,1,sps)
!!$         do iii = 2,il
!!$            do jjj = 2,jl
!!$               do kkk = 2,kl
!!$                  istart = -3
!!$                  jstart = -3
!!$                  kstart = -3
!!$                  iend2 = 2
!!$                  jend2 = 2
!!$                  kend2 = 2
!!$                  if(iii==2) istart=-2
!!$                  if(jjj==2) jstart=-2
!!$                  if(kkk==2) kstart=-2
!!$                  if(iii==il) iend2=1
!!$                  if(jjj==jl) jend2=1
!!$                  if(kkk==kl) kend2=1
!!$                  do iiii = istart,iend2
!!$                     do jjjj = jstart,jend2
!!$                        do kkkk = kstart,kend2
!!$                           do n = 1,3
!!$                              i = iii+iiii
!!$                              j = jjj+jjjj
!!$                              k = kkk+kkkk
!!$                              write(unitx,10) i,j,k,n,x(i,j,k,n)
!!$10                            format(1x,'res',4I8,f20.14)
!!$                           enddo
!!$                        enddo
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$            enddo
!!$         enddo
!!$         call setPointersAdj(nn,1,sps)
!!$      enddo
!!$   endif   
 end subroutine xhalo
