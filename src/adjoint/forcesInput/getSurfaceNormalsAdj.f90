!
!      ******************************************************************
!      *                                                                *
!      * File:          getSurfaceNormalsAdj.f90                        *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi                                    *
!      * Starting date: 12-18-2007                                      *
!      * Last modified: 12-18-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine getSurfaceNormalsAdj(xAdj,siAdj,sjAdj,skAdj,normAdj, &
           iiBeg,iiEnd,jjBeg,jjEnd,mm,level,nn,sps,righthanded)
!
!
!      ******************************************************************
!      *                                                                *
!      * boundarySurfaceNormals computes the outward normals of the     *
!      * boundary subfaces at the given time instance and the given     *
!      * multi grid level.                                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockpointers, only: ie,je,ke,il,jl,kl,zero,one,half,bcfaceid
       use communication
!       use section
       implicit none
!
!      Subroutine arguments.
!

       integer(kind=intType), intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
       integer(kind=intType), intent(in) :: level,nn,sps,mm

       real(kind=realType), dimension(0:ie,0:je,0:ke,3),intent(in) :: xAdj

!       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3), intent(out) :: siAdj
!       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3), intent(out) :: sjAdj
!       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3), intent(out) :: skAdj
!       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3), intent(out) :: normAdj

       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3) :: siAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3) :: sjAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3) :: skAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: normAdj

!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: i,j,k,ii,jj,kk,l,m,n
       real(kind=realType)   :: mult,xp,yp,zp,fact
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: ss
       real(kind=realType), dimension(3)   :: v1, v2!,v12,v22
       
       logical ,intent(in)::righthanded
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************

!
!      **************************************************************
!      *                                                            *
!      * Computation of the face normals in i-, j- and k-direction. *
!      * Formula's are valid for a right handed block; for a left   *
!      * handed block the correct orientation is obtained via fact. *
!      * The normals point in the direction of increasing index.    *
!      * The absolute value of fact is 0.5, because the cross       *
!      * product of the two diagonals is twice the normal vector.   *
!      *                                                            *
!      * Note that also the normals of the first level halo cells   *
!      * are computed. These are needed for the viscous fluxes.     *
!      *                                                            *
!      **************************************************************
!

      ! Set the factor in the surface normals computation. For a
      ! left handed block this factor is negative, such that the
      ! normals still point in the direction of increasing index.
      ! The formulae used later on assume a right handed block
      ! and fact is used to correct this for a left handed block,
      ! as well as the scaling factor of 0.5

       if( rightHanded ) then
          fact =  half
       else
          fact = -half
       endif

       select case (BCFaceID(mm))

       case(iMin)
          do k=jjBeg,jjEnd
             n = k -1
             do j=iiBeg,iiEnd
                m = j -1
                do i=1,1
   
                   ! to get siAdj(1,:,:,3)
                   if(i==1)  ii=1
                   
                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,n,1) - xAdj(i,m,k,1)
                   v1(2) = xAdj(i,j,n,2) - xAdj(i,m,k,2)
                   v1(3) = xAdj(i,j,n,3) - xAdj(i,m,k,3)
                   
                   v2(1) = xAdj(i,j,k,1) - xAdj(i,m,n,1)
                   v2(2) = xAdj(i,j,k,2) - xAdj(i,m,n,2)
                   v2(3) = xAdj(i,j,k,3) - xAdj(i,m,n,3)
!!$
!!$                   v12(1) = x(i,j,n,1) - x(i,m,k,1)
!!$                   v12(2) = x(i,j,n,2) - x(i,m,k,2)
!!$                   v12(3) = x(i,j,n,3) - x(i,m,k,3)
!!$                   
!!$                   v22(1) = x(i,j,k,1) - x(i,m,n,1)
!!$                   v22(2) = x(i,j,k,2) - x(i,m,n,2)
!!$                   v22(3) = x(i,j,k,3) - x(i,m,n,3)

                   
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   siAdj(ii,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   siAdj(ii,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   siAdj(ii,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                   
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', ii,j,k,siAdj(ii,j,k,1)-si(i,j,k,1),siAdj(ii,j,k,2)-si(i,j,k,2),siAdj(ii,j,k,3)-si(i,j,k,3)                
                enddo
             enddo
          enddo
          

       case(iMax)
          do k=jjBeg,jjEnd
             n = k -1
             do j=iiBeg,iiEnd
                m = j -1
                do i=il,il
                   
                   ! to get siAdj(2,:,:,3)
                   if(i==il) ii=2 
                   
                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,n,1) - xAdj(i,m,k,1)
                   v1(2) = xAdj(i,j,n,2) - xAdj(i,m,k,2)
                   v1(3) = xAdj(i,j,n,3) - xAdj(i,m,k,3)
                   
                   v2(1) = xAdj(i,j,k,1) - xAdj(i,m,n,1)
                   v2(2) = xAdj(i,j,k,2) - xAdj(i,m,n,2)
                   v2(3) = xAdj(i,j,k,3) - xAdj(i,m,n,3)
                   
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   siAdj(ii,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   siAdj(ii,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   siAdj(ii,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                   
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', ii,j,k,siAdj(ii,j,k,1)-si(i,j,k,1),siAdj(ii,j,k,2)-si(i,j,k,2),siAdj(ii,j,k,3)-si(i,j,k,3)                
                enddo
             enddo
          enddo
          
       case(jMin)
          ! Projected areas of cell faces in the j direction.
          do k=jjBeg,jjEnd
             n = k -1
             do j=1,1
                ! to get sjAdj(:,1,:,:)
                if(j==1)  jj=1
                do i=iiBeg,iiEnd
                   l = i -1

                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,n,1) - xAdj(l,j,k,1)
                   v1(2) = xAdj(i,j,n,2) - xAdj(l,j,k,2)
                   v1(3) = xAdj(i,j,n,3) - xAdj(l,j,k,3)
                   
                   v2(1) = xAdj(l,j,n,1) - xAdj(i,j,k,1)
                   v2(2) = xAdj(l,j,n,2) - xAdj(i,j,k,2)
                   v2(3) = xAdj(l,j,n,3) - xAdj(i,j,k,3)
                   
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   sjAdj(i,jj,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   sjAdj(i,jj,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   sjAdj(i,jj,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
               
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', &
                   !i,jj,k,sjAdj(i,jj,k,1)-sj(i,j,k,1),sjAdj(i,jj,k,2)-sj(i,j,k,2),sjAdj(i,jj,k,3)-sj(i,j,k,3)                 
                enddo
             enddo
          enddo

       case(jMax)
       ! Projected areas of cell faces in the j direction.
          do k=jjBeg,jjEnd
             n = k -1
             do j=jl,jl
                ! to get sjAdj(:,1:2:,:,3)
                if(j==jl) jj=2 

                do i=iiBeg,iiEnd
                   l = i -1
                   
                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,n,1) - xAdj(l,j,k,1)
                   v1(2) = xAdj(i,j,n,2) - xAdj(l,j,k,2)
                   v1(3) = xAdj(i,j,n,3) - xAdj(l,j,k,3)
                   
                   v2(1) = xAdj(l,j,n,1) - xAdj(i,j,k,1)
                   v2(2) = xAdj(l,j,n,2) - xAdj(i,j,k,2)
                   v2(3) = xAdj(l,j,n,3) - xAdj(i,j,k,3)
                   
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   sjAdj(i,jj,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   sjAdj(i,jj,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   sjAdj(i,jj,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', &
                   !i,jj,k,sjAdj(i,jj,k,1)-sj(i,j,k,1),sjAdj(i,jj,k,2)-sj(i,j,k,2),sjAdj(i,jj,k,3)-sj(i,j,k,3)                 
                enddo
             enddo
          enddo
          


       case(kMin)
          ! Projected areas of cell faces in the k direction.
          do k=1,1
             ! to get skAdj(:,::,1:2,3)
             if(k==1)  kk=1
             do j=jjBeg,jjEnd
                m = j -1
                do i=iiBeg,iiEnd
                   l = i -1
                   
                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,k,1) - xAdj(l,m,k,1)
                   v1(2) = xAdj(i,j,k,2) - xAdj(l,m,k,2)
                   v1(3) = xAdj(i,j,k,3) - xAdj(l,m,k,3)
                   
                   v2(1) = xAdj(l,j,k,1) - xAdj(i,m,k,1)
                   v2(2) = xAdj(l,j,k,2) - xAdj(i,m,k,2)
                   v2(3) = xAdj(l,j,k,3) - xAdj(i,m,k,3)
                   !print *,'kmin vectors',v1,'v2',v2,'fact',fact
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   skAdj(i,j,kk,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   skAdj(i,j,kk,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   skAdj(i,j,kk,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', &
                   !i,j,kk,skAdj(i,j,kk,1)-sk(i,j,k,1),skAdj(i,j,kk,2)-sk(i,j,k,2),skAdj(i,j,kk,3)-sk(i,j,k,3)                
                enddo
             enddo
          enddo



       case(kMax)
          ! Projected areas of cell faces in the k direction.
          do k=kl,kl
             ! to get skAdj(:,::,1:2,3)
             if(k==kl) kk=2 
             
             do j=jjBeg,jjEnd
                m = j -1
                do i=iiBeg,iiEnd
                   l = i -1
                   
                   ! Determine the two diagonal vectors of the face.
                   
                   v1(1) = xAdj(i,j,k,1) - xAdj(l,m,k,1)
                   v1(2) = xAdj(i,j,k,2) - xAdj(l,m,k,2)
                   v1(3) = xAdj(i,j,k,3) - xAdj(l,m,k,3)
                   
                   v2(1) = xAdj(l,j,k,1) - xAdj(i,m,k,1)
                   v2(2) = xAdj(l,j,k,2) - xAdj(i,m,k,2)
                   v2(3) = xAdj(l,j,k,3) - xAdj(i,m,k,3)
                   
                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.
                   
                   skAdj(i,j,kk,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   skAdj(i,j,kk,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   skAdj(i,j,kk,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
                   !if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi ', &
                   !i,j,kk,skAdj(i,j,kk,1)-sk(i,j,k,1),skAdj(i,j,kk,2)-sk(i,j,k,2),skAdj(i,j,kk,3)-sk(i,j,k,3)                
                enddo
             enddo
          enddo
          
       end select

       ! Determine the block face on which this subface is located
       ! and set ss and mult accordingly.
       
       select case (BCFaceID(mm))
          
       case (iMin)
          mult = -one; ss(:,:,:) = siAdj(1,:,:,:)
          !print *,'imin',sum(ss)
       case (iMax)
          mult = one;  ss(:,:,:) = siAdj(2,:,:,:) ! which was si(il,:,:,:)
          !print *,'imax',sum(ss)
       case (jMin)
          mult = -one; ss(:,:,:) = sjAdj(:,1,:,:)
          !print *,'jmin',sum(ss)
       case (jMax)
          mult = one;  ss(:,:,:) = sjAdj(:,2,:,:) ! which was sj(:,jl,:,:)
          !print *,'jmax',sum(ss)
       case (kMin)
          mult = -one; ss(:,:,:) = skAdj(:,:,1,:)
          !print *,'kmin',sum(ss)
       case (kMax)
          mult = one;  ss(:,:,:) = skAdj(:,:,2,:) ! which was sk(:,:,kl,:)
          !print *,'kmax',sum(ss)
       end select
       

       
!       do k=1,ke
!          do j=1,je
!             do i=1,2

                ! to get siAdj(1:2,:,:,3)
!                if(i==1)  ii=1
!                if(i==il) ii=2
                
                
!if(myID==0.and.sps==1) write(*,'(a5,3i5,3e)')'hi2 ',ii,j,k,siAdj(ii,j,k,1),siAdj(ii,j,k,2),siAdj(ii,j,k,3)                
!             enddo
!          enddo
!       enddo


       ! Loop over the boundary faces of the subface.
       
       do j=jjBeg,jjEnd
          do i=iiBeg,iiEnd
             
             ! Compute the inverse of the length of the normal vector
             ! and possibly correct for inward pointing.

             
             xp = ss(i,j,1);  yp = ss(i,j,2);  zp = ss(i,j,3)
!!$             fact = sqrt(xp*xp + yp*yp + zp*zp)
!!$             if(fact > zero) fact = mult/fact
!!$             
!!$             ! Compute the unit normal.
!!$             
!!$             normAdj(i,j,1) = fact*xp
!!$             normAdj(i,j,2) = fact*yp
!!$             normAdj(i,j,3) = fact*zp
             
             !alternate form to allow inclusion of degenrate halos???
             if( xp**2>zero .or. yp**2>zero .or. zp**2>zero)then
                !if (fact > zero)then
                !compute length
                fact = sqrt(xp*xp + yp*yp + zp*zp)
                !set factor to 1/length
                fact = mult/fact
                !compute unit normal...
                normAdj(i,j,1) = fact*xp
                normAdj(i,j,2) = fact*yp
                normAdj(i,j,3) = fact*zp
             else
                !Length is zero
                normAdj(i,j,:) = zero
             endif
             !print *,'normInd',normAdj,i,j,fact,xp,yp,zp
          enddo
       enddo
       
       
     end subroutine getSurfaceNormalsAdj
