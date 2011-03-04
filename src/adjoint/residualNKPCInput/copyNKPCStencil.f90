
subroutine copyNKPCStencil(iCell, jCell, kCell, nn, level, sps, wAdj, &
     siAdj, sjAdj, skAdj, sAdj, sfaceIAdj, sfaceJAdj, sfaceKAdj, rotRateAdj,&
     voladj)

  use blockPointers  
  use flowVarRefState  
  use inputPhysics
  use inputTimeSpectral 
  use cgnsgrid   
  use inputMotion    
  implicit none

  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: iCell, jCell, kCell, nn, level, sps
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,&
       nTimeIntervalsSpectral), intent(out) :: wAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,&
       nTimeIntervalsSpectral),intent(out) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,&
       nTimeIntervalsSpectral),intent(out) :: sAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,&
       nTimeIntervalsSpectral),intent(out) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType),dimension(3),intent(out) :: rotRateAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral),intent(out) :: volAdj

  !     Local variables.
  !
  integer(kind=intType) :: sps2,i,j,k,l,ii,jj,kk
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  ! 
  ! Copy the wAdj from w
  wAdj = 0.0
  siAdj = 0.0
  sjAdj = 0.0
  skAdj = 0.0
  sfaceiadj = 0.0
  sfacejadj = 0.0
  sfacekadj = 0.0

  do sps2 = 1,nTimeIntervalsSpectral
     do l=1,nw
        do kk=-2,2
           do jj=-2,2
              do ii=-2,2
                 wAdj(ii,jj,kk,l,sps2) = &
                      flowdoms(nn,level,sps2)%w(iCell+ii, jCell+jj, kCell+kk,l)
              enddo
           enddo
        enddo
     enddo
     voladj(sps2) = flowdoms(nn,level,sps2)%vol(icell,jcell,kcell)


     ! Projected areas of cell faces in the i direction.

     kStart=-2; kEnd=2
     jStart=-2; jEnd=2
     iStart=-3; iEnd=2

     if(iCell==il) iEnd=1
     if(iCell==2) iStart=-2

     if(jCell==2)  jStart=-1
     if(jCell==jl) jEnd=1 

     if(kCell==2) kStart=-1
     if(kCell==kl) kEnd=1

     do kk=kStart,kEnd !-2,2
        do jj=jStart,jEnd !-2,2
           do ii=iStart,iEnd !-2,2
              siAdj(ii,jj,kk,:,sps2) = flowdoms(nn,level,sps2)% &
                   si(icell+ii,jcell+jj,kcell+kk,:)
           end do
        end do
     end do

     ! Projected areas of cell faces in the j direction.
     kStart=-2; kEnd=2
     jStart=-3; jEnd=2
     iStart=-2; iEnd=2

     if(iCell==2)  iStart=-1
     if(iCell==il) iEnd=1
     
     if(jCell==jl) jEnd=1 
     if(jCell==2) jStart=-2
     
     if(kCell==2) kStart=-1
     if(kCell==kl) kEnd=1

     do kk=kStart,kEnd !-2,2
        do jj=jStart,jEnd !-2,2
           do ii=iStart,iEnd !-2,2
              sjAdj(ii,jj,kk,:,sps2) = flowdoms(nn,level,sps2)% &
                   sj(icell+ii,jcell+jj,kcell+kk,:)
           end do
        end do
     end do

     ! Projected areas of cell faces in the k direction.
           
     kStart=-3; kEnd=2
     jStart=-2; jEnd=2
     iStart=-2; iEnd=2
     
     if(iCell==2)  iStart=-1
     if(iCell==il) iEnd=1
     
     if(jCell==2)  jStart=-1
     if(jCell==jl) jEnd=1 
     
     if(kCell==kl) kEnd=1
     if(kCell==2)  kStart=-2
     
     do kk=kStart,kEnd
        do jj=jStart,jEnd 
           do ii=iStart,iEnd 
              skAdj(ii,jj,kk,:,sps2) = flowdoms(nn,level,sps2)% &
                   sk(icell+ii,jcell+jj,kcell+kk,:)
           end do
        end do
     end do

     kStart=-2; kEnd=2
     jStart=-2; jEnd=2
     iStart=-2; iEnd=2
     
     if(iCell==2)  iStart=-1
     if(iCell==il) iEnd=1
     
     if(jCell==2)  jStart=-1
     if(jCell==jl) jEnd=1 
     
     if(kCell==kl) kEnd=1
     if(kCell==2)  kStart=-1
     testMoving: if( blockIsMoving ) then
        do kk=kstart,kend
           do jj=jstart,jend
              do ii=istart,iend
                 sadj(ii,jj,kk,:,sps2) = flowdoms(nn,level,sps)% &
                      s(icell+ii,jcell+jj,kcell+kk,:)
              end do
           end do
        end do
        
        do kk=kstart,kend
           do jj=jstart,jend
              do ii=istart,iend
                 sfaceiadj(ii,jj,kk,sps2) = flowdoms(nn,level,sps2)% &
                      sfacei(icell+ii,jcell+jj,kcell+kk)
                 sfacejadj(ii,jj,kk,sps2) = flowdoms(nn,level,sps2)% &
                      sfacej(icell+ii,jcell+jj,kcell+kk)
                 sfacekadj(ii,jj,kk,sps2) = flowdoms(nn,level,sps2)% &
                      sfacek(icell+ii,jcell+jj,kcell+kk)
              end do
           end do
        end do
     else
        sadj = 0.0
        sfaceiadj = 0.0
        sfacejadj = 0.0
        sfacekadj = 0.0
     end if testMoving
  end do

  j = nbkGlobal
  rotRateAdj   = timeRef*cgnsDoms(j)%rotRate

end subroutine copyNKPCStencil



