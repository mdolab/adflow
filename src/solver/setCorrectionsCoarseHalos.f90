!
!      ******************************************************************
!      *                                                                *
!      * File:          setCorrectionsCoarseHalos.f90                   *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 05-14-2003                                      *
!      * Last modified: 08-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setCorrectionsCoarseHalos(sps, nn, coarseLevel, &
                                            fact, nVarInt)
!
!      ******************************************************************
!      *                                                                *
!      * setCorrectionsCoarseHalos sets the values of the coarse        *
!      * grid boundary halo corrections. For all boundaries, either a   *
!      * homogeneous Dirichlet condition (fact = 0.0) or a Neumann      *
!      * condition (fact = 1.0) is used. Exception are symmetry planes, *
!      * where a mirroring takes place.                                 *
!      *                                                                *
!      ******************************************************************
!
       use block
       use BCTypes
       use constants
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps, nn, coarseLevel
       integer(kind=intType), intent(in) :: nVarInt
       real(kind=realType), intent(in)   :: fact
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, mm
       integer(kind=intType) :: il, jl, kl, ie, je, ke

       real(kind=realType) :: nnx, nny, nnz, vn

       real(kind=realType), dimension(:,:,:,:), pointer :: ww
       real(kind=realType), dimension(:,:,:),   pointer :: ww1, ww2

       type(BCDataType), dimension(:), pointer :: BCData
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointer ww to the coarse grid variables. At the moment
       ! when this routine is called, these contain the corrections in a
       ! normal grid cycle and the true variables when used to
       ! interpolate the solution to the next finer mesh level in full
       ! multigrid mode. Also set the pointer for BCData, such that the
       ! unit normals are accessed easier.

       ww     => flowDoms(nn,coarseLevel,sps)%w
       BCData => flowDoms(nn,coarseLevel,sps)%BCData

       ! Easier storage of the upper coarse block range.

       il = flowDoms(nn,coarseLevel,sps)%il
       jl = flowDoms(nn,coarseLevel,sps)%jl
       kl = flowDoms(nn,coarseLevel,sps)%kl

       ie = flowDoms(nn,coarseLevel,sps)%ie
       je = flowDoms(nn,coarseLevel,sps)%je
       ke = flowDoms(nn,coarseLevel,sps)%ke

       ! Loop over the number of boundary subfaces to correct the
       ! boundary halo values.

       subfacesCoarse: do mm=1,flowDoms(nn,coarseLevel,sps)%nBocos

         ! Set the pointers for ww1 and ww2, depending on the block face
         ! on which this subface is located. Note that bcFaceID is the
         ! same for all spectral modes and only the 1st is allocated.
         ! Therefore the value of the 1st spectral mode is used here.

         select case (flowDoms(nn,coarseLevel,1)%BCFaceID(mm))

           case (iMin)
             ww1 => ww(1, 1:,1:,:); ww2 => ww(2, 1:,1:,:)
           case (iMax)
             ww1 => ww(ie,1:,1:,:); ww2 => ww(il,1:,1:,:)
           case (jMin)
             ww1 => ww(1:,1 ,1:,:); ww2 => ww(1:,2 ,1:,:)
           case (jMax)
             ww1 => ww(1:,je,1:,:); ww2 => ww(1:,jl,1:,:)
           case (kMin)
             ww1 => ww(1:,1:,1 ,:); ww2 => ww(1:,1:,2 ,:)
           case (kMax)
             ww1 => ww(1:,1:,ke,:); ww2 => ww(1:,1:,kl,:)

         end select

         ! Choose what to do based on the BC type. Note that BCType is
         ! the same for all spectral modes and only the 1st is allocated.
         ! Therefore the value of the 1st spectral mode is used here.

         select case (flowDoms(nn,coarseLevel,1)%BCType(mm))

           case (SlidingInterface,   OversetOuterBound,     &
                 DomainInterfaceAll, DomainInterfaceRhoUVW, &
                 DomainInterfaceP,   DomainInterfaceRho,    &
                 DomainInterfaceTotal)

             ! None of these are physical boundary conditions and thus
             ! nothing needs to be done. The halos already contain the
             ! corrections or solution.

           case (Symm)

             ! This is a symmetry plane. Loop over the faces of this 
             ! subface.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
               do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ! Compute twice the normal velocity component of the
                 ! internal cell.

                 nnx = BCData(mm)%norm(i,j,1)
                 nny = BCData(mm)%norm(i,j,2)
                 nnz = BCData(mm)%norm(i,j,3)

                 vn = two*(ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny &
                    +      ww2(i,j,ivz)*nnz)

                 ! Compute the halo state. Make sure that the average
                 ! normal velocity component is zero.

                 ww1(i,j,irho)  = ww2(i,j,irho)
                 ww1(i,j,ivx)   = ww2(i,j,ivx) - vn*nnx
                 ww1(i,j,ivy)   = ww2(i,j,ivy) - vn*nny
                 ww1(i,j,ivz)   = ww2(i,j,ivz) - vn*nnz
                 ww1(i,j,irhoE) = ww2(i,j,irhoE)

                 do l=nt1,nVarInt
                   ww1(i,j,l) = ww2(i,j,l)
                 enddo

               enddo
             enddo

           case default

             ! Other type of boundary subface. Set the boundary halo's
             ! as fact times the internal cell value.

             do l=1,nVarInt
               do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                   ww1(i,j,l) = fact*ww2(i,j,l)
                 enddo
               enddo
             enddo

         end select
       enddo subfacesCoarse

       end subroutine setCorrectionsCoarseHalos
