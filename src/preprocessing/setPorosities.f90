!
!      ******************************************************************
!      *                                                                *
!      * File:          setPorosities.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-06-2003                                      *
!      * Last modified: 11-06-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setPorosities(level)
!
!      ******************************************************************
!      *                                                                *
!      * setPorosities sets the porosities for the faces to a certain   *
!      * flag. Default is normalFlux. The two other possibilities are   *
!      * boundFlux, used for solid wall boundaries, and noFlux for a    *
!      * conservative treatment of non matching block boundaries. In    *
!      * the latter case the flux is constructed differently and the    *
!      * flux computation in the block must be neglected.               *
!      * Note that only the 1st spectral solution is treated, because   *
!      * this informations is the same for all of them.                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use constants
       use inputDiscretization
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i, j, k

       integer(kind=intType), dimension(2) :: ri, rj, rk

       integer(kind=porType) :: por
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of domains.

       domains: do nn=1,nDom

         ! Store the number of nodes in this block a bit easier.

         il = flowDoms(nn,level,1)%il
         jl = flowDoms(nn,level,1)%jl
         kl = flowDoms(nn,level,1)%kl

         ! Allocate the memory for the porosities.

         allocate(flowDoms(nn,level,1)%porI(1:il,2:jl,2:kl), &
                  flowDoms(nn,level,1)%porJ(2:il,1:jl,2:kl), &
                  flowDoms(nn,level,1)%porK(2:il,2:jl,1:kl), stat=ierr)
         if(ierr /= 0)                     &
           call terminate("setPorosities", &
                          "Memory allocation failure for porosities")

         ! Set the pointers for this block to make the source
         ! more readable.

         call setPointers(nn, level, 1_intType)

         ! Initialize the porosities to normalFlux.

         porI = normalFlux
         porJ = normalFlux
         porK = normalFlux

         ! Loop over the subfaces to alter the porosities.

         subface: do mm=1,nsubface

           ! Set the porosity for this subface or continue with the
           ! next if the porosity should not be changed.

           if(BCType(mm) == NSWallAdiabatic  .or. &
              BCType(mm) == NSWallIsothermal .or. &
              BCType(mm) == EulerWall        .or. &
              BCType(mm) == Extrap)    then
             por = boundFlux
           else if(BCType(mm)        == B2BMismatch .and. &
                   nonMatchTreatment == Conservative) then
             por = noFlux
           else
             cycle
           endif

           ! Set the range for the faces on this subface.

           ri(1) = min(inBeg(mm), inEnd(mm)) +1
           ri(2) = max(inBeg(mm), inEnd(mm))

           rj(1) = min(jnBeg(mm), jnEnd(mm)) +1
           rj(2) = max(jnBeg(mm), jnEnd(mm))

           rk(1) = min(knBeg(mm), knEnd(mm)) +1
           rk(2) = max(knBeg(mm), knEnd(mm))

           ! Determine the block face this subface is located on and
           ! set the corresponding porosities correctly.

           select case( BCFaceID(mm) )

             case (iMin)
               do k=rk(1),rk(2)
                 do j=rj(1),rj(2)
                   porI(1,j,k) = por
                 enddo
               enddo

             !===========================================================

             case (iMax)
               do k=rk(1),rk(2)
                 do j=rj(1),rj(2)
                   porI(il,j,k) = por
                 enddo
               enddo

             !===========================================================

             case (jMin)
               do k=rk(1),rk(2)
                 do i=ri(1),ri(2)
                   porJ(i,1,k) = por
                 enddo
               enddo

             !===========================================================

             case (jMax)
               do k=rk(1),rk(2)
                 do i=ri(1),ri(2)
                   porJ(i,jl,k) = por
                 enddo
               enddo

             !===========================================================

             case (kMin)
               do j=rj(1),rj(2)
                 do i=ri(1),ri(2)
                   porK(i,j,1) = por
                 enddo
               enddo

             !===========================================================

             case (kMax)
               do j=rj(1),rj(2)
                 do i=ri(1),ri(2)
                   porK(i,j,kl) = por
                 enddo
               enddo

           end select

         enddo subface

       enddo domains

       end subroutine setPorosities
