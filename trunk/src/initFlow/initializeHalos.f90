!
!      ******************************************************************
!      *                                                                *
!      * File:          initializeHalos.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-13-2007                                      *
!      * Last modified: 09-14-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initializeHalos(halosRead)
!
!      ******************************************************************
!      *                                                                *
!      * initializeHalos sets the flow variables in the halo cells      *
!      * using a constant extrapolation. If the halos are read only the *
!      * second halos are initialized, otherwise both.                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: halosRead
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j, k, l
       integer(kind=intType) :: jj, kk
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and blocks.

       spectralLoop: do mm=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn,mgStartlevel,mm)

           ! Determine the situation we are dealing with.

           testHalosRead: if( halosRead ) then

             ! The first layer of halo cells have been read. Initialize
             ! the second layer from the first layer.

             ! Halo cells in the i-direction.

             do k=0,kb
               kk = max(1_intType,min(k,ke))
               do j=0,jb
                 jj = max(1_intType,min(j,je))

                 do l=1,nw
                   w(0, j,k,l) = w(1, jj,kk,l)
                   w(ib,j,k,l) = w(ie,jj,kk,l)
                 enddo

                 p(0, j,k) = p(1, jj,kk)
                 p(ib,j,k) = p(ie,jj,kk)

               enddo
             enddo

             ! Halo cells in j-direction. Note that the i-halo's have
             ! already been set.

             do k=0,kb
               kk = max(1_intType,min(k,ke))
               do i=1,ie

                 do l=1,nw
                   w(i,0, k,l) = w(i,1, kk,l)
                   w(i,jb,k,l) = w(i,je,kk,l)
                 enddo

                 p(i,0 ,k) = p(i,1, kk)
                 p(i,jb,k) = p(i,je,kk)

               enddo
             enddo

             ! Halo cells in k-direction. Note that the halo's in both
             ! i and j direction have already been set.

             do j=1,je
               do i=1,ie

                 do l=1,nw
                   w(i,j,0, l) = w(i,j,1, l)
                   w(i,j,kb,l) = w(i,j,ke,l)
                 enddo

                 p(i,j,0)  = p(i,j,1)
                 p(i,j,kb) = p(i,j,ke)

               enddo
             enddo

           else testHalosRead

             ! No halo cells have been read. Initialize both layers
             ! using the internal value.

             ! Halo cells in the i-direction.

             do k=0,kb
               kk = max(2_intType,min(k,kl))
               do j=0,jb
                 jj = max(2_intType,min(j,jl))

                 do l=1,nw
                   w(0, j,k,l) = w(2, jj,kk,l)
                   w(1, j,k,l) = w(2, jj,kk,l)
                   w(ie,j,k,l) = w(il,jj,kk,l)
                   w(ib,j,k,l) = w(il,jj,kk,l)
                 enddo

                 p(0, j,k) = p(2, jj,kk)
                 p(1, j,k) = p(2, jj,kk)
                 p(ie,j,k) = p(il,jj,kk)
                 p(ib,j,k) = p(il,jj,kk)

               enddo
             enddo

             ! Halo cells in j-direction. Note that the i-halo's have
             ! already been set.

             do k=0,kb
               kk = max(2_intType,min(k,kl))
               do i=2,il

                 do l=1,nw
                   w(i,0, k,l) = w(i,2, kk,l)
                   w(i,1, k,l) = w(i,2, kk,l)
                   w(i,je,k,l) = w(i,jl,kk,l)
                   w(i,jb,k,l) = w(i,jl,kk,l)
                 enddo

                 p(i,0 ,k) = p(i,2, kk)
                 p(i,1 ,k) = p(i,2, kk)
                 p(i,je,k) = p(i,jl,kk)
                 p(i,jb,k) = p(i,jl,kk)

               enddo
             enddo

             ! Halo cells in k-direction. Note that the halo's in both
             ! i and j direction have already been set.

             do j=2,jl
               do i=2,il

                 do l=1,nw
                   w(i,j,0, l) = w(i,j,2, l)
                   w(i,j,1, l) = w(i,j,2, l)
                   w(i,j,ke,l) = w(i,j,kl,l)
                   w(i,j,kb,l) = w(i,j,kl,l)
                 enddo

                 p(i,j,0)  = p(i,j,2)
                 p(i,j,1)  = p(i,j,2)
                 p(i,j,ke) = p(i,j,kl)
                 p(i,j,kb) = p(i,j,kl)

               enddo
             enddo

           endif testHalosRead

           ! Initialize the laminar and eddy viscosity, if appropriate,
           ! such that no uninitialized memory is present.
           ! As the viscosities are dependent variables their values
           ! are not read during a restart.

           if( viscous )   rlv = muInf
           if( eddyModel ) rev = eddyVisInfRatio*muInf

         enddo domains
       enddo spectralLoop

       end subroutine initializeHalos
