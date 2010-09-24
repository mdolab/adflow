!
!      ******************************************************************
!      *                                                                *
!      * File:          nondimBoundData.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-21-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine nondimBoundData
!
!      ******************************************************************
!      *                                                                *
!      * nondimBoundData nondimensionalizes the boundary data           *
!      * specified in the cgns file.                                    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, sps, nn, mm, i
       real(kind=realType)   :: hRef, uRef

       type(BCDataType), dimension(:), pointer :: BCData
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the number of multigrid levels and the values of hRef
       ! and uRef.

       nLevels = ubound(flowDoms,2)
       hRef    = pRef/rhoRef
       uRef    = sqrt(hRef)

       ! Loop over the number of multigrid levels, spectral solutions
       ! and local blocks.

       do mm=1,nLevels
         do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom

             ! Set the pointer for BCData to make the code more readable.

             BCData => flowDoms(nn,mm,sps)%BCData

             ! Loop over the number of boundary faces.

             do i=1,flowDoms(nn,mm,sps)%nBocos

               ! Nondimensionalize the data if the pointer is associated
               ! with data.

               if( associated(BCData(i)%TNS_Wall) ) &
                 BCData(i)%TNS_Wall = BCData(i)%TNS_Wall/TRef

               if( associated(BCData(i)%ptInlet) ) &
                 BCData(i)%ptInlet = BCData(i)%ptInlet/pRef

               if( associated(BCData(i)%ttInlet) ) &
                 BCData(i)%ttInlet = BCData(i)%ttInlet/TRef

               if( associated(BCData(i)%htInlet) ) &
                 BCData(i)%htInlet = BCData(i)%htInlet/HRef

               if( associated(BCData(i)%turbInlet) ) &
                 call nondimTurb(BCData(i)%turbInlet)

               if( associated(BCData(i)%rho) ) &
                 BCData(i)%rho = BCData(i)%rho/rhoRef

               if( associated(BCData(i)%velx) ) &
                 BCData(i)%velx = BCData(i)%velx/uRef

               if( associated(BCData(i)%vely) ) &
                 BCData(i)%vely = BCData(i)%vely/uRef

               if( associated(BCData(i)%velz) ) &
                 BCData(i)%velz = BCData(i)%velz/uRef

               if( associated(BCData(i)%ps) ) &
                 BCData(i)%ps = BCData(i)%ps/pRef

             enddo
           enddo
         enddo
       enddo

       !=================================================================

       contains

         !===============================================================

         subroutine nondimTurb(turbInlet)
!
!        ****************************************************************
!        *                                                              *
!        * NondimTurb nondimensionalizes the turbulent data for inlet  *
!        * boundary conditions.                                         *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(:,:,:), pointer :: turbInlet
!
!        Local variables.
!
         integer(kind=intType) :: nn
         real(kind=realType)   :: nuRef, tmp

         real(kind=realType), dimension(nt1:nt2) :: ref
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Set the reference values depending on the turbulence model.
 
         nuRef = muRef/rhoRef
         select case (turbModel)

           case (spalartAllmaras, spalartAllmarasEdwards)
             ref(itu1) = nuRef

           case (komegaWilcox, komegaModified, menterSST)
             ref(itu1) = pRef/rhoRef
             ref(itu2) = ref(itu1)/nuRef

           case (ktau)
             ref(itu1) = pRef/rhoRef
             ref(itu2) = nuRef/ref(itu1)

           case (v2f)
             ref(itu1) = pRef/rhoRef
             ref(itu4) = ref(itu1)/nuRef
             ref(itu2) = ref(itu1)*ref(itu4)
             ref(itu3) = ref(itu1)

         end select

         ! Loop over the number of turbulence variables and make
         ! them nondimensional.

         do nn=nt1,nt2
           tmp = one/ref(nn)
           turbInlet(:,:,nn) = turbInlet(:,:,nn)*tmp
         enddo

         end subroutine nondimTurb

       end subroutine nondimBoundData
