!
!      ******************************************************************
!      *                                                                *
!      * File:          setBCTurb.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBCVarNamesTurb(offset)
!
!      ******************************************************************
!      *                                                                *
!      * setBCVarNamesTurb sets the names for the turbulence            *
!      * variables to be determined. This depends on the turbulence     *
!      * model. If not the RANS equations are solved an immediate       *
!      * return is made.                                                *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use inputPhysics
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: offset
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if not the RANS equations are solved.

       if(equations /= RANSEquations) return

       ! Determine the turbulence model and set the names accordingly.

       select case (turbModel)
         case (spalartAllmaras, spalartAllmarasEdwards)
           bcVarNames(offset+1) = cgnsTurbSaNu

         case (komegaWilcox, komegaModified, menterSST)
           bcVarNames(offset+1) = cgnsTurbK
           bcVarNames(offset+2) = cgnsTurbOmega

         case (ktau)
           bcVarNames(offset+1) = cgnsTurbK
           bcVarNames(offset+2) = cgnsTurbTau

         case (v2f)
           bcVarNames(offset+1) = cgnsTurbK
           bcVarNames(offset+2) = cgnsTurbEpsilon
           bcVarNames(offset+3) = cgnsTurbV2
           bcVarNames(offset+4) = cgnsTurbF

       end select

       end subroutine setBCVarNamesTurb

       !=================================================================

       logical function setBCVarTurb(offset, boco, turbInlet)
!
!      ******************************************************************
!      *                                                                *
!      * SetBCVarTurb sets the array for the turbulent halo data        *
!      * for inlet boundaries. This function returns .true. If all      *
!      * turbulence variables could be interpolated and .false.         *
!      * otherwise.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: offset, boco
       real(kind=realType), dimension(:,:,:), pointer :: turbInlet
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j
       real(kind=realType)   :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize setBCVarTurb to .true. And return immediately
       ! if not the rans equations are solved.

       setBCVarTurb = .true.
       if(equations /= RANSEquations) return

       ! Loop over the number of turbulent variables. mm is the counter
       ! in the arrays bcVarArray and bcVarPresent.

       mm = offset
       turbLoop: do nn=nt1,nt2
         mm = mm + 1

         ! Check if the variable is present. If so, use the
         ! interpolated data.

         if( bcVarPresent(mm) ) then

           ! Conversion to SI units if possible.

           call siTurb(mass(mm), length(mm), time(mm), temp(mm), &
                       bcVarNames(mm), mult, trans)

           ! Set the turbulent variables.

           do j=jBeg,jEnd
             do i=iBeg,iEnd
                turbInlet(i,j,nn) = mult*bcVarArray(i,j,mm) + trans
             enddo
           enddo

         else

           ! Turbulent variable not present.
           ! Set setBCVarTurb to .false.

           setBCVarTurb = .false.

         endif
       enddo turbLoop

       ! If not all the turbulence variables could be interpolated
       ! reset them to zero and store this face in the array
       ! turbFreestreamSubfaces, which stores these subfaces.
       ! They will be set to free stream values later on.

       if(.not. setBCVarTurb) then

         call storeTurbFreestreamSubface

         do nn=nt1,nt2
           do j=jBeg,jEnd
             do i=iBeg,iEnd
                turbInlet(i,j,nn) = zero
             enddo
           enddo
         enddo
       endif

       !=================================================================

       contains

         !===============================================================

         subroutine storeTurbFreestreamSubface
!
!        ****************************************************************
!        *                                                              *
!        * storeTurbFreestreamSubface stores the currently active       *
!        * subface in the array turbFreestreamSubfaces, such that the   *
!        * turbulence variables can be set to the free stream values    *
!        * later on in setInletFreestreamTurb.                          *
!        *                                                              *
!        ****************************************************************
!
         use blockPointers
         implicit none
!
!        Local variables.
!
         integer :: ierr

         integer(kind=intType) :: nn

         integer(kind=intType), dimension(:,:), allocatable :: tmp
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the situation we are dealing with here.

         testAllocated: if( allocated(turbFreestreamSubfaces) ) then

           ! TurbFreestreamSubfaces has already been allocated and
           ! thus contains information. It must be reallocated and the
           ! current subface should be added.

           ! Allocate the memory for tmp and copy the data from
           ! turbFreestreamSubfaces.

           allocate(tmp(nTurbFreestreamSubfaces,3), stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("storeTurbFreestreamSubface", &
                            "Memory allocation failure for tmp")

           tmp = turbFreestreamSubfaces

           ! Release turbFreestreamSubfaces and allocate it again
           ! with an increased dimension.

           deallocate(turbFreestreamSubfaces, stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("storeTurbFreestreamSubface", &
                            "Deallocation failure for &
                            &turbFreestreamSubfaces")

           nTurbFreestreamSubfaces = nTurbFreestreamSubfaces + 1

           allocate(turbFreestreamSubfaces(nTurbFreestreamSubfaces,3), &
                    stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("storeTurbFreestreamSubface", &
                            "Memory allocation failure for &
                            &turbFreestreamSubfaces")

           ! Copy the data back from tmp into turbFreestreamSubfaces
           ! and release tmp again.

           do nn=1,(nTurbFreestreamSubfaces-1)
             turbFreestreamSubfaces(nn,1) = tmp(nn,1)
             turbFreestreamSubfaces(nn,2) = tmp(nn,2)
             turbFreestreamSubfaces(nn,3) = tmp(nn,3)
           enddo

           deallocate(tmp, stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("storeTurbFreestreamSubface", &
                            "Deallocation failure for tmp")

         else testAllocated

           ! TurbFreestreamSubfaces has not been allocated yet. This
           ! is the first subface to store in this array. Allocate the
           ! array and set nTurbFreestreamSubfaces to 1.

           nTurbFreestreamSubfaces = 1
           allocate(turbFreestreamSubfaces(nTurbFreestreamSubfaces,3), &
                    stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("storeTurbFreestreamSubface", &
                            "Memory allocation failure for &
                            &turbFreestreamSubfaces")

         endif testAllocated

         ! Store the current subface in turbFreestreamSubfaces.

         nn = nTurbFreestreamSubfaces
         turbFreestreamSubfaces(nn,1) = nbkLocal
         turbFreestreamSubfaces(nn,2) = boco
         turbFreestreamSubfaces(nn,3) = spectralSol

         end subroutine storeTurbFreestreamSubface

       end function setBCVarTurb
