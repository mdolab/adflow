!
!      ******************************************************************
!      *                                                                *
!      * File:          searchMod.f90                                   *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-12-2005                                      *
!      * Last modified: 08-17-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module searchMod
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains several parameters and variables to *
!      * assist in the donor search process, including the set of       *
!      * possible results returned by the stencil search.               *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

!      ******************************************************************
!      *                                                                *
!      * The definition of the parameters for the possible stencil      *
!      * search results. See the stencilSearch.f90 routine for the      *
!      * the details. Note it is importants that they are less than 0.  *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: Success    = -1_intType
       integer(kind=intType), parameter :: BadDonor   = -2_intType
       integer(kind=intType), parameter :: StopAtBoco = -3_intType
       integer(kind=intType), parameter :: HitMaxIter = -4_intType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! interpolationType, nInterp: self-explanatory - these are stored
       !                             for a given level to make coding
       !                             a bit easier.
       ! ignoreCutoffQuality:        Whether or not to require perfect
       !                             donor quality during a search.

       integer(kind=intType) :: nInterp, interpolationType

       logical :: ignoreCutoffQuality

       end module searchMod
