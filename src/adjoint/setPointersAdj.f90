!
!      ******************************************************************
!      *                                                                *
!      * File:          setPointersAdj.f90                              *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setPointersAdj(nn,mm,ll)
!
!      ******************************************************************
!      *                                                                *
!      * setPointers makes the variables in blockPointers point to      *
!      * block nn for grid level mm and spectral solution ll.           *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: nn, mm, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       globalCell =>flowDoms(nn,mm,ll)%globalCell
       globalNode =>flowDoms(nn,mm,ll)%globalNode
       
       call setPointers(nn,mm,ll)

     end subroutine setPointersAdj
       
