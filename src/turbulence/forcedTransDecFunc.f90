!
!      ******************************************************************
!      *                                                                *
!      * File:          forcedTransDecFunc                              *
!      * Author:        Eran Arad                                       *
!      * Starting date: 08-04-2009                                      *
!      * Last modified:                                                 *
!      *                                                                *
!      ******************************************************************
!
       subroutine  forcedTransDecFunc (Xs,Xe,pp,transC,Xpoint,ft2,modCb1)
!
!      ******************************************************************
!      *                                                                *
!      * Computes the decay function of the SA model in both RANS and   *
!      * DES modes. First 5 arguments are input. Last 2: output         *
!      *                                                                *
!      ******************************************************************
!
       use paramTurb
       use inputPhysics
       use constants

       implicit none

!      subroutine arguments.
!
       real(kind=realType), intent(in) :: Xs,Xe,transC,Xpoint
       integer(kind=intType), intent(in) :: pp
       real(kind=realType), intent(out) :: ft2,modCb1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
       if (TransHLength <= zero)then
!
!-- ft2 modification transition
!
          if(Xpoint < xTransition)then

             ft2     = one
                      
          end if

       else
!
! ---Cb1 modification transition (Roy-Blottner)
!                                          
          ft2            = zero

          if(Xpoint <= Xe)then

             modCb1 = zero

          else if(Xpoint >= Xe)then

             modCb1 = rsaCb1
		
          else

             modCb1 = rsaCb1*transC*(Xpoint-Xs)**pp 

          end if

       end if

     end subroutine forcedTransDecFunc
