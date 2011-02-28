
!
!      ******************************************************************
!      *                                                                *
!      * File:          riemannFluxAdj.f90                              *
!      * Author:        Edwin van der Weide, C.A.(Sandy) Mader          *
!      * Starting date: 04-25-2008                                      *
!      * Last modified: 04-25-2008                                      *
!      *                                                                *
!      ******************************************************************
!


         subroutine riemannFluxNKPC(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)
!
!        ****************************************************************
!        *                                                              *
!        * riemannFlux computes the flux for the given face and left    *
!        * and right states.                                            *
!        *                                                              *
!        ****************************************************************
!
         use iteration
         use constants
         use inputDiscretization

         implicit none
!
!        Subroutine arguments.
!
         integer(kind=porType) :: por
         real(kind=realType) :: sx, sy, sz, gammaFace
         real(kind=realType) :: sFace
         real(kind=realType), dimension(*), intent(in)  :: left, right
         real(kind=realType), dimension(*), intent(out) :: flux
         logical, intent(in) :: fineGrid,correctfork
!
!        Local variables.
!
         real(kind=realType) :: porFlux, rFace
         real(kind=realType) :: Etl, Etr, z1l, z1r, tmp
         real(kind=realType) :: dr, dru, drv, drw, drE, drk
         real(kind=realType) :: rAvg, uAvg, vAvg, wAvg, hAvg, kAvg
         real(kind=realType) :: alphaAvg, a2Avg, aAvg, unAvg
         real(kind=realType) :: ovaAvg, ova2Avg, area, Eta
         real(kind=realType) :: gm1, gm53
         real(kind=realType) :: lam1, lam2, lam3
         real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7

         integer(kind=intType) :: limUsed, riemannUsed


         real(kind=realType), dimension(2) :: rhotmp, utmp, vtmp, wtmp
         real(kind=realType), dimension(2) :: ptmp, ktmp, Etmp
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Set the porosity for the flux. The default value, 0.5*rFil, is
         ! a scaling factor where an rFil != 1 is taken into account.
         !print *,'adjflux2*****************************',half,rfil
         !stop
         porFlux = half*rFil
         if(por == noFlux .or. por == boundFlux)    porFlux = zero

         ! Determine the limiter scheme to be used. On the fine grid the
         ! user specified scheme is used; on the coarse grid a first order
         ! scheme is computed.
         !print *,'selecting limiter',limiter,limused,firstorder,finegrid
         limUsed = firstOrder
         if( fineGrid ) limUsed = limiter
         
         ! Determine the riemann solver which must be used.
         !print *,'selecting riemann',riemann,riemannused,riemannCoarse,finegrid
         riemannUsed = riemannCoarse
         if( fineGrid ) riemannUsed = riemann
         
         ! Abbreviate some expressions in which gamma occurs.
         !print *,'some constantts',gammaface
         gm1  = gammaFace - one
         gm53 = gammaFace - five*third

         ! Determine which riemann solver must be solved.

         select case (riemannUsed)

           case (Roe)
              !print *,'using roe solver'
             ! Determine the preconditioner used.

             select case (precond)

               case (noPrecond)

                 ! No preconditioner used. Use the Roe scheme of the
                 ! standard equations.

                 ! Compute the square root of the left and right densities
                 ! and the inverse of the sum.

                 z1l = sqrt(left(irho))
                 z1r = sqrt(right(irho))
                 tmp = one/(z1l + z1r)
                 !print *,'tmp',one,z1l , z1r

                 ! Compute some variables depending whether or not a
                 ! k-equation is present.

                 if( correctForK ) then
                    !print *,'correcting for k'
                   ! Store the left and right kinetic energy in ktmp,
                   ! which is needed to compute the total energy.

                   ktmp(1) = left(itu1)
                   ktmp(2) = right(itu1)

                   ! Store the difference of the turbulent kinetic energy
                   ! per unit volume, i.e. the conserved variable.

                   drk = right(irho)*right(itu1) - left(irho)*left(itu1)

                   ! Compute the average turbulent energy per unit mass
                   ! using Roe averages.

                   kAvg = tmp*(z1l*left(itu1) + z1r*right(itu1))

                 else

                    !print *,'not correcting for k'
                   ! Set the difference of the turbulent kinetic energy
                   ! per unit volume and the averaged kinetic energy per
                   ! unit mass to zero.

                   drk  = 0.0
                   kAvg = 0.0

                 endif

                 ! Compute the total energy of the left and right state.
                 !print *,'calculating temp vaars'
                 rhotmp(1) = left(irho)
                 rhotmp(2) = right(irho)

                 utmp(1) = left(ivx)
                 utmp(2) = right(ivx)

                 vtmp(1) = left(ivy)
                 vtmp(2) = right(ivy)

                 wtmp(1) = left(ivz)
                 wtmp(2) = right(ivz)

                 ptmp(1) = left(irhoE)
                 ptmp(2) = right(irhoE)

!                 call etotArrayAdj(rhotmp, utmp, vtmp, wtmp, ptmp, ktmp, &
!                                Etmp, correctForK, 2_intType)
                 call etotArrayNKPC(rhotmp, utmp, vtmp, wtmp, ptmp, ktmp, &
                                Etmp, correctForK, 2)

                 Etl = Etmp(1)
                 Etr = Etmp(2)

                 ! Compute the difference of the conservative mean
                 ! flow variables.

                 dr  = right(irho) - left(irho)
                 dru = right(irho)*right(ivx) - left(irho)*left(ivx)
                 drv = right(irho)*right(ivy) - left(irho)*left(ivy)
                 drw = right(irho)*right(ivz) - left(irho)*left(ivz)
                 drE = Etr - Etl

                 ! Compute the Roe average variables, which can be
                 ! computed directly from the average Roe vector.

                 !rAvg = fourth*(z1r + z1l)**2 deadend code!
                 uAvg = tmp*(z1l*left(ivx) + z1r*right(ivx))
                 vAvg = tmp*(z1l*left(ivy) + z1r*right(ivy))
                 wAvg = tmp*(z1l*left(ivz) + z1r*right(ivz))
                 hAvg = tmp*((Etl+left(irhoE)) /z1l &
                      +      (Etr+right(irhoE))/z1r)
                 !print *,'uavg',tmp,z1l,left(ivx),z1r,right(ivx)
                 ! Compute the unit vector and store the area of the
                 ! normal. Also compute the unit normal velocity of the face.

                 area  = sqrt(sx**2 + sy**2 + sz**2)
                 tmp   = one/max(1.e-25_realType,area)
                 sx    = sx*tmp
                 sy    = sy*tmp
                 sz    = sz*tmp
                 rFace = sFace*tmp

                 ! Compute some dependent variables at the Roe
                 ! average state.

                 alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
                 a2Avg    = abs(gm1*(hAvg - alphaAvg) - gm53*kAvg)
                 aAvg     = sqrt(a2Avg)
                 unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
                 !print *,'unavg',uAvg,sx,vAvg,sy,wAvg,sz

                 ovaAvg  = one/aAvg
                 ova2Avg = one/a2Avg

                 ! Set for a boundary the normal velocity to rFace, the
                 ! normal velocity of the boundary.

                 if(por == boundFlux) unAvg = rFace

                 ! Compute the coefficient eta for the entropy correction.
                 ! At the moment a 1D entropy correction is used, which
                 ! removes expansion shocks. Although it also reduces the
                 ! carbuncle phenomenon, it does not remove it completely.
                 ! In other to do that a multi-dimensional entropy fix is
                 ! needed, see Sanders et. al, JCP, vol. 145, 1998,
                 ! pp. 511 - 537. Although relatively easy to implement,
                 ! an efficient implementation requires the storage of
                 ! all the left and right states, which is rather
                 ! expensive in terms of memory.

                 eta = half*(abs((left(ivx) - right(ivx))*sx        &
                     +           (left(ivy) - right(ivy))*sy        &
                     +           (left(ivz) - right(ivz))*sz)       &
                     +       abs(sqrt(gammaFace*left(irhoE)/left(irho)) &
                     -           sqrt(gammaFace*right(irhoE)/right(irho))))

                 ! Compute the absolute values of the three eigenvalues.

                 lam1 = abs(unAvg - rFace + aAvg)
                 lam2 = abs(unAvg - rFace - aAvg)
                 lam3 = abs(unAvg - rFace)
                 
                 !print *,'lam3',unAvg, rFace
                 ! Apply the entropy correction to the eigenvalues.

                 tmp = two*eta
                 if(lam1 < tmp) lam1 = eta + fourth*lam1*lam1/eta
                 if(lam2 < tmp) lam2 = eta + fourth*lam2*lam2/eta
                 if(lam3 < tmp) lam3 = eta + fourth*lam3*lam3/eta

                 ! Multiply the eigenvalues by the area to obtain
                 ! the correct values for the dissipation term.

                 lam1 = lam1*area
                 lam2 = lam2*area
                 lam3 = lam3*area

                 ! Some abbreviations, which occur quite often in the
                 ! dissipation terms.

                 abv1 = half*(lam1 + lam2)
                 abv2 = half*(lam1 - lam2)
                 abv3 = abv1 - lam3

                 abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                      -      wAvg*drw + drE) - gm53*drk
                 abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

                 abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
                 !print *,'abv6',abv3,abv4,ova2Avg,abv2,abv5,ovaAvg
                 abv7 = abv2*abv4*ovaAvg  + abv3*abv5

                 ! Compute the dissipation term, -|a| (wr - wl), which is
                 ! multiplied by porFlux. Note that porFlux is either
                 ! 0.0 or 0.5*rFil.

                 !print *,'inriemannend',-porFlux,lam3,dr, abv6

                 flux(irho)  = -porFlux*(lam3*dr  + abv6)
                 flux(imx)   = -porFlux*(lam3*dru + uAvg*abv6 &
                                                  + sx*abv7)
                 flux(imy)   = -porFlux*(lam3*drv + vAvg*abv6 &
                                                  + sy*abv7)
                 flux(imz)   = -porFlux*(lam3*drw + wAvg*abv6 &
                                                  + sz*abv7)
                 flux(irhoE) = -porFlux*(lam3*drE + hAvg*abv6 &
                                                  + unAvg*abv7)

      !          tmp = max(lam1,lam2,lam3)

      !          flux(irho)  = -porFlux*(tmp*dr)
      !          flux(imx)   = -porFlux*(tmp*dru)
      !          flux(imy)   = -porFlux*(tmp*drv)
      !          flux(imz)   = -porFlux*(tmp*drw)
      !          flux(irhoE) = -porFlux*(tmp*drE)

               case (Turkel)
                 call terminate("riemannFlux", "Turkel preconditioner not implemented yet")

               case (ChoiMerkle)
                 call terminate("riemannFlux", "choi merkle preconditioner not implemented yet")

             end select

           case (vanLeer)
             call terminate("riemannFlux", "van leer fvs not implemented yet")

           case (ausmdv)
             call terminate("riemannFlux", "ausmdv fvs not implemented yet")

         end select
         !print *,'riemanflux complete'
       end subroutine riemannFluxNKPC
