!        ================================================================

         subroutine riemannFluxAdj(left,right,flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
!
!        ****************************************************************
!        *                                                              *
!        * riemannFlux computes the flux for the given face and left    *
!        * and right states.                                            *
!        *                                                              *
!        ****************************************************************
!
         use precision
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
         real(kind=realType), dimension(*) :: flux
         logical :: correctForK
!
!        Local variables.
!
         real(kind=realType) :: porFlux, rFace
         real(kind=realType) :: Etl, Etr, z1l, z1r, tmp
         real(kind=realType) :: dr, dru, drv, drw, drE, drk
         real(kind=realType) :: rAvg, uAvg, vAvg, wAvg, hAvg, kAvg
         real(kind=realType) :: alphaAvg, a2Avg, aAvg, unAvg
         real(kind=realType) :: ovaAvg, ova2Avg, area, eta
         real(kind=realType) :: gm1, gm53, ovgm1, factK
         real(kind=realType) :: lam1, lam2, lam3
         real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Set the porosity for the flux. The default value, 0.5*, is
         ! a scaling factor.

         porFlux = half
         if(por == noFlux) porFlux = zero

         ! Abbreviate some expressions in which gamma occurs.

         gm1   =  gammaFace - one
         ovgm1 =  one/gm1
         gm53  =  gammaFace - five*third
         factK = -ovgm1*gm53

         ! Determine which riemann solver must be solved.


         ! No preconditioner used. Use the Roe scheme of the
         ! standard equations.

         ! Compute the total energy. Calorically perfect gas
         ! for now.

         Etl = ovgm1*left(irhoE)                            &
             + half*left(irho)*(left(ivx)**2 + left(ivy)**2 &
             +                  left(ivz)**2)
         Etr = ovgm1*right(irhoE)                              &
             + half*right(irho)*(right(ivx)**2 + right(ivy)**2 &
             +                  right(ivz)**2)

         ! Compute the square root of the left and right densities
         ! and the inverse of the sum.

         z1l = sqrt(left(irho))
         z1r = sqrt(right(irho))
         tmp = one/(z1l + z1r)

         ! Compute some variables depending whether or not a
         ! k-equation is present.

         if( correctForK ) then

            ! Store the difference of the turbulent kinetic energy
            ! per unit volume, i.e. the conserved variable.
            
            drk = right(irho)*right(itu1) - left(irho)*left(itu1)
            
            ! Compute the average turbulent energy per unit mass
            ! using Roe averages.

            kAvg = tmp*(z1l*left(itu1) + z1r*right(itu1))
            
            ! Correct the total energy for the present of the
            ! turbulent kinetic energy.
            
            Etl = Etl - factK* left(irho)* left(itu1)
            Etr = Etr - factK*right(irho)*right(itu1)
            
         else
            
            ! Set the difference of the turbulent kinetic energy
            ! per unit volume and the averaged kinetic energy per
            ! unit mass to zero.
            
            drk  = 0.0
            kAvg = 0.0
            
         endif
         
         ! Compute the difference of the conservative mean
         ! flow variables.

         dr  = right(irho) - left(irho)
         dru = right(irho)*right(ivx) - left(irho)*left(ivx)
         drv = right(irho)*right(ivy) - left(irho)*left(ivy)
         drw = right(irho)*right(ivz) - left(irho)*left(ivz)
         drE = Etr - Etl
         
         ! Compute the Roe average variables, which can be
         ! computed directly from the average Roe vector.
         
         rAvg = fourth*(z1r + z1l)**2
         uAvg = tmp*(z1l*left(ivx) + z1r*right(ivx))
         vAvg = tmp*(z1l*left(ivy) + z1r*right(ivy))
         wAvg = tmp*(z1l*left(ivz) + z1r*right(ivz))
         hAvg = tmp*((Etl+left(irhoE)) /z1l &
              +      (Etr+right(irhoE))/z1r)

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
         abv7 = abv2*abv4*ovaAvg  + abv3*abv5

         ! Compute the dissipation term, -|a| (wr - wl), which is
         ! multiplied by porFlux. Note that porFlux is either
         ! 0.0 or 0.5.

         flux(irho)  = -porFlux*(lam3*dr  + abv6)
         flux(imx)   = -porFlux*(lam3*dru + uAvg*abv6 &
                                          + sx*abv7)
         flux(imy)   = -porFlux*(lam3*drv + vAvg*abv6 &
                                          + sy*abv7)
         flux(imz)   = -porFlux*(lam3*drw + wAvg*abv6 &
                                          + sz*abv7)
         flux(irhoE) = -porFlux*(lam3*drE + hAvg*abv6 &
                                          + unAvg*abv7)
         
       end subroutine riemannFluxAdj
