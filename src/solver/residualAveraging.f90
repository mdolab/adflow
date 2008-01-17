!
!      ******************************************************************
!      *                                                                *
!      * File:          residualAveraging.f90                           *
!      * Author:        Juan J. Alonso, Steve Repsher (blanking)        *
!      * Starting date: 09-17-2004                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine residualAveraging
!
!      ******************************************************************
!      *                                                                *
!      * Implicit residual smoothing is a simple procedure that         *
!      * replaces the residual at each point by a weighted sum of all   *
!      * of the residuals in the block (although the residuals that are *
!      * closer to the cell under consideration are weighted more       *
!      * heavily). This smoothing can be applied explicitly, but may    *
!      * result in zero smoothed residuals for non-zero initial         *
!      * residual modes.  For this reason, the smoothing is applied     *
!      * implicitly in the following form:                              *
!      *                                                                *
!      * -epz R{i+1} + (1 + 2 epz) R{i} -epz R{i-1} = r{i}              *
!      *                                                                *
!      * Where r{i} is the original residual at point i, and R{i} is    *
!      * the implicitly smoothed residual at point i.  The analysis for *
!      * the 1-D scalar convection-diffusion equation shows that if     *
!      *                                                                *
!      * Epz >= (1/4)*((lambda/lambda*)^2 -1),                          *
!      *                                                                *
!      * where lambda is the cfl number desired to be used, and         *
!      * lambda* is the CFL limit of the unsmoothed scheme, the scheme  *
!      * can be made unconditionally stable (arbitrarily large lambda). *
!      * In practice, lambda = 6-8 is common for Euler solutions.  For  *
!      * RANS solutions lambda = 3-4 is what we have used in practice   *
!      * resulting in a slight improvement in convergence rate, but,    *
!      * more importantly, an increase in the robustness of the solver. *
!      *                                                                *
!      * Note that this theoretical result can be shown for infinite    *
!      * 1-D problems, but also for finite-periodic 1-D problems and    *
!      * finite-size 1-D problems (i.e. with block boundaries).  Such   *
!      * theoretical results are not available for 3-D cases, where the *
!      * scheme is applied in an ADI fashion:                           *
!      *                                                                *
!      * (1 -epzI d_ii)(1 -epzJ d_jj)(1 -epzK d_kk) r = r               *
!      *                                                                *
!      * Where d_ii, d_jj, d_kk are second difference operators in each *
!      * of the mesh coordinate directions.                             *
!      *                                                                *
!      * For each of the coordinate direction solves, the initial       *
!      * matrix problem is of the form:                                 *
!      *                                                                *
!      *   -                                             - - -   - -    *
!      *   | (1+2 epz)   -epz                            | |r|   |r|    *
!      *   |   -epz    (1 +2 epz)    -epz                | |r|   |r|    *
!      *   |             -epz       (1 + 2 epz)   -epz   | |r| = |r|    *
!      *   |                 .           .           .   | |.|   |.|    *
!      *   |                   .            .          . | |.|   |.|    *
!      *   -                                             - - -   - -    *
!      *                                                                *
!      * And after the forward elimination phase a normalization is     *
!      * applied so the result looks like:                              *
!      *                                                                *
!      *   -                   - - -   - -                              *
!      *   |  1  -d            | |r|   |r|                              *
!      *   |  0   1  -d        | |r|   |r|                              *
!      *   |      0   1  -d    | |r| = |r|                              *
!      *   |       .   .   .   | |.|   |.|                              *
!      *   |          .  .   . | |.|   |.|                              *
!      *   -                   - - -   - -                              *
!      *                                                                *
!      * Which can then be used with a straightforward backsolve to     *
!      * obtain the answer.                                             *
!      *                                                                *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputIteration
       use flowVarRefState
       use iteration

       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l

       real(kind=realType), parameter :: b = 2.0_realType

       real(kind=realType) :: cflim, currentCfl, rfl0, plim
       real(kind=realType) :: dpi, dpj, dpk, r
       real(kind=realType), dimension(il,max(jl,kl)) :: epz, d, t, rfl

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      rfl0 is a measure of the ratio lambda/lambda*
!
!      Hardwire cflim for the time being (=lambda*).  It should be
!       defaulted to something else later on.
!
       cflim = 3.0

       if(currentLevel <= groundLevel) then
         currentCfl = cfl
       else
         currentCfl = cflCoarse
       endif

       rfl0  = half*currentCfl/cflim

       plim  = 0.001_realType*pInfCorr

!      ******************************************************************
!      *                                                                *
!      * Smoothing in the i-direction. Only done when enough cells are  *
!      * present in the i-direction.                                    *
!      *                                                                *
!      ******************************************************************
!
       if(nx > 1) then

         do k=2,kl

           ! Compute smoothing coeficients

           do j=2,jl
             do i=2,il
               dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                   /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
               dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                   /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
               dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                   /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
               rfl(i,j) = one/(one + b*(dpi  +dpj  +dpk))
             end do
           end do

           do j=2,jl
             do i=2,nx
                r = rfl0*(rfl(i,j) +rfl(i+1,j))
                epz(i,j) = fourth*smoop*dim(r*r,one)*iblank(i,j,k)
             end do
           end do

           ! Zero out coefficients for boundary condition treatment

           do j=2,jl
             epz(1,j)  = zero
             epz(il,j) = zero
             d(1,j)    = zero
           end do

           ! Compute coefficients for forward elimination process

           do i=2,il
             do j=2,jl
               t(i,j) = one &
                      / (one +epz(i,j) +epz(i-1,j) -epz(i-1,j)*d(i-1,j))
               d(i,j) = t(i,j)*epz(i,j)
             end do
           end do

           ! Apply same transformation to the rhs vector of residuals

           do i=2,il
             do j=2,jl
               do l=1,nMGVar
                 dw(i,j,k,l) = t(i,j) &
                             * (dw(i,j,k,l) +epz(i-1,j)*dw(i-1,j,k,l))
               end do
             end do
           end do

           ! Backsolve operation.  Smoothed residuals are left
           !  in dw(i,j,k,l)

           do i=nx,2,-1
             do j=2,jl
               do l=1,nMGVar
                 dw(i,j,k,l) = dw(i,j,k,l) +d(i,j)*dw(i+1,j,k,l)
               end do
             end do
           end do

         enddo
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Smoothing in the j-direction. Only done when enough cells are  *
!      * present in the j-direction.                                    *
!      *                                                                *
!      ******************************************************************
!
       if(ny > 1) then

         do k=2,kl

           ! Compute smoothing coeficients

           do j=2,jl
             do i=2,il
               dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                   /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
               dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                   /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
               dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                   /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
               rfl(i,j) = one/(one + b*(dpi  +dpj  +dpk))
             end do
           end do

           do j=2,ny
             do i=2,il
               r = rfl0*(rfl(i,j) +rfl(i,j+1))
               epz(i,j) = fourth*smoop*dim(r*r,one)*iblank(i,j,k)
             end do
           end do

           ! Zero out coefficients for boundary condition treatment

           do i=2,il
             epz(i,1)  = zero
             epz(i,jl) = zero
             d(i,1)    = zero
           end do

           ! Compute coefficients for forward eliMination process

           do j=2,jl
             do i=2,il
               t(i,j) = one &
                      / (one +epz(i,j) +epz(i,j-1) -epz(i,j-1)*d(i,j-1))
               d(i,j) = t(i,j)*epz(i,j)
             end do
           end do

           ! Apply same transformation to the rhs vector of residuals

           do j=2,jl
             do i=2,il
               do l=1,nMGVar
                 dw(i,j,k,l) = t(i,j) &
                             * (dw(i,j,k,l) +epz(i,j-1)*dw(i,j-1,k,l))
               end do
             end do
           end do

           ! Backsolve operation.  Smoothed residuals are left
           !  in dw(i,j,k,l)

           do j=ny,2,-1
             do i=2,il
               do l=1,nMGVar
                 dw(i,j,k,l) = dw(i,j,k,l) +d(i,j)*dw(i,j+1,k,l)
               end do
             end do
           end do

         enddo
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Smoothing in the k-direction. Only done when enough cells are  *
!      * present in the k-direction.                                    *
!      *                                                                *
!      ******************************************************************
!
       if(nz > 1) then

         do j=2,jl

           ! Compute smoothing coeficients

           do k=2,kl
             do i=2,il
               dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                   /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
               dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                   /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
               dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                   /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
               rfl(i,k) = one/(one + b*(dpi  +dpj  +dpk))
             end do
           end do

           do k=2,nz
             do i=2,il
               r = rfl0*(rfl(i,k) +rfl(i,k+1))
               epz(i,k) = fourth*smoop*dim(r*r,one)*iblank(i,j,k)
             end do
           end do

           ! Zero out coefficients for boundary condition treatment

           do i=2,il
             epz(i,1)  = zero
             epz(i,kl) = zero
             d(i,1)    = zero
           end do

           ! Compute coefficients for forward eliMination process

           do k=2,kl
             do i=2,il
               t(i,k) = one &
                      / (one +epz(i,k) +epz(i,k-1) -epz(i,k-1)*d(i,k-1))
               d(i,k) = t(i,k)*epz(i,k)
             end do
           end do

           ! Apply same transformation to the rhs vector of residuals

           do k=2,kl
             do i=2,il
               do l=1,nMGVar
                 dw(i,j,k,l) = t(i,k) &
                             * (dw(i,j,k,l) +epz(i,k-1)*dw(i,j,k-1,l))
               end do
             end do
           end do

           ! Backsolve operation.  Smoothed residuals are left
           !  in dw(i,j,k,l)

           do k=nz,2,-1
             do i=2,il
               do l=1,nMGVar
                 dw(i,j,k,l) = dw(i,j,k,l) +d(i,k)*dw(i,j,k+1,l)
               end do
             end do
           end do

         enddo
       endif

       end subroutine residualAveraging
