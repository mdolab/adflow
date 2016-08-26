!
!       Finds quadratic uvw weights [-1:1] for interpolant point       
!       xSearch by solving for its co-ord in donor element defined by  
!       xElem using newton-raphson iteration.                          
!       Fractions uvwQuadratic come with initial guess from ADT search 
!       used for linear interpolation.                                 

subroutine computeQuadraticWeights(xSearch,xElem,uvwQuadratic)

  use precision
  use constants

  implicit none

  ! Input variables
  real(kind=realType), intent(in)    :: xSearch(3), xElem(3, -1:1, -1:1, -1:1)
  real(kind=realType), intent(inout) :: uvwQuadratic(3)
 
  ! Working variables
  ! -----------------------------------------------------------------
  ! newton related
  integer(kind=intType) :: n, niter
  real(kind=realType)   :: resid,residtol
  real(kind=realType)   :: B(3,4)
  
  ! others
  integer(kind=intType) :: j, l, iii, jjj, kkk
  real(kind=realType)   :: ff(27), shp(3,3), psi(3)
  real(kind=realType)   :: dff(27,3),dshp(3,3,3) !-> differentials wrt 3 psi dirs
  logical               :: ok_flag
  ! -----------------------------------------------------------------

  ! Initialize newton parameters
  niter = 0 
  resid = 1.0e10
  residtol = 1.0e-15
  

  ! Begin newton iterations to solve for uvw weights wrt quadratic element


  newton_loop: do while (niter < 10 .and. resid > residtol)

     !step 1: find weights
     !--------------------

     ! Initialize weights 
     psi(1:3) = uvwQuadratic(1:3)
     
     ! Precopute the FE shape functions for each j-direction
     do j=1,3
        shp(1, j) = half*psi(j)*(psi(j) - one)
        shp(2, j) = -(psi(j)**2-1)
        shp(3, j) = half*psi(j)*(psi(j) + one)
     end do

     ! These are the 27 quadratic weights
     ff(1 )   = shp(1, 1)*shp(1, 2)*shp(1, 3)
     ff(2 )   = shp(2, 1)*shp(1, 2)*shp(1, 3)
     ff(3 )   = shp(3, 1)*shp(1, 2)*shp(1, 3)
     
     ff(4 )   = shp(1, 1)*shp(2, 2)*shp(1, 3)
     ff(5 )   = shp(2, 1)*shp(2, 2)*shp(1, 3)
     ff(6 )   = shp(3, 1)*shp(2, 2)*shp(1, 3)
     
     ff(7 )   = shp(1, 1)*shp(3, 2)*shp(1, 3)
     ff(8 )   = shp(2, 1)*shp(3, 2)*shp(1, 3)
     ff(9 )   = shp(3, 1)*shp(3, 2)*shp(1, 3)
     
     ff(10)   = shp(1, 1)*shp(1, 2)*shp(2, 3)
     ff(11)   = shp(2, 1)*shp(1, 2)*shp(2, 3)
     ff(12)   = shp(3, 1)*shp(1, 2)*shp(2, 3)
     
     ff(13)   = shp(1, 1)*shp(2, 2)*shp(2, 3)
     ff(14)   = shp(2, 1)*shp(2, 2)*shp(2, 3)
     ff(15)   = shp(3, 1)*shp(2, 2)*shp(2, 3)

     ff(16)   = shp(1, 1)*shp(3, 2)*shp(2, 3)
     ff(17)   = shp(2, 1)*shp(3, 2)*shp(2, 3)
     ff(18)   = shp(3, 1)*shp(3, 2)*shp(2, 3)
     
     ff(19)   = shp(1, 1)*shp(1, 2)*shp(3, 3)
     ff(20)   = shp(2, 1)*shp(1, 2)*shp(3, 3)
     ff(21)   = shp(3, 1)*shp(1, 2)*shp(3, 3)

     ff(22)   = shp(1, 1)*shp(2, 2)*shp(3, 3)
     ff(23)   = shp(2, 1)*shp(2, 2)*shp(3, 3)
     ff(24)   = shp(3, 1)*shp(2, 2)*shp(3, 3)
     
     ff(25)   = shp(1, 1)*shp(3, 2)*shp(3, 3)
     ff(26)   = shp(2, 1)*shp(3, 2)*shp(3, 3)
     ff(27)   = shp(3, 1)*shp(3, 2)*shp(3, 3)

     !step 2: find differentials of weights
     !-------------------------------------

     ! Linearize the FE shape functions wrt each psi(j) direction
     ! Note: only derivatives wrt psi(j) for any j are non-zero, rest are zero

     dshp(:, :, :) = 0.d0
     do j=1,3
        dshp(1, j, j) = half*(psi(j) - one) + half*psi(j)
        dshp(2, j, j) = -2.d0*psi(j)
        dshp(3, j, j) = half*(psi(j) + one) + half*psi(j)
     end do
   
     ! Linearize 27 quadratic weights wrt each psi(j) dir, build from dshp

     loop_psi_j: do j=1,3
        dff(1, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(1, 3   ) &
                   +  shp(1, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
                   +  shp(1, 1   )* shp(1, 2   )*dshp(1, 3, j)

        dff(2, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(1, 3   ) &
                   +  shp(2, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
                   +  shp(2, 1   )* shp(1, 2   )*dshp(1, 3, j)

        dff(3, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(1, 3   ) &
                   +  shp(3, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
                   +  shp(3, 1   )* shp(1, 2   )*dshp(1, 3, j)

        dff(4, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(1, 3   ) &
                   +  shp(1, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
                   +  shp(1, 1   )* shp(2, 2   )*dshp(1, 3, j)

        dff(5, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(1, 3   ) &
                   +  shp(2, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
                   +  shp(2, 1   )* shp(2, 2   )*dshp(1, 3, j)

        dff(6, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(1, 3   ) &
                   +  shp(3, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
                   +  shp(3, 1   )* shp(2, 2   )*dshp(1, 3, j)

        dff(7, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(1, 3   ) &
                   +  shp(1, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
                   +  shp(1, 1   )* shp(3, 2   )*dshp(1, 3, j)

        dff(8, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(1, 3   ) &
                   +  shp(2, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
                   +  shp(2, 1   )* shp(3, 2   )*dshp(1, 3, j)

        dff(9, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(1, 3   ) &
                   +  shp(3, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
                   +  shp(3, 1   )* shp(3, 2   )*dshp(1, 3, j)

        dff(10, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(2, 3   ) &
                    +  shp(1, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
                    +  shp(1, 1   )* shp(1, 2   )*dshp(2, 3, j)

        dff(11, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(2, 3   ) &
                    +  shp(2, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
                    +  shp(2, 1   )* shp(1, 2   )*dshp(2, 3, j)

        dff(12, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(2, 3   ) &
                    +  shp(3, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
                    +  shp(3, 1   )* shp(1, 2   )*dshp(2, 3, j)

        dff(13, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(2, 3   ) &
                    +  shp(1, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
                    +  shp(1, 1   )* shp(2, 2   )*dshp(2, 3, j)

        dff(14, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(2, 3   ) &
                    +  shp(2, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
                    +  shp(2, 1   )* shp(2, 2   )*dshp(2, 3, j)

        dff(15, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(2, 3   ) &
                    +  shp(3, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
                    +  shp(3, 1   )* shp(2, 2   )*dshp(2, 3, j)

        dff(16, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(2, 3   ) &
                    +  shp(1, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
                    +  shp(1, 1   )* shp(3, 2   )*dshp(2, 3, j)

        dff(17, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(2, 3   ) &
                    +  shp(2, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
                    +  shp(2, 1   )* shp(3, 2   )*dshp(2, 3, j)

        dff(18, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(2, 3   ) &
                    +  shp(3, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
                    +  shp(3, 1   )* shp(3, 2   )*dshp(2, 3, j)

        dff(19, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(3, 3   ) &
                    +  shp(1, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
                    +  shp(1, 1   )* shp(1, 2   )*dshp(3, 3, j)

        dff(20, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(3, 3   ) &
                    +  shp(2, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
                    +  shp(2, 1   )* shp(1, 2   )*dshp(3, 3, j)

        dff(21, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(3, 3   ) &
                    +  shp(3, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
                    +  shp(3, 1   )* shp(1, 2   )*dshp(3, 3, j)

        dff(22, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(3, 3   ) &
                    +  shp(1, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
                    +  shp(1, 1   )* shp(2, 2   )*dshp(3, 3, j)

        dff(23, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(3, 3   ) &
                    +  shp(2, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
                    +  shp(2, 1   )* shp(2, 2   )*dshp(3, 3, j)

        dff(24, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(3, 3   ) &
                    +  shp(3, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
                    +  shp(3, 1   )* shp(2, 2   )*dshp(3, 3, j)

        dff(25, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(3, 3   ) &
                    +  shp(1, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
                    +  shp(1, 1   )* shp(3, 2   )*dshp(3, 3, j)

        dff(26, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(3, 3   ) &
                    +  shp(2, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
                    +  shp(2, 1   )* shp(3, 2   )*dshp(3, 3, j)

        dff(27, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(3, 3   ) &
                    +  shp(3, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
                    +  shp(3, 1   )* shp(3, 2   )*dshp(3, 3, j)

     end do loop_psi_j


     ! Step 3: construct Jacobian d(R(psi(:))/d(psi(:)) and residue R(psi(:))
     ! ----------------------------------------------------------------------

     ! loop over x,y,z dirs, stored row-wise
     loop_xyz: do n=1,3 

        ! construct LHS -d(R(psi(:))/d(psi(:))

        ! loop over each psi dir j
        do j=1,3

           ! loop over nodes
           l = 0
           b(n, j) = 0.d0

           do kkk=-1,1
              do jjj=-1,1
                 do iii=-1,1
                    l = l + 1
                    b(n,j) = b(n,j) + dff(l,j) * xElem(n, iii, jjj, kkk)
                 end do
              end do
           end do
        end do !j
       
        ! construct RHS (Xp - sum_i(ff_i * X_i), i =1,27) for nth row (x,y or z)
        l = 0
        b(n, 4) = xSearch(n) 

        ! loop over nodes
        do kkk=-1,1
           do jjj=-1,1
              do iii=-1,1
                 l = l + 1
                 b(n,4) = b(n,4) - ff(l) * xElem(n, iii, jjj, kkk)
              end do
           end do
        end do
       
     end do loop_xyz

     ! Get d(uvwQuadratic) weights
     ! invert 3x3 matrix returns solution in b(:,4) 
     call matrixinv3by3(b, ok_flag)
     if (.not. ok_flag) stop 'Can not invert B in computeQuadraticWeights'
     
     ! update uvwQuadratic weights
     do n=1,3
        uvwQuadratic(n) = uvwQuadratic(n) + b(n,4)
        
        ! sanity check: if weights lie outside [-1:1] reset to 0.5
        if( (uvwQuadratic(n)-1)*(uvwQuadratic(n)+1) > 0.d0) uvwQuadratic(n) = half
     end do

     ! L2-norm of residue
     resid = 0.d0
     do n=1,3
        resid = resid + b(n,4)**2
     end do
  
     niter = niter +1

  end do newton_loop
end subroutine computeQuadraticWeights
