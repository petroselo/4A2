      subroutine super_flow_guess

      use common_block

      implicit none

! Local stuff
      real :: mflow,machlim,tlim,mach_num,Tdown,e,dy,dx,dxy, vdown, rodown
      real, dimension(i_max) :: Tstatic,Pstatic
      integer :: i, j

! In this subroutine we make an initial guess of the primary variables
! i.e.  ro, rovx, rovy and roe.
! The guess does not need to be very accurate but the better
! it is the faster your program will converge.
! You should assign values to ro(i,j), rovx(i,j), rovy(i,j) and roe(i,j)
! at every grid point in this subroutine.

! Work out the length of each "i" line between grid points "i,1" snd "i,nj"
! and call it  "aflow(i)" .

      aflow = sqrt( (xhigh-xlow)**2 + (yhigh-ylow)**2 )

! Make an initial guess of the density and velocity at the exit by
! assuming isentropic flow between the inlet stagnation pressure pstagin
! and temperature tstagin and the exit static pressure pdown.
! Use these together with "aflow(ni)" to estimate the mass flow rate.
! call this "mflow".

      ! Isentropic perfect gas relation from inlet to exit.
      Tdown = tstagin * (pdown / pstagin)**fga

      ! Ideal gas law at exit
      rodown = pdown / (rgas * Tdown)

      ! Isentropic => T01 == T02 => exit velocity:
      vdown = sqrt(2*cp*(tstagin - Tdown))

      ! Mass flow rate mflow = roAV at exit:
      mflow = rodown * aflow(ni) * vdown
      

! Set a limit to the maximum allowable mach number in the initial
! guess. call this "machlim". calculate the corresponding temperature.

      machlim = 10 !### Updated for supersonic flow
      tlim = tstagin/(1.0 + 0.5*(gamma-1.0)*machlim*machlim)

      v_guess = vdown
      ro_guess = rodown

! Direct the velocity found above along the j= constant grid lines to find
! the velocity vx(i,j) in the  x  direction and vy(i,j) in the y.
! Use these and ro_guess(i) to set rovx(i,j), rovy(i,j) and roe(i,j).
! Also set ro(i,j).
! Note that roe(i,j) includes the kinetic energy component of the
! internal energy.

      ! Set velocity components based on upstream j vector.
      
      do j = 1,nj
        ! Set inlet based on downstream as no upstream avilable.
        
        dx = dljy(1,j)
        dy = -dljx(1,j)
        dxy = sqrt(dljy(1,j)**2 + dljx(1,j)**2)
        vx(1, j) = v_guess(1) * dx / dxy
        vy(1, j) = v_guess(1) * dy / dxy

        do i = 2,ni
          dx = dljy(i-1,j)
          dy = -dljx(i-1,j)
          dxy = sqrt(dljy(i-1,j)**2 + dljx(i-1,j)**2)
          vx(i, j) = v_guess(i) * dx / dxy
          vy(i, j) = v_guess(i) * dy / dxy
        end do

      end do

      ! Set variables
      do i=1,ni
        do j=1,nj
          ro(i,j) = ro_guess(i)
          rovx(i,j) = ro(i,j) * vx(i,j)
          rovy(i,j) = ro(i,j) * vy(i,j)
        end do
      end do


! Store the "old" values of the variables for use in the first
! convergence check in subroutine "check_conv"

      do i=1,ni
        do j=1,nj
          ro_old(i,j)   = ro(i,j)
          rovx_old(i,j) = rovx(i,j)
          rovy_old(i,j) = rovy(i,j)
        end do
      end do

      end
