      subroutine set_others

      use common_block

      implicit none

! Local stuff
      integer  ::  i, j
!      real, dimension(ni,nj) :: v, tstatic
      real, dimension(i_max,j_max) :: halfvsq, tstatic
! This routine calculates secondary flow variables from the primary ones
! at every grid point.

! The primary variables are ro, rovx, rovy and roe

! The secondary variables are the velocity components vx(i,j) and vy(i,j),
! the static pressure p(i,j) and the stagnation enthalpy hstag(i,j).
! Note:  "hstag"  not  "ho".


      ! Maybe try subslicing for greater efficiency?
      vx = rovx / ro
      vy = rovy / ro
      
      tstatic = tstagin - oneover2Cp * (vx**2 + vy**2)
      p = ro*rgas*tstatic
      ! hstag = cp * tstagin set only once before loop to save time: 
      end
