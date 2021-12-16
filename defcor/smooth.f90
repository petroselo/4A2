      subroutine smooth(prop, corr_prop)

      use common_block

      implicit none

! This subroutine smooths the variable "prop" (i.e. it adds the
! artificial viscosity) by taking (1-sf) * the calculated value of
! "prop" + sf * (the average of the surrounding values of "prop").
! where sf is the smoothing factor.
! It is modified if the deferred correction extension is adopted.

      real, dimension(i_max,j_max) ::  prop, corr_prop

! Local stuff
      real, dimension(i_max,j_max) ::  store
      real  ::  sf, sfm1, avg, avg1, avgnj
      integer  ::  i, j, ip1, im1

! Deferred corrections
      real :: corrnew, corrnew1, corrnewnj

! To avoid using already smoothed values for smoothing other values
! the smoothed values are initially stored in an array "store".
! Wall nodes do not have 4 neighbours, so the formula for smoothing them
! is different.
! Inlet and exit values are constrained by boundary conditions and extrapolation from
! previously smoothed nodes, so how they are smoothed has little material effect
! on the solution.  Remember that they do not have 4 neighbours.

      sf   = smooth_fac
      sfm1 = 1.0 - sf

      do i=1,ni
        ip1 = i+1
        if( i==ni ) ip1 = ni
        im1 = i-1
        if( i==1  ) im1 = 1

         do j=2,nj-1
           avg = 0.25*(prop(ip1,j)+prop(im1,j)+prop(i,j-1)+prop(i,j+1))
           
           corrnew = fcorr * (prop(i,j) - avg)
           corr_prop(i,j) = .99 * corr_prop(i,j) + 0.01 * corrnew


           store(i,j) = sfm1*prop(i,j) + sf*( avg + corr_prop(i,j) )
         end do

! On the surfaces j=1 and j=nj take the average as shown below.

         avg1  = (prop(im1,1)  + prop(ip1,1)  + 2.*prop(i,2)    - prop(i,3)    ) / 3.0
         corrnew1 = fcorr * ( prop(i,1) - avg1 )
         corr_prop(i,1) = .99 * corr_prop( i,1 ) + 0.01 * corrnew1

         avgnj = (prop(im1,nj) + prop(ip1,nj) + 2.*prop(i,nj-1) - prop(i,nj-2) ) / 3.0
         corrnewnj = fcorr * ( prop(i,nj) - avgnj )
         corr_prop(i,nj) = .99 * corr_prop( i,nj ) + 0.01 * corrnewnj

         store(i,1) = sfm1*prop(i,1) + sf*( avg1 + corr_prop(i,1) )
         store(i,nj) = sfm1*prop(i,nj) + sf*( avgnj + corr_prop(i,nj) )

      end do

! Set prop to the smoothed value before returning to the main program.

      do i=1,ni
        do j=1,nj
          prop(i,j) = store(i,j)
        end do
      end do

      end