      subroutine read_data

      use common_block

      implicit none

! Local stuff
      real ::   conlim_in, smooth_fac_in
      integer :: i, jmid

! Assign unit 1 to the file 'geom'
! Assign unit 2 to the file 'flow'

      open(unit=1,file='geom')
      open(unit=2,file='flow')

! Read in the title and ni and nj from unit 1.
! Before you do this, check the format of the data file "geom" and read the
! appropriate values for ni & nj.

      read(1,*) title
      read(1,*) ni, nj

! Check that ni and nj are less than the array dimensions i_max and j_max.
! Then read in xlow, ylow, xhigh, yhigh for each i between 1 and ni.

      if (ni > i_max .or. nj > j_max) then
        write(6,'(4(A,I7))') 'Geometry file ni=',ni,' , nj=',nj,' exceeded i_max=',i_max,' ,j_max=',j_max
        stop
      end if
      
      do i = 1, ni
        ! These are initialised in the common_block
        read(1,*) xlow(i), ylow(i), xhigh(i), yhigh(i)
      end do

! Now read in the flow data from unit 2.
! You should read in the following variables sequentially:

!       rgas, gamma
!       pstagin, tstagin, alpha1, pdown
!       cfl, smooth_fac_in
!       nsteps, conlim_in

      read(2,*) rgas, gamma
      read(2,*) pstagin, tstagin, alpha1, pdown
      read(2,*) cfl, smooth_fac_in
      read(2,*) nsteps, conlim_in
      

! Set some other variables that will be used throughout the
! calculation. Change alpha1 to radians. Scale the smoothing factor
! and convergence limits by cfl: changes over a timestep should be proportional
! to deltat (which is proportional to cfl).

      emax       = 1000000.0
      eavg       = emax
      cp         = rgas*gamma/(gamma-1.0)
      cv         = cp/gamma
      fga        = (gamma - 1.0)/gamma
      smooth_fac = smooth_fac_in*cfl
      conlim     = conlim_in*cfl
      alpha1     = alpha1*3.14159/180.0

! Close the flow and geom files

      close(1)
      close(2)

! Calculate the reference values which are used to check convergence

      ncells    =  ni * nj
      jmid      =  (1 + nj)/2
      roin      =  pstagin/rgas/tstagin
      ref_ro    =  (pstagin-pdown)/rgas/tstagin
      ref_t     =  tstagin*(pdown/pstagin)**fga
      ref_v     =  sqrt(2*cp*(tstagin-ref_t))
      ref_rovx  =  roin*ref_v
      ref_rovy  =  ref_rovx
      ref_roe   =  roin*cv*(tstagin-ref_t)

      end
