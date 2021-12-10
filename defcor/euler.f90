      program euler

      use common_block

      implicit none

! Local stuff
      integer i,j

! Open files to store the convergence history. Plotting is done via a separate
! program. "euler.log" is for use by pltconv. "pltcv.csv" is for use by paraview.

      open(unit=3,file='euler.log')
      open(unit=31,file='pltcv.csv')

      call read_data

      call generate_grid

      call check_grid

     !call crude_guess
     call flow_guess

! "output" here to plot out your initial guess of the flow field.
      !call output_hg(ro,1)
      !call output_mat(1)
      !call output(1)

! "set_timestep": to set the length of the timestep.
! initially this is a constant time step based on a conservative guess
! of the mach number.

      call set_timestep

      conv_check_steps = 5

!************************************************************************
!     start the time stepping do loop for "nsteps" loops.
!************************************************************************

      ! Constant stagnation enthalpy:
      hstag = cp * tstagin

      nrkut_max = 4

      do nstep = 1, nsteps

        do i=1,ni
          do j=1,nj
            ro_start(i,j) = ro(i,j)
            rovx_start(i,j) = rovx(i,j)
            rovy_start(i,j) = rovy(i,j)
          end do
        end do

        ! Runge--Kutta Loop
        
        do nrkut = 1,nrkut_max
          
          frkut = 1./(1.+nrkut_max - nrkut)
          
          call set_others

          call apply_bconds
  
          call set_fluxes

          !               i-flux      j-flux      delprop  prop_inc  delta-t
          call sum_fluxes(fluxi_mass, fluxj_mass, delro  , ro_inc  , frkut)
          call sum_fluxes(fluxi_xmom, fluxj_xmom, delrovx, rovx_inc, frkut)
          call sum_fluxes(fluxi_ymom, fluxj_ymom, delrovy, rovy_inc, frkut)
 
        


! Update solution

        do i=1,ni
          do j=1,nj
            ro  (i,j) = ro_start  (i,j) + ro_inc  (i,j)
            rovx(i,j) = rovx_start(i,j) + rovx_inc(i,j)
            rovy(i,j) = rovy_start(i,j) + rovy_inc(i,j)
          end do
        end do

! Smooth the problem to ensure it remains stable.

        call smooth(ro, corr_ro)
        call smooth(rovx, corr_rovx)
        call smooth(rovy, corr_rovy)

      end do

! Check convergence and write out summary every 5 steps

        if(mod(nstep, conv_check_steps) == 0) then
          call check_conv
        end if

! Stop looping if converged to the input tolerance "conlim"

        if( emax < conlim .and. eavg < (0.5*conlim) ) then
          write(6,*) ' Calculation converged in ',nstep,' iterations'
          write(6,*) ' To a convergence limit of ', conlim
          exit
        endif

      end do

!************************************************************************
!  end of time stepping do loop for "nsteps" loops.
!************************************************************************

! Calculation finished. call "output" to write the plotting file.
! N.B. Solution hasn't necessarily converged.
      !write(6,*) 'Write output'
      call output(ro, 1)   ! Paraview
      call output_hg(ro,1) ! Eulplt
      call output_mat(1)   ! Matlab

      close(3)
      close(31)

      end