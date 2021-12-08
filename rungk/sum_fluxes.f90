      subroutine sum_fluxes(iflux, jflux, delprop, prop_inc, tfrac)

      use common_block

      implicit none

! This subroutine sums the fluxes for each element, calculates the changes
! in the variable "prop" appropriate to a cell (delprop) and distributes them to the
! four corners of the element (stored in prop_inc).

      real, dimension(i_max,j_max) ::  iflux, jflux, prop_inc, delprop!,  &
                                       !previous, store
      real :: tfrac

! Local stuff
      integer ::    i, j

! Find the change in the variable "prop" in each cell over the
! time step "deltat" and save it in "delprop(i,j)".

      delprop(1:ni-1, 1:nj-1) = ( tfrac * deltat / area(1:ni-1, 1:nj-1) ) * &
                                ( iflux(1:ni-1, 1:nj-1) &
                                + jflux(1:ni-1, 1:nj-1) &
                                - iflux(2:ni, 1:nj-1)  &
                                - jflux(1:ni-1, 2:nj) )


! Now distribute the changes equally to the four corners of each
! cell (nodes). Each interior grid point receives one quarter of the change
! from each of the four cells adjacent to it.  Edge nodes do not have four adjacent
! cells.

      ! Interior Nodes
      do i=2,ni-1
        do j=2,nj-1
          prop_inc(i,j) = 0.25 * ( delprop(i-1,j-1) + delprop(i,j) + &
                                   delprop(i,j-1) + delprop(i-1,j) )
        enddo
      enddo

! Now deal with the changes to the upper and lower boundaries.
! These receive half the change from each of the two cells adjacent to them.

      do i=2,ni-1

        prop_inc(i,1) = 0.5 * ( delprop(i-1,1) + delprop(i,1) )

        prop_inc(i,nj) = 0.5 * ( delprop(i-1,nj-1) + delprop(i,nj-1) )

      enddo

! Now deal with changes to the inlet & outlet boundary points.
! These receive half the change from each of the two cells adjacent to them.

      do j=2,nj-1

        prop_inc(ni,j) = 0.5 * ( delprop(ni-1,j-1) + delprop(ni-1,j) )

        prop_inc(1,j) = 0.5 * ( delprop(1,j-1) + delprop(1,j) )

      enddo

! Finally find the changes to be added to the four corner points.
! These receive the full change from the single cell of which they form one corner.

     prop_inc(1,1) = delprop(1,1)

     prop_inc(1,nj) = delprop(1, nj-1)

     prop_inc(ni,1) = delprop(ni-1,1)

     prop_inc(ni,nj)  = delprop(ni-1,nj-1)

! If the values in delprop have been hi-jacked by the second order timestep extension,
! restore them here to the true residuals. (The second order extension should have
! loaded them into the variable "store" before hi-jacking them)
! These are used in the convergence check and also in future improvements to the scheme.

      end
