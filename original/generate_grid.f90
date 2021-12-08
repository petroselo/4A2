      subroutine generate_grid

      use common_block

      implicit none

! Local variables - maybe needed for looping?
      integer :: i, j

! Calculate x and y values  x(i,j),y(i,j) of the grid nodes.

! For each value of "i" the i-grid line joins (xlow(i),ylow(i)) to
! (xhigh(i),yhigh(i)). for each value of "i" grid points (nodes) should be
! linearly interpolated between these values for j between 1 and nj.
! i.e.  x(i,1) should be xlow(i), x(i,nj) should be xhigh(i), etc.

      ! Could probably vectorise this, might have to reshape some arrays.
      ! OR at least move common calculations out of the inner loop
      do i = 1, ni
        do j = 1, nj
          ! Careful with integer division here. Order of Operations should make this ok.
          x(i,j) = xlow(i) + ( xhigh(i)-xlow(i) ) * (j-1) / (nj-1)
          y(i,j) = ylow(i) + ( yhigh(i)-ylow(i) ) * (j-1) / (nj-1)    
        end do
      end do

! Calculate the areas of the cells area(i,j)
! (N.B. there are (ni-1) x (nj-1) cells.

! The area of a quadrilateral (regular or irregular) can be shown to be
! half of the cross product of the vectors forming the diagonals.
! see Hirsch volume 1, section 6.2.1. (or lecture).
! Make sure that the area comes out positive!

      do i = 1, (ni-1)
        do j = 1, (nj-1)
          area(i,j) = 0.5 * ( (x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i+1,j)) - (y(i+1,j+1)-y(i,j))*(x(i,j+1)-x(i+1,j)) )
        end do
      end do

! Calculate the x and y components of the length vector of the i-faces
! (i.e. those corresponding to i = constant).
! The length vector of a face is a vector normal to the face wi
! magnitude equal to the length of the face.
! It is positive in the direction of an inward normal to the cell i,j .
! Call these lengths dlix(i,j) and dliy(i,j)

    ! set initial dmin to largest possible real
     dmin = huge(0.0)

      do i = 1, ni
        do j = 1, (nj-1)
          dlix(i,j) = y(i,j+1) - y(i,j)
          dliy(i,j) = x(i,j) - x(i,j+1)
          dmin = min(dmin, dlix(i,j)**2 + dliy(i,j)**2)
        end do
      end do

! Now calculate the x and y components of the length vector of the j-faces. (i.e. those corresponding to j = constant)
! Call these lengths dljx(i,j) and dljy(i,j)

      do i = 1, (ni-1)
        do j = 1, nj
          dljx(i,j) = y(i,j) - y(i+1,j)
          dljy(i,j) = x(i+1,j) - x(i,j)
          dmin = min(dmin, dljx(i,j)**2 + dljy(i,j)**2)
        end do
      end do

! Now find "dmin" the minimum length scale of any element. This is
! defined as the length of the shortest side of the element.
! Call this minimum "dmin". it is used to set the time step from the cfl no.

      ! Above we found the square of dmin. 
      ! Take the root here to avoid repeating expensive sqrt.
      dmin = sqrt(dmin)

      write(6,*)  ' overall minimum element size = ', dmin

      end
