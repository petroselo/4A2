      subroutine check_grid

      use common_block

      implicit none

! Local stuff
      integer :: i, j,test
      real :: xSum, ySum,small

! Check your grid and areas for both the "bump" and the "bend"
! test data.

! First check that all areas are numbers then that they are positive (by program or writing out)

      if ( any(isnan(area(1:ni-1, 1:nj-1))) ) then
        write(6,*) 'A grid area was NaN'
        stop
      end if

      if ( any(area(1:ni-1, 1:nj-1) .le. 0) ) then
        write(6,*) 'A grid area was less than or equal to 0'
        stop
      end if

! Next check that the sum of the length vectors of the 4 faces
! of every element is very nearly zero in each coordinate direction.
! It is absolutely essential that this is correct !
! If not go back and check your subroutine "generate_grid".

! Careful with a test of the form
!          if( a == 0.0 ) then .....
! This will probably never be true.  Computers work to a finite number of
! Significant figures and "a" will probably be +0.0000001 or -0.0000001.
! Test for something like
!          if( abs(a) <= small_number ) then ...
      small = dmin * 0.000001

      do i = 1, (ni-1)
        do j = 1, (nj - 1)
          xSum = dlix(i,j) + dljx(i,j) - dlix(i+1,j) - dljx(i,j+1)
          ySum = dliy(i,j) + dljy(i,j) - dliy(i+1,j) - dljy(i,j+1)
          if (xSum .gt. small .or. ySum .gt. small) then
            write(6,*) 'Sum of edges was greater than 0'
            stop
          end if
        end do
      end do

      write(6,*) 'Max area value = ', maxval(area)
      ! should be around 30/(19x15) for test 0

! Any other tests that you can think of. For example you could plot
! contours of  "area(i,j)" by using  --  call output_hg(area) .

      !call output_hg(area(1:ni-1, 1:nj-1),0) 
      !-can't really get this to work.

      end
