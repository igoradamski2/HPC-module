module hw3
  use omp_lib
  implicit none
  contains

    subroutine rwnet(h,dim,n,display,x)
      !Simulate random network model
      !Input variables
      !h: Height at which nodes are initially introduced
      !n: Total number of nodes
      !dim: dimension of node coordinates
      !display: Coordinate of every 100th node added to network is printed to screen when display=1
      !Output variables
      !x: Final node coordinates

      implicit none
      integer, intent(in) :: h,dim,n,display
      integer, dimension(dim,n), intent(out) ::  x
      integer :: i1,j1,deltaX,dist
      integer, dimension(dim) :: r,xx

      !initialize
      x = 0

      do i1 = 2,n
        xx(1:dim-1) = 0 !introduce new node
        xx(dim) = h
        dist=h
        do j1=1,i1-1
          deltaX = maxval(abs(x(:,j1)-xx))
          dist = min(dist,deltaX)
        end do

        do while ((xx(dim)>0) .and. (dist>1))
          call random_int(dim,r)
          r(1:dim-1) = r(1:dim-1)*2-1
          r(dim) = -r(dim)
          xx = xx + r
          j1=i1-1
          do while (j1>0 .and. dist>1)
            deltaX = maxval(abs(x(:,j1)-xx))
            dist = min(dist,deltaX)
            j1=j1-1
          end do
        end do
        if (display==1  .and. mod(i1,100)==0) print *, 'i1,xx=',i1,xx
        x(:,i1) = xx
      end do
    end subroutine rwnet
!------------------------------------------------------------
subroutine simulate_m(h,dim,n,m,numthreads,x,hmax,hmin,hmean)
  !Simulate m networks
  !Input variables -- same as in rwnet except:
  !m: Number of simulations
  !numthreads: number of threads used in parallel regions
  !Output variables
  !x: Coordinates for all m networks
  !hmax,hmin,hmean: max,min, and average final network heights

  !MY APPROACH TO PARALLELIZATION:
  !In the below function I have parallelized the do loop which
  !loops through m number of times each time calling the network
  !generation algorithm rwnet. My approach to parallelization is
  !such that each of the numthreads is going to simultaneously
  !generate a network with rwnet. Therefore one core will do one
  !network generation at a time so using multiple threads will
  !speed up the simulation because more than one network can be
  !generated at a time. In the parallel loop I do three reductions
  !one for the max height, one for the min height and one for the hmean. This
  !way the threads can keep track of the actual minimum and maximum
  !of the network heights across all simulations. The hmean gets summed every
  !iteration and at the end i divide it by m to get the mean height. This
  !is faster than for example keeping the means in one vector and then outside
  !the paralell region sum the components (for example) because this way
  !the summation is distributed over numthread threads rather than being computed
  !by one core.
  !So to summarize, basically instead of doing m loops with one core
  !which would take some time t, I now basically parallelize the loop
  !so that numthreads cores are dividing the m iterations among
  !themselves, which theoretically should take time t/m (assuming that
  !t/m is not lower than a single network generation).

  implicit none
  integer, intent(in) :: h,dim,n,m,numthreads
  integer, dimension(dim,n,m), intent(out) ::  x
  integer, intent(out) :: hmax,hmin
  integer :: curr_height
  real(kind=8), intent(out) :: hmean
  integer i1
  integer, dimension(dim,n) :: dummy

  hmin = 100
  !$ call omp_set_num_threads(numthreads)
  !$OMP parallel do reduction(max:hmax), reduction(min:hmin), reduction(+:hmean)
  do i1 = 1,m
    call rwnet(h,dim,n,0,dummy)
    x(:, :, i1) = dummy
    curr_height = maxval(dummy(dim,:))
    hmean = hmean + curr_height
    hmax = max(hmax, curr_height)
    hmin = min(hmin, curr_height)
  end do
  !$OMP end parallel do
  hmean = hmean/m

end subroutine simulate_m


!---------------------------
subroutine random_int(m,r)
  !Generate m random numbers each of which is 0 or 1
  !with equal probability and return numbers in array, r
  implicit none
  integer, intent(in) :: m
  integer, dimension(m), intent(out) :: r
  integer :: i1
  real(kind=8),dimension(m) :: rtemp

  call random_number(rtemp)
  r = 0
  do i1 = 1,m
    if (rtemp(i1)>=0.5) r(i1) = 1
  end do

end subroutine random_int

end module hw3
