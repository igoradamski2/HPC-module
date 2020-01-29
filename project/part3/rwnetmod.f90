module rwnetmod
	implicit none
	contains

	subroutine rwnet_adjacency(x,n,a)
		!Generate adjacency matrix a corresponding to n node coordinates
		!stored in x
			implicit none
			integer, dimension(:,:), intent(in) :: x
			integer, intent(in) :: n
			integer, dimension(n,n), intent(out) :: a
			integer :: i1, j1,dim,nx,deltaX

			dim = size(x,1)
			nx = size(x,2)
			a = 0

			!Build adjacency matrix
			do i1=1,nx
				do j1=i1+1,nx
					deltaX = maxval(abs(x(:,i1)-x(:,j1)))
					if (deltaX<=1) then
						a(i1,j1) = 1
						a(j1,i1) = 1
					end if
				end do
			end do

	end subroutine rwnet_adjacency

	subroutine rwnet_edgelist(x,l,e)
		!Generate edge list, e, corresponding to l links between coordinates x
			implicit none
			integer, dimension(:,:), intent(in) :: x
			integer, intent(in) :: l
			integer, dimension(l,2), intent(out) :: e
			integer :: i1, j1,k1,dim,n,deltaX

			dim = size(x,1)
			n = size(x,2)
			e = 0
			k1 = 0

			!Build edge list
			do i1=1,n
				do j1=i1+1,n
					deltaX = maxval(abs(x(:,i1)-x(:,j1)))
					if (deltaX<=1) then
						k1 = k1 + 1
						e(k1,1) = i1
						e(k1,2) = j1
					end if
				end do
			end do

	end subroutine rwnet_edgelist

subroutine remove_isolated(x,nnodes,xnew)
	!Remove zero-link nodes from x and store new coordinates in xnew which contains
	!nnodes nodes; nnodes must be pre-computed
	implicit none
	integer, dimension(:,:), intent(in) :: x
	integer, intent(in) :: nnodes
	integer, dimension(size(x,1),nnodes),intent(out) :: xnew
	integer :: i1,j1,count,n,deltaX

	n = size(x,2)
	count = 0

	do i1=1,n
		j1 = 1
		deltaX=2
		do while (deltaX>1 .and. j1<=n)
			if (j1 .ne. i1) deltaX = maxval(abs(x(:,i1)-x(:,j1)))
			j1 = j1 + 1
		end do
		if (deltaX<=1) then
				count = count + 1
				xnew(:,count) = x(:,i1)
			end if
	end do

end subroutine remove_isolated


subroutine count_links(x,l)
	!count number of links in network based on node coordinates x
	!l is needed to compute edge list
	implicit none
	integer, dimension(:,:), intent(in) :: x
	integer, intent(out) :: l
	integer :: i1,j1,n,deltaX

	n = size(x,2)
	l=0
	do i1=1,n
		do j1=i1+1,n
			deltaX = maxval(abs(x(:,i1)-x(:,j1)))
			if (deltaX<=1) then
				l = l + 1
			end if
		end do
	end do

end subroutine count_links

subroutine count_nodes(x,n)
	!count number of zero-node nodes in network based on node coordinates x
	!n is needed by remove_isolated
	implicit none
	integer, dimension(:,:), intent(in) :: x
	integer, intent(out) :: n
	integer :: i1,j1,n0,deltaX

	n0 = size(x,2)
	n = 0

	do i1=1,n0
		j1 = 1
		deltaX=2
		do while (deltaX>1 .and. j1<=n0)
			if (j1 .ne. i1) deltaX = maxval(abs(x(:,i1)-x(:,j1)))
			j1 = j1 + 1
		end do
		if (deltaX<=1) then
				n = n + 1
			end if
	end do

end subroutine count_nodes


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


end module rwnetmod
