!Igor Adamski, CID 01069152
module part3
  use rwnetmod
  implicit none
  integer, allocatable, dimension(:,:) :: a !adjacency matrix
  integer, allocatable, dimension(:,:) :: e !edge list
  integer, allocatable, dimension(:,:) :: xnew !network node coordinates after isolated-node removal
  real(kind=8), allocatable, dimension(:,:) :: u !wave amplitude
  real(kind=8) :: c=1 !wave speed
contains

subroutine wave_init(h,dim,n,flag_init,nnodes,nlinks)
  !simulate random network with height h, dim d, and n total nodes
  !flag_init = 1: adjacency matrix is generated
  !flag_init = 2: edge list is generated
  !flag_init = 3: both are generated
  !output: nnodes and nlinks: total number of nodes and links in graph
  ! (with no zero-link nodes) corresponding to random network
  implicit none
  integer, intent(in) :: h,dim,n,flag_init
  integer, dimension(dim,n) :: x !network coordinates
  integer,intent(out) :: nnodes,nlinks

  call rwnet(h,dim,n,0,x) !generate network
  call count_nodes(x,nnodes)
  call count_links(x,nlinks)
  if (allocated(xnew)) deallocate(xnew)
  allocate(xnew(dim,nnodes))
  call remove_isolated(x,nnodes,xnew) !remove isolated nodes

!Generate adjacency matrix and/or edge list based on xnew (with isolated nodes removed)
  if (flag_init==1 .or. flag_init==3) then
      if (allocated(a)) deallocate(a)
      allocate(a(nnodes,nnodes))
      call rwnet_adjacency(xnew,nnodes,a)
  end if
  if (flag_init==2 .or. flag_init==3) then
      if (allocated(e)) deallocate(e)
      allocate(e(nlinks,2))
      call rwnet_edgelist(xnew,nlinks,e)
  end if

end subroutine wave_init


subroutine wave(h,dim,nt,dt,ntotal,energy)
  implicit none
  integer, intent(in) :: h,dim,nt,ntotal
  real(kind=8), intent(in) :: dt
  real(kind=8), dimension(nt+1), intent(out) :: energy
  integer :: nnodes,nlinks, i1, t1
  real(kind=8), allocatable :: y0(:), y(:,:)
  real(kind=8) :: dummy

  !"generate" graph
  call wave_init(h,dim,ntotal,3,nnodes,nlinks) !set flag_init as appropriate

  allocate(y0(nnodes*2))
  allocate(y(nnodes*2,nt+1))
  y0 = 0.d0
  y0(nnodes) = 1
  call rk4(0.d0,y0,dt,nt,y)
  u = y(1:nnodes,:)

  do t1=1,nt+1
    dummy=0
    do i1=1,nnodes
      dummy = dummy + u(i1, t1)**2
    end do
    energy(t1) = dummy
  end do


end subroutine wave



subroutine rk4(t0,y0,dt,nt,y)
    !4th order RK method
    !input:
    !t0: initial time
    !y0: initial condition (array)
    !dt: time step
    !nt: number of time steps
    !output:
    !y: solution at each time step (array)
    !order: order parameter
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
    integer, intent (in) :: nt
    real(kind=8), dimension(size(y0),nt+1), intent(out) :: y
    real(kind=8), dimension(size(y0)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac
    integer:: k
        halfdt = 0.5d0*dt
        fac = 1.d0/6.d0

        y(:,1) = y0 !initial condition
        t = t0 !initial time
        do k = 1, nt !advance nt time steps

	         f1 = dt*RHS(t, y(:,k))

	         f2 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f1)

           f3 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f2)

           f4 = dt*RHS(t + dt, y(:,k) + f3)

	         y(:,k+1) = y(:,k) + (f1 + 2*f2  + 2*f3 + f4)*fac

           t = t + dt

        end do


end subroutine rk4

!---------------------------
function RHS(t,f)
    !RHS for rk4 correspoinding to ODEs,
    !df/dt = RHS(f)
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)) :: RHS
    integer :: i1,j1,n
    real(kind=8) :: dummy

    n = size(f)/2

    do i1=1,size(f)
      if (i1 <= n) then
        RHS(i1) = f(i1+n)
      else
        dummy = 0
        do j1=1,n
          dummy = dummy + a(i1-n,j1)*(f(j1) - f(i1-n))
        end do
        RHS(i1) = (c**2)*dummy
      end if
    end do

end function RHS

!---------------------------


end module part3
