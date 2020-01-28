!Simple driver program to run routines in hw3 module
!To compile: gfortran -fopenmp -O3 -o main.exe hw3_dev.f90 hw3_main.f90
program test
  use hw3
  implicit none
  integer :: n,h,i1,dim
  integer, allocatable, dimension(:,:) :: x

  dim = 2
  h = 100
  n = 400
  allocate(x(dim,n))
  call rwnet(h,dim,n,1,x)

  open(unit=12,file='x.dat')
  do i1=1,n
    write(12,*) x(:,i1)
  end do
  close(12)

!To load data in python: np.loadtxt('x.dat')

end program test
