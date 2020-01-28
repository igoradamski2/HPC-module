!Simple main program that can be used with the cost module
!and developed/modified as needed
!Should not be submitted with assignment
!To compile: $ gfortran -o main.exe cost.f90 hw2_main.f90
!To run: $ ./main.exe
program hw2_main
  use cost
  implicit none
  integer, parameter :: ndim=5
  real(kind=8) :: xguess(ndim), j

  xguess = (/5.d0,4.d0,3.d0,4.d0,2.d0/)

  call costj(xguess,j)
  print *, 'x guess=', xguess
  print *, 'cost=', j

end program
