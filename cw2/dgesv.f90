module mydgesv
  implicit none

contains

  subroutine lineq_sol(H, g, x_change)

    implicit none
      real(kind=8), dimension(:,:), intent(in) :: H
      real(kind=8), dimension(:,:), intent(in) :: g
      real(kind=8), allocatable, dimension(:,:), intent(out) :: x_change
      integer :: N,NRHS,LDA,LDB,INFO
      integer, allocatable, dimension(:) :: IPIV


      !set dimension for H and g
      print *, 'we are here'
      N = size(g)
      NRHS = 1
      LDA = N
      LDB = N
      allocate(IPIV(N))

      call dgesv(N, NRHS, H, LDA, IPIV, g, LDB, INFO)

      !extract solution
      allocate(x_change(N,NRHS))
      x_change = g(1:N,:)

  end subroutine lineq_sol

end module mydgesv
