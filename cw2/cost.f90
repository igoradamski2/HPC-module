!Module containing cost function and related routines for use with optimization routines
!IGOR ADAMSKI: 01069152
module cost
  implicit none
  logical :: c_noise = .False. !when true, white noise is added to cost function using add_noise
  real(kind=8) :: c_noise_amp=0.d0 !amplitude of noise

contains

  subroutine costj(x,j)
    !evaluate cost function at x, store cost in j
    implicit none
      real(kind=8), dimension(:), intent(in) :: x !do not need to explicitly specify dimension of input variable when subroutine is within a module
      real(kind=8), intent(out) :: j
      integer i1

      j = 0.d0
      do i1=1,size(x)-1
        j = j + (1.d0 - x(i1))**2 + 20.d0*(x(i1+1) - x(i1)**2)**2
      end do

      if (c_noise) call add_noise(x,j)

  end subroutine costj

  subroutine costj_grad2d(x,j_grad)
    !evaluate gradient of cost function at x, store grad in j_grad
    !does not account for noise, only correct for 2d case.
    implicit none
      real(kind=8), dimension(:), intent(in) :: x
      real(kind=8), intent(out) :: j_grad(size(x))

    j_grad=0.d0
    j_grad(1) = 2*(x(1)-1) - 4*20*x(1)*(x(2) - x(1)**2)
    j_grad(2) = 2*20*(x(2) - x(1)**2)

  end subroutine costj_grad2d

  subroutine costj_hess2d(x,j_hess)
    !evaluate Hessian of cost function at x, store Hessian in j_hess
    !does not account for noise, only correct for 2d case
    implicit none
      real(kind=8), dimension(:), intent(in) :: x
      real(kind=8), intent(out) :: j_hess(size(x),size(x))

      j_hess=0.d0
      j_hess(1,1) = 12*20*x(1)**2-4*20*x(2)+2
      j_hess(1,2) = -4*20*x(1)
      j_hess(2,1) = -4*20*x(1)
      j_hess(2,2) = 2*20

  end subroutine costj_hess2d


  subroutine add_noise(x,j)
    !adds noise with amplitude c_noise_amp to j at x
    implicit none
    real(kind=8), dimension(:), intent(in) :: x
    real(kind=8), intent(inout) :: j

    integer :: i1,k
    integer,allocatable,dimension(:) :: seed
    real(kind=8) :: noise

    call random_seed(size=k)
    allocate(seed(k))
    do i1=1,size(x)
      seed(:) = nint(x(i1)*1000000)
      call random_seed(put=seed)
      call random_number(noise)
      j = j + c_noise_amp*(noise-0.5d0)
    end do

  end subroutine add_noise


end module cost
