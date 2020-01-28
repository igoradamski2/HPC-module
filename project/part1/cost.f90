!Module containing cost function and related routines for use with optimization routines
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

      j = -x(1)*sin(sqrt(abs(x(1)-(x(2)+47)))) - (x(2)+47)*sin(sqrt(abs(x(1)/2+(x(2)+47))))

      if (c_noise) call add_noise(x,j)

  end subroutine costj


  subroutine costjold(x,j)
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

  end subroutine costjold


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
