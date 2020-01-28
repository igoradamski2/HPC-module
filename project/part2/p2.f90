!Dev code for project, part 2: Simulate Pedestrian motion model
!Igor Adamski, CID 01069152
module part2
  implicit none
  ! integer :: isample
  ! real(kind=8), allocatable :: x(:,:),y(:,:)
contains

  subroutine simulate(m,nt,s0,l,d,a,x,y,psi)
    !Pedestrian motion model
    !input variables:
    !n =  m^2 = number of students
    !nt: number of time steps
    !s0: student speed
    !l: initial spacing between students
    !d: student motion influence by all other students within distance <= d
    !a: noise amplitude
    !output variables:
    !x,y: all m^2 student paths from i=1 to nt+1
    !psi: synchronization parameter, store at all nt+1 times (including initial condition)
    implicit none
    integer, intent(in) :: m,nt
    real(kind=8), intent(in) :: s0,l,d,a
    real(kind=8), dimension(m*m,nt+1), intent(out) :: x,y
    real(kind=8), dimension(nt+1), intent(out) :: psi
    integer :: i1,j1,k1,m1,n1
    real(kind=8) :: init_dir(2, m*m), u(m*m,2:nt+1), v(m*m,2:nt+1), pi, diff(m*m,2)
    real(kind=8) :: theta, matrix_positions(m*m, 2), vector_norm(m*m)
    real(kind=8) :: u_bar(m*m), v_bar(m*m), R(m*m), dummy(2)


    !initialize pi
    pi = acos(-1.d0)
    ! initialize student positions
    !x = 0.d0
    !y = 0.d0
    k1 = 0
    do i1 = 1,m
      do j1=1,m
        k1 = k1 + 1
        x(k1,1) = (j1-1)*l/2 - (m-1)*l/4
        y(k1,1) = (i1-1)*l/2 - (m-1)*l/4
      end do
    end do
    x(:,1) = x(:,1)/(m-1)
    y(:,1) = y(:,1)/(m-1)

    psi(1) = 0d0

    !initial direction
    call random_number(init_dir)
    init_dir = 2*pi*init_dir - pi
    u(:,2) = init_dir(1,:)
    v(:,2) = init_dir(2,:)
    x(:,2) = x(:,1) + u(:,2)
    y(:,2) = y(:,1) + v(:,2)
    dummy = 0
    do k1=1,m*m
      !dummy = dummy + [u(k1,2),v(k1,2)]
      dummy(1) = dummy(1) + u(k1,2)
      dummy(2) = dummy(2) + v(k1,2)
    end do
    psi(2) = sqrt((dummy(1)**2)+(dummy(2)**2))/(m*m*s0)

    do i1=3,nt+1
      matrix_positions(:,1) = x(:,i1-1)
      matrix_positions(:,2) = y(:,i1-1)
      call random_number(R)
      R = 2*a*pi*R - a*R
      do j1=1, m*m
        diff(:,1) = matrix_positions(:,1) - matrix_positions(j1,1)
        diff(:,2) = matrix_positions(:,2) - matrix_positions(j1,2)
        vector_norm = diff(:,1)*diff(:,1) + diff(:,2)*diff(:,2)
        k1 = 0
        u_bar(j1) = 0
        v_bar(j1) = 0
        do m1=1, m*m
          if (sqrt(vector_norm(m1)) < d) then
            u_bar(j1) = u_bar(j1) + u(m1, i1-1)
            v_bar(j1) = v_bar(j1) + v(m1, i1-1)
            k1 = k1 + 1
          end if
        end do
        u_bar(j1) = u_bar(j1)/k1
        v_bar(j1) = v_bar(j1)/k1
        theta = atan(v_bar(j1)/u_bar(j1)) + R(j1)
        u(j1,i1) = s0*cos(theta)
        v(j1,i1) = s0*sin(theta)
      end do
      x(:,i1) = x(:,i1-1) + u(:,i1)
      y(:,i1) = y(:,i1-1) + v(:,i1)
      do n1=1, m*m
        if ((abs(x(n1,i1)) - L) > 0) then
          if (x(n1,i1) > 0) then
            x(n1,i1) = x(n1,i1) - 2*L
          else
            x(n1,i1) = x(n1,i1) + 2*L
          end if
        end if
        if ((abs(y(n1,i1)) - L) > 0) then
          if (y(n1,i1) > 0) then
            y(n1,i1) = y(n1,i1) - 2*L
          else
            y(n1,i1) = y(n1,i1) + 2*L
          end if
        end if
      end do

      dummy = 0
      do k1=1,m*m
        dummy(1) = dummy(1) + u(k1,i1)
        dummy(2) = dummy(2) + v(k1,i1)
      end do
      psi(i1) = sqrt((dummy(1)**2)+(dummy(2)**2))/(m*m*s0)
    end do


  end subroutine simulate


  subroutine simulate_omp(m,nt,s0,l,d,a,numthreads,x,y,psi)
    !Pedestrian motion model
    !input variables:
    !n = m^2 = number of students
    !nt: number of time steps
    !s0: student speed
    !l: initial spacing between students
    !d: student motion influence by all other students within distance <= d
    !a: noise amplitude
    !numthreads: number of threads to use in parallel regions
    !output variables:
    !x,y: all m^2 student paths from i=1 to nt+1
    !psi: synchronization parameter, store at all nt+1 times (including initial condition)

    !MY APPROACH TO PARALLELIZATION:

    !My approach is as follows. Every loop before the main loop, so some of the initialization
    !loops, are parallelized with omp, with paralell do's. So for example the double loop in which
    !we initialize the positions of walkers, is paralellized with omp so that one of the threads takes
    !half of the iterations. When we get to the main loop, so the marching in time, we cannot paralellize
    !this one, because timestep t uses information from timestep t-1 and so on. Therefore I paralellize as
    !many do loops as possible in the main time loop. The first one, is the distance calculations, where I calculate
    !all the needed variables needed to update the position of a walker.
    !I use a paralell do with several private variables, so that each node computes the relevant portion of the loop.
    !The next loop that I paralellize is the loop in which we check that our steps did not go beyond the box
    !in which we want to keep our walkers.
    !So basically every do loop that is possible to paralellize with omp (so ones that do not require
    !knowledge of all the previous iterations like with the time loop) I paralellize with omp.
    use omp_lib
    implicit none
    integer, intent(in) :: m,nt,numthreads
    real(kind=8), intent(in) :: s0,l,d,a
    real(kind=8), dimension(m*m,nt+1), intent(out) :: x,y
    real(kind=8), dimension(nt+1), intent(out) :: psi
    integer :: i1,j1,k1,m1,n1
    real(kind=8) :: init_dir(2, m*m), u(m*m,2:nt+1), v(m*m,2:nt+1), pi, diff(m*m,2)
    real(kind=8) :: theta, matrix_positions(m*m, 2), vector_norm(m*m)
    real(kind=8) :: u_bar(m*m), v_bar(m*m), R(m*m), dummy(2)


    !setup omp
    !$ call omp_set_num_threads(numthreads)

    !initialize pi
    pi = acos(-1.d0)
    ! initialize student positions
    !x = 0.d0
    !y = 0.d0

    k1 = 0
    !$OMP parallel do private(j1)
    do i1 = 1,m
      do j1=1,m
        k1 = k1 + 1
        x(k1,1) = (j1-1)*l/2 - (m-1)*l/4
        y(k1,1) = (i1-1)*l/2 - (m-1)*l/4
      end do
    end do
    !$OMP end parallel do

    x(:,1) = x(:,1)/(m-1)
    y(:,1) = y(:,1)/(m-1)

    psi(1) = 0

    !initial direction
    call random_number(init_dir)
    init_dir = 2*pi*init_dir - pi
    u(:,2) = init_dir(1,:)
    v(:,2) = init_dir(2,:)
    x(:,2) = x(:,1) + u(:,2)
    y(:,2) = y(:,1) + v(:,2)
    dummy = 0
    !$OMP parallel do reduction(+:dummy)
    do k1=1,m*m
      dummy = dummy + [u(k1,2),v(k1,2)]
    end do
    !$OMP end parallel do
    psi(2) = sqrt((dummy(1)**2)+(dummy(2)**2))/(m*m*s0)

    do i1=3,nt+1
      matrix_positions(:,1) = x(:,i1-1)
      matrix_positions(:,2) = y(:,i1-1)
      call random_number(R)
      R = 2*a*pi*R - a*R
      !$OMP parallel do private(diff, vector_norm, theta, k1, m1)
      do j1=1, m*m
        diff(:,1) = matrix_positions(:,1) - matrix_positions(j1,1)
        diff(:,2) = matrix_positions(:,2) - matrix_positions(j1,2)
        vector_norm = diff(:,1)*diff(:,1) + diff(:,2)*diff(:,2)
        k1 = 0
        u_bar(j1) = 0
        v_bar(j1) = 0
        do m1=1, m*m
          if (sqrt(vector_norm(m1)) < d) then
            u_bar(j1) = u_bar(j1) + u(m1, i1-1)
            v_bar(j1) = v_bar(j1) + v(m1, i1-1)
            k1 = k1 + 1
          end if
        end do
        u_bar(j1) = u_bar(j1)/k1
        v_bar(j1) = v_bar(j1)/k1
        theta = atan(v_bar(j1)/u_bar(j1)) + R(j1)
        u(j1,i1) = s0*cos(theta)
        v(j1,i1) = s0*sin(theta)
      end do
      !$OMP end parallel do
      x(:,i1) = x(:,i1-1) + u(:,i1)
      y(:,i1) = y(:,i1-1) + v(:,i1)
      !$OMP parallel do
      do n1=1, m*m
        if ((abs(x(n1,i1)) - L) > 0) then
          if (x(n1,i1) > 0) then
            x(n1,i1) = x(n1,i1) - 2*L
          else
            x(n1,i1) = x(n1,i1) + 2*L
          end if
        end if
        if ((abs(y(n1,i1)) - L) > 0) then
          if (y(n1,i1) > 0) then
            y(n1,i1) = y(n1,i1) - 2*L
          else
            y(n1,i1) = y(n1,i1) + 2*L
          end if
        end if
      end do
      !$OMP end parallel do

      dummy = 0
      do k1=1,m*m
        dummy(1) = dummy(1) + u(k1,i1)
        dummy(2) = dummy(2) + v(k1,i1)
      end do
      psi(i1) = sqrt((dummy(1)**2)+(dummy(2)**2))/(m*m*s0)
    end do



  end subroutine simulate_omp





end module part2
