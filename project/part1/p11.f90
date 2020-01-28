!Igor Adamski, CID 01069152.
module part1
  use cost
  implicit none
  integer :: itermax = 10000 !maximum number of iterations used by an optimizer
  real(kind=8) :: tol=1.0e-6 !stopping criteria for optimizer
  real(kind=8), allocatable :: jpath(:), xpath(:,:) !allocatable 1-d and 2-d arrays, should contain cost and location values used during optimization iterations
  real(kind=8), allocatable :: jpathf(:,:), xpathf(:,:,:)
  integer, allocatable, dimension(:) :: niter
  contains

  subroutine bracket_descent(xguess,xf,jf)
    !Use bracket descent method to minimize cost function, costj
    !input: xguess -- initial guess for loaction of minimum
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    !Assumes size(xguess) = 2
    implicit none
    real(kind=8), dimension(2), intent(in) :: xguess
    real(kind=8), intent(out) :: xf(2),jf !location of minimum, minimum cost
    integer :: i1,k1
    real(kind=8) :: j3(3),x3(3,2),j3temp(3),x3temp(3,2),xtemp(itermax+1,2),jtemp(itermax+1)
    real(kind=8) :: xm(2),dx(2),xr(2),jr,j_old
    logical :: convergence_flag

    !initialize triangle, costs
    call bd_initialize(xguess,x3,j3)
    call bd_sort(x3,j3,x3temp,j3temp)
    j3 = j3temp
    x3 = x3temp
    i1=0
    xtemp(i1+1,:) = sum(x3,dim=1)/3
    call costj(xtemp(i1+1,:),j_old)
    jtemp(i1+1) = j_old
    convergence_flag = .False.


    !Iterate using b-d scheme
    do while(i1<itermax .and. (.not. convergence_flag))
      !reflect
      xm = 0.5d0*(x3(1,:)+x3(2,:))
      dx = x3(3,:)-xm
      xf = xm-dx

      call costj(xf,jf)

      if (jf<j3(3)) then
        if (jf>j3(1)) then !accept new vertex
            j3(3) = jf
            x3(3,:) = xf
        else !extend

          xr = xm-2*dx
          call costj(xr,jr)
          if (jr<jf) then !accept extension
            j3(3) = jr
            x3(3,:) = xr
          else !accept original reflection
            j3(3) = jf
            x3(3,:) = xf
          end if
        end if
      else
        xf = xm - dx/2 !contract
        call costj(xf,jf)
        if (jf<j3(3)) then !accept contraction
          j3(3) = jf
          x3(3,:) = xf
        else !contract further
          xf = xm + dx/2
          call costj(xf,jf)
          if (jf<j3(3)) then !accept 2nd contraction
            j3(3) = jf
            x3(3,:) = xf
          else !shrink
            do k1=2,3
              x3(k1,:) = 0.5d0*(x3(1,:)+x3(k1,:))
              xf = x3(k1,:)
              call costj(xf,jf)
              j3(k1) = jf
            end do
          end if
        end if
      end if

      !compute centroid, sort vertices, check for convergence
      i1 = i1 + 1
      xtemp(i1+1,:) = sum(x3,dim=1)/3
      xf = xtemp(i1+1,:)
      call costj(xf,jf)
      jtemp(i1+1) = jf
      call bd_sort(x3,j3,x3temp,j3temp)
      x3 = x3temp
      j3 = j3temp
        call convergence_check(j_old,jf,convergence_flag)
      j_old = jf
    end do

    !set xpath, jpath
    if (allocated(xpath)) deallocate(xpath)
    if (allocated(jpath)) deallocate(jpath)
    allocate(xpath(i1,2),jpath(i1))
    xpath = xtemp(1:i1,:)
    jpath = jtemp(1:i1)
    xf = x3(1,:)
    jf = j3(1)
  end subroutine bracket_descent


  subroutine onestepbd(positions, x3, jf)
    implicit none
    real(kind=8), dimension(3,2), intent(in) :: positions
    real(kind=8), intent(out) :: jf, x3(3,2) !location of minimum, minimum cost
    real(kind=8) :: xf(2), j3(3)
    real(kind=8) :: j3temp(3),x3temp(3,2),xtemp(itermax+1,2),jtemp(itermax+1)
    real(kind=8) :: xm(2),dx(2),xr(2),jr,j_old
    integer :: r1, k1

    !initialize triangle, costs
    !call bd_initialize(xguess,x3,j3)
    x3 = positions
    call costj(positions(1,:), j3(1))
    call costj(positions(2,:), j3(2))
    call costj(positions(3,:), j3(3))
    call bd_sort(x3,j3,x3temp,j3temp)
    j3 = j3temp
    x3 = x3temp
    r1=0
    xtemp(r1+1,:) = sum(x3,dim=1)/3
    call costj(xtemp(r1+1,:),j_old)
    jtemp(r1+1) = j_old

    xm = 0.5d0*(x3(1,:)+x3(2,:))
    dx = x3(3,:)-xm
    xf = xm-dx

    call costj(xf,jf)

    if (jf<j3(3)) then
      if (jf>j3(1)) then !accept new vertex
          j3(3) = jf
          x3(3,:) = xf
      else !extend

        xr = xm-2*dx
        call costj(xr,jr)
        if (jr<jf) then !accept extension
          j3(3) = jr
          x3(3,:) = xr
        else !accept original reflection
          j3(3) = jf
          x3(3,:) = xf
        end if
      end if
    else
      xf = xm - dx/2 !contract
      call costj(xf,jf)
      if (jf<j3(3)) then !accept contraction
        j3(3) = jf
        x3(3,:) = xf
      else !contract further
        xf = xm + dx/2
        call costj(xf,jf)
        if (jf<j3(3)) then !accept 2nd contraction
          j3(3) = jf
          x3(3,:) = xf
        else !shrink
          do k1=2,3
            x3(k1,:) = 0.5d0*(x3(1,:)+x3(k1,:))
            xf = x3(k1,:)
            call costj(xf,jf)
            j3(k1) = jf
          end do
        end if
      end if
    end if

    xf = sum(x3,dim=1)/3
    call costj(xf,jf)

  end subroutine onestepbd


  subroutine bdglobal(xguess,nb,d,alpha,xf,jf)
    !Use bracket descent method to minimize cost function, costj
    !input: xguess -- initial guess for location of minimum
    !nb: use nb * nb triangles
    !d: grid spacing for initial nb x nb grid of triangles
    !alpha: agressiveness in centroid shift
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    !Assumes size(xguess) = 2
    implicit none
    real(kind=8), dimension(2), intent(in) :: xguess
    integer, intent(in) :: nb
    real(kind=8), intent(in) :: d, alpha
    real(kind=8), intent(out) :: xf(2),jf !location of minimum, minimum cost
    integer :: i1,j1,k1,m1,n1
    real(kind=8), dimension(2,nb*nb) :: xg, matrix_centroids, new_centroids
    logical, dimension(nb*nb) :: flags
    real(kind=8), dimension(nb*nb, 3, 2) :: matrix_positions, old_matrix_positions
    real(kind=8), dimension(nb*nb) :: vector_costs
    real(kind=8) :: curr_jstar, curr_xstar(2), curr_maxdenom, change_factor
    real(kind=8) :: epsilon, dummy(3,2), change(2), j_old, dummy2(3)

    !Set up initial grid of triangles
    k1 = 0
    do i1 = 1,nb
      do j1=1,nb
        k1 = k1 + 1
        xg(1,k1) = (j1-1)*d - (nb-1)*d/2 + xguess(1)
        xg(2,k1) = (i1-1)*d - (nb-1)*d/2 + xguess(2)
      end do
    end do

    !initialize variables
    epsilon = 1e-12
    matrix_centroids = xg

    do i1 = 1,nb*nb
      call bd_initialize(xg(:,i1), old_matrix_positions(i1,:,:), dummy2)
    end do

    i1 = 0

    if (allocated(xpath)) deallocate(xpath)
    if (allocated(jpath)) deallocate(jpath)

    allocate(xpath(itermax+1,2))
    allocate(jpath(itermax+1))

    xpath(1,:) = xguess
    call costj(xguess, jpath(1))

    !iterate all triangles
    do while (i1<itermax)
      !setup a matrix with positions and costs
      i1 = i1 + 1
      do m1 = 1, nb*nb
        if ((flags(m1) .eqv. .False.) .or. i1 == 1) then
          call onestepbd(old_matrix_positions(m1, :, :), matrix_positions(m1, :, :), vector_costs(m1))
        end if
      end do

      curr_jstar = minval(vector_costs)
      dummy = matrix_positions(minloc(vector_costs, dim=1), :, :)
      curr_xstar = sum(dummy, dim=1)/3
      curr_maxdenom = maxval(vector_costs - curr_jstar)

      !shift the vertices
      do n1 = 1, nb*nb
        if ((flags(n1) .eqv. .false.) .or. i1 == 1) then
          change_factor = (vector_costs(n1) - curr_jstar)/(curr_maxdenom + epsilon)
          change = change_factor * alpha * (sum(matrix_positions(n1,:,:),dim=1)/3 - curr_xstar)
          matrix_positions(n1,1,:) = matrix_positions(n1, 1, :) - change
          matrix_positions(n1,2,:) = matrix_positions(n1, 2, :) - change
          matrix_positions(n1,3,:) = matrix_positions(n1, 3, :) - change
          new_centroids(:,n1) = sum(matrix_positions(n1,:,:),dim=1)/3
          call costj(matrix_centroids(:,n1), j_old)
          call costj(new_centroids(:,n1), vector_costs(n1))
          call convergence_check(j_old, vector_costs(n1), flags(n1))
        end if
      end do
      matrix_centroids = new_centroids
      old_matrix_positions = matrix_positions
      xpath(i1+1,:) = new_centroids(:, minloc(vector_costs, dim=1))
      jpath(i1+1) = minval(vector_costs)
      if (all(flags) .eqv. .true.) exit
    end do
    print *, i1
    jf = jpath(i1+1)
    xf = xpath(i1+1,:)
    xpath = xpath(1:i1+1,:)
    jpath = jpath(1:i1+1)

  end subroutine bdglobal


  subroutine bd_initialize(xguess,x3,j3)
    !given xguess, generates vertices (x3) and corresponding costs (j3) for initial
    !bracket descent step
    implicit none
    real(kind=8), intent(in) :: xguess(2)
    real(kind=8), intent(out) :: j3(3),x3(3,2) !location of minimum
    integer :: i1
    real(kind=8), parameter :: l=1.d0

    x3(1,1) = xguess(1)
    x3(2,1) = xguess(1)+l*sqrt(3.d0)/2
    x3(3,1) = xguess(1)-l*sqrt(3.d0)/2
    x3(1,2) = xguess(2)+l
    x3(2,2) = xguess(2)-l/2
    x3(3,2) = xguess(2)-l/2

    do i1=1,3
      call costj(x3(i1,:),j3(i1))
    end do
  end subroutine bd_initialize

  subroutine bd_sort(x,c,xs,cs)
    !Three element insertion sort
    implicit none
    real(kind=8), intent(in) :: c(3),x(3,2)
    real(kind=8), intent(out) :: cs(3),xs(3,2)
    real(kind=8), dimension(2) :: temp

    cs = c
    xs = x
    if (c(2)<c(1)) then
      cs(1) = c(2)
      cs(2) = c(1)
      xs(1,:) = x(2,:)
      xs(2,:) = x(1,:)
    end if

    if (c(3)<cs(2)) then
      cs(3) = cs(2)
      cs(2) = c(3)
      xs(3,:) = xs(2,:)
      xs(2,:) = x(3,:)

      if (c(3)<cs(1)) then
        temp(1) = cs(1)
        cs(1) = cs(2)
        cs(2) = temp(1)
        temp = xs(1,:)
        xs(1,:) = xs(2,:)
        xs(2,:) = temp
      end if
    end if

  end subroutine bd_sort


  subroutine convergence_check(j1,j2,flag_converged)
    !check if costs j1 and j2 satisfy convergence criteria
    implicit none
    real(kind=8), intent(in) :: j1,j2
    real(kind=8) :: test
    logical, intent(out) :: flag_converged

    test = abs(j1-j2)/max(abs(j1),abs(j2),1.d0)
    if (test .le. tol) then
      flag_converged = .True.
    else
      flag_converged = .False.
    end if
  end subroutine convergence_check


end module part1
