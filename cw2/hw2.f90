!Module for solving n-d optimization problems with Newton's method and 2-d problems
!with bracket descent. Necessary cost function details are provided in separate cost
!module.
!IGOR ADAMSKI: 01069152
module hw2
  use cost
  implicit none
  integer :: itermax = 1000 !maximum number of iterations used by an optimizer
  real(kind=8) :: tol=1.0e-6 !stopping criteria for optimizer
  real(kind=8), allocatable :: jpath(:), xpath(:,:) !allocatable 1-d and 2-d arrays, should contain cost and location values used during optimization iterations

  contains

  !subroutine to solve linear equations
  subroutine lineq_sol(H, g, x_change)

      implicit none
      real(kind=8), dimension(:,:) :: H
      real(kind=8), dimension(:) :: g
      real(kind=8), dimension(:), intent(out) :: x_change
      integer :: N,NRHS,LDA,LDB,INFO
      integer, allocatable, dimension(:) :: IPIV
      !set dimension for H and g
      N = size(g)
      NRHS = 1
      LDA = N
      LDB = N
      allocate(IPIV(N))
      call dgesv(N, NRHS, H, LDA, IPIV, g, LDB, INFO)
      x_change = g(1:N)

  end subroutine lineq_sol

  subroutine newton(xguess,xf,jf)
    !Use Newton's method to minimize cost function, costj
    !input: xguess -- initial guess for loaction of minimum
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    implicit none
    real(kind=8), dimension(:), intent(in) :: xguess !do not need to explicitly specify dimension of input variable when subroutine is within a module
    real(kind=8), intent(out) :: xf(size(xguess)),jf !location of minimum, minimum cost
    real(kind=8) :: curr_hess(size(xguess),size(xguess)), curr_grad(size(xguess)), x_change(size(xguess))
    !real(kind=8), allocatable, dimension(:) :: x_change
    real(kind=8) :: j_last, j_new
    logical :: flag_converged
    integer i1

    !if xpath,jpath already taken, then discard them
    if (allocated(xpath)) deallocate(xpath)
    if (allocated(jpath)) deallocate(jpath)

    !allocate xpath,jpath
    allocate(xpath(itermax+1,2))
    allocate(jpath(itermax+1))
    xpath(1,:) = xguess
    call costj(xguess, j_last)
    jpath(1)=j_last
    do i1=1, itermax
      !get hess, grad and solve equation on current x
      call costj_hess2d(xpath(i1,:), curr_hess)

      call costj_grad2d(xpath(i1,:), curr_grad)

      call lineq_sol(curr_hess, -curr_grad, x_change)

      xpath(i1+1,:) = xpath(i1,:) + x_change
      call costj(xpath(i1+1,:), j_new)
      jpath(i1+1) = j_new
      !if below is true then stop as we converged
      call convergence_check(jpath(i1+1),jpath(i1),flag_converged)
      if (flag_converged) exit
    end do
    !get xpath to the right size (discard zeros)
    xpath = xpath(1:i1+1,:)
    jpath = jpath(1:i1+1)
    xf=xpath(i1+1,:)
    jf=j_new


  end subroutine newton

  !subroutine to arrange vertices in the traingle
  subroutine arrange_vertices(mat, v_a, v_b, v_c, Ja, Jb, Jc)
    implicit none
    real(kind=8), intent(in) :: mat(3,3)
    real(kind=8), intent(out) :: v_a(2), v_b(2), v_c(2), Ja, Jb, Jc
    real(kind=8) :: Jmax, Jmin
    integer i1
    integer, allocatable :: indeces(:,:)

    Jmax = max(mat(1,3), mat(2,3), mat(3,3))
    Jmin = min(mat(1,3), mat(2,3), mat(3,3))

    allocate(indeces(3,3))
    !initialize matrix of permutations
    indeces(1,:) = [1,2,3]
    indeces(2,:) = [2,3,1]
    indeces(3,:) = [3,1,2]
    !trust the developer that this works
    do i1=1, 3
      outer : if (mat(i1,3) .eq. Jmax) then
        v_a(1) = mat(i1,1)
        v_a(2) = mat(i1,2)
        Ja = mat(i1,3)
        inner1 : if (mat(indeces(2,i1),3) .eq. Jmin) then
          v_c(1) = mat(indeces(2,i1),1)
          v_c(2) = mat(indeces(2,i1),2)
          v_b(1) = mat(indeces(3,i1),1)
          v_b(2) = mat(indeces(3,i1),2)
          Jb = mat(indeces(3,i1),3)
          Jc = mat(indeces(2,i1),3)
        else inner1
          v_b(1) = mat(indeces(2,i1),1)
          v_b(2) = mat(indeces(2,i1),2)
          v_c(1) = mat(indeces(3,i1),1)
          v_c(2) = mat(indeces(3,i1),2)
          Jb = mat(indeces(2,i1),3)
          Jc = mat(indeces(3,i1),3)
        end if inner1

      end if outer
    end do
    deallocate(indeces)

  end subroutine arrange_vertices


  subroutine bracket_descent(xguess,xf,jf)
    !Use bracket descent method to minimize cost function, costj
    !input: xguess -- initial guess for loaction of minimum
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    !Assumes size(xguess) = 2
    implicit none
    real(kind=8), dimension(2), intent(in) :: xguess
    real(kind=8), intent(out) :: xf(2),jf !location of minimum, minimum cost
    real(kind=8) :: v_a(2), v_b(2), v_c(2)
    real(kind=8) :: xm(2), xstar(2), xdstar(2), Jstar, Jdstar, Ja_new, Jb_new, Jc_new
    real(kind=8) :: mat(3,3), dummy1, dummy2, Ja, Jb, Jc
    integer :: i1
    logical :: flag_converged

    !discard xpath,jpath if already allocated
    if (allocated(xpath)) deallocate(xpath)
    if (allocated(jpath)) deallocate(jpath)
    allocate(xpath(itermax+1,2))
    allocate(jpath(itermax+1))
    !finding our matrix at the first triangle
    mat(1,1) = xguess(1) - 0.5 * sqrt(3.d0)
    mat(1,2) = xguess(2) - 0.5
    mat(2,1) = xguess(1) + 0.5 * sqrt(3.d0)
    mat(2,2) = xguess(2) - 0.5
    mat(3,1) = xguess(1)
    mat(3,2) = xguess(2) + 1

    v_a = mat(1,1:2)
    v_b = mat(2,1:2)
    v_c = mat(3,1:2)
    xpath(1,:) = xguess
    call costj(xpath(1,:), jpath(1))
    !finalize the initial matrix
    call costj(v_a, mat(1,3))
    call costj(v_b, mat(2,3))
    call costj(v_c, mat(3,3))

    do i1=1, itermax
      !arrange vertices in matrix
      call arrange_vertices(mat, v_a, v_b, v_c, Ja, Jb, Jc)
      xm = (1.0/2.0) * (v_b + v_c)
      xstar = v_a + 2.0 * (xm - v_a)
      !evaluate cost at new point and call it Jstar
      call costj(xstar, Jstar)
      if (Jstar < min(Ja, Jb, Jc)) then
        xdstar = v_a + 3.0 * (xm - v_a)
        call costj(xdstar, Jdstar)
        inner1 : if (Jdstar < Jstar) then
          v_a = xdstar
        else inner1
          v_a = xstar
        end if inner1

      else if (Jstar >= Jc .and. Jstar <= Ja)  then
        v_a = xstar

      else if (Jstar > Ja) then
        xdstar = v_a + (3.0/2.0) * (xm - v_a)
        call costj(xdstar, Jdstar)
        inner2 : if (Jdstar < Ja) then
          v_a = xdstar
        else inner2
          xdstar = v_a + (1.0/2.0) * (xm - v_a)
          call costj(xdstar, Jdstar)
          deeper : if (Jdstar < Ja) then
            v_a = xdstar
          else deeper
            v_a = (1.0/2.0) * (v_c + v_a)
            v_b = (1.0/2.0) * (v_c + v_b)
          end if deeper
        end if inner2
      end if

      call costj(v_a, Ja_new)
      call costj(v_b, Jb_new)
      call costj(v_c, Jc_new)

      !fill the new matrix
      mat(1,:) = [v_a, Ja_new]
      mat(2,:) = [v_b, Jb_new]
      mat(3,:) = [v_c, Jc_new]

      !allocate xpath with the centroid
      xpath(i1+1,:) = (v_a + v_b + v_c)/3.0
      call costj(xpath(i1+1,:), jpath(i1+1))

      dummy1 = abs(Ja_new) + abs(Jb_new) + abs(Jc_new)
      dummy2 = abs(Ja) + abs(Jb) + abs(Jc)
      !check if we converged
      call convergence_check(dummy1, dummy2, flag_converged)
      if (flag_converged) exit
    end do
    xpath = xpath(1:i1+1,:)
    jpath = jpath(1:i1+1)
    xf = xpath(i1+1,:)
    call costj(xf,jf)

  end subroutine bracket_descent



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


end module hw2
