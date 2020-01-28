!Igor Adamski, CID 01069152
module part1_mpi
  use cost
  implicit none
  integer :: itermax = 10000 !maximum number of iterations used by an optimizer
  real(kind=8) :: tol=1.0e-6 !stopping criteria for optimizer
contains

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

    r1 = 1
    xf = sum(x3,dim=1)/3
    call costj(xf,jf)

  end subroutine onestepbd

  subroutine bdglobal_mpi(comm,numprocs,xguess,nb,d,alpha,xf,jf)
    !Use global bracket descent method to minimize cost function, costj
    !input: xguess -- initial guess for location of minimum
    !nb: use nb * nb triangles
    !d: grid spacing for initial nb x nb grid of triangles
    !alpha: agressiveness in centroid shift
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Assumes size(xguess) = 2

    !MY APPROACH TO PARALELLIZATION:

    !This code is paralellized with MPI, so a scheme in which the individual threads
    !do not share common memory. My approach to distributing this code over numprocs
    !threads starts by distributing the centroids (triangles) among the processes.
    !So if we have 2 processes, half of the centroid will belong to one process and
    !the other half to the other process. The decomposition is carried out by the
    !function MPE_DECOMP1D. Once each process knows which centroids it will iterate
    !I initialize the triangles (using bd_initialize) and all the processes are ready
    !to move to the main loop. The main loop is a while loop that terminates when the
    !number of iterations reaches maxiter(and when all flags converge but more below).
    !In the loop, each process moves each centroid assigned to it with one bracket
    !descent step, using the subroutine onestepbd. Then after computing the bd step
    !each process computes the minimum cost and puts it in the proc_jstar variable.
    !Then I use MPI_ALLREDUCE with MPI_MIN on all the proc_jstar's to find the
    !glob_jstar and all_reduce sends this result to all the processes. Then the
    !processes proceed to calculate the max (in the denominator for the change fraction)
    !and again Im using MPI_ALLREDUCE to get the max out of all processes. Then if
    !the glob_jstar == proc_jstar xstar is gathered and broadcasted to all the other
    !processes (tis if statement ensures that the correct xstar will be sent). Then
    !each process computes the shift for all the centroids (all which are not frozen)
    !towards the one with the lowest cost and each centroid calls convergence_check
    !at the end. Then after the shift I make a if statement: if all the flags (so
    !all the centroids have converged) are true then set a variable im_finished=1.
    !Then I reduce the im_finished variable with MPI_ALLREDUCE using MPI_SUM to get
    !the variable all_finished. And now if the variable all_finished is equal to
    !the number of processes we are using then all processes should quit the loop.
    !The reason why I cant quit the loop when a process converged immediately is because
    !the other processes will wait on the MPI_ALLREDUCE schemes for the process which
    !exited. This way when all the centroids in a process converge, this process will
    !stay in the loop untill all processes converge, but it will not change the position
    !of its triangles (because of the if statements in the do loops).
    !This approach is bound to be faster, but up to a point, because there is a lot of
    !communication between nodes and this could slow down the process.

    use mpi
    implicit none
    real(kind=8), dimension(2), intent(in) :: xguess
    integer, intent(in) :: comm,numprocs,nb
    real(kind=8), intent(in) :: d, alpha
    real(kind=8), intent(out) :: xf(2),jf !location of minimum, minimum cost
    integer :: myid, ierr
    integer :: istart, iend, nlocal
    integer :: i1,j1,k1,l1,m1
    real(kind=8), allocatable, dimension(:,:) :: xg, matrix_centroids, new_centroids
    real(kind=8), allocatable, dimension(:,:,:) :: matrix_positions, old_matrix_positions
    real(kind=8), allocatable, dimension(:) :: vector_costs
    logical, allocatable, dimension(:) :: flags
    integer :: im_finished, all_finished, send_id, glob_send_id, final_sender, glob_final_sender
    real(kind=8) :: proc_jstar, glob_jstar, curr_xstar(2), proc_maxdenom, glob_maxdenom, xf_proc(2)
    real(kind=8) :: jf_proc, change_factor, epsilon, dummy(3,2), change(2), j_old, j_new, dummy2(3)
    integer, dimension(MPI_STATUS_SIZE) :: status


    call MPI_COMM_RANK(comm, myid, ierr) !use comm instead of MPI_COMM_WORLD
    im_finished = 0
    print *, 'im process', myid, 'out of', numprocs, 'and im starting'

    !Create decomposition and initial condition
    call mpe_decomp1d(nb*nb,numprocs,myid,istart,iend)
    nlocal = iend - istart + 1
    !allocate variables
    allocate(xg(2,nlocal))
    allocate(matrix_positions(nlocal,3,2))
    allocate(old_matrix_positions(nlocal,3,2))
    allocate(matrix_centroids(2,nlocal))
    allocate(flags(nlocal))
    allocate(vector_costs(nlocal))
    allocate(new_centroids(2,nlocal))

    l1 = 0
    k1 = 0
    do i1 = 1,nb
      do j1=1,nb
        k1 = k1 + 1
        if (k1>=istart .and. k1<=iend) then !check if triangle i1,j1 "belongs" to this process
          l1 = l1 +1
          xg(1,l1) = (j1-1)*d - (nb-1)*d/2 + xguess(1)
          xg(2,l1) = (i1-1)*d - (nb-1)*d/2 + xguess(2)
        end if
      end do
    end do

    epsilon = 1e-12
    i1 = 0
    l1 = 0
    k1 = 0
    matrix_centroids = xg
    do i1 = 1,nlocal
      call bd_initialize(xg(:,i1), old_matrix_positions(i1,:,:), dummy2)
    end do

    i1 = 0


    do while (i1<itermax)
      i1 = i1 + 1
      k1 = 0
      l1 = 0
      do m1=istart,iend
        k1 = k1 + 1
        if ((i1 == 1) .or. (flags(k1) .eqv. .false.)) then
          call onestepbd(old_matrix_positions(k1,:,:), matrix_positions(k1, :, :), vector_costs(k1))
        end if
      end do

      proc_jstar = minval(vector_costs)
      call MPI_ALLREDUCE(proc_jstar, glob_jstar, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)

      proc_maxdenom = maxval(vector_costs - glob_jstar)
      call MPI_ALLREDUCE(proc_maxdenom, glob_maxdenom, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      if (proc_jstar == glob_jstar) then
        send_id = myid
      else
        send_id = 0
      end if

      call MPI_ALLREDUCE(send_id, glob_send_id, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
      curr_xstar = sum(matrix_positions(minloc(vector_costs, dim=1), :, :), dim=1)/3
      call MPI_BCAST(curr_xstar, 2, MPI_DOUBLE_PRECISION, glob_send_id, comm, ierr)

      do m1=istart,iend
        l1 = l1 + 1
        if ((flags(l1) .eqv. .false.) .or. i1 == 1) then
          change_factor = (vector_costs(l1) - glob_jstar)/(glob_maxdenom + epsilon)
          change = change_factor * alpha * (sum(matrix_positions(l1,:,:),dim=1)/3 - curr_xstar)
          matrix_positions(l1,1,:) = matrix_positions(l1, 1, :) - change
          matrix_positions(l1,2,:) = matrix_positions(l1, 2, :) - change
          matrix_positions(l1,3,:) = matrix_positions(l1, 3, :) - change
          new_centroids(:,l1) = sum(matrix_positions(l1,:,:),dim=1)/3
          call costj(matrix_centroids(:,l1), j_old)
          call costj(new_centroids(:,l1), vector_costs(l1))
          call convergence_check(j_old, vector_costs(l1), flags(l1))
        end if
      end do
      jf_proc = minval(vector_costs)
      xf_proc = new_centroids(:, minloc(vector_costs,dim=1))
      matrix_centroids = new_centroids
      old_matrix_positions = matrix_positions

      if (all(flags) .eqv. .true.) im_finished = 1
      call MPI_ALLREDUCE(im_finished, all_finished, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      if (all_finished == numprocs) exit

    end do
    print *, 'it took i1=', i1

    call MPI_REDUCE(jf_proc, jf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    if (jf == jf_proc) then
      final_sender = myid
    else
      final_sender = 0
    end if
    call MPI_ALLREDUCE(final_sender, glob_final_sender, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (myid == glob_final_sender) call MPI_SEND(xf_proc, 2, MPI_DOUBLE_PRECISION, 0, 1, comm, ierr)
    if (myid == 0) call MPI_RECV(xf,2,MPI_DOUBLE_PRECISION,glob_final_sender,1,comm,status,ierr)

    call MPI_BARRIER(comm, ierr)
  end subroutine bdglobal_mpi

!--------------------------------------------------------------------
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in online MPE documentation.
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
  subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
    implicit none
    integer :: n, numprocs, myid, s, e
    integer :: nlocal
    integer :: deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
        nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

  end subroutine MPE_DECOMP1D

end module part1_mpi

program p12
  !main program can be used to call bdglobal
  use mpi
  use part1_mpi
  implicit none
  integer :: nb
  integer :: myid,ierr,numprocs
  real(kind=8) :: xguess(2),d,alpha
  real(kind=8) :: xf(2),jf

  ! Initialize MPI
     call MPI_INIT(ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

     open(unit=11,file='data.in')
     read(11,*) xguess(1)
     read(11,*) xguess(2)
     read(11,*) nb
     read(11,*) d
     read(11,*) alpha
     close(11)

     call bdglobal_mpi(MPI_COMM_WORLD,numprocs,xguess,nb,d,alpha,xf,jf)

     !output xf jf to file
     call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
     if (myid==0) then
       open(unit=12,file='p12.dat')
       write(12,*) xf(1)
       write(12,*) xf(2)
       write(12,*) jf
       close(12)
     end if

     call MPI_FINALIZE(ierr)


end program p12
