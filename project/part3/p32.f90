!Igor Adamski, CID 01069152
module wavemod
  use rwnetmod
  implicit none
  integer, allocatable, dimension(:,:) :: a !adjacency matrix
  integer, allocatable, dimension(:,:) :: e !edge list
  integer, allocatable, dimension(:,:) :: xnew !network node coordinates after isolated-node removal
  real(kind=8), allocatable, dimension(:) :: u !wave amplitude
  real(kind=8) :: c=1 !wave speed
	save

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

end module wavemod
!-------------------------------

program part3_mpi
    use mpi
    use wavemod
    implicit none
    integer :: i1,j1
    integer :: nt !number of time steps, number of oscillators
    real(kind=8) :: dt !time step
    integer :: myid, numprocs, ierr
    integer :: h,dim,nnodes,nlinks, ntotal !graph properties
    real(kind=8), allocatable, dimension(:) :: f0,f ! initial condition, frequencies, solution
    real(kind=8), allocatable, dimension(:) :: energy !order parameter
		integer :: nrand
		integer, allocatable :: seed(:)
		real(kind=8) :: test_rand

 ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!gather input
    open(unit=10,file='data.in')
        read(10,*) h
        read(10,*) dim
        read(10,*) ntotal
        read(10,*) nt
        read(10,*) dt
        read(10,*) c
   close(10)


		call random_seed(size = nrand)
		allocate(seed(nrand))
		seed = 0
		call random_seed(put = seed )
		call random_number(test_rand)
		print *, 'test_rand=',test_rand

    call wave_init(h,dim,ntotal,3,nnodes,nlinks) !set flag_init as appropriate
    allocate(f0(2*nnodes),f(2*nnodes),energy(nt+1))

		!set initial condition
		f0 = 0.d0
	  f0(nnodes) = 1
		!compute initial energy, energy(1) =
		print *, nnodes

		!compute solution
    call euler_mpi(MPI_COMM_WORLD,numprocs,nnodes,0.d0,f0,dt,nt,f,energy)


!output solution (after completion of gather in euler_mpi)
     call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      if (myid==0) then
        open(unit=11,file='f.dat')
        do i1=1,2*nnodes
            write(11,*) f(i1)
        end do
        close(11)

        open(unit=12,file='energy.dat')
        do i1=1,nt+1
	       write(12,*) energy(i1)
	      end do
      close(12)
      end if
    !can be loaded in python with: f=np.loadtxt('theta.dat')

    call MPI_FINALIZE(ierr)
end program part3_mpi



subroutine euler_mpi(comm,numprocs,n,t0,y0,dt,nt,y,energy)
    !explicit Euler method, parallelized with mpi
    !input:
    !comm: MPI communicator
    !numprocs: total number of processes
    !n: number of vertices
    !t0: initial time
    !y0: initial condition
    !dt: time step
    !nt: number of time steps
    !output: y, final solution
    !energy: energy at each time step
    use mpi
    use wavemod
    implicit none
    integer, intent (in) :: n,nt
    real(kind=8), dimension(2*n), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
    real(kind=8), dimension(2*n), intent(out) :: y
		real(kind=8), intent(out) :: energy(nt+1)
    real(kind=8) :: t
    integer :: i1,k,istart,iend,en
    integer :: comm,myid,ierr,numprocs,u_comm
		real(kind=8), allocatable :: ylocal(:),Rpart(:), u_part(:), u_total(:)
    integer, allocatable :: Nper_proc(:), disps(:)
		real(kind=8) :: dummy, y_temp(2*n)
		integer :: receiver, my_sender, my_receiver, sender, color, s
		integer, dimension(MPI_STATUS_SIZE) :: status


    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y)/2,numprocs,myid,istart,iend)
    print *, 'istart,iend,threadID=',istart,iend,myid, 'n=', n

    s = iend-istart+1
    allocate(ylocal(2*s), Rpart(2*s), u_part(s), u_total(n))
    allocate(disps(numprocs+1), Nper_proc(numprocs))
		ylocal(1:s) = y(istart:iend)
    ylocal(s+1:2*s) = y(istart+n:iend+n)

		!initial energy
		if (myid==0) energy(1) = 1
    !time marching
    do k = 1,nt


      !call RHS_mpi(istart,iend,t,y,Rpart)
      call RHS_mpi(comm, istart, iend, n, t, ylocal, Rpart)

      !march solution at nodes istart to iend forward in time

      ylocal= ylocal + dt*Rpart !ylocal must be declared and defined, Rpart must be declared, and
                                  !should be returned by RHS_mpi
			t = t + dt

      dummy = 0
      do i1=1,s
        dummy = dummy + ylocal(i1) * ylocal(i1)
      end do

      call MPI_REDUCE(dummy, energy(k+1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    end do

		if (myid == 0) print *, 'energy is=', energy

    print *, 'before collection',myid, maxval(abs(ylocal))


    call MPI_GATHER(size(ylocal)/2, 1, MPI_INT, Nper_proc, 1, MPI_INT, 0, comm, ierr)
      !collect ylocal from each processor onto myid=0

    if (myid==0) then
      disps(1) = 1
      do i1=1, numprocs
        disps(i1+1) = 1 + sum(Nper_proc(1:i1))
      end do
    end if

    !collect ylocal from each processor onto myid=0
	  call MPI_GATHERV(ylocal(1:s),size(ylocal)/2,MPI_DOUBLE_PRECISION,y_u,Nper_proc,disps,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_GATHERV(ylocal(s+1:2*s),size(ylocal)/2,MPI_DOUBLE_PRECISION,y_v,Nper_proc,disps,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (myid == 0) then
      y(1:s) = y_u
      y(s+1:2*s) = y_v
    end if

    if (myid==0) print *, 'finished', maxval(abs(y))


end subroutine euler_mpi
!-------------------------

subroutine RHS_mpi(comm,start,end,n,t,f,rhs)
    !called by euler_mpi
    !rhs = df/dt
    use wavemod
    use mpi
    implicit none
    integer, intent(in) :: start, end, n, comm
    real(kind=8), intent(in) :: t
!dimensions of variables below must be added
    real(kind=8), dimension(2*(end-start+1)), intent(in) :: f
    real(kind=8), dimension(size(f)), intent(out) :: rhs
		integer :: i1, j1, k1, m1, sender, glob_sender, s,ierr,myid
		real(kind=8) :: u_part(end-start+1), dummy
    integer :: curr_array(2)

    call MPI_COMM_RANK(comm, myid, ierr)

    s = end-start+1
    u_part = f(1:s)
    rhs(1:s) = f(s+1:2*s)
    j1 = 0
    m1 = 0
    do i1=1,n
      if (i1>= start .and. i1<=end) then
        j1 = j1 + 1
        curr_array(1) = f(j1)
        curr_array(2) = start+j1-1
        sender = myid
      else
        sender = 0
      end if
      call MPI_ALLREDUCE(sender, glob_sender, 1, MPI_INTEGER, MPI_MAX,comm, ierr)
      call MPI_BCAST(curr_array, 2, MPI_INTEGER, glob_sender, comm, ierr)
      !now everyone has the current ui
      dummy = 0
      m1 = 0
      do k1=start,end
        m1 = m1 + 1
        dummy = dummy + a(curr_array(2), k1)*(u_part(m1) - curr_array(1))
      end do
      dummy = (c**2)*dummy
      call MPI_REDUCE(dummy, rhs(j1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, glob_sender, comm, ierr)
      !print *, 'after reduction on iteration', i1, 'out of', n
    end do

end subroutine RHS_mpi


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
