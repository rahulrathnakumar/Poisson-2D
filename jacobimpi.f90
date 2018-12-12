!	MPI Implementation of Solution to 2D Poisson Equation through 
!	Finite Differencing and Jacobi Iterations
!	PROGRAMMER: RAHUL RATHNAKUMAR 
!	E-MAIL: RRATHNAK@ASU.EDU
!	REVISION:
!	v1.0 - 11/28/18 : Current Version
!	DESCRIPTION: This is a parallel implementation of the program
!	assuming a unit square grid using MPI. The communication topology is bi-directional with a 4x4 
!	process arrangement that can be extended to 8x8, 16x16, etc by making a few modifications to the source code.
!---------------------------------------------------------------------------------------------
program main
	use, intrinsic :: iso_fortran_env
	use mpi
	use jacobi
	implicit none
	!GLOBAL COMPUTATION VARIABLES
	integer, parameter ::m = 255 !GRID RESOLUTION
	real(DOUBLE), target, allocatable :: u(:,:)
	real(DOUBLE), allocatable:: f(:,:), uanalytic(:,:), x(:), y(:)
	real(DOUBLE) :: delta, bcn, bcs, bcw, bce, uinit, n_error, validation_error, globmax,globvalerror
	!LOCAL VARIABLES (COMPUTATION)
	integer :: j,k, flag, itercount,q
	!MPI & SHELL VARIABLES
	integer :: rank
	integer:: npes  ! Number of Processing Elements (Has to be >0 and powers of 4)
	!MPI JOB MANAGEMENT VARIABLES
	integer :: extent
	integer :: left, right, up, down!Variables to let each process know the process numbers of its neighbour
	!---------------------------------------------------------------------------------------------
	!Initialize MPI
	call mpi_init(ierr)  ! should be first executable statement
	call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)  ! my rank
	call mpi_comm_size(MPI_COMM_WORLD, npes, ierr) ! number of processes 
	if((npes .NE. 16)) then !MODIFY THIS TO 4,16,64, ETC DEPENDING ON YOUR REQUIREMENT
		print*, "Error in number of processes!" 
		goto 911
	endif
	extent=floor(m/sqrt(real(npes))) !LOCAL EXTENTS FOR 16 PROCESSES
	!Allocate U, UANALYTIC, F
	allocate(u(0:extent,0:extent))
	allocate(uanalytic(0:extent, 0:extent))
	allocate(f(0:extent, 0:extent))
	allocate(x(0:extent))
	allocate(y(0:extent))
	!INITIALIZE OTHER COMPUTATION VARIABLES
	itercount=1
	flag=0
	delta = real(1/m)
	!READ BOUNDARY CONDITIONS AND INITIAL GUESS VALUE 
	! BCN - NORTH END BOUNDARY
	! BCS - SOUTH END BOUNDARY
	! BCE - EAST END BOUNDARY
	! BCW - WEST END BOUNDARY
	! q - CONVERGENCE CHECK FREQUENCY
	! uinit - INITIAL GUESS
	bcn=0
	bcs=0
	bce=0
	bcw=0
	q=10
	uinit=0
	!CALCULATE LOCAL EXTENTS FOR EACH PROCESS 
	call domaindef(x, y, extent, rank, npes, delta)	!DEFINES THE SPATIAL COORDINATES FOR EACH PROCESS
	call init_source_function(f, x, y, delta, extent, rank, npes, m)	!SOURCE FUNCTION IS GENERATED FOR EACH PROCESS LOCALLY
	call init_problem(u, x, y, uinit, bce, bcn, bcs, bcw, extent, rank, npes, m, delta) !U IS GENERATED FOR EACH PROCESS LOCALLY
	call defneighbour(x, y, left, right, up, down, extent, rank, npes) 
	do while((flag .EQ. 0) .AND. (itercount .NE. 100000))
		if(rank .EQ. ROOT) then
			print*, "Iteration: ", itercount
		endif
		call numeric(m, extent, rank, npes, itercount, q, left, right, up, down, flag, n_error, f, u, delta)
		call mpi_allreduce(n_error, globmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
		if((itercount .GT. q) .AND. (globmax .LT. 0.0001)) then
			flag=1
		else
			flag=0
		endif
		if(rank .EQ. ROOT) then
			print*, "Iteration " , itercount, "completed."
		endif
		itercount=itercount+1
	enddo
	if(rank .EQ. ROOT) then
		if(flag .EQ. 1) then
			print*, "Solution converged in ", itercount, "iterations."
			print*, "Error:", globmax
		else
			print*, "Solution failed to converge."
			print*, "Error:", globmax

		endif
	endif
	if((bce .EQ. 0) .AND. (bcw .EQ. 0) .AND. (bcs .EQ. 0) .AND. (bcn .EQ. 0)) then
		call validation(m, extent, uanalytic, validation_error, u, x, y)
		call mpi_allreduce(validation_error, globvalerror, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
		if(rank .EQ. ROOT) then
			print*, "Error wrt Analytical Solution:", globvalerror
		endif
	endif
	911 continue
	call mpi_finalize(ierr)  
	stop
end
