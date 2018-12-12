!	Co-Array Implementation of Solution to 2D Poisson Equation through 
!	Finite Differencing and Jacobi Iterations
!	PROGRAMMER: RAHUL RATHNAKUMAR 
!	E-MAIL: RRATHNAK@ASU.EDU
!	DESCRIPTION: This is a parallel implementation of the program
!	assuming a unit square grid using Co-Array Fortran.
!---------------------------------------------------------------------------------------------
program main
	use, intrinsic :: iso_fortran_env
	use jacobiCAF
	implicit none
	!GLOBAL COMPUTATION VARIABLES
	integer, parameter ::m = 255 !DIMENSION OF THE GRID
	real(DOUBLE), target, allocatable :: u(:,:)
	real(DOUBLE), allocatable :: f(:,:), uanalytic(:,:), x(:), y(:), ghost_left(:)[:], ghost_right(:)[:], ghost_up(:)[:], ghost_down(:)[:], unew(:,:)
	real(DOUBLE) :: delta, bcn, bcs, bcw, bce, uinit, n_error[*], globmax, globvalerror, validation_error[*]
	!LOCAL VARIABLES (COMPUTATION)
	integer :: j,k, flag, itercount,q
	!CO-ARRAY VARIABLES
	integer :: image, nimages, npes
	!JOB MANAGEMENT VARIABLES
	integer :: extent
	integer :: left, right, up, down!Variables to let each process know the image numbers of its neighbour
	!---------------------------------------------------------------------------------------------
	!INITIALIZE IMAGES AND EQUIVALENT RANK:
	image=this_image()
	npes=16
	extent=floor(m/sqrt(real(npes))) !DEBUG: extent = 63 for 16 images
	!Allocate U, UANALYTIC, F
	allocate(u(0:extent,0:extent))
	allocate(uanalytic(0:extent, 0:extent))
	allocate(f(0:extent, 0:extent))
	allocate(x(0:extent))
	allocate(y(0:extent))
	allocate(ghost_down(0:extent)[*])
	allocate(ghost_up(0:extent)[*])
	allocate(ghost_left(0:extent)[*])
	allocate(ghost_right(0:extent)[*])
	allocate(unew(0:extent,0:extent))
	!INITIALIZE OTHER COMPUTATION VARIABLES
	itercount=1
	flag=0
	delta = 0.00390625
	!READ BOUNDARY CONDITIONS AND INITIAL GUESS VALUE
	bcn=0
	bcs=0
	bce=0
	bcw=0
	q=10
	uinit=0
	!CALCULATE LOCAL EXTENTS FOR EACH IMAGE 
	call domaindef(x, y, extent, image, npes, delta)
	call init_source_function(f, x, y, delta, extent, image, npes, m)	
	call init_problem(u, x, y, uinit, bce, bcn, bcs, bcw, extent, image, npes, m, delta)
	call defneighbour(x, y, left, right, up, down, extent, image, npes) 
	do while((flag .EQ. 0) .AND. (itercount .NE. 100000))
		call numeric(m, extent, image, npes, itercount, q, left, right, up, down, flag, n_error, f, u, unew, delta, ghost_up, ghost_down, ghost_left, ghost_right)
		sync all
		do j=1, npes-1
			globmax=max(n_error[j], n_error[j+1])
		enddo
		if((itercount .GT. q) .AND. (globmax .LT. 0.0001)) then
			flag=1
		else
			flag=0
		endif
		if(image .EQ. 1) then
			print*, "Iteration " , itercount, "completed."
		endif
		itercount=itercount+1
	enddo
	if(image .EQ. 1) then
		if(flag .EQ. 1) then
			print*, "Solution converged in ", itercount, "iterations."
			print*, "Error:", globmax
		else
			print*, "Solution failed to converge."
			print*, "Error:", globmax

		endif
	endif
	if((bce .EQ. 0) .AND. (bcw .EQ. 0) .AND. (bcs .EQ. 0) .AND. (bcn .EQ. 0)) then
		call validation(m, extent, uanalytic, validation_error, u, x, y, image)
		sync all
		do j=1, npes-1
			globvalerror=max(validation_error[j], validation_error[j+1])
		enddo
		if(image .EQ. 1) then
			print*, "Error wrt Analytical Solution:", globvalerror
		endif
	endif
end