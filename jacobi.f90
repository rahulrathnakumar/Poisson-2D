module jacobi
!THIS MODULE CONTAINS ALL THE SUBROUTINES NEEDED FOR THE JACOBI COMPUTATION.
!DESCRIPTION FOR EACH SUBROUTINE IS GIVEN ALONGISDE THE SUBROUTINE.
use mpi
implicit none
integer, parameter :: DOUBLE=kind(1.0d0)
integer, parameter :: ROOT=0
integer, parameter :: INTERNAL_ERROR=1  ! shell error code; must be nonzero
integer, parameter :: NONE=-1
integer :: ierr
integer status(MPI_STATUS_SIZE)
contains
subroutine defneighbour(x, y, left, right, up, down, extent, rank, npes) 
!DESCRIPTION:
!Subroutine takes in the rank and uses it to identify the processors' neighbours.
!------------------------------------------------------------------------------------
	integer, intent(inout) :: left, right, up, down
	integer, intent(in) :: rank, extent, npes
	real(DOUBLE), intent(in) :: x(0:extent), y(0:extent)
	integer :: j
	if(mod(rank,4) .EQ. 0) then 
		left=NONE
	else
		left=rank-1
	endif
	if(mod(rank,4) .EQ. int(sqrt(real(npes))-1)) then
		right=NONE
	else
		right=rank+1
	endif
	if((rank .LT. (sqrt(real(npes))-1))) then
		down=NONE
	else
		down=rank - sqrt(real(npes))
	endif
	if((rank .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then
		up=NONE
	else
		up=rank + sqrt(real(npes))
	endif
	return
end subroutine defneighbour
		
subroutine domaindef(x, y, extent, rank, npes, delta)
!	DESCRIPTION:
!	Defines the local subdomains for each process
!-----------------------------------------------------------------------
	integer, intent(in) :: extent, rank, npes
	real(DOUBLE), intent(in) :: delta
	real(DOUBLE), intent(inout) :: x(0:extent), y(0:extent)
	integer :: rowcount, colcount !VARIABLES THAT KNOW WHAT ROW & COLUMN THE PROCESSOR IS AT IN RELATION TO THE DECOMPOSED SYSTEM OF BLOCKS
	integer :: j
	rowcount = rank/sqrt(real(npes))
	colcount = mod(rank,int(sqrt(real(npes))))
	do j=0, extent
		x(j)=(delta*colcount*(extent+1)) + (delta*j)
		y(j)=(delta*rowcount*(extent+1)) + (delta*j)
!		print*, x(j), y(j)
	enddo	
	return
end subroutine domaindef

subroutine init_source_function(f, x, y, delta, extent, rank, npes, m)
!	DESCRIPTION:
!	Defines the forcing function for each sub-domain, process-wise.
!-------------------------------------------------------------------------
	integer, intent(in) :: extent, rank, npes, m
	real(DOUBLE), intent(in) :: delta, x(0:extent), y(0:extent)
	real(DOUBLE), intent(inout) :: f(0:extent, 0:extent)
	!LOCAL VARIABLES
	integer :: j,k,l
	!SOURCE FUNCTION IS DEFINED ONLY IN THE GLOBAL INTERIOR POINTS
	do j=0, extent
		do k=0, extent
			if((x(j) .NE. 0) .AND. (y(j) .NE. 0) .AND. (x(j) .NE. (m*delta)) .AND. (y(j) .NE. (m*delta))) then !CONDITION CHECKS IF THE COORDINATE IS A PROCESSOR BOUNDARY
				f(j,k)=exp(x(j)+y(k))*((((x(j)**2)+(3*x(j)))*((y(k)**2)-y(k)))+ &
				& (((y(k)**2)+(3*y(k)))*((x(j)**2)-x(j))))
			endif
		enddo
	enddo
	return
end subroutine init_source_function

subroutine init_problem(u, x, y, uinit, bce, bcn, bcs, bcw, extent, rank, npes, m, delta)
! DESCRIPTION:
! This subroutine defines the boundary conditions and initial guess value for
! the problem. The matrix 'u' is of size [(extent+1)x(extent+1)].
!--------------------------------------------------------------------------
	integer, intent(in) :: rank, extent, npes, m
	real(DOUBLE), intent(inout) :: u(0:extent,0:extent)
	real(DOUBLE), intent(in) :: x(0:extent), y(0:extent)
	real(DOUBLE), intent(in) :: bce, bcn, bcs, bcw, uinit, delta
	integer :: j,k
	!IDENTIFY AND FILL UP GLOBAL BOUNDARIES
	if((rank .LT. (sqrt(real(npes))-1))) then !CONDITION CHECKS FOR THE BOTTOM PROCESSES
		u(:,0) = bcs
	else if((mod(rank,4) .EQ. 0)) then !CONDITION CHECKS FOR ALL THE PROCESSES AT THE LEFT SIDE
		u(0,:) = bcw
	else if((mod(rank,4) .EQ. (sqrt(real(npes))-1))) then !CONDITION CHECKS FOR ALL THE PROCESSES ON THE RIGHT SIDE 
		u(extent,:) = bce
	else if((rank .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then!CONDITION CHECKS FOR ALL THE PROCESSES ON THE TOP SIDE
		u(:,extent) = bcn
	endif
	!FILL OUT THE INTERIOR POINTS FOR EACH PROCESS WITH THE GUESS VALUE, IE USE THE GUESS VALUE FOR ALL 
	!THE COORDINATES THAT ARE NOT PROCESSOR BOUNDARIES
	do j=0, extent
		do k=0, extent
			if((x(j) .NE. 0) .AND. (y(j) .NE. 0) .AND. (x(j) .NE. (m*delta)) .AND. (y(j) .NE. (m*delta))) then !CONDITION CHECKS IF THE COORDINATE IS A GLOBAL BOUNDARY
				u(j,k)=uinit
			endif
		enddo
	enddo		
	return	
end subroutine init_problem

subroutine numeric(m, extent, rank, npes, itercount, q, left, right, up, down, flag, n_error, f, u, delta)
!	DESCRIPTION:
!	This subroutine performs the jacobi iteration in parallel for each process separately. Each process has 
!	its own interior and boundary points. Boundary points of every process needs to exchange data
!	with its neighbours. The data to be exchanged is the data from the previous iteration. Each process, therefore
!	has to be in sync at the end of each iteration, which is taken care of by the MPI_BARRIER call in the main
!	routine. 
!	Each processor gets assigned set of "ghost" points from its neighbor so that it can readily utilize it for the 
!	jacobi iterations. 
!----------------------------------------------------------------------------------------------------------------
	use mpi
    	implicit none
	integer, parameter :: DOUBLE = kind(1.0d0)
	integer, intent(in) :: m, extent, rank, npes, itercount ,q, left, right, up, down
	integer, intent(inout) :: flag
	real(DOUBLE), intent(inout) :: n_error
	real(DOUBLE), intent(in) :: f(0:extent, 0:extent), delta
	real(DOUBLE), target, intent(inout) :: u(0:extent,0:extent)
	real(DOUBLE),pointer:: north(:,:),east(:,:),west(:,:), south(:,:)
	real(DOUBLE) :: ghost_left(0:extent), ghost_right(0:extent), ghost_up(0:extent), ghost_down(0:extent)
	real(DOUBLE) :: unew(0:extent,0:extent)
	integer :: j, k, tag_right, tag_left, tag_up, tag_down, input
!	SET UP TAGS FOR MESSAGE PASSING
	tag_right=1
	tag_left=2
	tag_up=3
	tag_down=4
!	DEFINE THE LOCAL BOUNDARIES OF EACH PROCESS AS THE GHOST VECTORS
	ghost_left=u(0,:)
	ghost_right=u(extent,:)
	ghost_up=u(:,extent)
	ghost_down=u(:,0)
!	MPI COMMUNICATIONS BEGIN-USING SENDRECV EXHCHANGES
	if(mod(rank,4) .EQ. 0) then !PROCESESS IN THE LEFT BDRY
		call mpi_sendrecv(ghost_right, extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_right, ghost_left,extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_left, MPI_COMM_WORLD, status, ierr)
		if(rank .EQ. 0) then !PROCESS NUMBER 0, IE, PROC AT THE BOTTOM-LEFT CORNER
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
		else if(rank .EQ. (sqrt(real(npes))*(sqrt(real(npes))-1))) then !PROC AT THE TOP-LEFT CORNER
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
		else !OTHER PROCESSES IN THE LEFT BDRY
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
		endif
	else if(mod(rank,4) .EQ. (sqrt(real(npes))-1)) then !PROCESESS IN THE RIGHT BDRY
		call mpi_sendrecv(ghost_left, extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_left, ghost_right,extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_right, MPI_COMM_WORLD, status, ierr)
		if(rank .EQ. (sqrt(real(npes))-1)) then !PROCESS ON THE BOTTOM-RIGHT CORNER
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
		else if(rank .EQ. (npes-1)) then!PROCESS ON THE TOP-RIGHT CORNER
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
		else !OTHER PROCESESS IN THE RIGHT BDRY
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
		endif
	else !PROCESSES NOT IN THE LEFT AND RIGHT BDRY
		if(rank .LT. sqrt(real(npes)) -1) then !PROCESSES AT THE BOTTOM BDRY EXCL LEFT AND RIGHT END
			call mpi_sendrecv(ghost_left, extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_left, ghost_right,extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_right, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_right, extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_right, ghost_left,extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_left, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
		else if((rank .GT. (sqrt(real(npes))*(sqrt(real(npes))-1))) .AND. (rank .LT. (npes-1))) then !PROCESSES AT THE TOP BDRY EXCL LEFT AND RIGHT END
			call mpi_sendrecv(ghost_left, extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_left, ghost_right,extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_right, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_right, extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_right, ghost_left,extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_left, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
		else !ALL OTHER INTERIOR PROCESSES
			call mpi_sendrecv(ghost_left, extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_left, ghost_right,extent+1, MPI_DOUBLE_PRECISION, (rank-1), tag_right, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_right, extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_right, ghost_left,extent+1, MPI_DOUBLE_PRECISION, (rank+1), tag_left, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_down, extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_down, ghost_up,extent+1, MPI_DOUBLE_PRECISION, (rank-4), tag_up, MPI_COMM_WORLD, status, ierr)
			call mpi_sendrecv(ghost_up, extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_up, ghost_down,extent+1, MPI_DOUBLE_PRECISION, (rank+4), tag_down, MPI_COMM_WORLD, status, ierr)
		endif
	endif
!	CONSIDER WAIT FOR ALL PROCESSES TO COMPLETE SENDING AND RECIEVING THEIR GHOST BOUNDARIES
!	call mpi_barrier(MPI_COMM_WORLD, ierr) 
!*************************************************************************************************
!	JACOBI COMPUTATION SECTION
	north => u(1:extent-1,2:extent) !north(1,1)=u(1,2) 
	south => u(1:extent-1,0:extent-2) !south(1,1)=u(1,0) 
	west => u(0:extent-2,1:extent-1) !west(1,1)=u(0,1) 
	east => u(2:extent,1:extent-1) !east(1,1)=u(2,1) 
	if(rank .LT. sqrt(real(npes))) then !CONDITION TO CHECK FOR BOTTOM PROCESSES
		if(mod(rank,4) .EQ. 0) then !BOTTOM LEFT PROCESS
			do j=1, extent
				do k=1, extent
					if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT CORNER
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else if((k .EQ. extent) .AND. (j .LT. extent)) then !TOP BDRY EXCL TOP RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .LT. extent)) then! RIGHT BDRY EXCL TOP-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))						
					else !INTERIOR PTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		else if(mod(rank,4) .EQ. (sqrt(real(npes))-1)) then !BOTTOM-RIGHT PROCESS
			do j=0, (extent-1)
				do k=1, extent
					if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .LT. extent)) then !LEFT END EXCL TOP LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((k .EQ. extent) .AND. (j .GT. 0)) then !TOP END EXCL TOP LEFT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else !INTERIOR PTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		else !ALL OTHER PROCESSES IN THE BOTTOM ROW
			do j=0, extent
				do k=1, extent
					if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .LT. extent)) then!LEFT BOUNDARY EXCL TOP LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT POINT
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .LT. extent)) then !RIGHT BOUNDARY EXCL TOP RIGHT POINT
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. extent)) then !TOP BOUNDARY EXCL TOP-LEFT AND TOP-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
					else !INTERIOR POINTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		endif
	else if((rank .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !CONDITION CHECKS FOR TOP PROCESSES
		if(mod(rank,4) .EQ. 0) then !TOP-LEFT PROCESS
			do j=1, extent
				do k=0, (extent-1)
					if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+ghost_left(k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BOUNDARY EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .GT. 0)) then !RIGHT BDRY PTS EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+ghost_left(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo					
		else if(rank .EQ. (npes-1)) then !TOP-RIGHT PROCESS
			do j=0, (extent-1)
				do k=0, (extent-1)
					if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (k .EQ. 0)) then !BOTTOM BOUNDARY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .GT. 0)) then !LEFT BNDRY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo	
		else !OTHER PROCESESS ON TOP SIDE
			do j=0, extent
				do k= 0, (extent-1)
					if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. 0)) then!BOTTOM BDRY EXCL BOTTOM-LEFT & BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+ghost_up(k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .GT. 0)) then !LEFT BDRY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .GT. 0)) then !RIGHT BDRY EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		endif
	else if((rank .GT. (sqrt(real(npes))-1)) .AND. (mod(rank, 4) .EQ. 0) .AND. (rank .LT. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !LEFT BOUNDARY PROCESSES EXCL TOP & BOTTOM PROCS
		do j=1, extent
			do k=0, extent
				if((j .EQ.  extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BDRY EXCL BOTTOM-RIGHT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .LT. extent) .AND. (k .GT. 0)) then !RIGHT BDRY EXCL TOP-RIGHT PT AND BOTTOM-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .LT. extent) .AND. (k .EQ. extent)) then !TOP BDRY EXCL TOP-RIGHT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else !INTERIOR PTS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				endif
			enddo
		enddo
	else if((rank .GT. (sqrt(real(npes))-1)) .AND. (mod(rank, 4) .EQ. (sqrt(real(npes))-1)) .AND. (rank .LT. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !RIGHT BDRY PROCS EXCL TOP & BOTTOM PROCS
		do j=0, extent - 1
			do k=0, extent
				if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !LEFT BDRY PTS EXCL BOTTOM-LEFT PT AND TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (k .EQ. extent)) then !TOP BDRY PTS EXCL TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (k .EQ. 0))then !BOTTOM BDRY PTS EXCL BOTTOM-LEFT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else !INTERIOR PTS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				endif
			enddo
		enddo
	else !INTERIOR PROCESSES
		do j=0,extent
			do k=0, extent
				if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT CORNER
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT CORNER
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BDRY EXCL CORNERS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT CORNER
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT CORNER
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. extent)) then !TOP BDRY EXCL CORNERS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !LEFT BDRY EXL CORNERS
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !RIGHT BDRY EXCL CORNERS
					unew(j,k)=0.25*(ghost_left(k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				endif
			enddo
		enddo
	endif
!	UPDATE BOUNDARIES OF UNEW TO U FOR COMPARISON PURPOSES
	if(left .EQ. NONE) then
		unew(0,:) = u(0,:)
		if(down .EQ. NONE) then
			unew(:,0) = u(:,0)
		else if(up .EQ. NONE) then
			unew(:,extent)=u(:,extent)
		endif
	else if(right .EQ. NONE) then 
		unew(extent, :) = u(extent, :)
		if(down .EQ. NONE) then
			unew(:,0)=u(:,0)
		else if(up .EQ. NONE) then 
			unew(:,extent)=u(:,extent)
		endif
	else if(down .EQ. NONE) then 
		unew(:,0) = u(:,0)
	else if(up .EQ. NONE) then 
		unew(:,extent) = u(:,extent)
	endif
!	END OF JACOBI COMPUTATION	
!*************************************************************************************************
!	CONVERGENCE CHECK - LOCAL MAXIMA
	if(mod(itercount,q) .EQ. 0) then
		n_error=maxval(abs(unew-u))
	endif
	u=unew
	return
end subroutine numeric
subroutine validation(m, extent, uanalytic, validation_error, u, x, y)
integer, intent(in) :: m, extent
real(DOUBLE), intent(inout) :: uanalytic(0:extent,0:extent), validation_error
real(DOUBLE), intent (in) :: u(0:extent,0:extent), x(0:extent),y(0:extent)
integer :: j,k
do j=0 , extent 
	do k=0 , extent
 		uanalytic(j,k)=exp(x(j)+y(k))*(((x(j)**2)-x(j))*((y(k)**2)-y(k)))
	enddo
enddo
validation_error = maxval(abs(u-uanalytic))
return
end subroutine validation
end module jacobi
