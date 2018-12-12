module jacobiCAF
!THIS MODULE CONTAINS ALL THE SUBROUTINES NEEDED FOR THE JACOBI COMPUTATION.
!DESCRIPTION FOR EACH SUBROUTINE IS GIVEN ALONGISDE THE SUBROUTINE.
implicit none
integer, parameter :: DOUBLE=kind(1.0d0)
integer, parameter :: ROOT=0
integer, parameter :: INTERNAL_ERROR=1  ! shell error code; must be nonzero
integer, parameter :: NONE=-1
contains
subroutine defneighbour(x, y, left, right, up, down, extent, image, npes) 
!DESCRIPTION:
!Subroutine takes in the image number and uses it to identify the image neighbours.
!------------------------------------------------------------------------------------
	integer, intent(inout) :: left, right, up, down, image
	integer, intent(in) :: extent, npes
	real(DOUBLE), intent(in) :: x(0:extent), y(0:extent)
	integer :: j
	if(mod((image-1),4) .EQ. 0) then 
		left=NONE
	else
		left=(image-1)
	endif
	if(mod((image-1),4) .EQ. int(sqrt(real(npes))-1)) then
		right=NONE
	else
		right=(image+1)
	endif
	if((((image-1)) .LT. (sqrt(real(npes))-1))) then
		down=NONE
	else
		down=image - sqrt(real(npes))
	endif
	if((((image-1)) .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then
		up=NONE
	else
		up=image + sqrt(real(npes))
	endif
	return
end subroutine defneighbour

subroutine domaindef(x, y, extent, image,npes, delta)
!	DESCRIPTION:
!	Each image gets its own portion of the 255x255 grid in this subroutine
!-------------------------------------------------------------------------------
	integer, intent(in) :: extent, npes
	integer, intent(inout) :: image
	real(DOUBLE), intent(in) :: delta
	real(DOUBLE), intent(inout) :: x(0:extent), y(0:extent)
	integer :: rowcount, colcount !VARIABLES THAT KNOW WHAT ROW & COLUMN THE IMAGEOR IS AT IN RELATION TO THE DECOMPOSED SYSTEM OF BLOCKS
	integer :: j
	rowcount = ((image-1))/sqrt(real(npes))
	colcount = mod((image-1),int(sqrt(real(npes))))
	do j=0, extent
		x(j)=(delta*colcount*(extent+1)) + (delta*j)
		y(j)=(delta*rowcount*(extent+1)) + (delta*j)
!		print*, x(j), y(j)
	enddo	
	return
end subroutine domaindef

subroutine init_source_function(f, x, y, delta, extent, image, npes, m)
!	DESCRIPTION:
!	Source function is created image-wise
!--------------------------------------------------------------------
	integer, intent(in) :: extent, npes, m
	integer, intent(inout) :: image
	real(DOUBLE), intent(in) :: delta, x(0:extent), y(0:extent)
	real(DOUBLE), intent(inout) :: f(0:extent, 0:extent)
	!LOCAL VARIABLES
	integer :: j,k,l
	!SOURCE FUNCTION IS DEFINED ONLY IN THE GLOBAL INTERIOR POINTS
	do j=0, extent
		do k=0, extent
			if((x(j) .NE. 0) .AND. (y(j) .NE. 0) .AND. (x(j) .NE. (m*delta)) .AND. (y(j) .NE. (m*delta))) then !CONDITION CHECKS IF THE COORDINATE IS A IMAGEOR BOUNDARY
				f(j,k)=exp(x(j)+y(k))*((((x(j)**2)+(3*x(j)))*((y(k)**2)-y(k)))+ &
				& (((y(k)**2)+(3*y(k)))*((x(j)**2)-x(j))))
			endif
		enddo
	enddo
	return
end subroutine init_source_function

subroutine init_problem(u, x, y, uinit, bce, bcn, bcs, bcw, extent, image, npes, m, delta)
! DESCRIPTION:
! This subroutine defines the boundary conditions and initial guess value for
! the problem. The matrix 'u' is of size [(extent+1)x(extent+1)].
!--------------------------------------------------------------------------
	integer, intent(in) :: npes, m, extent
	integer, intent(inout) :: image
	real(DOUBLE), intent(inout) :: u(0:extent,0:extent)
	real(DOUBLE), intent(in) :: x(0:extent), y(0:extent)
	real(DOUBLE), intent(in) :: bce, bcn, bcs, bcw, uinit, delta
	integer :: j,k
	!IDENTIFY AND FILL UP GLOBAL BOUNDARIES
	if((((image-1)) .LT. (sqrt(real(npes))-1))) then !CONDITION CHECKS FOR THE BOTTOM IMAGEES
		u(:,0) = bcs
	else if((mod((image-1),4) .EQ. 0)) then !CONDITION CHECKS FOR ALL THE IMAGEES AT THE LEFT SIDE
		u(0,:) = bcw
	else if((mod((image-1),4) .EQ. (sqrt(real(npes))-1))) then !CONDITION CHECKS FOR ALL THE IMAGEES ON THE RIGHT SIDE 
		u(extent,:) = bce
	else if((((image-1)) .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then!CONDITION CHECKS FOR ALL THE IMAGEES ON THE TOP SIDE
		u(:,extent) = bcn
	endif
	!FILL OUT THE INTERIOR POINTS FOR EACH IMAGE WITH THE GUESS VALUE, IE USE THE GUESS VALUE FOR ALL 
	!THE COORDINATES THAT ARE NOT IMAGEOR BOUNDARIES
	do j=0, extent
		do k=0, extent
			if((x(j) .NE. 0) .AND. (y(j) .NE. 0) .AND. (x(j) .NE. (m*delta)) .AND. (y(j) .NE. (m*delta))) then !CONDITION CHECKS IF THE COORDINATE IS A GLOBAL BOUNDARY
				u(j,k)=uinit
			endif
		enddo
	enddo		
	return	
end subroutine init_problem

subroutine numeric(m, extent, image, npes, itercount, q, left, right, up, down, flag, n_error, f, u, unew, delta, ghost_up, ghost_down, ghost_left, ghost_right)
!	DESCRIPTION:
!	This subroutine performs the jacobi iteration in parallel for each image separately. Each image has 
!	its own interior and boundary points. Boundary points of every image needs to exchange data
!	with its neighbours. The data to be exchanged is the data from the previous iteration.
!	Ghost points for each image is defined before a jacobi computation takes place.
!---------------------------------------------------------------------------------------------------------------------
    implicit none
	integer, parameter :: DOUBLE = kind(1.0d0)
	integer, intent(in) :: m, extent, npes, itercount ,q, left, right, up, down
	integer, intent(inout) :: flag, image
	real(DOUBLE), intent(inout) :: n_error[*]
	real(DOUBLE), intent(in) :: f(0:extent, 0:extent), delta
	real(DOUBLE), target, intent(inout) :: u(0:extent,0:extent)
	real(DOUBLE), intent(inout) :: unew(0:extent,0:extent)
	real(DOUBLE), intent(inout) :: ghost_left(0:extent)[*], ghost_right(0:extent)[*], ghost_up(0:extent)[*], ghost_down(0:extent)[*]
	real(DOUBLE),pointer :: north(:,:),east(:,:),west(:,:), south(:,:)
	integer :: j, k
!	DEFINE THE LOCAL BOUNDARIES OF EACH IMAGE AS THE GHOST VECTORS
	ghost_left(:)[image]=u(0,:)
	ghost_right(:)[image]=u(extent,:)
	ghost_up(:)[image]=u(:,extent)
	ghost_down(:)[image]=u(:,0)
	sync all
!*************************************************************************************************
!	JACOBI COMPUTATION SECTION
	north => u(1:extent-1,2:extent) !north(1,1)=u(1,2) 
	south => u(1:extent-1,0:extent-2) !south(1,1)=u(1,0) 
	west => u(0:extent-2,1:extent-1) !west(1,1)=u(0,1) 
	east => u(2:extent,1:extent-1) !east(1,1)=u(2,1) 
	sync all
	if(((image-1)) .LT. sqrt(real(npes))) then !CONDITION TO CHECK FOR BOTTOM IMAGES
		if(mod((image-1),4) .EQ. 0) then !BOTTOM LEFT IMAGE
			do j=1, extent
				do k=1, extent
					if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT CORNER
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
					else if((k .EQ. extent) .AND. (j .LT. extent)) then !TOP BDRY EXCL TOP RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .LT. extent)) then! RIGHT BDRY EXCL TOP-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))						
					else !INTERIOR PTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		else if(mod((image-1),4) .EQ. (sqrt(real(npes))-1)) then !BOTTOM-RIGHT IMAGE
			do j=0, (extent-1)
				do k=1, extent
					if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(3)]+south(j,k)+ghost_down(j)[(8)]-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .LT. extent)) then !LEFT END EXCL TOP LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(3)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((k .EQ. extent) .AND. (j .GT. 0)) then !TOP END EXCL TOP LEFT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(8)]-((delta**2)*f(j,k)))
					else !INTERIOR PTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		else !ALL OTHER IMAGEES IN THE BOTTOM ROW
			do j=0, extent
				do k=1, extent
					if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .LT. extent)) then!LEFT BOUNDARY EXCL TOP LEFT POINT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT POINT
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .LT. extent)) then !RIGHT BOUNDARY EXCL TOP RIGHT POINT
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. extent)) then !TOP BOUNDARY EXCL TOP-LEFT AND TOP-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
					else !INTERIOR POINTS
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		endif
	else if(((image-1) .GE. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !CONDITION CHECKS FOR TOP IMAGEES
		if(mod((image-1),4) .EQ. 0) then !TOP-LEFT IMAGE
			do j=1, extent
				do k=0, (extent-1)
					if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+ghost_left(k)[(image+1)]+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BOUNDARY EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .GT. 0)) then !RIGHT BDRY PTS EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+ghost_left(k)[(image+1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo					
		else if((image-1) .EQ. (npes-1)) then !TOP-RIGHT IMAGE
			do j=0, (extent-1)
				do k=0, (extent-1)
					if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (k .EQ. 0)) then !BOTTOM BOUNDARY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .GT. 0)) then !LEFT BNDRY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo	
		else !OTHER IMAGES ON TOP SIDE
			do j=0, extent
				do k= 0, (extent-1)
					if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. 0)) then!BOTTOM BDRY EXCL BOTTOM-LEFT & BOTTOM-RIGHT PT
						unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+ghost_up(k)[image]+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. 0) .AND. (k .GT. 0)) then !LEFT BDRY EXCL BOTTOM-LEFT PT
						unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else if((j .EQ. extent) .AND. (k .GT. 0)) then !RIGHT BDRY EXCL BOTTOM-RIGHT PT
						unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					else
						unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
					endif
				enddo
			enddo
		endif
	else if(((image-1) .GT. (sqrt(real(npes))-1)) .AND. (mod((image-1), 4) .EQ. 0) .AND. ((image-1) .LT. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !LEFT BOUNDARY IMAGEES EXCL TOP & BOTTOM PROCS
		do j=1, extent
			do k=0, extent
				if((j .EQ.  extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BDRY EXCL BOTTOM-RIGHT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .LT. extent) .AND. (k .GT. 0)) then !RIGHT BDRY EXCL TOP-RIGHT PT AND BOTTOM-RIGHT PT
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .LT. extent) .AND. (k .EQ. extent)) then !TOP BDRY EXCL TOP-RIGHT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else !INTERIOR PTS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				endif
			enddo
		enddo
	else if(((image-1) .GT. (sqrt(real(npes))-1)) .AND. (mod((image-1), 4) .EQ. (sqrt(real(npes))-1)) .AND. ((image-1) .LT. (sqrt(real(npes))*(sqrt(real(npes))-1)))) then !RIGHT BDRY PROCS EXCL TOP & BOTTOM PROCS
		do j=0, extent - 1
			do k=0, extent
				if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !LEFT BDRY PTS EXCL BOTTOM-LEFT PT AND TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (k .EQ. extent)) then !TOP BDRY PTS EXCL TOP-LEFT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (k .EQ. 0))then !BOTTOM BDRY PTS EXCL BOTTOM-LEFT PT
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else !INTERIOR PTS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				endif
			enddo
		enddo
	else !INTERIOR IMAGEES
		do j=0,extent
			do k=0, extent
				if((j .EQ. 0) .AND. (k .EQ. 0)) then !BOTTOM-LEFT CORNER
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. 0)) then !BOTTOM-RIGHT CORNER
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. 0)) then !BOTTOM BDRY EXCL CORNERS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+ghost_up(j)[(image-4)]+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .EQ. extent)) then !TOP-LEFT CORNER
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .EQ. extent)) then !TOP-RIGHT CORNER
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .GT. 0) .AND. (j .LT. extent) .AND. (k .EQ. extent)) then !TOP BDRY EXCL CORNERS
					unew(j,k)=0.25*(east(j,k)+west(j,k)+south(j,k)+ghost_down(j)[(image+4)]-((delta**2)*f(j,k)))
				else if((j .EQ. 0) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !LEFT BDRY EXL CORNERS
					unew(j,k)=0.25*(east(j,k)+ghost_right(k)[(image-1)]+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
				else if((j .EQ. extent) .AND. (k .GT. 0) .AND. (k .LT. extent)) then !RIGHT BDRY EXCL CORNERS
					unew(j,k)=0.25*(ghost_left(k)[(image+1)]+west(j,k)+south(j,k)+north(j,k)-((delta**2)*f(j,k)))
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
	sync all
!	END OF JACOBI COMPUTATION	
!*************************************************************************************************
!	CONVERGENCE CHECK - LOCAL MAXIMA
	if(mod(itercount,q) .EQ. 0) then
!		print*, itercount, "Inside convergence check"
		n_error=maxval(abs(unew-u))
!		print*, "ProcID:", rank
!		print*, n_error
	endif
	u=unew
!	print*, "ProcID:", rank
!	do j=0, extent
!		do k=0,extent
!			print*, unew(j,k)
!		enddo
!	enddo
	return
end subroutine numeric

subroutine validation(m, extent, uanalytic, validation_error, u, x, y, image)
!DESCRIPTION:
!
integer, intent(in) :: m, extent
integer, intent(inout) :: image
real(DOUBLE), intent(inout) :: uanalytic(0:extent,0:extent), validation_error[*]
real(DOUBLE), intent (in) :: u(0:extent,0:extent), x(0:extent),y(0:extent)
integer :: j,k
do j=0 , extent 
	do k=0 , extent
 		uanalytic(j,k)=exp(x(j)+y(k))*(((x(j)**2)-x(j))*((y(k)**2)-y(k)))
	enddo
enddo
validation_error[image] = maxval(abs(u-uanalytic))
return
end subroutine validation
end module jacobiCAF
