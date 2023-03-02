subroutine init_scc(npart,dimn,l_box,coord)

	!==================================================================
	!Performs the calculation of the initial positions of the particles
	!------------------------------------------------------------------
	! ... INPUT ...
	!npart - integer : number of particles of the system 
	!dimn - integer : dimension of the system
	!l_box - real*8 : length simulation box
	! ... OUTPUT ...
	!coord - real*8, dimension(npart,dmin) : matrix 
	!				of xyz coordinates of the npart particles
	!==================================================================
	implicit none
	integer, intent(in) :: npart, dimn
	real*8, intent (out) :: coord
	integer :: i,j,k,ncont,cells
	real*8 :: rho, l_box, lat
	dimension :: coord(npart,dimn)
	
	!initialization code
	cells = nint(dble(npart)**(1d0/dble(dimn))) !unit cells in each dimension
	lat = l_box/dble(cells) !lattice spacing
	ncont=1 !matrix positions counter 
	
	do i = 1, cells
		do j = 1, cells
			do k = 1, cells
				coord(ncont,1)= i*lat
				coord(ncont,2)= j*lat
				coord(ncont,3)= k*lat
				ncont = ncont + 1
			end do
		end do
	end do
	return
	
end subroutine init_scc