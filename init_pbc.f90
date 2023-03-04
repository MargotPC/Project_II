module init_pbc

contains
subroutine init_scc(npart,dimn,l_box,pos)

	!==================================================================
	!Performs the calculation of the initial positions of the particles
	!------------------------------------------------------------------
	! ... INPUT ...
	!npart - integer : number of particles of the system 
	!dimn - integer : dimension of the system
	!l_box - real*8 : length simulation box
	! ... OUTPUT ...
	!pos - real*8, dimension(npart,dmin) : matrix 
	!				of xyz positions of npart particles
	!==================================================================
	implicit none
	integer, intent(in) :: npart, dimn
	real*8, intent (out) :: pos
	integer :: i,j,k,ncont,cells
	real*8 :: rho, l_box, lat
	dimension :: pos(npart,dimn)
	
	!initialization code
	cells = nint(dble(npart)**(1d0/dble(dimn))) !unit cells in each dimension
	lat = l_box/dble(cells) !lattice spacing
	ncont=1 !matrix position counter 
	
	do i = 1, cells
		do j = 1, cells
			do k = 1, cells
				pos(ncont,1)= i*lat
				pos(ncont,2)= j*lat
				pos(ncont,3)= k*lat
				ncont = ncont + 1
			end do
		end do
	end do
	
	!Save xyz coordinates file 
	open(unit=10, file='input.xyz')
	write(10,"(i4)") npart
	write(10,*) ""
	do i = 1,npart
		write(10,"(a3,f14.8,f14.8,f14.8)") 'A',pos(i,1), pos(i,2),pos(i,3)
	end do
	close(10)
	return
	
end subroutine init_scc

subroutine pbc(pos,l_box)

	!==============================================================
	!Assign coordinates of particles with periodic box
	!--------------------------------------------------------------
	! ... INPUT ...
	!pos - real*8, dimension(3) : distance between particles 
	!l_box - real*8 			: length simulation box
	! ... OUTPUT ...
	!pos - real*8, dimension(3) : return new distance between 
	!							  particles 
	!==============================================================

	implicit none
	real*8, intent(in) :: l_box
	integer :: i
	real*8,dimension (3) :: pos
	
	!pbc code
	do i = 1, size(pos)
		pos(i) = pos(i) - nint(pos(i)/l_box)*l_box
	end do
	return
	
end subroutine pbc

end module init_pbc
