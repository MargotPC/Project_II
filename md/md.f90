program md

	implicit none
	integer :: npart, dimn
	real*8 :: rho, l_box
	
	!Initial Conditions
	npart = 125 !number of particles
	dimn = 3 ! dimension of the system
	rho = 0.8d0 !density
	l_box = (npart/rho)**(1d0/dble(dimn)) !length simulation box
	
	!output - xyz file of initial coordinates of the particles
	call output(npart,dimn,l_box)
	
end program md
