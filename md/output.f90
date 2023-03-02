subroutine output(npart,dimn,l_box)

	!==========================================================
	!performs the creation of the output files -
	!initial coordinates of the particles 
	!----------------------------------------------------------
	! ... INPUT ...
	!npart - integer : number of particles of the system 
	!dimn - integer : dimension of the system
	!l_box - real*8 : length simulation box
	! ... OUTPUT ...
	!Save xyz file of the initial positions of the particles
	!==========================================================
	
	implicit none
	integer, intent (in) :: npart, dimn
	real*8, intent (in) :: l_box
	integer :: i
	real*8 :: coord
	dimension :: coord(npart,dimn)
	
	!initialization 
	call init_scc(npart,dimn,l_box,coord)
	
	!Save xyz coordinates file 
	open(unit=10, file='coordinates.xyz')
	write(10,"(i4)") npart
	write(10,*) ""
	do i = 1,npart
		write(10,"(a3,f14.8,f14.8,f14.8)") 'A',coord(i,1), coord(i,2),coord(i,3)
	end do
	close(10)

end subroutine output