function pbc(dist,l_box)

	!==============================================================
	!Assign coordinates of particles with periodic box
	!--------------------------------------------------------------
	! ... INPUT ...
	!dist - real*8 : distance between two particles 
	!l_bodist - real*8 : length simulation box
	! ... OUTPUT ...
	!coord - real*8, dimension(npart,dmin) : return xyz coordinates 
	!				matrix of the npart particles with periodic box
	!==============================================================

	implicit none
	real*8, intent(in) :: l_box
	real*8 :: pbc, dist
	
	!pbc code
	if ((dist.le.(l_box/2D0)) .and. (dist.ge.(-l_box/2D0))) then
	pbc=dist
	else if (dist.gt.(l_box/2D0)) then
	pbc=dist-l_box
	else if (dist.lt.(l_box/2D0)) then
	pbc=dist+l_box
	endif
	return
	
end function pbc