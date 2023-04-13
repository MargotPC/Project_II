module init_pbc
use mpi
contains
subroutine init_scc(rank,nprocs,particles_division,counts,displs,pos_index,npart,l_box,pos,filename)

    !===========================================================================
	!Performs the calculation of the initial positions of the particles
	!---------------------------------------------------------------------------
	! *** INPUT ***
	!rank      - integer : Processors IDs.
	!nprocs    - integer : Number of processors.
	!particles_division - integer,dimension(:,:) : Initial and final particle
	!                      for each processor.
	!counts    - integer,dimension(nprocs) : number of particles for each 
	!                      processors.
	!pos_index - integer,dimenssion(:,:) : Index table of positions particles.
	!npart     - integer : number of particles of the system. 
	!l_box     - real*8  : length simulation box.
	! *** OUTPUT ***
	!pos       - real*8, dimension(npart,dmin) : matrix of xyz coordinates of
	!				       the npart particles.
	!===========================================================================
	
	implicit none
	integer, intent(in) :: npart, particles_division,pos_index
	real*8,intent (out) :: pos(npart,3)
	integer 	:: i,j,k,ii,ncont,cells, initial, final
	integer 	:: ierror, rank, nprocs, local_size
	real*8 		:: l_box, lat
	dimension 	:: particles_division(0:nprocs-1,2), pos_index(npart,3)
	integer 	:: counts(nprocs), displs(nprocs)
	character(64) :: filename

	initial = particles_division(rank,1)
	final = particles_division(rank,2)
	local_size = final-initial+1
    !print*, 'init= ', initial,'final= ',final,'rank= ',rank        

	cells = nint(dble(npart)**(1d0/dble(3)))
	lat = l_box/dble(cells) !lattice spacing
	
	!fill local position array
	ncont=initial !matrix positions counter 
	do i = initial,final
		pos(i,1)= pos_index(ncont,1)*lat
		pos(i,2)= pos_index(ncont,2)*lat
		pos(i,3)= pos_index(ncont,3)*lat
		ncont = ncont + 1
	end do
	
	do ii=1,3
		call MPI_ALLGATHERV(pos(initial:final,ii), local_size, MPI_REAL8,&
		pos(:,ii), counts, displs, MPI_REAL8, MPI_COMM_WORLD,ierror)
	end do
	
	!Save xyz coordinates file 
	if (rank==0) then
	open(unit=10, file=trim(filename)//'.xyz')
	write(10,"(i4)") npart
	write(10,*) ""
	do i = 1,npart
	write(10,"(a3,f14.8,f14.8,f14.8)") 'A',pos(i,1), pos(i,2), pos(i,3)
	end do
	close(10)
	end if
	return
    
end subroutine init_scc

subroutine pbc(pos,l_box)

    !==============================================================
    !Assign coordinates of particles with periodic box
    !--------------------------------------------------------------
    ! ... INPUT ...
    !pos - real*8, dimension(3) : distance between particles 
    !l_box - real*8             : length simulation box
    ! ... OUTPUT ...
    !pos - real*8, dimension(3) : return new distance between 
    !                             particles 
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

subroutine pbc2(pos,L,npar)
    !===========================================================
    ! Applies the periodic boundary conditions to any position
    ! of the npar particles.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! L - real*8 : dimension of the simmulation box. 
    ! npar - integer : number of particles of the system
    !
    ! ... OUTPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the new x,y,z 
    !       coordinates of the npar particles once pdb have been applied.
    !===========================================================

    implicit none

    real*8 :: L 
    integer :: npar
    real*8, dimension(npar,3) :: pos
    integer :: i,j


    do i=1,npar
        do j=1,3

            if (pos(i,j).gt.L/2.d0) then
                pos(i,j)=pos(i,j)-L
            endif

            if (pos(i,j).lt.(-L/2.d0)) then
                pos(i,j)=pos(i,j)+L
            endif

        enddo
    enddo  

end subroutine


subroutine read_xyz(filename,npar,pos)

    !===========================================================
    ! Reads a .xyz file and extracts the positions and number of
    ! particles of the system.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! filename - character(64) : name of the file where the initial
    !            condfiguation is written.
    !
    ! ... OUTPUT ...
    ! npar - integer : number of particles of the system.
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    !===========================================================

    implicit none
    integer :: npar,i
    real*8, dimension(npar,3) :: pos 
    real*8 :: x,y,z
    character(64) :: filename
    character(1) :: part

    open(100, file=trim(filename)//'.xyz',status='old')

            read(100,*)
            read(100,*)

            do i=1,npar

                read(100,*) part,x,y,z
                ! print*, part,x,y,z

                pos(i,:)=(/x,y,z/)

            enddo

    close(100)

end subroutine

end module init_pbc
