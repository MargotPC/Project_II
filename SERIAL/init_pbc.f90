module init_pbc

contains
subroutine init_scc(npart,dimn,l_box,pos,filename)

    !==================================================================
    !Performs the calculation of the initial positions of the particles
    !------------------------------------------------------------------
    ! ... INPUT ...
    !npart - integer : number of particles of the system 
    !dimn - integer : dimension of the system
    !l_box - real*8 : length simulation box
    !filename - character(64) : name of the initial configuation file
    ! ... OUTPUT ...
    !pos - real*8, dimension(npart,dmin) : matrix 
    !               of xyz positions of npart particles
    !==================================================================
    implicit none
    integer, intent(in) :: npart, dimn
    real*8, intent (out) :: pos
    integer :: i,j,k,ncont,cells
    real*8 :: rho, l_box, lat
    dimension :: pos(npart,dimn)
    character(64) :: filename
    
    !initialization code
    cells = nint(dble(npart)**(1.d0/dble(dimn))) !unit cells in each dimension
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
    open(unit=10, file=trim(filename)//'.xyz')
    write(10,"(i4)") npart
    write(10,*) ""
    do i = 1,npart
        write(10,"(a3,f14.8,f14.8,f14.8)") 'Xe',pos(i,1), pos(i,2),pos(i,3)
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
