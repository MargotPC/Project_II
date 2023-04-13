program Molecular_dynamics

use iso_c_binding
use init_pbc
use module_forces
use integrators
use properties
use mpi
use iso_fortran_env


implicit none

integer :: M,npar,i,n,j,k,icont,seed,steps,limit,nbins,h,min,sec
real*8 :: density,cutoff,upot,cutoff4,cutoff6,cutoff12,cutoff2,L,dt,E,E_tot,mom,eps,mass,sig,msd
real*8,allocatable, dimension(:,:) :: pos,lj_force,vel,initial_pos
real*8,allocatable, dimension(:) :: gdR, distances_gdR
character(64) :: filename
character(8) :: fmt,ext,ext2
integer,dimension(3) :: M_values
real*8 :: sigma,x1,x2,xout1,xout2,T,nu,t_i,P,time,T1,T2,time1,time2
integer :: myid,comm,ierr,nprocs,ii,part_per_worker,pair_per_worker,rest1,rest2
integer, allocatable, dimension(:,:) :: pairs,particles,pairindex,posindex
integer, allocatable, dimension(:) :: displs, counts, seeds
real*8 :: utime,utemp,upress,udens

!Get the input from the input file
call get_command_argument(1, filename)

open(101,file=trim(filename))
    read(101,*) ! INITIAL CONFIGURATION SC
    read(101,*) M !size of the system (npar=M**3)
    read(101,*) mass !in g/mol
    read(101,*) sig !in Armstrongs (lj parameter)
    read(101,*) eps !in J (lj parameter)
    read(101,*) density !in m/sig^3

    read(101,*) ! INITIAL CONDITIONS
    read(101,*) vel !initial velocity of all the particles
    read(101,*) nu !probability of accepting a change with thermal andersen
    read(101,*) dt !timestep
    read(101,*) T1 !initial temperature of the system so it is disordered
    read(101,*) T2 !temperature around which the system will be equilibrated
    read(101,*) steps !steps of the simulation

    read(101,*) ! GDR
    read(101,*) nbins ! Number of points calculated in the gdR function
close(101)

npar=M**3 ! particles


fmt='(f5.3)'

write(ext,fmt) density
write(ext2,'(I6)') npar
ext2 = adjustl(ext2)
ext = adjustl(ext)

call All_reduce(eps,mass,sig,utime,utemp,upress,udens)

!Reduced units for all inputs:

density = density/udens
T1 = T1/utemp
T2 = T2/utemp



! ###########################
! INTIALIZE MPI ENVIRONMENT
! ###########################

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

!initialisation of a random seed for each processor
allocate(seeds(nprocs))


! Start computing the CPU time
if (myid.eq.0) then
    call CPU_time(time1)
    do i=1,nprocs
        seeds(i)=271846+(i)*3596
    enddo
endif
call random_seed(put=seeds)

!----------------------------------------------------------------------
!Define vector of index of particles (for init_scc subroutine)
allocate(posindex(npar,3))
icont=0
do i = 1, M
	do j = 1, M
		do k = 1, M
		icont = icont + 1
		posindex(icont,:)= [i,j,k]
		end do
	end do
end do
!-----------------------------------------------------------------------
!Define vector of pairs of particles (for the force and gdr subroutines)

allocate(pairindex((npar*(npar-1))/2,2))
ii = 1
do i = 1,npar-1
    do j = i+1,npar
        pairindex(ii,:) = (/i,j/)
        ii = ii + 1
    enddo
enddo
!_-----------------------------------------------------------------------

!Define the division of information per worker. Pairs contains the pairs that each worker has to evaluate while particles contains the particles.
!if it is not divisible, the partition is done as equitative as possible:
!         EXAMPLE: nprocs = 4 while npar = 35             
!                  part_per_worker = 8 --> partition would be (9:9:9:8)
     

allocate(particles(0:nprocs-1,2)) !starting from 0 since the workers ids go from 0 to nprocs-1
allocate(pairs(0:nprocs-1,2))
allocate(displs(nprocs), counts(nprocs)) ! prepareing vectors for the integrator

part_per_worker = npar/nprocs
pair_per_worker = (npar*(npar-1)/2)/nprocs

rest1 = mod(npar,nprocs)
rest2 = mod(npar*(npar-1)/2,nprocs)


do i = 0, nprocs-1
    if (rest1.gt.0) then
        particles(i,1) = i*(part_per_worker+1)+1
        particles(i,2) = (i+1)*(part_per_worker+1)
        counts(i+1) = (part_per_worker+1)
        rest1 = rest1 - 1
    else if (mod(npar,nprocs).eq.0) then
        particles(i,1) = i*(part_per_worker)+1
        particles(i,2) = (i+1)*(part_per_worker)
        counts(i+1) = (part_per_worker)
    else
        particles(i,1) = particles(i-1,2) + 1
        particles(i,2) = particles(i,1) + (part_per_worker-1)
        counts(i+1) = (part_per_worker)
    endif


    if (rest2.gt.0) then
        pairs(i,1) = i*(pair_per_worker+1)+1
        pairs(i,2) = (i+1)*(pair_per_worker+1)
        rest2 = rest2 - 1
    else if (mod(npar,nprocs).eq.0) then
        pairs(i,1) = i*(pair_per_worker)+1
        pairs(i,2) = (i+1)*(pair_per_worker)
    else
        pairs(i,1) = pairs(i-1,2) + 1
        pairs(i,2) = pairs(i,1) + (pair_per_worker-1)
    endif
enddo

displs(1) = 0
do i = 2, nprocs
    displs(i) = displs(i-1)+counts(i-1)
enddo


! ####################################
!         START MAIN PROGRAM
! ####################################


! ####################################################
!                INITIAL CONFIGURATION SC
! ####################################################


filename='results/initial_conf_'//trim(ext)//'_sc'
filename=trim(filename)

allocate(pos(npar,3),lj_force(npar,3),vel(npar,3),initial_pos(npar,3))

! ####################################################
!                  INITIAL CONDITIONS
! ####################################################

fmt='(f6.4)'
sigma=sqrt(T1)

! call bimodal(vel,sigma,npar)!set the initial velocities as a bimodal distribution

limit=steps/5000 !write 5000 data
L=(npar/density)**(1.d0/3.d0) !box length
cutoff=L/3.d0 !cutoff is set to a third of the box lenght
allocate(gdR(nbins), distances_gdR(nbins))

call init_scc(myid,nprocs,particles,counts,displs,posindex,npar,L,pos,filename) !initial configuration could be generated if needed
! call read_xyz(filename,npar,pos) !read the initial configuration from xyz file
initial_pos = pos ! assign the initial positions to save it for the msd calculation


if (myid.eq.0) then

    open(100, file='results/'//trim(ext2)//'_dynamics_'//trim(ext)//'.xyz')
    open(101, file='results/'//trim(ext2)//'_energy_'//trim(ext)//'.dat')
    open(103, file='results/'//trim(ext2)//'_gdr_'//trim(ext)//'.dat')
    open(102, file='get_input.dat')

! ####################################################
!                       DINÀMICA
! ####################################################

    write(*,*)''
    write(*,'(x,A,x,I3)') 'Number of particles:', npar
    write(*,'(x,A,x,f7.3)') 'Density of the system:',density*udens
    write(*,'(x,A,x,f9.3,A,x,f8.3)') 'Thermostat temperature:',T1*utemp,' --->',T2*utemp 
    write(*,'(x,A,x,f7.3)') 'Box length:',L
    write(*,*) 'Andersen Thermostat: ON'
    write(*,'(x,A,x,I3)') 'Number of processors:',nprocs
    print*,''

    write(*,*) 'Output files generated into /data:'
    write(*,'(8x,A)') trim(ext2)//'_dynamics_'//trim(ext)//'.xyz'
    write(*,'(8x,A)') trim(ext2)//'_energy_'//trim(ext)//'.dat'
    write(*,'(8x,A)') 'get_input.dat'
    write(102,'(A)') 'INPUT FILE'
    write(102,'(A)') 'results/'//trim(ext2)//'_energy_'//trim(ext)//'.dat'
    ! write(*,'(8x,A)') 'initial_vel.dat'
    ! write(*,'(8x,A)') 'final_vel.dat'
    write(*,'(8x,A)') 'initial_conf_'//trim(ext)//'_sc.xyz'
    write(*,'(8x,A)') trim(ext2)//'_gdr_'//trim(ext)//'.dat'
    write(102,'(A)') 'results/'//trim(ext2)//'_gdr_'//trim(ext)//'.dat'
    write(102,'(A)') '1000'
    write(102,'(A)') ''
    write(102,'(A)') 'THIS INPUT FILE READS THE NAME OF THE OUTPUT FILES AND THE NUMBER OF STEPS NEEDED TO EQUILIBRATE THE SYSTEM'

    print*, ''
    write(*,*) '~~~~ STARTING MOLECULAR DYNAMICS ~~~~'

    write(100,*) npar
    write(100,*) ''

    do j=1,npar!writes the initial configuration.

        write(100,*) 'Xe',pos(j,:)

    enddo

    write(100,*) ''
    write(100,*) ''

    write(101,*) '# Density in g/cm^3:',density*udens

endif

call calc_radial_dist_func(nprocs, myid, pairs, pairindex, L, npar, density, nbins, pos, gdR, distances_gdR, 1)

T=T1
do i=0,steps

    if (i.eq.1d4) then!condition for, when then system is disordered, 
                      !starting the real simulation.
        T=T2 
        sigma=sqrt(T)
    endif

    if ((mod(i,limit)==0).and.(myid.eq.0)) then

        do j=1,npar

            write(100,*) 'Xe',pos(j,:)

        enddo

        write(100,*) ''
        write(100,*) ''
    endif

    ! if (i.eq.10000) then
    !     sigma=sqrt(2.d0)
    ! endif
    
    call time_step_v_verlet(pos,vel,L,dt,npar,cutoff,nprocs,myid,displs,counts,particles,pairs,pairindex,nu,sigma) !activate use verlet steps
    
    if (mod(i,limit)==0) then

        call calc_kinetic_energy_momentum(nprocs,myid,particles,npar,vel,E,mom)
        call calc_Tinstant(myid,npar,E,t_i)
        call LENNARD_JONNES_FORCES(nprocs,myid,pairs,pairindex,pos, L, cutoff, npar, lj_force, Upot)
        call calc_pressure(nprocs,myid,particles,npar,density,t_i,pos,lj_force,L,P)
        call calc_msd(nprocs,myid,particles,L,npar,density,initial_pos,pos,msd)
        call calc_radial_dist_func(nprocs, myid, pairs, pairindex, L, npar, density, nbins, pos, gdR, distances_gdR, 2)

        
        E_tot=(E+Upot)*eps !so it is in kJ/mol
        t_i=utemp*t_i !so it is in K
        P=upress*P !so it is in Pa
        time=utime*i*dt !so it is in ps
        msd=msd*sig*sig    !so it is in Armstrongs

        if ((mod(i,limit*500)==0).and.(myid.eq.0)) then
            write(*,*) time,Upot*eps,E*eps,E_tot,t_i,mom,P,msd
            print*, ''
        endif


        if (myid.eq.0) then
            write(101,*) time,Upot*eps,E*eps,E_tot,t_i,mom,P,msd
        endif

    endif

enddo

if (myid.eq.0) then

    call CPU_time(time2)

    h=int((time2-time1)/3600.d0)
    min=int((time2-time1)/60.d0-60.d0*dble(h))
    sec=(time2-time1)-h*3600-min*60

    write(*,*) 'TIME: ',h,'h ',min,'m ',sec,'s'
endif

call calc_radial_dist_func(nprocs, myid, pairs, pairindex, L, npar, density, nbins, pos, gdR, distances_gdR, 3)

if (myid.eq.0) then
    write(103,*) "#gdR", "#distances"
    do i=0,nbins
        write(103,*) gdR(i), distances_gdR(i)*sig*sig
    enddo

    
    write(*,*) ""

    close(100)
    close(101)
    close(103)

endif

! Uncomment to write a file with the final velocities in order to do a distribution
! open(100, file='final_vel.dat')
! do i=1,npar

!     write(100,*) vel(i,:)

! enddo

! close(100)


!SPEED RUN
if (myid.eq.0) then
    write(*,*) time2-time1
endif

! Finalize MPI
call MPI_Finalize(ierr)

end program



! ###########################################
!               SUBROUTINE
! ###########################################

subroutine All_reduce(uener,mass,sigma,utime,utemp,upress,udens)

    implicit none
    real*8 :: mass,sigma,utime,utemp,uener,upress,udens

    upress = uener*1d3/(6.022d0*1d23*(sigma*1d-10)**3)   !pressure in atm (J/mol)/A^3=10^30/6.022E23 Pa * 1atm/101325Pa
    utemp  = uener*1d3/(6.022d0*1.380649d0)              ! temperature in Kelvin units (uener-->J/mol, 8.31---> J/(K·mol)
    udens  = mass/(6.022d0*1d23*(sigma*1d-8)**3)         ! 0.6022: factor (10^-8)^3*(6.022^10^(-23))
    utime  = sqrt(mass*1d-3/(uener*1d3))*sigma*1d2       ! unit time in picoseconds

end subroutine
