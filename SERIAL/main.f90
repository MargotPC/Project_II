program Molecular_dynamics

use iso_c_binding
use init_pbc
use module_forces
use integrators
use properties


implicit none

integer M,npar,i,n,j,seed,steps,limit,nbins
real*8 :: density,cutoff,upot,cutoff4,cutoff6,cutoff12,cutoff2,L,dt,E,E_tot,mom,eps,mass,sig,msd
real*8,allocatable, dimension(:,:) :: pos,lj_force,vel,initial_pos
real*8,allocatable, dimension(:) :: gdR, distances_gdR
character(64) :: filename
character(8) :: fmt,ext,ext2
integer,dimension(3) :: M_values
real*8 :: sigma,x1,x2,xout1,xout2,T,nu,t_i,P,time,T1,T2
real*8 :: utime,utemp,upress,udens

seed=165432156
call srand(seed)

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


! ####################################################
!                INITIAL CONFIGURATION SC
! ####################################################

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

call init_scc(npar,3,L,pos,filename) !initial configuration could be generated if needed
! call read_xyz(filename,npar,pos) !read the initial configuration from xyz file
initial_pos = pos ! assign the initial positions to save it for the msd calculation

open(99, file='results/initial_conf_'//trim(ext)//'_sc.xyz')
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
write(*,'(x,A,x,f7.3,A,x,f7.3)') 'Thermostat temperature:',T1*utemp,' --->',T2*utemp 
write(*,'(x,A,x,f7.3)') 'Box length:',L
write(*,*) 'Andersen Thermostat: ON'
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
close(102)
print*, ''
write(*,*) '~~~~ STARTING MOLECULAR DYNAMICS ~~~~'

write(100,*) npar
write(100,*) ''

do j=1,npar!writes the initial configuration.

    write(99,*) 'Xe',pos(j,:)

enddo

close(99)

write(100,*) ''
write(100,*) ''

write(101,*) '# Density in g/cm^3:',density*udens


call calc_radial_dist_func(npar, density, nbins, pos, gdR,distances_gdR, 1)

T=T1
do i=0,steps

    if (i.eq.1d4) then!condition for, when then system is disordered, 
                      !starting the real simulation.
        T=T2 
        sigma=sqrt(T)
    endif

    if (mod(i,limit)==0) then

        do j=1,npar

            write(100,*) 'Xe',pos(j,:)

        enddo

        write(100,*) ''
        write(100,*) ''
    endif

    ! if (i.eq.10000) then
    !     sigma=sqrt(2.d0)
    ! endif
    
    call time_step_v_verlet(pos,vel,L,dt,npar,cutoff) !activate use verlet steps
    ! call time_step_Euler_pbc(pos,L,dt,npar,cutoff,vel) !activate to use euler steps
    call therm_Andersen(vel,nu,sigma,npar) !activate to add a thermostat

    call pbc(pos,L)

    if (mod(i,limit)==0) then

        call calc_kinetic_energy_momentum(npar,vel,E,mom)
        ! call calc_kinetic_energy(npar, vel, E)
        call calc_Tinstant(npar,E,t_i)
        ! call lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)
        call LENNARD_JONNES_FORCES(pos, L, cutoff, npar, Upot, lj_force)
        call calc_pressure(npar,density,t_i,pos,lj_force,L,P)
        call calc_msd(npar,density,initial_pos,pos,msd)
        call calc_radial_dist_func(npar, density, nbins, pos, gdR, distances_gdR, 2)
        
        E_tot=(E+Upot)*eps !so it is in kJ/mol
        t_i=utemp*t_i !so it is in K
        P=upress*P !so it is in Pa
        time=utime*i*dt !so it is in ps
        msd=msd*sig*sig    !so it is in Armstrongs

        write(101,*) time,Upot*eps,E*eps,E_tot,t_i,mom,P,msd

    endif

    if (mod(i,limit)==0) then

        ! write(*,'(I6,2x,A,1x,I6)',advance='no') i,'out of',steps
        ! write(*,'(A)',advance='no') repeat(c_backspace,100)

        call progress_bar(i,steps)
        
    endif

enddo

call calc_radial_dist_func(npar, density, nbins, pos, gdR,distances_gdR, 3)
write(103,*) "#gdR", "#distances"
do i=0,nbins
    write(103,*) gdR(i), distances_gdR(i)*sig*sig
enddo
    
write(*,*) ""

close(100)
close(101)
close(103)

! Uncomment to write a file with the final velocities in order to do a distribution
! open(100, file='final_vel.dat')
! do i=1,npar

!     write(100,*) vel(i,:)

! enddo

! close(100)


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

