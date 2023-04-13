module integrators
use mpi
use init_pbc
use module_forces
contains
subroutine time_step_v_verlet(pos,vel,L,dt,npar,cutoff,nprocs,myid,displs,counts,particles, pairs,pairindex,nu,sigma)
    
    !=============================================================================
    ! Subroutine: time_step_v_verlet
    ! Description: Performs one Verlet time-step for a molecular dynamics simulation.
    !
    !----------------------------- INPUT -----------------------------------------
    ! pos - real8, dimension(npar,3): Matrix containing the x, y, z coordinates of the npar particles.
    ! pos_old - real8, dimension(npar,3): Matrix containing the two previous x, y, z coordinates of the npar particles.
    ! pos_aux - real8, dimension(npar,3): Matrix containing the previous x, y, z coordinates of the npar particles.
    ! L - real8: Dimension of the simulation box.
    ! dt - real8: Size of the time step.
    ! npar - integer: Number of particles in the system.
    ! cutoff - real8: Cutoff applied to the system.
    ! nprocs - integer: Total number of processors in the MPI communicator.
    ! myid - integer: Rank of the current processor in the MPI communicator.
    ! displs - integer, dimension(nprocs): Displacements for the MPI_ALLGATHERV operation.
    ! counts - integer, dimension(nprocs): Counts for the MPI_ALLGATHERV operation.
    ! particles - integer, dimension(0:nprocs-1,2): Particle indices for each processor.
    ! pairs - integer, dimension((npar*(npar-1))/2,2): Pairs of particles for computing forces.
    ! nu - real8: Parameter for Andersen thermostat.
    ! sigma - real8: Parameter for Andersen thermostat.
    !
    !----------------------------- OUTPUT ----------------------------------------
    ! pos - real8, dimension(npar,3): Matrix containing the new x, y, z coordinates of the npar particles after applying PBC.
    ! pos_aux - real8, dimension(npar,3): Previous pos.
    !=============================================================================

    integer :: npar, nprocs, myid, ierr,comm,i,ii,dim
    real*8 :: L,cutoff2,cutoff6,cutoff4,cutoff12,Upot,dt,cutoff 
    real*8,dimension(npar,3) :: pos,F,lj_force,vel
    integer, dimension(0:nprocs-1,2) :: particles,pairs
    integer, dimension((npar*(npar-1))/2,2) :: pairindex
    integer :: npar_local, npar_per_proc
    real*8,allocatable, dimension(:,:) :: pos_local, vel_local, lj_force_local
    integer, dimension(nprocs) :: counts, displs
    real*8 :: nu,sigma

    npar_local = counts(myid+1)
    
    allocate(pos_local(npar_local,3),vel_local(npar_local,3))
    allocate(lj_force_local(npar_local,3))

    ! Compute LJ forces for local particles
    call LENNARD_JONNES_FORCES(nprocs,myid,pairs,pairindex,pos, L, cutoff, npar, lj_force,Upot)

    !distribute the particles among the workers
    ii = 1
    do i = particles(myid,1), particles(myid,2)
        pos_local(ii,:) = pos(i,:)
        vel_local(ii,:) = vel(i,:)
        lj_force_local(ii,:) = lj_force(i,:)
        ii = ii + 1
    enddo

    ! if (myid.eq.0) then
    !     print*, 'after'
    !     do i=1,npar_local
    !         print*, vel_local(i,:)
    !     enddo 
    ! endif

    ! Update positions and velocities for local particles
    pos_local=pos_local+vel_local*dt+0.5d0*lj_force_local*dt*dt
    vel_local=vel_local+lj_force_local*0.5d0*dt 

    call pbc2(pos_local,L,npar_local)

    do dim=1,3

        call MPI_ALLGATHERV (pos_local(:,dim), npar_local, MPI_REAL8,  &
                         pos(:,dim), counts, displs, MPI_REAL8,   &
                         MPI_COMM_WORLD, ierr)


        call MPI_ALLGATHERV (vel_local(:,dim), npar_local, MPI_REAL8,  &
                         vel(:,dim), counts, displs, MPI_REAL8,   &
                         MPI_COMM_WORLD, ierr)

    enddo


    call LENNARD_JONNES_FORCES(nprocs,myid,pairs,pairindex,pos, L, cutoff, npar, lj_force,Upot)

    !distribute the particles among the workers
    ii = 1
    do i = particles(myid,1), particles(myid,2)
        lj_force_local(ii,:) = lj_force(i,:)
        ii = ii + 1
    enddo


    vel_local=vel_local+lj_force_local*0.5d0*dt 

    call therm_Andersen(vel_local,nu,sigma,npar_local)

    do dim=1,3
        call MPI_ALLGATHERV (vel_local(:,dim), npar_local, MPI_REAL8,  &
                         vel(:,dim), counts, displs, MPI_REAL8,   &
                         MPI_COMM_WORLD, ierr)
    enddo

end subroutine

subroutine box_muller(sigma,x1,x2,xout1,xout2)
    !===========================================================
    ! Gives two random numbers dstributed along a normal distribution.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! sigma - real*8 : width of the normal distribution
    ! x1,x2 - real*8 : random numbers between [0,1).
    !
    ! ... OUTPUT ...
    ! x1,x2 - real*8 : random numbers distributed along the 
    !                  correspondin normal distribution.
    !===========================================================


    implicit none
    real*8 sigma,x1,x2,xout1,xout2,pi

    pi=4.d0*atan(1.d0)

    xout1=sigma*sqrt(-2.d0*(log(1.d0-x1)))*cos(2.d0*pi*x2)
    xout2=sigma*sqrt(-2.d0*(log(1.d0-x1)))*sin(2.d0*pi*x2)

end subroutine


subroutine therm_Andersen(vel,nu,sigma,npar)

    !===========================================================
    ! Changes the velocities of the npar particles of the system
    ! acording to a thermostat of temperature sigma^2.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! vel - real*8,dimension(npar,3) : matrix containning the  
    !            x,y,z coordinates of the velocity of the npar particles.
    ! sigma - real*8 : square root of the temperature of the thermostat
    ! nu - real*8 : provides a selection criteria for applying or not
    !               the thermostat effects ([0,1)).
    ! npar - integer : number of particles of the system.
    !
    ! ... OUTPUT ...
    ! vel - real*8,dimension(npar,3) : matrix containning the new 
    !            x,y,z coordinates of the velocity of the npar particles.
    !===========================================================

    implicit none
    real*8 :: nu,sigma,x1,x2,xout1,xout2,x3,xout3,xout4,x4,x5
    real*8,dimension(npar,3) :: vel
    integer :: i,npar

    do i=1,npar
        call random_number(x1)
        call random_number(x2)
        call random_number(x3)
        call random_number(x4)
        call random_number(x5)

        if (x5.lt.nu) then
            call box_muller(sigma,x1,x2,xout1,xout2)
            call box_muller(sigma,x3,x4,xout3,xout4)            

            vel(i,:)=(/xout1,xout2,xout3/)
        endif

    enddo

end subroutine
end module

