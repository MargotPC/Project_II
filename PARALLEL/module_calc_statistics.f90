module properties
use init_pbc
use mpi
contains

!This module contains the thermodynamics properties calculations.

subroutine calc_kinetic_energy_momentum(numproc, taskid, particles_division, N, vel, fkin, total_momentum)
!====================================
! Calculate the kinetic energy and momentum of the system.
!-----------------------------------
!...INPUT...
! numproc - integer : number of processors
! taskid - integer : task ID of each processor
! particles_division - integer, dimension (0:numproc-1, 2) : Table indicating the initial and final number of particles each processor need to calculate
! N - integer : number of particles of the system
! vel -  real*8, dimension(N,3) : matrix containing the velocity of each particle on each direction(vx,vy,vz)
!
!...OUTPUT...
! fkin - real*8 : kinetic energy of the system
! total_momentum - real*8 : total momentum of the system
!===================================
implicit none
integer, intent(in) :: N, numproc, taskid
integer, dimension(0:numproc-1,2), intent(in) :: particles_division
real*8, dimension(:,:), intent(in) :: vel
real*8, intent(out) :: fkin, total_momentum
integer :: ierror, Nini, Nfin
real*8 :: kin

Nini = particles_division(taskid,1)
Nfin = particles_division(taskid,2)

kin = 0.0

!Calculating the kinetic energy and momentum of the system

kin = sum(0.5d0*vel(Nini:Nfin,:)*vel(Nini:Nfin,:))

call MPI_REDUCE(kin, fkin, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

if ( taskid == 0 ) then
    total_momentum = dsqrt(2.d0*fkin)
end if
end subroutine



subroutine calc_Tinstant(taskid, N, kin, Tempins)
!====================================
! Calculate the instant temperature of the system
!-----------------------------------
!...INPUT...
! taskid - integer : task ID of each processor
! N - integer : number of particles of the system
! kin -  real*8, dimension(N,3) : kinetic energy
!
!...OUTPUT...
! Tempins - real*8 : instant temperature of the system
!===================================
implicit none
integer, intent(in) :: N, taskid
real*8, intent(in) :: kin
real*8, intent(out) :: Tempins

!Calculate the instant temperature only the master
if ( taskid == 0 ) then
   Tempins = kin*2.d0/(3.d0*N)
end if
end subroutine        


subroutine calc_pressure(numproc, taskid, particles_division, N, density, Tempins, pos, lj_forces,L, P)
!====================================
! Calculate the pressure of the system
!-----------------------------------
!...INPUT...
! numproc - integer : number of processors
! taskid - integer : task ID of each processor
! particles_division - integer, dimension (0:numproc-1, 2) : Table indicating the initial and final number of particles each processor need to calculate
! N - integer : number of particles of the system
! density - real*8 : density of the system
! Tempins - real*8 : instant temperature of the system
! pos - real*8, dimension(N,3) : matrix containing the position of each particle in cartesian coordinates
! lj_forces - real*8, dimension(N,3) : matrix containing the forces of each particle on each direction
! L - real*8 : dimension of the box
!                                            
!...OUTPUT...                                
! P - real*8 : Pressure of the system
!===================================    
implicit none
integer, intent(in) :: N, taskid, numproc
integer, dimension(0:numproc-1,2), intent(in) :: particles_division
real*8, intent(in) :: density, Tempins,L
real*8, dimension(:,:), intent(in) :: pos, lj_forces
real*8, intent(out) :: P
integer :: ierror, Nini, Nfin 
real*8 :: accumulate, faccumulate

Nini = particles_division(taskid,1)
Nfin = particles_division(taskid,2)

! Calculates the pressure of the system.
accumulate = 0
accumulate = sum(abs(lj_forces(Nini:Nfin,:)*pos(Nini:Nfin,:)))

call MPI_REDUCE(accumulate, faccumulate, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

if ( taskid == 0 ) then
    faccumulate = faccumulate/N/(3.d0*(L)**3)
    P = density*Tempins+faccumulate
end if
end subroutine



subroutine calc_msd(numproc, taskid, particles_division, L, N, density, init_pos, pos, fmsd)
!====================================
! Calculate the mean square displacement        
!-----------------------------------         
!...INPUT...
! numproc - integer : number of processors
! taskid - integer : task ID of each processor
! particles_division - integer, dimension (0:numproc-1, 2) : Table indicating the initial and final number of particles each processor need to calculate
! L - real*8 : length of the box
! N - integer : number of particles of the system
! density - real*8 : density of the system
! init_pos - real*8, dimension(N,3) : matrix that contains the initial positions of the system
! pos -  real*8, dimension(N,3) : matrix that contains the actual positions of the system        
!                                            
!...OUTPUT...                                
! fmsd - real*8 : mean square displacement of the system      
!=================================== 
implicit none
integer, intent(in) :: N, taskid, numproc
real*8, intent(in) :: density, L
real*8, dimension(:,:), intent(in) :: pos, init_pos
integer, dimension(0:numproc-1,2), intent(in) :: particles_division
real*8, intent(out) :: fmsd
integer :: i, ierror, Nini, Nfin
real*8, dimension(3) :: rij
real*8 :: d2, msd 

Nini = particles_division(taskid,1)
Nfin = particles_division(taskid,2)

msd=0.d0
!Computes the distance between the position of the particle in the initial position and the actual position
do i=Nini,Nfin
   rij(:) = pos(i,:) - init_pos(i,:)
   call pbc(rij,L)
   d2 = sum(rij(:)*rij(:))

   ! Calculates the mean square displacement
   msd = msd + d2
end do

call MPI_REDUCE(msd, fmsd, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

! Normalizes the mean square displacement
if ( taskid == 0 ) then
     fmsd = fmsd/N
end if
end subroutine



subroutine calc_radial_dist_func(numproc, taskid, pairs_division, pair_index, L, N, density, nbins, pos, gdR, distances_gdR, switch_case)
!====================================
! Calculate the radial distribution function
!-----------------------------------
!...INPUT...
! numproc - integer : number of processors
! taskid - integer :: task ID of each processor
! pairs_division - integer,dimension(0:numproc-1,2) :: Table indicating the initial and final number of pairs of particles each
! processor needs to calculate
! pair_index - integer,dimension(N*(N-1)/2,2) :: Table indicating each pair of particles with the first column the first particle of
! the pair and the second one, the second particle of the pair
! L - real*8 : Length of the box
! N - integer : number of particles of the system
! density - real*8 : density of the system
! nbins - integer : number of points that calculate the gdR function
! pos -  real*8, dimension(N,3) : matrix that contains the actual positions of the system
!
!...OUTPUT...
! gdR - real, dimension(250) : array with the gdR values for each distance contained in the "distances_gdR"
! distances_gdR - real, dimension(250) : array with the distances of the gdR
!===================================
! ********** HOW TO IMPLEMENT THIS SUBROUTINE ************
! Before the loop of time write the following line for initialize the gdr
!
!    calc_radial_dist_func(N, density, nbins, pos, gdR,distances_gdR, 1)
!
! In the loop of time the gdr has to be calculated for different times at the trajectory. It can be done in the
! conditional used for save data into an output file. The line will be
!
!    calc_radial_dist_func(N, density, nbins, pos, gdR,distances_gdR, 2)
!
! Once the time loop is over, the last part of the subroutine has to be implemented. Also in that part of the subroutine
! the writing into a file can be added (it is not implemented yet).
!
!     calc_radial_dist_func(N, density, nbins, pos, gdR,distances_gdR, 3)
!
! In and out variables

integer, intent(in) :: N, switch_case, nbins, taskid, numproc
integer, dimension(:,:), intent(in) :: pair_index
integer, dimension(0:numproc-1,2), intent(in) :: pairs_division
real*8, intent(in) :: density, L
real*8, dimension(:,:), intent(in) :: pos
real*8, dimension(:), intent(inout) :: gdR
real*8, dimension(:), intent(inout) :: distances_gdR
! parameters
real*8, parameter :: pi = 4.d0*atan(1.d0)
! inner varables
real*8, dimension(nbins) :: fgdR
real*8 :: r, d, v
integer :: i, j, A, ind, numpair, ierror
real*8, dimension(3) :: rij
! variables than need to be save in the memory
real*8, save :: dr
integer*8, save :: ngdr_calcul, Pini, Pfin

select case (switch_case)
   case(1)
      ! CASE 1 : INICIALIZATION
      ! The radial distribution function has to be averaged over all the times step;
      ! therefore, before time starts the number of calculations of gdr is 0
      ngdr_calcul = 0

      dr= L/(2.d0*nbins) ! Increments of the distance in the radial distribution function,
                              ! taking into account the number of bins (paramenter), dr,
                              ! and that only half of the box needs to be calculated (grater distances than L/2 would imply duplicity due to pbc).

      gdR=0 ! The array of gdR is initialize now
      distances_gdR=0

      ! Says the pair of particles that each worker calculates
      Pini = pairs_division(taskid,1)
      Pfin = pairs_division(taskid,2)

   case(2)
      ! CASE 2 : CALCULATION OF GDR EACH TIMESTEP
      ngdr_calcul = ngdr_calcul + 1 ! Every time than the subroutine is executed for switch_case=2 a calculation of gdr will procced.
                                    ! This way it is known how much calculations have been done, in order to promediate over time.
      
      ! Computes the distance between the particle j and A.
      do numpairs = Pini,Pfin
         A = pair_index(numpairs,1)
         j = pair_index(numpairs,2)
         rij(:) = pos(A,:)-pos(j,:)
         call pbc(rij,L)
         d = dsqrt(sum(rij(:)*rij(:)))
         ! Counts the amount of particles that are between each distance r + dr.
         if (d .lt. L/2.d0) then
            ind = int(d/dr) + 1
            gdR(ind) = gdR(ind) + 2.d0
         end if
      end do

   case(3)
      ! CASE 3 : NORMALIZATION OF GDR.
      
      call MPI_REDUCE(gdR, fgdR, nbins, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

      if (taskid == 0) then
      gdR = fgdR
      do i = 1,nbins
         ! Calculate the new distance (r). It's at the centre of the bin
         r = dr*(i + 0.5)
         ! Append the distance on the distances_gdR array.
         distances_gdR(i)=r
         ! Calculate the volume of the spherical shell with inner radius r-0.5dr=dr*i and outer radius r+0.5dr=dr(i+1)
         v=4.d0/3.d0*pi*((dr*(i+1))**3-(dr*i)**3)
         ! Calculates the radial distribution value for a certain distance.
         gdR(i)=gdR(i)/(ngdr_calcul*N*density*v)
      end do
      end if
   end select
end subroutine

end module properties
