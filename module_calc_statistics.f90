module properties
use init_pbc
contains

!This module contains the thermodynamics properties calculations.

subroutine calc_kinetic_energy_momentum(N, vel, kin, total_momentum)
!====================================
! Calculate the kinetic energy and momentum of the system.
!-----------------------------------
!...INPUT...
! N - integer : number of particles of the system
! vel -  real*8, dimension(N,3) : matrix containing the velocity of each particle on each direction(vx,vy,vz)
!
!...OUTPUT...
! kin - real*8 : kinetic energy of the system
! total_momentum - real*8 : total momentum of the system
!===================================
implicit none
integer, intent(in) :: N
real*8, dimension(:,:), intent(in) :: vel
real*8, dimension(N,3) :: vel2
real*8, intent(out) :: kin, total_momentum
integer :: i
kin = 0.0

!Calculating the kinetic energy and momentum of the system
vel2 = 0.5*vel*vel
kin = sum(vel2)
total_momentum = sqrt(2*kin)

end subroutine



subroutine calc_Tinstant(N, kin, Tempins)
!====================================
! Calculate the instant temperature of the system
!-----------------------------------
!...INPUT...
! N - integer : number of particles of the system
! vel -  real*8, dimension(N,3) : matrix containing the velocity of each particle on each direction(vx,vy,vz)
!
!...OUTPUT...
! Tempins - real*8 : instant temperature of the system
!===================================
implicit none
integer, intent(in) :: N
real*8, intent(in) :: kin
real*8, intent(out) :: Tempins

!Calculate the instant temperature
Tempins = kin*2/(3*N)
end subroutine        



subroutine calc_pressure(N, density, Tempins, pos, lj_forces,L, P)
!====================================
! Calculate the pressure of the system
!-----------------------------------
!...INPUT...
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
integer, intent(in) :: N
real*8, intent(in) :: density, Tempins,L
real*8, dimension(:,:), intent(in) :: pos, lj_forces
real*8, intent(out) :: P

! Calculates the pressure of the system.
P = density*Tempins+sum(abs(lj_forces*pos))/dble(N)/(3.d0*(L)**3)

end subroutine



subroutine calc_msd(N, density, init_pos, pos, msd)
!====================================
! Calculate the mean square displacement        
!-----------------------------------         
!...INPUT...                                 
! N - integer : number of particles of the system
! density - real*8 : density of the system
! init_pos - real*8, dimension(N,3) : matrix that contains the initial positions of the system
! pos -  real*8, dimension(N,3) : matrix that contains the actual positions of the system        
!                                            
!...OUTPUT...                                
! msd - real*8 : mean square displacement of the system      
!=================================== 
integer, intent(in) :: N
real*8, intent(in) :: density
real*8, dimension(:,:), intent(in) :: pos
real*8, dimension(:,:), intent(in) :: init_pos
real*8, intent(out) :: msd
integer :: i
real*8, dimension(3) :: rij
real*8 :: d2, L

!Calculates the length of the system
L = (N/density)**(1.0/3.0)

msd=0.0
!Computes the distance between the position of the particle in the initial position and the actual position
do i=1,N
   rij(1) = (pos(i,1)-init_pos(i,1)); rij(2) = (pos(i,2)-init_pos(i,2)); rij(3) = (pos(i,3)-init_pos(i,3))
   call pbc(rij,L)
   d2 = (rij(1)**2)+(rij(2)**2)+(rij(3)**2)

   ! Calculates the mean square displacement
   msd = msd + d2
end do

! Normalizes the mean square displacement
msd = msd/N
end subroutine



subroutine calc_radial_dist_func(N, density, pos, gdR, distances_gdR, switch_case)
!====================================
! Calculate the radial distribution function
!-----------------------------------
!...INPUT...
! N - integer : number of particles of the system
! density - real*8 : density of the system
! pos -  real*8, dimension(N,3) : matrix that contains the actual positions of the system
!
!...OUTPUT...
! gdR - real, dimension(250) : array with the gdR values for each distance contained in the "distances_gdR"
! distances_gdR - real, dimension(250) : array with the distances of the gdR
!===================================

   ! ********** HOW TO IMPLEMENT THIS SUBROUTINE ************
   ! Before the loop of time write the following line for initialize the gdr
   !
   !    calc_radial_dist_func(N, density, pos, gdR,distances_gdR, 1)
   !
   ! In the loop of time the gdr has to be calculated for different times at the trajectory. It can be done in the 
   ! conditional used for save data into an output file. The line will be
   !
   !    calc_radial_dist_func(N, density, pos, gdR,distances_gdR, 2)
   !
   ! Once the time loop is over, the last part of the subroutine has to be implemented. Also in that part of the subroutine 
   ! the writing into a file can be added (it is not implemented yet).
   !
   !     calc_radial_dist_func(N, density, pos, gdR,distances_gdR, 3)
   !


! In and out variables
integer, intent(in) :: N, switch_case
real*8, intent(in) :: density
real*8, dimension(:,:), intent(in) :: pos
real*8, dimension(250), intent(inout) :: gdR
real*8, dimension(250), intent(out) :: distances_gdR
! parameters
integer, parameter :: nbins = 250
real*8, parameter :: pi = 4.d0*atan(1.d0)
! inner varables
real*8 :: r, d, L, v
integer :: i, j, A
real*8, dimension(3) :: rij
! variables than need to be save in the memory
real*8, save :: dr
integer, save :: ngdr_calcul

select case (switch_case)
   case(1)
      ! CASE 1 : INICIALIZATION
      ! The radial distribution function has to be averaged over all the times step; 
      ! therefore, before time starts the number of calculations of gdr is 0
      ngdr_calcul=0 

      ! Calculates the length of the system.
      L = (N/density)**(1.0/3.0)


      dr= int(L/(2.d0*nbins)) ! Increments of the distance in the radial distribution function, 
                              ! taking into account the number of bins (paramenter), dr,
                              ! and that only half of the box needs to be calculated (grater distances than L/2 would imply duplicity due to pbc).

      gdR=0 ! The array of gdR is initialize now

   case(2)
      ! CASE 2 : CALCULATION OF GDR EACH TIMESTEP
      ngdr_calcul=ngdr_calcul+1 ! Every time than the subroutine is executed for switch_case=2 a calculation of gdr will procced.
                                ! This way it is known how much calculations have been done, in order to promediate over time.

      
      ! Computes the distance between the particle j and A.
      do A = 1,N-1
         do j = A,N
            rij(:)=(pos(j,:)-pos(A,:))
            call pbc(rij,L)
            d = dsqrt(sum(rij(:)**2))
            
            ! Counts the amount of particles that are between each distance r + dr.
            if (d .lt. L/2.d0) then
               index = int(d/dr)
               gdR(index) = gdR(index) + 2
            end if
         end do
      end do

   case(3)
      ! CASE 3 : NORMALIZATION OF GDR. WRITING INTO A FILE CAN BE ADDED
      do i=1, nbins
         ! Calculate the new distance (r). It's at the centre of the bin
         r = dr*(i+0.5d0)
         ! Append the distance on the distances_gdR array.
         distances_gdR(i)=r
         ! Calculate the volume of the spherical shell with inner radius r-0.5dr=dr*i and outer radius r+0.5dr=dr(i+1)
         v=4.d0/3.d0*pi*((dr*(i+1))**3-(dr*i)**3)
         ! Calculates the radial distribution value for a certain distance.
         gdR(i)=gdR(i)/(ngr*N*density*v)
      end do

   end select
end subroutine

end module properties
