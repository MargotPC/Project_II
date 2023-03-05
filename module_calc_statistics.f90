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



subroutine calc_radial_dist_func(N, density, pos, gdR, distances_gdR)
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
integer, intent(in) :: N
real*8, intent(in) :: density
real*8, dimension(:,:), intent(in) :: pos
real*8, dimension(250), intent(out) :: gdR, distances_gdR
real*8 :: pi, dr, r, d, L
integer :: Ncount, i, j, A
real*8, dimension(3) :: rij

! Calculates the length of the system.
L = (N/density)**(1.0/3.0)


pi = 4*atan(1.0)
dr = 0.01 ! Increments of the distance in the radial distribution function.
r=0.1 ! Initial distance to start computing the radial distribution function.

do i = 1, size(gdR)
   Ncount = 0
   ! Computes the distance between the particle j and A.
   do A = 1,N
      do j = 1,N
         rij(1)=(pos(j,1)-pos(A,1)); rij(2)=abs(pos(j,2)-pos(A,2)); rij(3)=abs(pos(j,3)-pos(A,3))
         call pbc(rij,L)
         d = (rij(1)**2+rij(2)**2+rij(3)**2)**(0.5)
         
         ! Counts the amount of particles that are separeted with a distance between r and dr.
         if (d .gt. r .and. d .lt. r+dr) then
            Ncount = Ncount+1
         end if
      end do
   end do
   ! Calculates the radial distribution value for a certain distance.
   gdR(i)=Ncount/(density*4*pi*r*r*dr)
   ! Append the distance on the distances_gdR array.
   distances_gdR(i)=r
   ! Calculate the new distance (r).
   r = r+dr
end do
end subroutine

end module properties
