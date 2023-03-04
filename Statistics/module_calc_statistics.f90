module calcthermodynamicproperties
use pbc, only : pbc1
contains

!This module contains the thermodynamics properties calculations.

subroutine kinetic_energy(N, vel, kin)
!====================================
! Calculate the kinetic energy of the system.
!-----------------------------------
!...INPUT...
! N - integer : number of particles of the system
! vel -  real*8, dimension(N,3) : matrix containing the velocity on each direction(vx,vy,vz)
!
!...OUTPUT...
! kin - real*8 : kinetic energy of the system
!===================================
implicit none
integer, intent(in) :: N
real*8, dimension(:,:), intent(in) :: vel
real*8, intent(out) :: kin
integer :: i
kin = 0.0

!Calculating the kinetic energy of the system
do i = 1,N
   kin = kin + 0.5*(vel(i,1)**2 + vel(i,2)**2 + vel(i,3)**2)
end do
end subroutine



subroutine calc_momentum(N, vel, total_momentum)
!====================================
! Calculate the toal momentum of the system
!-----------------------------------
!...INPUT...
! N - integer : number of particles of the system
! vel -  real*8, dimension(N,3) : matrix containing the velocity on each direction(vx,vy,vz)
!
!...OUTPUT...
! total_momentum - real*8 : total momentum of the system
!===================================
integer, intent(in) :: N
real*8, dimension(:,:), intent(in) :: vel
real*8, intent(out) :: total_momentum
integer :: i
real*8, dimension(3) :: momentum
momentum = 0

! Calculate the momentum (sum of the velocities) on each direction
do i = 1,N
   momentum(1) = momentum(1) + vel(i,1)
   momentum(2) = momentum(2) + vel(i,2)
   momentum(3) = momentum(3) + vel(i,3)
end do

! Computes the module of the momentum, the total momentum
total_momentum = (momentum(1)**2.0+momentum(2)**2.0+momentum(3)**2.0)**(0.5)
end subroutine



subroutine calc_PT(N, kin, density, cutoff, pos, P, Tempins)
!====================================
! Calculate the pressure and the instant temperature of the system       
!-----------------------------------         
!...INPUT...                                 
! N - integer : number of particles of the system        
! kin -  real*8 : kinetic energy of the system
! density - real*8 : density of the system
! cutoff - real*8 : cutoff applied to the system
! pos - real*8, dimension(N,3) : matrix containing the position of each particle in cartesian coordinates
!                                            
!...OUTPUT...                                
! P - real*8 : Pressure of the system
! Temp_ins - real*8 : Instant temperature of the system
!===================================                 
integer, intent(in) :: N
real*8, intent(in) :: kin, density, cutoff
real*8, dimension(:,:), intent(in) :: pos
real*8, intent(out) :: P, Tempins
real*8, dimension(3) :: rij, Fij
real*8 :: V, L, virial, cutoff2, d2, d4, d6, d8, d12, d14
integer :: i, j

! Calculates the instant temperature (it needs the kinetic energy).
Tempins = kin*2/(3*N)

! Calculates the volume and length of the system.
V = N/density
L = V**(1./3.)

virial = 0.0
cutoff2 = cutoff*cutoff
! Computes the distances between particles.
do i = 1,N-1
   do j = i+1,N
      rij(1) = (pos(i,1)-pos(j,1)); rij(2) = (pos(i,2)-pos(j,2)); rij(3) = (pos(i,3)-pos(j,3))
      call pbc1(L, rij)
      d2 = (rij(1)**2)+(rij(2)**2)+(rij(3)**2)

      ! Computes the forces of the system.
      if (d2 .lt. cutoff2) then
         d4 = d2*d2; d6 = d4*d2; d8 = d4*d4; d12 = d6*d6; d14 = d12*d2
         Fij(1)=(48/d14-24/d8)*rij(1)
         Fij(2)=(48/d14-24/d8)*rij(2)
         Fij(3)=(48/d14-24/d8)*rij(3)
         
         ! Computes the virial term needed to calculate the pressure.
         virial = virial + dot_product(rij, Fij)
      end if
   end do
end do

! Calculates the pressure of the system.
P = density*Tempins+(1/(3*V))*virial
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
   call pbc1(L, rij)
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

do i = 1, size(gdR):
   Ncount = 0
   ! Computes the distance between the particle j and A.
   do A = 1,N
      do j = 1,N
         rij(1)=(pos(j,1)-pos(A,1)); rij(2)=abs(pos(j,2)-pos(A,2)); rij(3)=abs(pos(j,3)-pos(A,3))
         call pbc1(L, rij)
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

end module calcthermodynamicproperties
