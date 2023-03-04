!===========================================================
! Performs the  Lennard-Jonnes forces.
!-----------------------------------------------------------
! ... INPUT ...
! position - real*8,dimension(npar,3) : matrix containning the x,y,z 
!       coordinates of the npar particles.
! L - real*8 : dimension of the simmulation box. 
! N - integer : number of particles of the system.
! cutoff - real*8 : cutoff applied to the system.
!
! ... OUTPUT ...
! upot_pbc - real*8 : potential energy of the system.
! force - real*8,dimension(npar,3) : matrix containning the new 
!            x,y,z coordinates of the forces of the npar particles.
!===========================================================

module module_forces
use pbcs,only:pbc1 !this will change with the new pbcs module
contains

!This module contains the force calculation of the system.

subroutine LENNARD_JONNES_FORCES(position, L, cutoff, N, upot_pbc, force)
implicit none
real*8, intent(in) :: cutoff,L
real*8, dimension(:,:), intent(in) :: position
integer, intent(in) :: N
real*8, intent(out) :: upot_pbc
real*8, dimension(:,:), intent(out) :: force
integer :: i, j
real*8 :: cutoff2, d2, d4, d6, d8, d12, d14, cf4, cf6, cf12
real*8, dimension(3) :: rij
cutoff2 = cutoff*cutoff
force = force*0.0
upot_pbc = 0.0

!Computing the distances between particles.
do i = 1,N-1
   do j = i+1, N
      rij(1) = (position(i,1)-position(j,1)); rij(2) = (position(i,2)-position(j,2)); rij(3) = (position(i,3)-position(j,3))
      call pbc1(L, rij)
      d2 = (rij(1)**2)+(rij(2)**2)+(rij(3)**2)

      !Computing the forces , potential energy of the system between particles with a distance lower than the cutoff
      if (d2 .lt. cutoff2) then
         d4 = d2*d2; d6 = d4*d2; d8 = d6*d2; d12 = d6*d6; d14 = d12*d2
         cf4 = cutoff2*cutoff2; cf6 = cutoff2*cf4; cf12 = cf6*cf6

         !Here the potential energy of the system is computed.
         upot_pbc = upot_pbc + 4.0*((1.0/d12) - (1.0/d6)) - 4.0*((1.0/cf12) - (1.0/cf6))

         !Here the forces  for each particle are calculated.
         force(i,1) = force(i,1) + (48/d14-24/d8)*rij(1)
         force(i,2) = force(i,2) + (48/d14-24/d8)*rij(2)
         force(i,3) = force(i,3) + (48/d14-24/d8)*rij(3)
         force(j,1) = force(j,1) - (48/d14-24/d8)*rij(1)
         force(j,2) = force(j,2) - (48/d14-24/d8)*rij(2)
         force(j,3) = force(j,3) - (48/d14-24/d8)*rij(3)
      end if
   end do
end do
end subroutine
end module
