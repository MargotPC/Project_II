module module_forces
use init_pbc
use mpi 
contains

!This module contains the force calculation of the system.

subroutine LENNARD_JONNES_FORCES(numproc,taskid,pairs_division,pair_index,position, L, cutoff, N,Ffinal,upot_fin)
!===========================================================
! Performs the  Lennard-Jonnes forces.
!-----------------------------------------------------------
! ... INPUT ...
! position - real*8,dimension(npar,3) : matrix containning the x,y,z 
!       coordinates of the npar particles.
! L - real*8 : dimension of the simmulation box. 
! N - integer : number of particles of the system.
! cutoff - real*8 : cutoff applied to the system.
! pairs_division - integer,dimension(:,:) : Table that shows the first pair of particles and the last pair that each worker needs to compute 
! taskid 
! numproc 
! pair_index  - interger,dimenssion(:,;) : table indication each pair if particles with the first column being the first particle of the pair ,and the the second component 
! of the last particle     
!
! ... OUTPUT ...
! Ffinal - real*8,dimension(N,3) : Matrix contanning the new x,y,z coordinates of the forces of the N particles
! upot_fin - real*8 : Potential energy of the system
!===========================================================
implicit none
integer,dimension((N*(N-1))/2,2), intent(in) :: pair_index
real*8, intent(in) :: cutoff,L
real*8, dimension(N,3), intent(in) :: position
integer, intent(in) :: N, numproc, taskid
integer,dimension(0:numproc-1,2), intent(in) :: pairs_division
real*8, intent(out) :: upot_fin
real*8, dimension(N,3) :: force
real*8, dimension(N,3), intent(out) :: Ffinal
integer :: i, j,iii,Ninici,Nfinal,ierror,request
real*8 :: cutoff2, d2, d4, d6, d8, d12, d14, cf4, cf6, cf12, upot_pbc
real*8, dimension(3) :: rij

cutoff2 = cutoff*cutoff
force = 0.d0
upot_pbc = 0.d0

Ninici=pairs_division(taskid,1) !First particle of the particle division
Nfinal=pairs_division(taskid,2) !last particle of the particle division


!Computing the distances between particles.
do iii = Ninici,Nfinal
      i=pair_index(iii,1)
      j=pair_index(iii,2)
      rij(:) = position(i,:)-position(j,:)
      call pbc(rij,L)
      d2 = sum(rij(:)*rij(:))

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


!calling this MPI FUNCITION TO SUM ALL THE MATRIX FORCES IN A ONE LAST MATRIX
call MPI_ALLREDUCE(force,Ffinal,3*N,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
!calling the sum over all the  potential energy values from all the workers 
call MPI_REDUCE(upot_pbc,upot_fin,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierror)
end subroutine

end module
