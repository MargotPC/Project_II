module magnitudes
contains

subroutine Ekin(vel,npar,E,mom)

    !===========================================================
    ! Computes the kinetic energy and the momentum of the system
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! vel - real*8,dimension(npar,3) : matrix containning the  
    !            x,y,z coordinates of the velocity of the npar particles.
    ! npar - integer : number of particles of the system.
    !
    ! ... OUTPUT ...
    ! E - real*8 : Kynetic energy of the system.
    ! mom - real*8 : Momentum of the system
    !===========================================================


    implicit none
    real*8,dimension(npar,3) :: vel,vel2
    integer :: npar
    real*8 :: E,mom

    vel2=0.5d0*vel*vel
    E=sum(vel2)
    mom=sqrt(2*E)

end subroutine


subroutine t_inst(ekin,npar,t)

    !===========================================================
    ! Computes the instant temperature of the system.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! ekin - real*8 : Kynetic energy of the system.
    ! npar - integer : number of particles of the system.
    !
    ! ... OUTPUT ...
    ! t - real*8 : instant temperature of the system.
    !===========================================================

    implicit none
    integer :: npar,i
    real*8 :: t,ekin

    t=2.d0/3.d0*ekin/dble(npar)

end subroutine

subroutine P_inst(t_i,density,lj_force,npar,pos,L,P)

    !===========================================================
    ! Computes the instant pressure of the system.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! t_i - real*8 : instant temperature of the system.
    ! density - real*8 : density of the system in reduced units.
    ! lj_force - real*8,dimension(npar,3) : matrix containing the x,y,z
    !            components of the force done to each particle.
    ! npar - integer : number of particles of the system.
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! L - real*8 : dimension of the simmulation box.
    !
    ! ... OUTPUT ...
    ! P - real*8 : Pressure of the system.
    !===========================================================

    implicit none
    integer :: npar
    real*8,dimension(npar,3) :: pos,lj_force
    real*8 :: t_i,P,L,density

    P=density*t_i+sum(abs(lj_force*pos))/dble(npar)/(3.d0*L**3)

end subroutine

end module magnitudes