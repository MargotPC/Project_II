module integrators
use init_pbc
contains

subroutine lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)
    !===========================================================
    ! Calculates the value of the lennard-jonnes forces of the 
    ! sistem of npar particles i pos positions. The function also
    ! provides the forces of the whole system.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! cutoff2 - real*8 : cutoff value (cf) squared
    ! cutoff6 - real*8 : cf^6
    ! cutoff12 - real*8 : cf^12
    ! npar - integer : number of particles of the system
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! L - real*8 : dimension of the simmulation box. 
    ! pb - logical : indicates if the periodic boundary conditions
    !                are used (.True.) or not (.False.).
    !
    ! ... OUTPUT ...
    ! Upot - real*8 : Potential energy of the system via lennard jonnes
    ! lj_force - real*8,dimension(npar,3) : matrix containing the x,y,z
    !            components of the force done to each particle.
    !===========================================================
    implicit none
    real*8 :: cutoff6,Upot,dist2,cutoff12,dist4,dist6,dist12,cutoff2,L
    real*8 :: dist14,dist8 
    real*8,dimension(npar,3) :: pos,lj_force
    integer :: npar,i,j 
    real*8,dimension(3) :: delta
    logical :: pb


    
    Upot=0.d0
    lj_force=0.d0 
    do i=1,npar
        do j=i+1,npar

            delta=pos(i,:)-pos(j,:)

            if (pb.eqv..True.) then
                call pbc(delta,L)
            endif
            dist2=delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3)

            if (dist2.lt.cutoff2) then
                dist4=dist2*dist2
                dist6=dist4*dist2
                dist8=dist6*dist2
                dist12=dist6*dist6
                dist14=dist12*dist2

                Upot=Upot+4.d0*(1/dist12-1/dist6)-4.d0*(1.d0/cutoff12-1.d0/cutoff6)
                lj_force(i,:)=lj_force(i,:)+(48.d0/(dist14)-24.d0/(dist8))*delta
                lj_force(j,:)=lj_force(j,:)-(48.d0/(dist14)-24.d0/(dist8))*delta
            endif

        enddo
    enddo

end subroutine

subroutine time_step_Euler_pbc(pos,L,dt,npar,cutoff,vel)

    !===========================================================
    ! Performs one Euler time-step. 
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! L - real*8 : dimension of the simmulation box. 
    ! dt - real*8 : size of the time step.
    ! npar - integer : number of particles of the system
    ! cutoff - real*8 : cutoff applied to the system
    !
    ! ... OUTPUT ...
    ! vel - real*8,dimension(npar,3) : matrix containning the new x,y,z 
    !       coordinates of the velocities for the npar particles.
    !===========================================================

    implicit none
    real*8 :: L,cutoff2,cutoff6,cutoff4,cutoff12,Upot,dt,cutoff 
    real*8,dimension(npar,3) :: pos,vel,lj_force
    integer :: npar
    logical :: pb


    cutoff2=cutoff*cutoff
    cutoff4=cutoff2*cutoff2
    cutoff6=cutoff4*cutoff2
    cutoff12=cutoff6*cutoff6

    pb=.True.

    call lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)

    pos=pos+vel*dt+0.5d0*lj_force*dt*dt

    vel=vel+lj_force*dt*0.5d0

end subroutine


subroutine time_step_verlet(pos,pos_old,pos_aux,L,dt,npar,cutoff)

    !===========================================================
    ! Performs one Verlet time-step.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! pos_old - real*8,dimension(npar,3) : matrix containning the  
    !           two previous x,y,z coordinates of the npar particles.
    ! pos_aux - real*8,dimension(npar,3) : matrix containning the  
    !           previous x,y,z coordinates of the npar particles.
    ! L - real*8 : dimension of the simmulation box. 
    ! dt - real*8 : size of the time step.
    ! npar - integer : number of particles of the system.
    ! cutoff - real*8 : cutoff applied to the system
    !
    ! ... OUTPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the new x,y,z 
    !       coordinates of the npar particles once pdb have been applied.
    ! pos_aux - real*8,dimension(npar,3) : previous pos.
    !===========================================================

    implicit none
    real*8 :: L,cutoff2,cutoff6,cutoff4,cutoff12,Upot,dt,cutoff 
    real*8,dimension(npar,3) :: pos,F,lj_force,pos_old,pos_aux
    integer :: npar
    logical :: pb

    cutoff2=cutoff*cutoff
    cutoff4=cutoff2*cutoff2
    cutoff6=cutoff4*cutoff2
    cutoff12=cutoff6*cutoff6

    pb=.True.

    call lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)

    pos_aux=pos

    pos=2.d0*pos-pos_old+lj_force*dt*dt

end subroutine


subroutine time_step_v_verlet(pos,vel,L,dt,npar,cutoff)

    !===========================================================
    ! Performs one velocity Verlet time-step.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the x,y,z 
    !       coordinates of the npar particles.
    ! vel - real*8,dimension(npar,3) : matrix containning the  
    !            x,y,z coordinates of the velocity of the npar particles.
    ! L - real*8 : dimension of the simmulation box. 
    ! dt - real*8 : size of the time step.
    ! npar - integer : number of particles of the system.
    ! cutoff - real*8 : cutoff applied to the system
    !
    ! ... OUTPUT ...
    ! pos - real*8,dimension(npar,3) : matrix containning the new x,y,z 
    !       coordinates of the npar particles once pdb have been applied.
    ! vel - real*8,dimension(npar,3) : matrix containning the new 
    !            x,y,z coordinates of the velocity of the npar particles.
    !===========================================================

    real*8 :: L,cutoff2,cutoff6,cutoff4,cutoff12,Upot,dt,cutoff 
    real*8,dimension(npar,3) :: pos,F,lj_force,vel
    integer :: npar
    logical :: pb

    cutoff2=cutoff*cutoff
    cutoff4=cutoff2*cutoff2
    cutoff6=cutoff4*cutoff2
    cutoff12=cutoff6*cutoff6

    pb=.True.

    call lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)

    pos=pos+vel*dt+0.5d0*lj_force*dt*dt
    vel=vel+lj_force*0.5d0*dt 

    call lj_f(cutoff2,cutoff6,cutoff12,npar,pos,L,pb,Upot,lj_force)

    vel=vel+lj_force*0.5d0*dt 

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
        x1=rand()
        x2=rand()
        x3=rand()
        x4=rand()
        x5=rand()

        if (x5.lt.nu) then
        call box_muller(sigma,x1,x2,xout1,xout2)
        call box_muller(sigma,x3,x4,xout3,xout4)
            vel(i,:)=(/xout1,xout2,xout3/)
        endif

    enddo

end subroutine

subroutine progress_bar(iteration,max_it)

    !===========================================================
    ! Prints a progress bar in the console that displays the
    ! progress of the iterations. Must be placed inside a loop.
    !-----------------------------------------------------------
    ! ... INPUT ...
    ! iteration - integer : current iteration.
    ! max_it - integer : maximimum of iterations of the loop.
    !
    ! ... OUTPUT ...
    ! Prints a progress bar in the console.
    !===========================================================

    use iso_c_binding
    implicit none
    integer :: iteration,max_it,prog,end,estimated_time,h,min,sec,bs
    character(3) :: percentage
    character(1) :: a
    character(2) :: h1,min1,sec1,h2,min2,sec2
    character(20) :: ET,T
    real*8 :: time

    a='#'

    call CPU_time(time)

    prog=nint(dble(iteration)/dble(max_it)*100.d0)

    ! if (prog.eq.5) then 
    !     estimated_time=nint(time*20)
    !     h=int(estimated_time/3600.d0)
    !     min=int(estimated_time/60.d0)
    !     sec=estimated_time-h*3600-min*60
    !     write(sec1,'(I2)') sec
    !     write(min1,'(I2)') min
    !     write(h1,'(I2)') h
    !     bs=80
    !     ET=h1//'h '//min1//'min '//sec1//'s'

    ! endif
        
    
    ! if (prog.lt.5) then

    !     ET=' ...'
    !     bs=100
    ! endif



    h=int(time/3600.d0)
    min=int(time/60.d0)
    sec=time-h*3600-min*60

    write(sec2,'(I2)') sec
    write(min2,'(I2)') min
    write(h2,'(I2)') h
    T=h2//'h '//min2//'m '//sec2//'s'



    write(percentage,'(I3)') prog




    write(*,'(A)',advance='no') percentage//'% |'//repeat(a,nint(3*prog/10.d0))//&
                                repeat('.',30-nint(3*prog/10.d0))//'| '//trim(T)
    

    end=int(aint(dble(iteration)/dble(max_it)*100.d0)) 

    if (end.ne.100) then
        write(*,'(A)',advance='no') repeat(c_backspace,60)
    endif

    if (end.eq.100) then 

        write(*,'(A)',advance='no') ' ---> DONE'
    endif

end subroutine

end module integrators