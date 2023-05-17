    !*************************************************************************
    ! A  new Euler-Lagrangian discontinuous Galerkin method for
    !      one-dimensional linear transport problem
    !
    !                 u_t + (a(x,t)u)_x = 0.
    !
    ! It is potential to be extend to solve
    !                         u_t + f(u)_x = 0,
    !
    !
    !                                  Authors:
    !                                                  Xiaofeng Cai
    !                                             University of Delaware
    !                                                 November 22th, 2019
    !*************************************************************************
    ! please set up the coefficients of the method:
    !    iexample == 1 is for u_tt + u_xx = 0;
    !    iexample == 2 is for u_t + (sin(x)u)_x = 0.
    !*************************************************************************
    program ELDGF1d
    use module_globals1d
    use module_data1d
    
    use module_examples_init_exact
    implicit none

    call parameters
    do j=1,1
        !cfl=0.1*(2.**j)
        !cfl=1.5
        !cfl=0.3*j
        cfl=0.1*j
        !cfl=0.18*j
        write(*,*)cfl
    do kkkk =4,4 ! loop for different meshes or different CFLs
        !call CPU_TIME(time_begin)  
        call setup   ! setup the coefficients of the method.
        !
        call allocate_var
        !
        call init    ! initial condition3
        !
        call boundary_eulerian(0) 
        ! set the boundary condition for the Eulerian mesh.
        ! write(*,*) kkkk
        time = 0.
        do while( time<time_final)

            call setdt
            !write(*,*) time
            call get_upstream_tn
            do iio = 1,irk
                if (iio >1 ) then
                    call boundary_eulerian(iio-1) 
                    
                    do io=0,iio-1
                       call projection(io,iio )
                    enddo
                endif
                
                !!SSPRK
                !if (irk==1) then
                !    io=iio-1
                !    call rhs_int_flux
                !elseif (irk==2) then
                !    if (iio==1) then
                !       io=iio-1
                !       call rhs_int_flux
                !    else
                !       io=0
                !       call rhs_int_flux 
                !       flux_LF01=flux_LF1
                !       flux_LF02=flux_LF2
                !       io=iio-1
                !       call rhs_int_flux 
                !    endif
                !elseif (irk==3) then
                !    if (iio==1) then
                !       io=iio-1
                !       call rhs_int_flux
                !    elseif (iio==2) then
                !       io=0
                !       call rhs_int_flux 
                !       flux_LF01=flux_LF1
                !       flux_LF02=flux_LF2
                !       io=iio-1
                !       call rhs_int_flux 
                !    else
                !       io=0
                !       call rhs_int_flux 
                !       flux_LF01=flux_LF1
                !       flux_LF02=flux_LF2
                !       io=1
                !       call rhs_int_flux 
                !       flux_LF001=flux_LF1
                !       flux_LF002=flux_LF2
                !       io=iio-1
                !       call rhs_int_flux 
                !    endif
                !elseif (irk==4) then
                !   if (iio<irk) then
                !   io=iio-1
                !   call rhs_int_flux
                !   else
                !   io=0
                !   call rhs_int_flux 
                !   flux_LF01=flux_LF1
                !   flux_LF02=flux_LF2
                !   io=1
                !   call rhs_int_flux 
                !   flux_LF001=flux_LF1
                !   flux_LF002=flux_LF2
                !   io=2
                !   call rhs_int_flux 
                !   flux_LF0001=flux_LF1
                !   flux_LF0002=flux_LF2
                !   io=iio-1
                !   call rhs_int_flux 
                !   endif
                !endif
                
                
                !!RK
                !if (irk==1 .or. irk==2) then
                !    io=iio-1
                !    call rhs_int_flux
                !elseif (irk==3) then
                !    if (iio<irk) then
                !       io=iio-1
                !       call rhs_int_flux
                !    else
                !       io=0
                !       call rhs_int_flux 
                !       flux_LF01=flux_LF1
                !       flux_LF02=flux_LF2
                !       io=iio-1
                !       call rhs_int_flux 
                !    endif
                !elseif (irk==4) then
                !   if (iio<irk) then
                !   io=iio-1
                !   call rhs_int_flux
                !   else
                !   io=0
                !   call rhs_int_flux 
                !   flux_LF01=flux_LF1
                !   flux_LF02=flux_LF2
                !   io=1
                !   call rhs_int_flux 
                !   flux_LF001=flux_LF1
                !   flux_LF002=flux_LF2
                !   io=2
                !   call rhs_int_flux 
                !   flux_LF0001=flux_LF1
                !   flux_LF0002=flux_LF2
                !   io=iio-1
                !   call rhs_int_flux 
                !   endif
                !endif
                !
                !!
                
                !if(io==0)then !the alpha of global Lax-Friedrich flux 
                !    call get_max_sl_wave
                !endif
 

                !call boundary_lag

                call ssp_runge_kutta
                !do i=1,nx
                !    write(*,*) element(i)%u1modal(1,iio),element(i)%u2modal(1,iio)
                !enddo
                !pause
            
                
            enddo !io
            
            !call computeLint
            !write(40,*)time,abs(uet-u210)
            ! update solution on background elements
            u1mean=0.
            u2mean=0.
            do i = 1,nx
                !element(i)%u1modal(1:nk+1,irk) = element(i)%R_j(1,1)*element_t1(i,irk,irk)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,irk,irk)%v2modal(1:nk+1)
                !element(i)%u2modal(1:nk+1,irk) = element(i)%R_j(2,1)*element_t1(i,irk,irk)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,irk,irk)%v2modal(1:nk+1)
                u1mean=u1mean+element(i)%u1modal(1,irk)*dx
                u2mean=u2mean+element(i)%u2modal(1,irk)*dx
                element(i)%u1modal(1:nk+1,0)=element(i)%u1modal(1:nk+1,irk)
                element(i)%u2modal(1:nk+1,0)=element(i)%u2modal(1:nk+1,irk)
                !write(*,*) element(i)%u1modal(1:nk+1)
            enddo
            !pause
            write(25,*) time,abs(u1mean-u10),abs(u2mean-u20)
            !write(*,*) time,u1mean,u2mean,u10,u20
            !call error_dg(time)
            !call Lerror_dg
            !call error_dg(time)

        enddo

 
        !do i = 1,nx
        !    element(i)%u1modal(1:nk+1) = element(i)%R_j(1,1)*element_t1(i,irk)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,irk)%v2modal(1:nk+1)
        !    element(i)%u2modal(1:nk+1) = element(i)%R_j(2,1)*element_t1(i,irk)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,irk)%v2modal(1:nk+1)
        !enddo
        !call CPU_TIME(time_end)
        !write(*,*) time_end-time_begin
        !call order
        call slowpost1
        call order_dg
        write(26,*) time,abs(u1meanpost-u10),abs(u2meanpost-u20)
        !call output_dg
        !pause

        call deallocate_var

    enddo
    enddo
    

    contains
    include "get_max_sl_wave.f90"
    !include "get_max_wave.f90"
    !include "weno_limiter.f90"

    include "boundary_eulerian.f90"
    include "boundary_lag.f90"
     
    include "parameters.f90"
    include "get_element_tmass.f90"
    include "setdt.f90"
    include "get_upstream_tn.f90"

    include "rhs_int_flux.f90"
    include "ssp_runge_kutta.f90"
    !include "error_distribution.f90"
    !*******************************
    include "wave_linear_j.f90"
    include "allocate_var.f90"
    include "init.f90"

    !include "order.f90"
    include "compute_num_flux_LF.f90"
    include "setup.f90"
    include "get_integral_pk.f90"
    include "ortho_poly1d.f90"
    include "order_dg.f90"

    !*******************************
    include "get_integral_f_phi_x.f90"
    include "projection.f90"
    include "slowpost.f90"

    include "compute_flux_in_upstream.f90"
    !*******************************
    !include "output_dg.f90"
    !*******************************

    
    !*******************************
    real function burgex2(xin,tin)
    implicit none
    real, intent(in) :: xin,tin
    real :: xx,tt
    ! q0 = 0.5 + sin(x*pi), domain is [0,2]
    integer :: nn,ii
    real :: q0,q1,qtemp
    real :: xxx(500),qq(500),xux(500)
    nn = 400

    xx = xin - 0.5*tin
    tt = tin
    if(xx<0.) xx = xx + 2!*pi
    ! domain
    do ii=1,nn
        xxx(ii) =  2./nn*(ii-0.5)
        qq(ii) = sin( xxx(ii)*pi )
    enddo
    !
    do ii = 1 , nn
        xux(ii) = xxx(ii) + qq(ii) * tt
    enddo
    ! we choose a initial point by upwind
    if(xx>=0. .and. xx<=1.)then
        ! upwind [0,1]
        do ii = 1 , nn-1
            if( xx-xux(ii) >=0. .and. xx-xux(ii+1) <=0. )then
                q0 = qq(ii)
                goto 10
            endif
        enddo
10      continue
        ! q1 = q0 - f(q0)/f'(q0)
        q1 = 100000.
        qtemp = q0
        do while(abs(q0-q1)>1.e-12)
            q0 = qtemp
            q1 = q0 - ( q0-sin( pi*(xx-q0*tt) )   )/ ( 1. + pi*tt*cos( pi*(xx-q0*tt) ))
            qtemp = q1
        enddo
        burgex2 = q1 +0.5
    else
        ! downwind [1,2]
        do ii = 1 , nn-1
            if( xx-xux(ii) >=0. .and. xx-xux(ii+1) <=0. )then
                q0 = qq(ii+1)
                goto 100
            endif
        enddo
100     continue
        ! q1 = q0 - f(q0)/f'(q0)
        q1 = 100000.
        qtemp = q0
        do while(abs(q0-q1)>1.e-12)
            q0 = qtemp
            q1 = q0 - ( q0-sin(pi*(xx-q0*tt) ) )/ ( 1. + pi*tt*cos(pi*( xx-q0*tt ) )    )
            qtemp = q1
        enddo
        burgex2 = q1 +0.5
    endif


    end function burgex2
    !**********************************************************************
    
    real function burgersriemann(xx,tt)
    implicit none
    real, intent(in) :: xx,tt
    ! u is Riemann solution, domain is [-10,10]
    
    
    if(xx<=0.5*tt)then
        burgersriemann=1.
    else
        burgersriemann=0.
    endif
    
    end function burgersriemann
    !**********************************************************************
 
    !**********************************************************************
    !subroutine st_speed( uleft,uright,x,speed )
    !implicit none
    !real,intent(in) :: uleft,uright,x
    !real,intent(out) :: speed
    !
    !if( abs(uright - uleft)<eps*10 )then
    !    speed = fp(uleft*0.5+uright*0.5,x)
    !else
    !    speed = ( f(uright,x) - f(uleft,x) )/(uright-uleft)
    !endif
    !
    !!!perturbation sin(x)dx for u_t+u_x=0
    !!speed=1.+sin(x)*dx
    !
    !
    !
    !end subroutine  st_speed
    
    subroutine st_speed1( uleft,uright,x,speed1 )
    implicit none
    real,intent(in) :: uleft,uright,x
    real,intent(out) :: speed1

    !!for u_tt=u_xx
    !speed1=1.
    !
    !add mesh a perturbation sin(x)dx
    !speed1=1.+sin(x)*dx
    !speed1=1.+0.5
    !speed1=1.+dx
    speed1=aa(x)+sin(x)*dx
    !speed1=aa(x)
    !speed1=0.
    

    end subroutine  st_speed1
    
    subroutine st_speed2( uleft,uright,x,speed2 )
    implicit none
    real,intent(in) :: uleft,uright,x
    real,intent(out) :: speed2

    !!for u_tt=u_xx
    !speed2=-1.
    !!
    !add mesh a perturbation sin(x)dx 
    !speed2=-1.-sin(x)*dx
    !speed2=-1.-0.5
    !speed2=-1.-dx
    speed2=-1.*aa(x)-sin(x)*dx
    !speed2=-1.*aa(x)
    !speed2=0.

    end subroutine  st_speed2
    
    real function A11(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    A11=0.

    end function A11
    
    real function A12(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    A12=-1.*(aa(x)**2)

    end function A12
    
    real function A21(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    A21=-1.

    end function A21
    
    real function A22(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    A22=0.

    end function A22
    
    real function aa(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    if (iexample == 1) then
        aa=1.
        !aa=2.
    elseif (iexample == 2) then
        if (x>=-1. .and. x<=1.) then 
            aa=3./5.
        else
            aa=1.
        endif
    elseif (iexample == 3) then
        if (x>=-1. .and. x<=1.) then 
            aa=3./5.
        else
            aa=1.
        endif
    elseif (iexample == 4) then
        aa=2.+sin(x)
        !aa=2.
    endif
     
    
    end function aa
    
    
    !****************************************************
    
    real function aa1(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    if (iexample == 1) then
        !aa1=1.+sin(x)*dx
        !aa1=1.+sin(x+1.)*dx
        !aa1=1.+dx*exp(sin(x+pi/2.))/exp(1.)
        !aa1=1.+x*dx/(2.*pi)
        !aa1=2.+sin(x)*dx
        !aa1=2.
        aa1=1.
    elseif (iexample == 2) then
        if (x>=-1. .and. x<=1.) then 
            aa1=3./5.
        else
            aa1=1.
        endif
    elseif (iexample == 3) then
        if (x>=-1. .and. x<=1.) then 
            aa1=3./5.
        else
            aa1=1.
        endif
    elseif (iexample == 4) then
        aa1=2.+sin(x)
        !aa=2.
    endif
     
    
    end function aa1
    
     real function RT11(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    RT11=-0.5/(aa1(x))

    end function RT11
    
    real function RT12(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    RT12=0.5
    end function RT12
    
    real function RT21(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    RT21=0.5/(aa1(x))

    end function RT21
    
    real function RT22(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    RT22=0.5

    end function RT22
    
    real function aaderi(x)
    implicit none
    real,intent(in) :: x
    !the derivative of aa1
    if (iexample == 1) then
        aaderi=0.
        !aaderi=cos(x)*dx
        !aaderi=cos(x+1.)*dx
        !aaderi=dx*exp(sin(x+pi/2.))*cos(x+pi/2.)/exp(1.)
        !aaderi=dx/(2.*pi)
    elseif (iexample == 2) then
        if (x>=-1. .and. x<=1.) then 
            aaderi=0.
        else
            aaderi=0.
        endif
    elseif (iexample == 3) then
        if (x>=-1. .and. x<=1.) then 
            aaderi=0.
        else
            aaderi=0.
        endif
    elseif (iexample == 4) then
        aaderi=cos(x)
        !aa=2.
    endif
    end function aaderi
    !****************************************************
    
     real function R11(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    R11=-(aa1(x))

    end function R11
    
    real function R12(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    R12=(aa1(x))
    end function R12
    
    real function R21(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    R21=1.

    end function R21
    
    real function R22(x)
    implicit none
    real,intent(in) :: x
    ! the element of A(x)
    R22=1.

    end function R22 
    
    ! !**********************************************************************
    !subroutine roe_speed( uleft,uright,x,speed )
    !implicit none
    !real,intent(in) :: uleft,uright,x
    !real,intent(out) :: speed
    !
    !if( abs(uright - uleft)<eps*10 )then
    !    speed = fp(uleft*0.5+uright*0.5,x)
    !else
    !    speed = ( f(uright,x) - f(uleft,x) )/(uright-uleft)
    !endif
    !
    !end subroutine  roe_speed
    !*******************************
    real function f(u,x)
    implicit none
    real,intent(in) :: u,x

    if( iexample == 1 )then
        f = u
    elseif( iexample == 2 )then
        f = sin(x)*u
    elseif( iexample == 3 )then
        f = 0.5*u**2
    elseif( iexample == 4 )then
        f = u**2/(  u**2 +(1-u)**2  )
    endif

    end function f
    !*******************************
    real function fp1(x)
    implicit none
    real,intent(in) :: x
    !the eigenvalue of A(x)
    fp1=aa(x)

    end function fp1
    
    real function fp2(x)
    implicit none
    real,intent(in) :: x
   !the eigenvalue of A(x)
    fp2=-1.*aa(x)

    end function fp2

    real function f_source(x,t)
    implicit none
    real,intent(in) :: x,t
   !the eigenvalue of f(x,t)
    if( iexample == 1 )then
        f_source = 0.
    elseif( iexample == 2 )then
        f_source = 0.
    elseif( iexample == 3 )then
        f_source = 0.
    elseif( iexample == 4 )then
        f_source=-4.*sin(x-2.*t)+sin(x-2.*t)*((2.+sin(x))**2)-2.*(2.+sin(x))*cos(x)*cos(x-2.*t)
        !f_source=3.*sin(x-t)
    endif
    
    end function f_source
    !!*******************************
    !real function f(u,x)
    !implicit none
    !real,intent(in) :: u,x
    !
    !if( iexample == 1 )then
    !    f = u
    !elseif( iexample == 2   )then
    !    f = sin(x)*u
    ! endif
    !
    !end function f
    !!*******************************
    !real function fp(u,x)
    !implicit none
    !real,intent(in) :: u,x
    !
    !if( iexample == 1 )then
    !    fp = 1
    !elseif( iexample == 2  )then
    !    fp = sin(x)
    ! endif
    !
    !end function fp


    end program ELDGF1d