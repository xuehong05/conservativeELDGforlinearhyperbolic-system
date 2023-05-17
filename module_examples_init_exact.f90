    !*************************************************************************
    ! This module contains
    !      the initials condition and the exact solutions of
    !              an example
    !*************************************************************************
    module module_examples_init_exact
    use module_globals1d , only : iexample
    implicit none

    contains
    !*******************************
    real function fun_init(x )
    implicit none
    real,intent(in) :: x

    real :: pi
    pi = 4.*atan(1.)
    if(iexample==1)then
        fun_init = sin( x  )! for u_tt+ u_xx=0
    elseif(iexample==2)then
        fun_init = 1.! u_t + (sin(x)u)_x = 0
    elseif(iexample==3)then
        !fun_init = sin( x *pi)+0.5
        if( x<=0. )then
            fun_init = 1.
        else
            fun_init = 0.
        endif
    elseif(iexample ==4)then
        if( x>=0. .and. x<=0.05 )then
            fun_init = 1.-20.*x
        elseif( x>=0.25 .and. x<=0.4 )then
            fun_init = 0.5
        else
            fun_init = 0.
        endif
    endif

    end function fun_init
    !*******************************
    real function exact(x,t)
    implicit none
    real,intent(in) :: x,t

    real :: pi
    pi = 4.*atan(1.)
    if( iexample==1 )then
        exact = sin( x+t) ! for wave equation system u_tt=u_xx
    elseif( iexample==2 )then
        exact = sin( 2.*atan(exp(-t)*tan(x/2.) ) )/sin(x)! u_t + (sin(x)u)_x = 0
    elseif(iexample==3)then
        !exact = sin( x *pi)+0.5 ! u_0 for u_t + (u^2/2)_x = 0
        if( x<=0. )then
            exact = 1.
        else
            exact = 0.
        endif
    elseif(iexample ==4)then
        if( x>=0. .and. x<=0.05 )then
            exact = 1.-20.*x
        elseif( x>=0.25 .and. x<=0.4 )then
            exact = 0.5
        else
            exact = 0.
        endif !u_0 for u_t + (u^2/(u^2+(1-u)^2))_x = 0
    endif

    end function exact
    
    real function exact1(x,t)
    implicit none
    real,intent(in) :: x,t

    real :: pi
    pi = 4.*atan(1.)
    if( iexample==1 )then
        !exact1 = cos( x+t) ! for wave equation system u_tt=u_xx
        !exact1 = 0.5*(sin(x+t)+cos(x+t)+sin(x-t)-cos(x-t))
        !exact1 = -2.*cos( x-2.*t)
        if(t<=0.05*pi) then
            if( x<0.95*pi+t .and. x>=0.95*pi-t) then
                exact1=0.75
            elseif( x<1.05*pi-t .and. x>=0.95*pi+t) then
                exact1=1.
            elseif( x<=1.05*pi+t .and. x>=1.05*pi-t) then
                exact1=0.75
            else
            exact1 = 0.5
            endif
        else
            if( x<1.05*pi-t .and. x>=0.95*pi-t) then
                exact1=0.75
            elseif( x<0.95*pi+t .and. x>=1.05*pi-t) then
                exact1=0.5
            elseif( x<=1.05*pi+t .and. x>=0.95*pi+t) then
                exact1=0.75
            else
            exact1 = 0.5
            endif
        endif
        !
    elseif( iexample==2 )then
        !exact1 = sin( 2.*atan(exp(-t)*tan(x/2.) ) )/sin(x)! u_t + (sin(x)u)_x = 0
         ! for wave equation system u_tt=u_xx
        if( x<=1. .and. x>=-1. )then
            exact1 = -3.*pi*cos(5.*pi*x)*sin(3.*pi*t)
        else
            exact1 = -3.*pi*cos(3.*pi*x)*sin(3.*pi*t)
        endif
    elseif(iexample==3)then
        !exact = sin( x *pi)+0.5 ! u_0 for u_t + (u^2/2)_x = 0
        if( x<=0. )then
            exact1 = 1.
        else
            exact1 = 0.
        endif
    elseif(iexample ==4)then
        !if( x>=0. .and. x<=0.05 )then
        !    exact1 = 1.-20.*x
        !elseif( x>=0.25 .and. x<=0.4 )then
        !    exact1 = 0.5
        !else
        !    exact1 = 0.
        !endif !u_0 for u_t + (u^2/(u^2+(1-u)^2))_x = 0
        exact1=-2.*cos(x-2.*t)
        !exact1=-cos(x-t)
    endif

    end function exact1
    
    real function exact2(x,t)
    implicit none
    real,intent(in) :: x,t

    real :: pi
    pi = 4.*atan(1.)
    if( iexample==1 )then
        !exact2 = cos( x+t) ! for wave equation system u_tt=u_xx
        !exact2 = 0.5*(sin(x+t)+cos(x+t)-sin(x-t)+cos(x-t))
        !exact2 = cos( x-2.*t)
         if(t<=0.05*pi) then
            if( x<0.95*pi+t .and. x>=0.95*pi-t) then
                exact2=0.75
            elseif( x<1.05*pi-t .and. x>=0.95*pi+t) then
                exact2=1.
            elseif( x<=1.05*pi+t .and. x>=1.05*pi-t) then
                exact2=1.25
            else
            exact2 = 1.
            endif
        else
            if( x<1.05*pi-t .and. x>=0.95*pi-t) then
                exact2=0.75
            elseif( x<0.95*pi+t .and. x>=1.05*pi-t) then
                exact2=1.
            elseif( x<=1.05*pi+t .and. x>=0.95*pi+t) then
                exact2=1.25
            else
            exact2 = 1.
            endif
        endif
        
    elseif( iexample==2 )then
        !exact2 = sin( 2.*atan(exp(-t)*tan(x/2.) ) )/sin(x)! u_t + (sin(x)u)_x = 0
         ! for wave equation system u_tt=u_xx
        if( x<=1. .and. x>=-1. )then
            exact2 = -5.*pi*cos(3.*pi*t)*sin(5.*pi*x)
        else
            exact2 = -3.*pi*cos(3.*pi*t)*sin(3.*pi*x)
        endif
    elseif(iexample==3)then
        !exact = sin( x *pi)+0.5 ! u_0 for u_t + (u^2/2)_x = 0
        if( x<=0. )then
            exact2 = 1.
        else
            exact2 = 0.
        endif
    elseif(iexample ==4)then
        !if( x>=0. .and. x<=0.05 )then
        !    exact2 = 1.-20.*x
        !elseif( x>=0.25 .and. x<=0.4 )then
        !    exact2 = 0.5
        !else
        !    exact2 = 0.
        !endif !u_0 for u_t + (u^2/(u^2+(1-u)^2))_x = 0
        exact2=cos(x-2.*t)
        !exact2=cos(x-t)
    endif

    end function exact2
    
        

    end module module_examples_init_exact