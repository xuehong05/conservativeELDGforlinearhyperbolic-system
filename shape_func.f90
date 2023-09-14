    real function shape_func(x ,y )
    implicit none

    real,intent(in) :: x ,y

    real :: rb,rc,rs
    real :: pi1

    pi1 = pi*2.

    rb = sqrt( ( x +0.25*pi1 )**2 + ( y )**2  )
    rc = sqrt( ( x )**2 + ( y +0.25*pi1 )**2 )
    rs = sqrt( ( x )**2 + ( y -0.25*pi1 )**2 )
    if(rb < 0.2*pi1 )then
        shape_func = 0.25*( 1. + cos(pi* rb/(0.2*pi1) ) )
    elseif(rc <0.15*pi1)then
        shape_func = 1 - rc/(0.15*pi1)
    elseif(rs < 0.15*pi1)then
        shape_func = 1.
    else
        shape_func = 0.
    endif

    if( x >-0.025*pi1 .and. x < 0.025*pi1 .and. y >0.1*pi1 .and. y <0.35*pi1)then
        shape_func = 0.
    endif


    end function shape_func