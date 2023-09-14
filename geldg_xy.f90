    subroutine geldg_xy(    ixy1, tnum,dt17,xy_split)
    
    implicit none
    integer,intent(in) :: ixy1
    real,intent(in) :: dt17,tnum
    !real,intent(in) :: euler_gl(1:n_gl)
    !real,intent(in) :: euler_gauss(1:n_g)
    real,intent(in) :: xy_split
    !real :: upstream_gl(1:n_gl), upstream_gauss(1:n_g)
    !real,intent(out) :: umod_temp(1:ng,1:nx_2d  )
    !real :: xl_star,xr_star
    !integer :: ia,ib
    !integer :: mx
    !integer :: kc,inter
    !integer :: isub , nsub
    !integer :: idx
    !real,allocatable :: zz(:)
    
    ! STEP1. Locate the foots of trajectory x_l^{n,star}, x_r^{n,star}.
    ! We numerically solve the following final-value problem (trajectory equation):
    !             d x(t) / dt = a( x(t) , t )
    ! with the final-value x( t^{n+1} ) = x_l^n, x_r^n by a high order numerical integrator such as a classical fourth-order
    ! numerical integrator such as a classical fourth-order Runge-Kutta method.
    
    
    ! call RK(vx(1:nx+1) ,vx_star(1:nx+1),nx+1,dt )
    
    ! STEP2. locate Gauss-lobatto points
    
    
    if( ixy1 == 1 )then
        ixy=1
        y_split=xy_split
        time=tnum
        dt=dt17
        call init_x_1D
        call RK_xy
        !call RK_x( dxy, tnum,dt17,xy_split,umod_temp )
    elseif(ixy1 == 2)then
        ixy=2
        x_split=xy_split
        time=tnum
        dt=dt17
        call init_y_1D
        call RK_xy
        !call RK_y( dxy, tnum,dt17,xy_split,umod_temp )
        !call RK_y(dxy, tnum,dt17,xy_split,umod_temp )
    endif
    
    
    end subroutine geldg_xy