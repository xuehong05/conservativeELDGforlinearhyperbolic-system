    !subroutine compute_flux_in_upstream(nk,u1mod,u2mod,RT_j1,RT_j2,xs,xc,dx,speed_l,speed_r,x_l,x_r,flux)
    !implicit none
    !integer,intent(in) :: nk
    !real,intent(in) :: u1mod(1:nk+1),u2mod(1:nk+1),RT_j1,RT_j2
    !real,intent(in) :: xc,dx
    !real,intent(in) :: xs,speed_l,speed_r,x_l,x_r
    !!real,intent(in) :: speed_l,speed_r,alpha,xl,xr
    !real,intent(out) :: flux
    !
    !real :: u1i,u2i,sp,flu11,flu12
    !real :: dtt,lamda
    !!**************************************
    !! this routine needs
    !!      rk_d(io) from global variables
    !!**************************************
    !
    !u1i = ortho_poly1d( u1mod(1:nk+1),xs ,xc,dx ,nk)
    !u2i = ortho_poly1d( u2mod(1:nk+1),xs ,xc,dx ,nk)
    !
    !!dtt = time-(time-dt + rk_d(io)*dt )
    !!sp = wave_linear_j(xs,speed_l,speed_r,xl,xr,(xr-xl),dtt)
    !!sp=element(i)%weno_t(1)
    !lamda=wave_linear_j(xs,speed_l,speed_r,x_l,x_r,(x_r-x_l))
    !!flux  = ((RT_j1*A11(xs)+RT_j(1,2)*A21(xs))*R_j(1,1)+(RT_j(1,1)*A12(xs)+RT_j(1,2)*A22(xs))*R_j(2,1)-lam1j)*v1i+((RT_j(1,1)*A11(xs)+RT_j(1,2)*A21(xs))*R_j(1,2)+(RT_j(1,1)*A12(xs)+RT_j(1,2)*A22(xs))*R_j(2,2))*v2i
    !flux  = RT_j1*((A11(xs)-lamda)*u1i+A12(xs)*u2i)+RT_j2*(A21(xs)*u1i+(A22(xs)-lamda)*u2i)
    !
    !!flux=0.
    !
    !end subroutine compute_flux_in_upstream
    !
    !
    !!subroutine compute_flux_in_upstream2(nk,v1mod,v2mod,xs,xc,dx,flux)
    !!implicit none
    !!integer,intent(in) :: nk
    !!real,intent(in) :: v1mod(1:nk+1),v2mod(1:nk+1)
    !!real,intent(in) :: xc,dx
    !!real,intent(in) :: xs
    !!!real,intent(in) :: speed_l,speed_r,alpha,xl,xr
    !!real,intent(out) :: flux
    !!
    !!real :: v1i,v2i,sp,flu11,flu12
    !!real :: dtt
    !!!**************************************
    !!! this routine needs
    !!!      rk_d(io) from global variables
    !!!**************************************
    !!
    !!!v1i = ortho_poly1d( v1mod(1:nk+1),xs ,xc,dx ,nk)
    !!!v2i = ortho_poly1d( v2mod(1:nk+1),xs ,xc,dx ,nk)
    !!!
    !!!dtt = time-(time-dt + rk_d(io)*dt )
    !!!!sp = wave_linear_j(xs,speed_l,speed_r,xl,xr,(xr-xl),dtt)
    !!!sp=element(i)%weno_t(1)
    !!!
    !!!flux  = f( ui,xs ) -  ui*sp
    !!flux=0.
    !!
    !!end subroutine compute_flux_in_upstream2
    !!
    !
    !!subroutine compute_flux_in_upstream(nk,umod, xs,xc,dx,speed_l,speed_r,alpha,xl,xr,flux)
    !!implicit none
    !!integer,intent(in) :: nk
    !!real,intent(in) :: umod(1:nk+1)
    !!real,intent(in) :: xc,dx
    !!real,intent(in) :: xs
    !!real,intent(in) :: speed_l,speed_r,alpha,xl,xr
    !!real,intent(out) :: flux
    !!
    !!real :: ui,sp
    !!real :: dtt
    !!!**************************************
    !!! this routine needs
    !!!      rk_d(io) from global variables
    !!!**************************************
    !!
    !!ui = ortho_poly1d( umod(1:nk+1),xs ,xc,dx ,nk)
    !!
    !!dtt = time-(time-dt + rk_d(io)*dt )
    !!sp = wave_linear_j(xs,speed_l,speed_r,xl,xr,(xr-xl),dtt)
    !!
    !!flux  = f( ui,xs ) -  ui*sp
    !!
    !!end subroutine compute_flux_in_upstream
    !
    !