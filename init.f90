
    !*************************************************************************
    !*************************************************************************
    subroutine init
    implicit none
    real :: u1temp,u2temp,ut1,ut2,u21temp,u22temp,ut21,ut22
    integer :: lx,kk

    ! x(i) ----> center

    dx = (xright-xleft)/nx

    DO I = 1 - nghost, nx + nghost
        x(I) = xleft + (I-0.5) * DX
    ENDDO

    DO I = 1 - nghost, nx + 1  + nghost
        Xgrid(I) = xleft + (I-1.) * DX
    ENDDO

    do i = 1,nx+1
        vertex(i)%coor = xleft + (i-1) * dx
    enddo

    do i = 1,nx
        element(i)%xl = x(i)-0.5*dx
        element(i)%xr = x(i)+0.5*dx
        !****************************
        do kk = 1,nk+1
            element(i)%xgl(kk) = x(i) + gau_lob(kk,1)*dx
        enddo
    enddo

!initial L2-projection
    ut1=0.
    ut2=0.
    ut21=0.
    ut22=0.
    do i = 1 , nx
        u21temp = 0.0
        u22temp = 0.0
        do kk = 0,nk
            u1temp = 0.0
            u2temp = 0.0
            do lx = 1 , 6
                u1temp=u1temp+exact1(x(i)+xg(lx)*dx,0.)*fle(kk,xg(lx))*wg(lx)
                u2temp=u2temp+exact2(x(i)+xg(lx)*dx,0.)*fle(kk,xg(lx))*wg(lx)
                !write(*,*)exact1(x(i)+xg(lx)*dx,0.)+2*exact2(x(i)+xg(lx)*dx,0.)
            enddo
            element(i)%u1modal(kk+1,0) = u1temp*ai(kk+1)
            element(i)%u2modal(kk+1,0) = u2temp*ai(kk+1)
            !write(*,*)u1temp,element(i)%u1modal(kk+1,0)
        enddo
        ut1=ut1+element(i)%u1modal(1,0)
        ut2=ut2+element(i)%u2modal(1,0)
        do lx = 1 , 6
            u21temp=u21temp+(exact1(x(i)+xg(lx)*dx,0.)**2.)*wg(lx)
            u22temp=u22temp+((aa(x(i)+xg(lx)*dx)*exact2(x(i)+xg(lx)*dx,0.))**2.)*wg(lx)
                !write(*,*)exact1(x(i)+xg(lx)*dx,0.)+2*exact2(x(i)+xg(lx)*dx,0.)
        enddo
        
        !element(i)%weno_t(1)=fp(element(i)%umodal(1),x(i))
        !!element(i)%weno_t(1)=1.+sin(x(i))*dx
        ut21=ut21+u21temp*dx
        ut22=ut22+u22temp*dx
        !write(*,*) x(i),element(i)%u1modal(1,0),element(i)%u2modal(1,0)
        !write(*,*) x(i),element(i)%u1modal(1,0),ut1!,element(i)%u2modal(1,0)
    enddo
    !write(*,*)ut1,ut2,dx
    u10=ut1*dx
    u20=ut2*dx
    u210=ut21+ut22
    !write(*,*)u10,u20,dx
    end subroutine init
    
 !********************************   
    subroutine computeLint
    implicit none
    real :: u1temp,u2temp,ut1,ut2,u21temp,u22temp,ut21,ut22
    integer :: lx,kk

  
!initial L2-energy
    ut21=0.
    ut22=0.
    do i = 1 , nx
        u21temp = 0.0
        u22temp = 0.0
        do lx = 1 , 6
            u21temp=u21temp+(ortho_poly1d(element(i)%u1modal(1:nk+1,irk),x(i)+xg(lx)*dx ,x(i) ,dx ,nk)**2.)*wg(lx)
            u22temp=u22temp+((aa(x(i)+xg(lx)*dx)*ortho_poly1d(element(i)%u2modal(1:nk+1,irk),x(i)+xg(lx)*dx ,x(i) ,dx ,nk))**2.)*wg(lx)
                !write(*,*)exact1(x(i)+xg(lx)*dx,0.)+2*exact2(x(i)+xg(lx)*dx,0.)
        enddo
        ut21=ut21+u21temp*dx
        ut22=ut22+u22temp*dx
        !write(*,*) x(i),element(i)%u1modal(1,0),element(i)%u2modal(1,0)
        !write(*,*) x(i),element(i)%u1modal(1,0),ut1!,element(i)%u2modal(1,0)
    enddo
    !write(*,*)ut1,ut2,dx
    uet=ut21+ut22
    !write(*,*)u10,u20,dx
    end subroutine computeLint

    !********************************
    real function fle(k,x)
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x
    ! 关于Legendre正交多项式的函数
    ! v_0=1.0, v_1=(x-x_j)/dx_j, v_2=( (x-x_j)/dx_j )^2-1.0/12.0,
    ! v_3=( (x-x_j)/dx_j )^3- 0.15*(x-x_j)/dx_j
    ! v_4=( ((x-x_j)/dx_j)^2-3.0/14.0 )*((x-x_j)/dx_j)^2 + 3.0/560.0

    if(k.eq.0) then
        fle=1.0
    elseif (k.eq.1) then
        fle= x
    elseif(k.eq.2) then
        fle= x*x-1./12.0
    elseif (k.eq.3) then
        fle=x*x*x-0.15*x
    elseif(k.eq.4) then
        fle=(x**2-3./14.0)*x*x+3.0/560.0
    endif

    return
    end function fle
