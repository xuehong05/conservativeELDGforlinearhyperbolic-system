    subroutine init_x_1D
    implicit none
    !real :: utemp
    !integer :: lx,kk
   
    xleft=xleft_2d
    xright=xright_2d
    nx=nx_2d
    dx=dx_2d
    x(:)=x_2d(:)
    Xgrid(:)=Xgrid_2d(:)
   
    do i = 1,nx+1
        vertex(i)%coor = xleft + (i-1) * dx
    enddo
    
    do i = 1,nx
        element(i)%xl = x(i)-0.5*dx
        element(i)%xr = x(i)+0.5*dx
        do kk = 1,nk+1
            element(i)%xgl(kk) = x(i) + gau_lob(kk,1)*dx
        enddo
        !element(i)%weno_t(1)=fp(ixy,element(i)%umodal(1),x(i),y_split,time)
    enddo
    end subroutine init_x_1D
   !********** 
    
    subroutine init_y_1D
    implicit none
    !real :: utemp
    !integer :: lx,kk
   
    xleft=yleft_2d
    xright=yright_2d
    nx=ny_2d
    dx=dy_2d
    x(:)=y_2d(:)
    Xgrid(:)=ygrid_2d(:)
   
    do i = 1,nx+1
        vertex(i)%coor = xleft + (i-1) * dx
    enddo
    
    do i = 1,nx
        element(i)%xl = x(i)-0.5*dx
        element(i)%xr = x(i)+0.5*dx
        do kk = 1,nk+1
            element(i)%xgl(kk) = x(i) + gau_lob(kk,1)*dx
        enddo
        !element(i)%weno_t(1)=fp(ixy,element(i)%umodal(1),x_split,x(i),time)
    enddo
    end subroutine init_y_1D
    