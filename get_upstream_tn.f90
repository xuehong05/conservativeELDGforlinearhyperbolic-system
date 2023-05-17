    !*************************************************************************
    ! if we let speed =0, that is irkdg == 1
    ! and set the solution as initial condition at background elements,
    ! the code reduces to the normal RKDG method.
    !*************************************************************************
    subroutine get_upstream_tn
    implicit none

    type(element1d_upstream), pointer :: p,p1
    type(element1d_dynamic), pointer :: pet1,pet2
    
    real :: xll,xrr
    real :: r

    !do i=1,nx
    !    !!write(*,*)i,nx
    !    !!element(i)%weno_t(1)=fp(element(i)%umodal(1),x(i))
    !    !!eigenvalue at x_j
    !    !element(i)%weno_t(1)=1
    !    !element(i)%weno_t(2)=-1
    !    !element(i)%R_j(1,1)=-1.
    !    !element(i)%R_j(1,2)=1.
    !    !element(i)%R_j(2,1)=1.
    !    !element(i)%R_j(2,2)=1.
    !    !element(i)%RT_j(1,1)=-0.5
    !    !element(i)%RT_j(1,2)=0.5
    !    !element(i)%RT_j(2,1)=0.5
    !    !element(i)%RT_j(2,2)=0.5
    !    !!
    !    
    !    !write(*,*)i,nx
    !    !element(i)%weno_t(1)=fp(element(i)%umodal(1),x(i))
    !    !eigenvalue at x_j
    !    !element(i)%weno_t(1)= aa(x(i))
    !    element(i)%weno_t(1)= aa(x(i))+dx
    !    element(i)%weno_t(2)=-1.*element(i)%weno_t(1)
    !    element(i)%R_j(1,1)=-1.*element(i)%weno_t(1)
    !    element(i)%R_j(1,2)=element(i)%weno_t(1)
    !    element(i)%R_j(2,1)=1.
    !    element(i)%R_j(2,2)=1.
    !    element(i)%RT_j(1,1)=-0.5/element(i)%weno_t(1)
    !    element(i)%RT_j(1,2)=0.5
    !    element(i)%RT_j(2,1)=0.5/element(i)%weno_t(1)
    !    element(i)%RT_j(2,2)=0.5
    !    
    !    ! !write(*,*)i,nx
    !    !!RK
    !    !!element(i)%weno_t(1)=fp(element(i)%umodal(1),x(i))
    !    !!eigenvalue at x_j
    !    !element(i)%weno_t(1)= 0.
    !    !element(i)%weno_t(2)= 0.
    !    !element(i)%R_j(1,1)= 1.
    !    !element(i)%R_j(1,2)= 0.
    !    !element(i)%R_j(2,1)= 0.
    !    !element(i)%R_j(2,2)= 1.
    !    !element(i)%RT_j(1,1)= 1.
    !    !element(i)%RT_j(1,2)= 0.
    !    !element(i)%RT_j(2,1)= 0.
    !    !element(i)%RT_j(2,2)= 1.
    !    
    !    !!!add eigenvalue to a perturbation:A_j,adjoint problem
    !    !!element(i)%weno_t(1)=1.+sin(x(i))*dx**2.
    !    !!element(i)%weno_t(2)=-1.-sin(x(i))*dx**2.
    !    !element(i)%weno_t(1)=1.+sin(x(i))*dx
    !    !element(i)%weno_t(2)=-1.-sin(x(i))*dx
    !    !!element(i)%weno_t(1)=1.+0.5*sin(x(i))
    !    !!element(i)%weno_t(2)=-1.-0.5*sin(x(i))
    !    !!element(i)%weno_t(1)=1.+0.5
    !    !!element(i)%weno_t(2)=-1.-0.5
    !    !!element(i)%weno_t(1)=1.+dx**2
    !    !!element(i)%weno_t(2)=-1.-dx**2
    !    !!element(i)%weno_t(1)=1.+dx
    !    !!element(i)%weno_t(2)=-1.-dx
    !    !element(i)%R_j(1,1)=element(i)%weno_t(2)
    !    !element(i)%R_j(1,2)=element(i)%weno_t(1)
    !    !element(i)%R_j(2,1)=1.
    !    !element(i)%R_j(2,2)=1.
    !    !element(i)%RT_j(1,1)=0.5/element(i)%weno_t(2)
    !    !element(i)%RT_j(1,2)=0.5
    !    !element(i)%RT_j(2,1)=0.5/element(i)%weno_t(1)
    !    !element(i)%RT_j(2,2)=0.5
    !enddo
    
    call boundary_eulerian(0)

    do i = 1,nx+1
        call st_speed1(element(i-1)%u1modal(1,0),element(i)%u1modal(1,0),vertex(i)%coor,speed1(i))
        call st_speed2(element(i-1)%u1modal(1,0),element(i)%u1modal(1,0),vertex(i)%coor,speed2(i))
        !if(irkdg==1)then
        !    speed(i) = 0.
        !endif
    enddo

    do iio = 1,irk
        do io = 0,iio
        do i = 1,nx
            !dtt = time-(time-dt+rk_d(io)*dt )
            dtt = (rk_d(iio)-rk_d(io))*dt
            pet1 => element_t1(i,io,iio)
            pet2 => element_t2(i,io,iio)
            p => element_star1(i,io,iio)
            p1 => element_star2(i,io,iio)
            xll = element(i)%xl
            xrr = element(i)%xr
            pet1%xl = vertex(i)%coor-speed1(i)*dtt
            pet1%xr = vertex(i+1)%coor-speed1(i+1)*dtt
            pet1%xc = (pet1%xl+pet1%xr)/2.
            pet1%dx = pet1%xr-pet1%xl
            pet1%midtime = time-dt+rk_d(io)*dt
            !pet1%weno_t(1)=x(i)-element(i)%weno_t(1)*dtt
            pet2%xl = vertex(i)%coor-speed2(i)*dtt
            pet2%xr = vertex(i+1)%coor-speed2(i+1)*dtt
            pet2%xc = (pet2%xl+pet2%xr)/2.
            pet2%dx = pet2%xr-pet2%xl
            pet2%midtime = time-dt+rk_d(io)*dt
            !pet2%weno_t(1)=x(i)-element(i)%weno_t(2)*dtt
            !p%point_origin%coor=pet1%xl
            element_star1(i,io,iio)%point_origin%coor=pet1%xl
            p%point_end%coor=pet1%xr
            !p1%point_origin%coor=pet2%xl
            element_star2(i,io,iio)%point_origin%coor=pet2%xl
            p1%point_end%coor=pet2%xr
            p%point_origin%id=ceiling( (pet1%xl-xleft)/dx )
            p%point_end%id=ceiling( (pet1%xr-xleft)/dx )
            p1%point_origin%id=ceiling( (pet2%xl-xleft)/dx )
            p1%point_end%id=ceiling( (pet2%xr-xleft)/dx )
            !call get_element_tmass(pet1,xll,xrr)
            !call get_element_tmass(pet2,xll,xrr)
        enddo
         !write(*,*) dtt
         !pause
        enddo
    enddo

    !getting the initial solution on the upstream at tn
    !io = 0
    !do i = 1,nx+1
    !    vertex_star1(i)%coor = vertex(i)%coor - speed1(i)*dt
    !    vertex_star2(i)%coor = vertex(i)%coor - speed2(i)*dt
    !    vertex_star1(i)%id = ceiling( (vertex_star1(i)%coor-xleft)/dx )
    !    vertex_star2(i)%id = ceiling( (vertex_star2(i)%coor-xleft)/dx )
    !enddo


    !alpha =0.
    !do i = 1,nx+1
    !    alpha = max( alpha, abs(fp( element(i-1)%umodal(1)*0.5+element(i)%umodal(1)*0.5,vertex(i)%coor ) - speed1(i)) )
    !    !alpha = max( alpha, abs(fp( element(i-1)%umodal(1)*0.5+element(i)%umodal(1)*0.5,vertex(i)%coor ) ) )
    !enddo
    !!alpha=0.1


    !write(11,*) alpha
    

    call projection(0,1 )
    !do i=1,nx
        
    !do i = 1,nx
    !    ! consider an upstream element
    !    p => element_star1(i,0,1)
    !    p1 => element_star2(i,0,1)
    !    !p%point_origin = vertex_star1(i)
    !    !p%point_end = vertex_star1(i+1)
    !    !p1%point_origin = vertex_star2(i)
    !    !p1%point_end = vertex_star2(i+1)
    !
    !    if( p%point_origin%coor > p%point_end%coor )then
    !        print *,'the first characteristic collide!', vertex_star1(i)%coor , vertex_star1(i+1)%coor
    !        pause
    !    endif
    !    
    !    if( p1%point_origin%coor > p1%point_end%coor )then
    !        print *,'the secoond characteristic collide!', vertex_star2(i)%coor , vertex_star2(i+1)%coor
    !        pause
    !    endif
    !
    !    call search_segment(p)
    !    call search_segment(p1)
    !    xll = element(i)%xl
    !    xrr = element(i)%xr
    !    !reconstruct u_h^* at T_n 
    !    call get_integral_pk(p,p1,0,1,xll,xrr)
    !    !write(*,*) element_t(i,0)%umodal(1:nk+1)
    !enddo
    !
    
    !write(*,*) element_t(i,0)%umodal(1:nk+1)
    !end of getting the initial solution on the upstream at tn
    end subroutine get_upstream_tn
    
    !*************************************************************
    subroutine search_segment(p)
    implicit none
    type(element1d_upstream), pointer :: p
    integer :: mx,kk
    integer :: inter

    mx = p%point_end%id - p%point_origin%id

    p%point_inter(0) = p%point_origin
    p%point_inter(1+mx) = p%point_end
    p%nsub = mx+1

    if(mx .ne. 0)then
        do kk = 1 , mx
            inter = p%point_origin%id + kk
            p%point_inter(kk)%coor = xgrid( inter )
            p%point_inter(kk)%id = inter
        enddo
    endif
    !
    do kk = 1 , 1 + mx
        p%segment(kk)%porigin = p%point_inter(kk-1)
        p%segment(kk)%pend = p%point_inter(kk)
        p%segment(kk)%id = p%point_inter(kk-1)%id
    enddo

    end subroutine search_segment
    !*******************************