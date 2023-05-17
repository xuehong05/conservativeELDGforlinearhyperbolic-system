    subroutine projection(io,iio )
    implicit none
    
    integer, intent(in) :: io,iio
    
    type(element1d_upstream), pointer :: p,p1
    type(element1d_dynamic), pointer :: pet,pet1
    
    real :: xll,xrr
    real :: u1l,u2l,u1r,u2r,u1bar,u2bar,u1jump,u2jump
    !integer :: i
    
    !real :: dx_t,dx
    do i = 1,nx
             ! consider an upstream element
        p => element_star1(i,io,iio)
        p1 => element_star2(i,io,iio)
        
        if( p%point_origin%coor > p%point_end%coor )then
            print *,'the first characteristic collide!', vertex_star1(i)%coor , vertex_star1(i+1)%coor
            pause
        endif
        
        if( p1%point_origin%coor > p1%point_end%coor )then
            print *,'the secoond characteristic collide!', vertex_star2(i)%coor , vertex_star2(i+1)%coor
            pause
        endif
        
        call search_segment(p)
        call search_segment(p1)
        xll = element(i)%xl
        xrr = element(i)%xr
        call get_umodal_dynamic(p,p1,io,iio,xll,xrr)
        
    enddo
    !write(*,*)io,iio
    call boundary_lag(io,iio)
    do i = 1,nx+1
             ! consider an upstream element
        alpha = 1.1*max( abs(fp1(element_t1(i-1,io,iio)%xr)-speed1(i)), &
            abs(fp2(element_t2(i,io,iio)%xl)-speed2(i) ) )
        pet=>element_t1(i-1,io,iio)
        u1l = ortho_poly1d(pet%u1modal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
        u2l = ortho_poly1d(pet%u2modal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
        pet=>element_t1(i,io,iio)
        u1r = ortho_poly1d(pet%u1modal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
        u2r = ortho_poly1d(pet%u2modal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
        
        u1bar=0.5*(u1l+u1r)
        u2bar=0.5*(u2l+u2r)
        u1jump=u1r-u1l
        u2jump=u2r-u2l
        flux_LF(i)=RT11(pet%xl)&
            *((A11(pet%xl)-speed1(i))*u1bar+A12(pet%xl)*u2bar-0.5*alpha*u1jump)&
            +RT12(pet%xl)*((A21(pet%xl))*u1bar+(A22(pet%xl)-speed1(i))*u2bar-0.5*alpha*u2jump)
        flux_LF11(i)=R11(pet%xl)*flux_LF(i)
        flux_LF21(i)=R21(pet%xl)*flux_LF(i)
        
        pet1=>element_t2(i-1,io,iio)
        u1l = ortho_poly1d(pet1%u1modal(1:nk+1),pet1%xr,pet1%xc,pet1%dx,nk)
        u2l = ortho_poly1d(pet1%u2modal(1:nk+1),pet1%xr,pet1%xc,pet1%dx,nk)
        pet1=>element_t2(i,io,iio)
        u1r = ortho_poly1d(pet1%u1modal(1:nk+1),pet1%xl,pet1%xc,pet1%dx,nk)
        u2r = ortho_poly1d(pet1%u2modal(1:nk+1),pet1%xl,pet1%xc,pet1%dx,nk)
        u1bar=0.5*(u1l+u1r)
        u2bar=0.5*(u2l+u2r)
        u1jump=u1r-u1l
        u2jump=u2r-u2l
        flux_LFp(i)=RT21(pet1%xl)&
            *((A11(pet1%xl)-speed2(i))*u1bar+A12(pet1%xl)*u2bar-0.5*alpha*u1jump)&
            +RT22(pet1%xl)*((A21(pet1%xl))*u1bar+(A22(pet1%xl)-speed2(i))*u2bar-0.5*alpha*u2jump)
        flux_LF12(i)=R12(pet1%xl)*flux_LFp(i)
        flux_LF22(i)=R22(pet1%xl)*flux_LFp(i)
        
        !if(abs(p%point_origin%coor-xgrid(p%point_origin%id))<10.**(-10.))then
        !flux_LF(i)=(RT11(p%point_origin%coor)&
        !    *((A11(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id-1)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)+A12(p%point_origin%coor)*ortho_poly1d(element(p%point_origin%id-1)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk))&
        !    +RT12(p%point_origin%coor)*((A21(p%point_origin%coor))*ortho_poly1d(element(p%point_origin%id-1)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)+(A22(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id-1)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)))
        !else
        !flux_LF(i)=(RT11(p%point_origin%coor)&
        !    *((A11(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+A12(p%point_origin%coor)*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk))&
        !    +RT12(p%point_origin%coor)*((A21(p%point_origin%coor))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+(A22(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)))
        !endif
        !flux_LF11(i)=R11(p%point_origin%coor)*flux_LF(i)
        !flux_LF21(i)=R21(p%point_origin%coor)*flux_LF(i)
        !write(*,*)flux_LF(i)
        !if(abs(p1%point_origin%coor-xgrid(p1%point_origin%id+1))<10.**(-10.))then
        !flux_LFp(i)=(RT21(p1%point_origin%coor)&
        !    *((A11(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id+1)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)+A12(p1%point_origin%coor)*ortho_poly1d(element(p1%point_origin%id+1)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk))&
        !    +RT22(p1%point_origin%coor)*((A21(p1%point_origin%coor))*ortho_poly1d(element(p1%point_origin%id+1)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)+(A22(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id+1)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)))
        !else
        !flux_LFp(i)=(RT21(p1%point_origin%coor)&
        !    *((A11(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+A12(p1%point_origin%coor)*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk))&
        !    +RT22(p1%point_origin%coor)*((A21(p1%point_origin%coor))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+(A22(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)))
        !endif
        !flux_LF12(i)=R12(p1%point_origin%coor)*flux_LFp(i)
        !flux_LF22(i)=R22(p1%point_origin%coor)*flux_LFp(i)
        !write(*,*)flux_LFp(i)
    enddo
   
    !flux_LF(nx+1)= flux_LF(1)
    !flux_LFp(nx+1)=flux_LFp(1)
    !flux_LF11(nx+1)=flux_LF11(1)
    !flux_LF21(nx+1)=flux_LF21(1)
    !flux_LF12(nx+1)=flux_LF12(1)
    !flux_LF22(nx+1)=flux_LF22(1)
       
    do i = 1,nx
             ! consider an upstream element
        p => element_star1(i,io,iio)
        p1 => element_star2(i,io,iio)
        xll = element(i)%xl
        xrr = element(i)%xr
        !reconstruct u_h^* at T_n 
        call get_integral_pk(p,p1,io,iio,xll,xrr)
        !write(*,*) i,io,iio,element_t1(i,io,iio)%u1modal(1:nk+1),element_t1(i,io,iio)%u2modal(1:nk+1)
        !pause
    enddo
    !write(*,*) element_t(i,0)%umodal(1:nk+1)
    !end of getting the initial solution on the upstream at tn
    
    end subroutine projection
