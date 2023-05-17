    !
    subroutine get_integral_pk(p,p1,io,iio,xl,xr)
    implicit none
    
    type(element1d_upstream), pointer :: p,p1
    integer,intent(in) :: io,iio
    real,intent(in) :: xl,xr
    
    integer :: kk
    real :: sum,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9
    integer :: idx
    real :: xi_end,xi_origin,term_1,term_2
    
    integer :: nm
    real ::  bb(1:nk+1)
    
    real :: amatrix(1:nk+1,1:nk+1),store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1)
    integer :: ig
    
    !real :: gau_t(1:nk+1)
    real :: gau_t(1:6)
    real :: gau_tt(1:3),st,st1,st2,st3,st4,st5,st6,st7,st8,st9
    
    real :: dx_t,dx,lamda
    
    integer :: ido,ide
    real :: xs_o,xs_e,xc_star,xi(5)
    integer :: igl
    
    
    !**************************
    ido = p%point_origin%id
    ide = p%point_end%id
    xs_o = p%point_origin%coor
    xs_e = p%point_end%coor
    !**************************
    p%aa(1,1) = 1.;  p%aa(2,1) = 0.;  p%aa(3,1) = 0.;
    p%aa(1,2) = 0.;  p%aa(2,2) = 1.;  p%aa(3,2) = 0.;
    p%aa(1,3) = 0.;  p%aa(2,3) = 0.;  p%aa(3,3) = 1.;
    
    p1%aa(1,1) = 1.;  p1%aa(2,1) = 0.;  p1%aa(3,1) = 0.;
    p1%aa(1,2) = 0.;  p1%aa(2,2) = 1.;  p1%aa(3,2) = 0.;
    p1%aa(1,3) = 0.;  p1%aa(2,3) = 0.;  p1%aa(3,3) = 1.;
    do nm = 1,nk+1
        
        !*********
    
        sum = 0.
        sum1= 0.
        sum2 = 0.
        sum3= 0.
        sum4= 0.
        sum5= 0.
        sum6= 0.
        sum7= 0.
        !sum8= 0.
        !sum9= 0.
        dx= xr-xl
        do kk = 1 ,p%nsub
            idx = p%segment(kk)%id
    
            gau_t(1:6) = (p%segment(kk)%pend%coor + p%segment(kk)%porigin%coor)/2. &
                + xg(1:6)* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            
            dx_t = p%point_end%coor-p%point_origin%coor
            xc_star = (p%point_origin%coor+p%point_end%coor)/2.
            st =0.
            st1=0.
            st2=0.
            st3=0.
            st4=0.
            st5=0.
            st6=0.
            st7=0.
            !st8=0.
            !st9=0.
            do ig = 1,6
                lamda=wave_linear_j(gau_t(ig),speed1(i),speed1(i+1),p%point_origin%coor,p%point_end%coor,dx_t)
                st = st +R11(gau_t(ig))* (RT11(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st1 = st1 + R21(gau_t(ig))* (RT11(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st2 = st2 +R11(gau_t(ig))* (RT11(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                !    +RT12(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                !    *diff_ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st3 = st3 +R21(gau_t(ig))* (RT11(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                !    +RT12(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                !    *diff_ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st4 = st4 +(-0.5*aaderi(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st5 = st5 +((0.5*aaderi(gau_t(ig))/(aa1(gau_t(ig))**2.))*((A11(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A12(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st6= st6 +R11(gau_t(ig))* RT11(gau_t(ig))*f_source(gau_t(ig),element_t1(i,io,iio)%midtime)&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st7= st7 +R21(gau_t(ig))* RT11(gau_t(ig))*f_source(gau_t(ig),element_t1(i,io,iio)%midtime)&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st8=st8 +ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                !    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st9=st9 +ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                !    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
            enddo
    
            sum = sum + st* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum1 = sum1 + st1* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            !sum2 = sum2 + st2* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            !sum3 = sum3 + st3* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum4 = sum4 + st4* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum5 = sum5 + st5* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum6 = sum6 + st6* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum7 = sum7 + st7* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            !sum8 = sum8 + st8* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            !sum9 = sum9 + st9* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
        enddo
        element_t1(i,io,iio)%sum_int(nm,1)=sum
        element_t1(i,io,iio)%sum_int(nm,2)=sum1
        !element_t1(i,io,iio)%phi_int(nm,1)=sum2
        !element_t1(i,io,iio)%phi_int(nm,2)=sum3
        element_t1(i,io,iio)%vector_int(nm,1)=sum4
        element_t1(i,io,iio)%vector_int(nm,2)=sum5
        element_t1(i,io,iio)%source_int(nm,1)=sum6
        element_t1(i,io,iio)%source_int(nm,2)=sum7
        !element_t1(i,io,iio)%u1modal(nm)=sum8*ai(nm)/(p%point_end%coor-p%point_origin%coor)
        !element_t1(i,io,iio)%u2modal(nm)=sum9*ai(nm)/(p%point_end%coor-p%point_origin%coor)
        !if(abs(p%point_end%coor-xgrid(p%point_end%id))<10.**(-10.))then
        !term_1=(RT11(p%point_end%coor)&
        !    *((A11(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id-1)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id-1) ,dx ,nk)+A12(p%point_end%coor)*ortho_poly1d(element(p%point_end%id-1)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id-1) ,dx ,nk))&
        !    +RT12(p%point_end%coor)*((A21(p%point_end%coor))*ortho_poly1d(element(p%point_end%id-1)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id-1) ,dx ,nk)+(A22(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id-1)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id-1) ,dx ,nk)))&
        !    *ortho_poly1d( p%aa(1:nk+1,nm),p%point_end%coor,xc_star ,dx_t,nk)
        !else
        !term_1=(RT11(p%point_end%coor)&
        !    *((A11(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)+A12(p%point_end%coor)*ortho_poly1d(element(p%point_end%id)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk))&
        !    +RT12(p%point_end%coor)*((A21(p%point_end%coor))*ortho_poly1d(element(p%point_end%id)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)+(A22(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)))&
        !    *ortho_poly1d( p%aa(1:nk+1,nm),p%point_end%coor,xc_star ,dx_t,nk)
        !endif
        !if(abs(p%point_origin%coor-xgrid(p%point_origin%id))<10.**(-10.))then
        !term_2=(RT11(p%point_origin%coor)&
        !    *((A11(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id-1)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)+A12(p%point_origin%coor)*ortho_poly1d(element(p%point_origin%id-1)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk))&
        !    +RT12(p%point_origin%coor)*((A21(p%point_origin%coor))*ortho_poly1d(element(p%point_origin%id-1)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)+(A22(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id-1)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id-1) ,dx ,nk)))&
        !    *ortho_poly1d( p%aa(1:nk+1,nm),p%point_origin%coor,xc_star ,dx_t,nk)
        !else
        !term_2=(RT11(p%point_origin%coor)&
        !    *((A11(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+A12(p%point_origin%coor)*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk))&
        !    +RT12(p%point_origin%coor)*((A21(p%point_origin%coor))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+(A22(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)))&
        !    *ortho_poly1d( p%aa(1:nk+1,nm),p%point_origin%coor,xc_star ,dx_t,nk)
        !endif
        !element_t1(i,io,iio)%flux1(nm,1)=R11(p%point_end%coor)*(term_1)
        !element_t1(i,io,iio)%flux0(nm,1)=R11(p%point_origin%coor)*(term_2)
        !element_t1(i,io,iio)%flux1(nm,2)=R21(p%point_end%coor)*(term_1)
        !element_t1(i,io,iio)%flux0(nm,2)=R21(p%point_origin%coor)*(term_2)
        element_t1(i,io,iio)%flux1(nm,1)=flux_LF11(i+1)*ortho_poly1d( p%aa(1:nk+1,nm),p%point_end%coor,xc_star ,dx_t,nk)
        element_t1(i,io,iio)%flux0(nm,1)=flux_LF11(i)*ortho_poly1d( p%aa(1:nk+1,nm),p%point_origin%coor,xc_star ,dx_t,nk)
        element_t1(i,io,iio)%flux1(nm,2)=flux_LF21(i+1)*ortho_poly1d( p%aa(1:nk+1,nm),p%point_end%coor,xc_star ,dx_t,nk)
        element_t1(i,io,iio)%flux0(nm,2)=flux_LF21(i)*ortho_poly1d( p%aa(1:nk+1,nm),p%point_origin%coor,xc_star ,dx_t,nk)
        !second charicteristic domain******************************************************************************************************************
        sum = 0.
        sum1= 0.
        sum2 = 0.
        sum3= 0.
        sum4= 0.
        sum5= 0.
        sum6= 0.
        sum7= 0.
        do kk = 1 ,p1%nsub
            idx = p1%segment(kk)%id
    
             gau_t(1:6) = (p1%segment(kk)%pend%coor + p1%segment(kk)%porigin%coor)/2. &
                + xg(1:6)* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            
            dx_t = p1%point_end%coor-p1%point_origin%coor
            xc_star = (p1%point_origin%coor+p1%point_end%coor)/2.
            st =0.
            st1=0.
            st2=0.
            st3=0.
            st4=0.
            st5=0.
            st6=0.
            st7=0.
            do ig = 1,6
                lamda=wave_linear_j(gau_t(ig),speed2(i),speed2(i+1),p1%point_origin%coor,p1%point_end%coor,dx_t)
                st = st +R12(gau_t(ig))* (RT21(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT22(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st1 = st1 + R22(gau_t(ig))* (RT21(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT22(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st2 = st2 +R12(gau_t(ig))* (RT21(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                !    +RT22(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                !    *diff_ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                !st3 = st3 +R22(gau_t(ig))* (RT21(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
                !    +RT22(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                !    *diff_ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st4 = st4 +(0.5*aaderi(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st5 = st5 +((-0.5*aaderi(gau_t(ig))/(aa1(gau_t(ig))**2.))*((A11(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A12(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st6= st6 +R12(gau_t(ig))* RT21(gau_t(ig))*f_source(gau_t(ig),element_t2(i,io,iio)%midtime)&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st7= st7 +R22(gau_t(ig))* RT21(gau_t(ig))*f_source(gau_t(ig),element_t2(i,io,iio)%midtime)&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
            enddo
    
            sum = sum + st* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum1 = sum1 + st1* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            !sum2 = sum2 + st2* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            !sum3 = sum3 + st3* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum4 = sum4 + st4* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum5 = sum5 + st5* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum6 = sum6 + st6* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum7 = sum7 + st7* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
        enddo
        element_t2(i,io,iio)%sum_int(nm,1)=sum
        element_t2(i,io,iio)%sum_int(nm,2)=sum1
        !element_t2(i,io,iio)%phi_int(nm,1)=sum2
        !element_t2(i,io,iio)%phi_int(nm,2)=sum3
        element_t2(i,io,iio)%vector_int(nm,1)=sum4
        element_t2(i,io,iio)%vector_int(nm,2)=sum5
        element_t2(i,io,iio)%source_int(nm,1)=sum6
        element_t2(i,io,iio)%source_int(nm,2)=sum7
        !if(abs(p1%point_end%coor-xgrid(p1%point_end%id+1))<10.**(-10.))then
        !term_1=(RT21(p1%point_end%coor)&
        !    *((A11(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id+1)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id+1) ,dx ,nk)+A12(p1%point_end%coor)*ortho_poly1d(element(p1%point_end%id+1)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id+1) ,dx ,nk))&
        !    +RT22(p1%point_end%coor)*((A21(p1%point_end%coor))*ortho_poly1d(element(p1%point_end%id+1)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id+1) ,dx ,nk)+(A22(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id+1)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id+1) ,dx ,nk)))&
        !    *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_end%coor,xc_star ,dx_t,nk)
        !else
        !term_1=(RT21(p1%point_end%coor)&
        !    *((A11(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)+A12(p1%point_end%coor)*ortho_poly1d(element(p1%point_end%id)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk))&
        !    +RT22(p1%point_end%coor)*((A21(p1%point_end%coor))*ortho_poly1d(element(p1%point_end%id)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)+(A22(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)))&
        !    *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_end%coor,xc_star ,dx_t,nk)
        !endif
        !if(abs(p1%point_origin%coor-xgrid(p1%point_origin%id+1))<10.**(-10.))then
        !term_2=(RT21(p1%point_origin%coor)&
        !    *((A11(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id+1)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)+A12(p1%point_origin%coor)*ortho_poly1d(element(p1%point_origin%id+1)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk))&
        !    +RT22(p1%point_origin%coor)*((A21(p1%point_origin%coor))*ortho_poly1d(element(p1%point_origin%id+1)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)+(A22(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id+1)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id+1) ,dx ,nk)))&
        !    *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_origin%coor,xc_star ,dx_t,nk)
        !else
        !term_2=(RT21(p1%point_origin%coor)&
        !    *((A11(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+A12(p1%point_origin%coor)*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk))&
        !    +RT22(p1%point_origin%coor)*((A21(p1%point_origin%coor))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+(A22(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)))&
        !    *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_origin%coor,xc_star ,dx_t,nk)
        !endif
        !
        !element_t2(i,io,iio)%flux1(nm,1)=R12(p1%point_end%coor)*(term_1)
        !element_t2(i,io,iio)%flux0(nm,1)=R12(p1%point_origin%coor)*(term_2)
        !element_t2(i,io,iio)%flux1(nm,2)=R22(p1%point_end%coor)*(term_1)
        !element_t2(i,io,iio)%flux0(nm,2)=R22(p1%point_origin%coor)*(term_2)
        
        element_t2(i,io,iio)%flux1(nm,1)=flux_LF12(i+1)*ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_end%coor,xc_star ,dx_t,nk)
        element_t2(i,io,iio)%flux0(nm,1)=flux_LF12(i)*ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_origin%coor,xc_star ,dx_t,nk)
        element_t2(i,io,iio)%flux1(nm,2)=flux_LF22(i+1)*ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_end%coor,xc_star ,dx_t,nk)
        element_t2(i,io,iio)%flux0(nm,2)=flux_LF22(i)*ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_origin%coor,xc_star ,dx_t,nk)
    enddo
    
    
    end subroutine get_integral_pk
    
    
    subroutine get_umodal_dynamic(p,p1,io,iio,xl,xr)
    implicit none
    
    type(element1d_upstream), pointer :: p,p1
    integer,intent(in) :: io,iio
    real,intent(in) :: xl,xr
    
    integer :: kk
    real :: sum,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9
    integer :: idx
    real :: xi_end,xi_origin,term_1,term_2
    
    integer :: nm
    real ::  bb(1:nk+1)
    
    real :: amatrix(1:nk+1,1:nk+1),store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1)
    integer :: ig
    
    !real :: gau_t(1:nk+1)
    real :: gau_t(1:6),gau_t1(1:6),xcc,xcc1
    real :: gau_tt(1:3),st,st1,st2,st3,st4,st5,st6,st7,st8,st9
    
    real :: dx_t,dx,lamda,dx1_t,lamda1
    
    integer :: ido,ide
    real :: xs_o,xs_e,xc_star,xi(5)
    integer :: igl
    
    
    !**************************
    ido = p%point_origin%id
    ide = p%point_end%id
    xs_o = p%point_origin%coor
    xs_e = p%point_end%coor
    !**************************
    p%aa(1,1) = 1.;  p%aa(2,1) = 0.;  p%aa(3,1) = 0.;
    p%aa(1,2) = 0.;  p%aa(2,2) = 1.;  p%aa(3,2) = 0.;
    p%aa(1,3) = 0.;  p%aa(2,3) = 0.;  p%aa(3,3) = 1.;
    
    p1%aa(1,1) = 1.;  p1%aa(2,1) = 0.;  p1%aa(3,1) = 0.;
    p1%aa(1,2) = 0.;  p1%aa(2,2) = 1.;  p1%aa(3,2) = 0.;
    p1%aa(1,3) = 0.;  p1%aa(2,3) = 0.;  p1%aa(3,3) = 1.;
    do nm = 1,nk+1
        
        !*********
    
       
        sum8= 0.
        sum9= 0.
        dx= xr-xl
        do kk = 1 ,p%nsub
            idx = p%segment(kk)%id
    
            gau_t(1:6) = (p%segment(kk)%pend%coor + p%segment(kk)%porigin%coor)/2. &
                + xg(1:6)* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            
            dx_t = p%point_end%coor-p%point_origin%coor
            xc_star = (p%point_origin%coor+p%point_end%coor)/2.
            
            st8=0.
            st9=0.
            do ig = 1,6
               
                st8=st8 +ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st9=st9 +ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
            enddo
    
            sum8 = sum8 + st8* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
            sum9 = sum9 + st9* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
        enddo
        
        element_t1(i,io,iio)%u1modal(nm)=sum8*ai(nm)/(p%point_end%coor-p%point_origin%coor)
        element_t1(i,io,iio)%u2modal(nm)=sum9*ai(nm)/(p%point_end%coor-p%point_origin%coor)
       
        !second charicteristic domain******************************************************************************************************************
        
        sum8= 0.
        sum9= 0.
        do kk = 1 ,p1%nsub
            idx = p1%segment(kk)%id
    
             gau_t(1:6) = (p1%segment(kk)%pend%coor + p1%segment(kk)%porigin%coor)/2. &
                + xg(1:6)* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            
            dx_t = p1%point_end%coor-p1%point_origin%coor
            xc_star = (p1%point_origin%coor+p1%point_end%coor)/2.
           
            st8=0.
            st9=0.
            do ig = 1,6
                st8=st8 +ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
                st9=st9 +ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)&
                    *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*wg(ig)
            enddo
            sum8 = sum8 + st8* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
            sum9 = sum9 + st9* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
        enddo
       
        element_t2(i,io,iio)%u1modal(nm)=sum8*ai(nm)/(p1%point_end%coor-p1%point_origin%coor)
        element_t2(i,io,iio)%u2modal(nm)=sum9*ai(nm)/(p1%point_end%coor-p1%point_origin%coor)
        
    enddo
    
    xcc=(p%point_end%coor + p%point_origin%coor)/2.
    dx_t = p%point_end%coor-p%point_origin%coor
    gau_t(1:6) = xcc + xg(1:6)* dx_t
    
    xcc1=(p1%point_end%coor + p1%point_origin%coor)/2.
    dx1_t = p1%point_end%coor-p1%point_origin%coor
    gau_t1(1:6) = xcc1 + xg(1:6)* dx1_t
   
    do nm = 1,nk+1
        st2=0.
        st3=0.
        st4=0.
        st5=0.
        do ig = 1,6
          lamda=wave_linear_j(gau_t(ig),speed1(i),speed1(i+1),p%point_origin%coor,p%point_end%coor,dx_t)
          st6=(RT11(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element_t1(i,io,iio)%u1modal(1:nk+1),gau_t(ig) ,xcc ,dx_t ,nk)+A12(gau_t(ig))*ortho_poly1d(element_t1(i,io,iio)%u2modal(1:nk+1),gau_t(ig) ,xcc ,dx_t ,nk))&
                +RT12(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element_t1(i,io,iio)%u1modal(1:nk+1),gau_t(ig) ,xcc ,dx_t ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element_t1(i,io,iio)%u2modal(1:nk+1),gau_t(ig) ,xcc ,dx_t ,nk)))&
                *diff_ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xcc ,dx_t,nk)*wg(ig)
          st2 = st2 +R11(gau_t(ig))*st6
          st3 = st3 +R21(gau_t(ig))*st6
          lamda1=wave_linear_j(gau_t1(ig),speed2(i),speed2(i+1),p1%point_origin%coor,p1%point_end%coor,dx1_t)
          st7=(RT21(gau_t1(ig))*((A11(gau_t1(ig))-lamda1)*ortho_poly1d(element_t2(i,io,iio)%u1modal(1:nk+1),gau_t1(ig) ,xcc1 ,dx1_t ,nk)+A12(gau_t1(ig))*ortho_poly1d(element_t2(i,io,iio)%u2modal(1:nk+1),gau_t1(ig) ,xcc1 ,dx1_t ,nk))&
                +RT22(gau_t1(ig))*((A21(gau_t1(ig)))*ortho_poly1d(element_t2(i,io,iio)%u1modal(1:nk+1),gau_t1(ig) ,xcc1 ,dx1_t ,nk)+(A22(gau_t1(ig))-lamda1)*ortho_poly1d(element_t2(i,io,iio)%u2modal(1:nk+1),gau_t1(ig) ,xcc1 ,dx1_t ,nk)))&
                *diff_ortho_poly1d( p1%aa(1:nk+1,nm),gau_t1(ig),xcc1 ,dx1_t,nk)*wg(ig)
          st4 = st4 +R12(gau_t1(ig))*st7
          st5 = st5 +R22(gau_t1(ig))*st7
        enddo
        element_t1(i,io,iio)%phi_int(nm,1)= st2* dx_t
        element_t1(i,io,iio)%phi_int(nm,2)= st3* dx_t
        element_t2(i,io,iio)%phi_int(nm,1)= st4* dx1_t
        element_t2(i,io,iio)%phi_int(nm,2)= st5* dx1_t
    enddo
    
    end subroutine get_umodal_dynamic
    ! subroutine get_integral_pk(p,p1,io,iio,xl,xr)
    !implicit none
    !
    !type(element1d_upstream), pointer :: p,p1
    !integer,intent(in) :: io,iio
    !real,intent(in) :: xl,xr
    !
    !integer :: kk
    !real :: sum,sum1,sum2,sum3,sum4,sum5,sum6,sum7
    !integer :: idx
    !real :: xi_end,xi_origin,term_1,term_2
    !
    !integer :: nm
    !real ::  bb(1:nk+1)
    !
    !real :: amatrix(1:nk+1,1:nk+1),store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1)
    !integer :: ig
    !
    !real :: gau_t(1:nk+1)
    !!real :: gau_t(1:6)
    !real :: gau_tt(1:3),st,st1,st2,st3,st4,st5,st6,st7
    !
    !real :: dx_t,dx,lamda
    !
    !integer :: ido,ide
    !real :: xs_o,xs_e,xc_star,xi(5)
    !integer :: igl
    !
    !
    !!**************************
    !ido = p%point_origin%id
    !ide = p%point_end%id
    !xs_o = p%point_origin%coor
    !xs_e = p%point_end%coor
    !!**************************
    !p%aa(1,1) = 1.;  p%aa(2,1) = 0.;  p%aa(3,1) = 0.;
    !p%aa(1,2) = 0.;  p%aa(2,2) = 1.;  p%aa(3,2) = 0.;
    !p%aa(1,3) = 0.;  p%aa(2,3) = 0.;  p%aa(3,3) = 1.;
    !
    !p1%aa(1,1) = 1.;  p1%aa(2,1) = 0.;  p1%aa(3,1) = 0.;
    !p1%aa(1,2) = 0.;  p1%aa(2,2) = 1.;  p1%aa(3,2) = 0.;
    !p1%aa(1,3) = 0.;  p1%aa(2,3) = 0.;  p1%aa(3,3) = 1.;
    !do nm = 1,nk+1
    !    
    !    !*********
    !
    !    sum = 0.
    !    sum1= 0.
    !    sum2 = 0.
    !    sum3= 0.
    !    sum4= 0.
    !    sum5= 0.
    !    sum6= 0.
    !    sum7= 0.
    !    dx= xr-xl
    !    do kk = 1 ,p%nsub
    !        idx = p%segment(kk)%id
    !
    !        gau_t(1:nk+1) = (p%segment(kk)%pend%coor + p%segment(kk)%porigin%coor)/2. &
    !            + gau(1:nk+1,1)* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        
    !        dx_t = p%point_end%coor-p%point_origin%coor
    !        xc_star = (p%point_origin%coor+p%point_end%coor)/2.
    !        st =0.
    !        st1=0.
    !        st2=0.
    !        st3=0.
    !        st4=0.
    !        st5=0.
    !        st6=0.
    !        st7=0.
    !        do ig = 1,nk+1
    !        
    !            lamda=wave_linear_j(gau_t(ig),speed1(i),speed1(i+1),p%point_origin%coor,p%point_end%coor,dx_t)
    !            st = st +R11(gau_t(ig))* (RT11(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st1 = st1 + R21(gau_t(ig))* (RT11(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st2 = st2 +R11(gau_t(ig))* (RT11(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                +RT12(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *diff_ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st3 = st3 +R21(gau_t(ig))* (RT11(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                +RT12(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *diff_ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st4 = st4 +(-0.5*aaderi(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st5 = st5 +((0.5*aaderi(gau_t(ig))/(aa1(gau_t(ig))**2.))*((A11(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A12(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st6= st6 +R11(gau_t(ig))* RT11(gau_t(ig))*f_source(gau_t(ig),element_t1(i,io,iio)%midtime)&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st7= st7 +R21(gau_t(ig))* RT11(gau_t(ig))*f_source(gau_t(ig),element_t1(i,io,iio)%midtime)&
    !                *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !        
    !        enddo
    !
    !        sum = sum + st* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum1 = sum1 + st1* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum2 = sum2 + st2* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum3 = sum3 + st3* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum4 = sum4 + st4* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum5 = sum5 + st5* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum6 = sum6 + st6* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !        sum7 = sum7 + st7* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)
    !    enddo
    !    !write(*,*)i
    !    element_t1(i,io,iio)%sum_int(nm,1)=sum
    !    element_t1(i,io,iio)%sum_int(nm,2)=sum1
    !    element_t1(i,io,iio)%phi_int(nm,1)=sum2
    !    element_t1(i,io,iio)%phi_int(nm,2)=sum3
    !    element_t1(i,io,iio)%vector_int(nm,1)=sum4
    !    element_t1(i,io,iio)%vector_int(nm,2)=sum5
    !    element_t1(i,io,iio)%source_int(nm,1)=sum6
    !    element_t1(i,io,iio)%source_int(nm,2)=sum7
    !    term_1=(RT11(p%point_end%coor)&
    !        *((A11(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)+A12(p%point_end%coor)*ortho_poly1d(element(p%point_end%id)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk))&
    !        +RT12(p%point_end%coor)*((A21(p%point_end%coor))*ortho_poly1d(element(p%point_end%id)%u1modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)+(A22(p%point_end%coor)-speed1(i+1))*ortho_poly1d(element(p%point_end%id)%u2modal(1:nk+1,io),p%point_end%coor ,x(p%point_end%id) ,dx ,nk)))&
    !        *ortho_poly1d( p%aa(1:nk+1,nm),p%point_end%coor,xc_star ,dx_t,nk)
    !    term_2=(RT11(p%point_origin%coor)&
    !        *((A11(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+A12(p%point_origin%coor)*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk))&
    !        +RT12(p%point_origin%coor)*((A21(p%point_origin%coor))*ortho_poly1d(element(p%point_origin%id)%u1modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)+(A22(p%point_origin%coor)-speed1(i))*ortho_poly1d(element(p%point_origin%id)%u2modal(1:nk+1,io),p%point_origin%coor ,x(p%point_origin%id) ,dx ,nk)))&
    !        *ortho_poly1d( p%aa(1:nk+1,nm),p%point_origin%coor,xc_star ,dx_t,nk)
    !    element_t1(i,io,iio)%flux1(nm,1)=R11(p%point_end%coor)*(term_1)
    !    element_t1(i,io,iio)%flux0(nm,1)=R11(p%point_origin%coor)*(term_2)
    !    element_t1(i,io,iio)%flux1(nm,2)=R21(p%point_end%coor)*(term_1)
    !    element_t1(i,io,iio)%flux0(nm,2)=R21(p%point_origin%coor)*(term_2)
    !    !write(*,*)i,io,iio,element_t1(i,io,iio)%u1modal(nm),element_t1(i,io,iio)%u2modal(nm)
    !    !pause
    !    !second charicteristic domain******************************************************************************************************************
    !    sum = 0.
    !    sum1= 0.
    !    sum2 = 0.
    !    sum3= 0.
    !    sum4= 0.
    !    sum5= 0.
    !    sum6= 0.
    !    sum7= 0.
    !    do kk = 1 ,p1%nsub
    !        idx = p1%segment(kk)%id
    !
    !        gau_t(1:nk+1) = (p1%segment(kk)%pend%coor + p1%segment(kk)%porigin%coor)/2. &
    !            + gau(1:nk+1,1)* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !         
    !        dx_t = p1%point_end%coor-p1%point_origin%coor
    !        xc_star = (p1%point_origin%coor+p1%point_end%coor)/2.
    !        st =0.
    !        st1=0.
    !        st2=0.
    !        st3=0.
    !        st4=0.
    !        st5=0.
    !        st6=0.
    !        st7=0.
    !        do ig = 1,nk+1
    !        
    !            lamda=wave_linear_j(gau_t(ig),speed2(i),speed2(i+1),p1%point_origin%coor,p1%point_end%coor,dx_t)
    !            st = st +R12(gau_t(ig))* (RT21(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT22(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st1 = st1 + R22(gau_t(ig))* (RT21(gau_t(ig))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+RT22(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st2 = st2 +R12(gau_t(ig))* (RT21(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                +RT22(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *diff_ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st3 = st3 +R22(gau_t(ig))* (RT21(gau_t(ig))*((A11(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+A12(gau_t(ig))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk))&
    !                +RT22(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig))-lamda)*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *diff_ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st4 = st4 +(0.5*aaderi(gau_t(ig))*((A21(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A22(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st5 = st5 +((-0.5*aaderi(gau_t(ig))/(aa1(gau_t(ig))**2.))*((A11(gau_t(ig)))*ortho_poly1d(element(idx)%u1modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)+(A12(gau_t(ig)))*ortho_poly1d(element(idx)%u2modal(1:nk+1,io),gau_t(ig) ,x(idx) ,dx ,nk)))&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st6= st6 +R12(gau_t(ig))* RT21(gau_t(ig))*f_source(gau_t(ig),element_t2(i,io,iio)%midtime)&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            st7= st7 +R22(gau_t(ig))* RT21(gau_t(ig))*f_source(gau_t(ig),element_t2(i,io,iio)%midtime)&
    !                *ortho_poly1d( p1%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
    !            
    !        enddo
    !
    !        sum = sum + st* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum1 = sum1 + st1* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum2 = sum2 + st2* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum3 = sum3 + st3* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum4 = sum4 + st4* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum5 = sum5 + st5* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum6 = sum6 + st6* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !        sum7 = sum7 + st7* (p1%segment(kk)%pend%coor - p1%segment(kk)%porigin%coor)
    !    enddo
    !    element_t2(i,io,iio)%sum_int(nm,1)=sum
    !    element_t2(i,io,iio)%sum_int(nm,2)=sum1
    !    element_t2(i,io,iio)%phi_int(nm,1)=sum2
    !    element_t2(i,io,iio)%phi_int(nm,2)=sum3
    !    element_t2(i,io,iio)%vector_int(nm,1)=sum4
    !    element_t2(i,io,iio)%vector_int(nm,2)=sum5
    !    element_t2(i,io,iio)%source_int(nm,1)=sum6
    !    element_t2(i,io,iio)%source_int(nm,2)=sum7
    !    term_1=(RT21(p1%point_end%coor)&
    !        *((A11(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)+A12(p1%point_end%coor)*ortho_poly1d(element(p1%point_end%id)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk))&
    !        +RT22(p1%point_end%coor)*((A21(p1%point_end%coor))*ortho_poly1d(element(p1%point_end%id)%u1modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)+(A22(p1%point_end%coor)-speed2(i+1))*ortho_poly1d(element(p1%point_end%id)%u2modal(1:nk+1,io),p1%point_end%coor ,x(p1%point_end%id) ,dx ,nk)))&
    !        *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_end%coor,xc_star ,dx_t,nk)
    !    term_2=(RT21(p1%point_origin%coor)&
    !        *((A11(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+A12(p1%point_origin%coor)*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk))&
    !        +RT22(p1%point_origin%coor)*((A21(p1%point_origin%coor))*ortho_poly1d(element(p1%point_origin%id)%u1modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)+(A22(p1%point_origin%coor)-speed2(i))*ortho_poly1d(element(p1%point_origin%id)%u2modal(1:nk+1,io),p1%point_origin%coor ,x(p1%point_origin%id) ,dx ,nk)))&
    !        *ortho_poly1d( p1%aa(1:nk+1,nm),p1%point_origin%coor,xc_star ,dx_t,nk)
    !    element_t2(i,io,iio)%flux1(nm,1)=R12(p1%point_end%coor)*(term_1)
    !    element_t2(i,io,iio)%flux0(nm,1)=R12(p1%point_origin%coor)*(term_2)
    !    element_t2(i,io,iio)%flux1(nm,2)=R22(p1%point_end%coor)*(term_1)
    !    element_t2(i,io,iio)%flux0(nm,2)=R22(p1%point_origin%coor)*(term_2)
    !enddo
    !
    !
    !end subroutine get_integral_pk
    