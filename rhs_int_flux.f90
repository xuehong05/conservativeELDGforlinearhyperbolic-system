!    subroutine rhs_int_flux
!    implicit none
!    type(element1d_dynamic), pointer :: pet,pet1
!    real :: xll,xrr
!    real :: u11l,u11r,u12l,u12r,u21l,u21r,u22l,u22r
!    real :: alpha,RL_j1(1:2,1:2),RL_j1p(1:2,1:2),RL_j2(1:2,1:2),RL_j2p(1:2,1:2)
!
!    do i = 1,nx
!        ! consider a dynamic element
!        !io=iio-1
!        pet => element_t1(i,io,iio)
!        pet1 => element_t2(i,io,iio)
!        xll = element(i)%xl
!        xrr = element(i)%xr
!!compute (F,phi_x)
!        call get_integral_f_phi_x(pet,pet1,element(i)%RT_j)
!        !call get_integral_f_phi_x(pet1,xll,xrr)
!
!    enddo
!
!    ! boundary condition of the dynamic elements
!    !element_t(1-1,io)%umodal(1:nk+1) = element_t(nx,io)%umodal(1:nk+1)
!    !element_t(nx+1,io)%umodal(1:nk+1) = element_t(1,io)%umodal(1:nk+1)
!    call boundary_lag
!    !element_t(1-1,io)%xr = element_t(1,io)%xl
!    !element_t(1-1,io)%xl = element_t(1,io)%xl-( element_t(nx,io)%xr-element_t(nx,io)%xl )
!    !element_t(1-1,io)%xc = ( element_t(1-1,io)%xl+element_t(1-1,io)%xr )/2.
!    !element_t(1-1,io)%dx = element_t(1-1,io)%xr-element_t(1-1,io)%xl
!    !element_t(1-1,io)%weno_t(1)=element_t(1,io)%weno_t(1)-(element(nx)%xr-element(nx)%xl)
!    !!*******************
!    !element_t(nx+1,io)%xl = element_t(nx,io)%xr
!    !element_t(nx+1,io)%xr = element_t(nx,io)%xr+ ( element_t(1,io)%xr-element_t(1,io)%xl )
!    !element_t(nx+1,io)%xc = ( element_t(nx+1,io)%xl+element_t(nx+1,io)%xr )/2.
!    !element_t(nx+1,io)%dx = element_t(nx+1,io)%xr-element_t(nx+1,io)%xl
!    !element_t(nx+1,io)%weno_t(1)=element_t(nx,io)%weno_t(1)+(element(nx)%xr-element(nx)%xl)
!
!    !***************************************************************************
!    ! compute alpha for global Lax-Friedrichs flux
!    !alpha = 0.
!    !do i = 1,nx+1
!    !    pet=>element_t(i-1,io)
!    !    ul = ortho_poly1d(pet%umodal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
!    !
!    !    pet=>element_t(i,io)
!    !    ur = ortho_poly1d(pet%umodal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
!    !    alpha1 = 1.0*max( abs(fp( ul,element_t(i-1,io)%xr  )-speed(i)), &
!    !        abs(fp( ur,element_t(i,io)%xl )-speed(i) ) )
!    !    if(iexample==4)then
!    !        if( (ul-0.5)*(ur-0.5)<0. )then
!    !            alpha1 = max( abs(  fp(0.5,element_t(i,io)%xl) -speed(i) ) ,&
!    !                abs(fp( ul,element_t(i-1,io)%xr  )-speed(i)), &
!    !        abs(fp( ur,element_t(i,io)%xl )-speed(i) )   )
!    !        endif
!    !    endif 
!    !    
!    !    alpha = max( alpha1,alpha) 
!    !enddo
!    !
!    !***************************************************************************
! 
!
!    do i =1,nx+1
!
!        pet=>element_t1(i-1,io,iio)
!        u11l = ortho_poly1d(pet%u1modal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
!        u12l = ortho_poly1d(pet%u2modal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
!        pet1=>element_t2(i-1,io,iio)
!        u21l = ortho_poly1d(pet1%u1modal(1:nk+1),pet1%xr,pet1%xc,pet1%dx,nk)
!        u22l = ortho_poly1d(pet1%u2modal(1:nk+1),pet1%xr,pet1%xc,pet1%dx,nk)
!
!        pet=>element_t1(i,io,iio)
!        u11r = ortho_poly1d(pet%u1modal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
!        u12r = ortho_poly1d(pet%u2modal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
!        pet1=>element_t2(i,io,iio)
!        u21r = ortho_poly1d(pet1%u1modal(1:nk+1),pet1%xl,pet1%xc,pet1%dx,nk)
!        u22r = ortho_poly1d(pet1%u2modal(1:nk+1),pet1%xl,pet1%xc,pet1%dx,nk)
!        
!        !local LF
!        alpha = 1.1*max( abs(fp1(element_t1(i-1,io,iio)%xr)-speed1(i)), &
!            abs(fp2(element_t2(i,io,iio)%xl)-speed2(i) ) )
!        
!        !alpha = 1. *max( abs(fp( ul,element_t(i-1,io)%xr  )-speed(i)), &
!        !    abs(fp( ur,element_t(i,io)%xl )-speed(i) ) )
!        !if(iexample==4)then
!        !    if( (ul-0.5)*(ur-0.5)<0. )then
!        !        alpha = 1. * max( abs( 1.1*fp(0.5,element_t(i,io)%xl) -speed(i) ) ,&
!        !            abs(fp( ul,element_t(i-1,io)%xr  )-speed(i)), &
!        !    abs(fp( ur,element_t(i,io)%xl )-speed(i) )   )
!        !    endif
!        !endif
!        !
!        !if (i==1) then
!        !    RL_j1(1,1)=element(i)%R_j(1,1)*element(i)%RT_j(1,1)
!        !    RL_j1(1,2)=element(i)%R_j(1,1)*element(i)%RT_j(1,2)
!        !    RL_j1(2,1)=element(i)%R_j(2,1)*element(i)%RT_j(1,1)
!        !    RL_j1(2,2)=element(i)%R_j(2,1)*element(i)%RT_j(1,2)
!        !    RL_j1p(1,1)=element(nx)%R_j(1,1)*element(nx)%RT_j(1,1)
!        !    RL_j1p(1,2)=element(nx)%R_j(1,1)*element(nx)%RT_j(1,2)
!        !    RL_j1p(2,1)=element(nx)%R_j(2,1)*element(nx)%RT_j(1,1)
!        !    RL_j1p(2,2)=element(nx)%R_j(2,1)*element(nx)%RT_j(1,2)
!        !    RL_j2(1,1)=element(i)%R_j(1,2)*element(i)%RT_j(2,1)
!        !    RL_j2(1,2)=element(i)%R_j(1,2)*element(i)%RT_j(2,2)
!        !    RL_j2(2,1)=element(i)%R_j(2,2)*element(i)%RT_j(2,1)
!        !    RL_j2(2,2)=element(i)%R_j(2,2)*element(i)%RT_j(2,2)
!        !    RL_j2p(1,1)=element(nx)%R_j(1,2)*element(nx)%RT_j(2,1)
!        !    RL_j2p(1,2)=element(nx)%R_j(1,2)*element(nx)%RT_j(2,2)
!        !    RL_j2p(2,1)=element(nx)%R_j(2,2)*element(nx)%RT_j(2,1)
!        !    RL_j2p(2,2)=element(nx)%R_j(2,2)*element(nx)%RT_j(2,2)
!        !elseif(i==nx+1) then
!        !    RL_j1(1,1)=element(1)%R_j(1,1)*element(1)%RT_j(1,1)
!        !    RL_j1(1,2)=element(1)%R_j(1,1)*element(1)%RT_j(1,2)
!        !    RL_j1(2,1)=element(1)%R_j(2,1)*element(1)%RT_j(1,1)
!        !    RL_j1(2,2)=element(1)%R_j(2,1)*element(1)%RT_j(1,2)
!        !    RL_j1p(1,1)=element(i-1)%R_j(1,1)*element(i-1)%RT_j(1,1)
!        !    RL_j1p(1,2)=element(i-1)%R_j(1,1)*element(i-1)%RT_j(1,2)
!        !    RL_j1p(2,1)=element(i-1)%R_j(2,1)*element(i-1)%RT_j(1,1)
!        !    RL_j1p(2,2)=element(i-1)%R_j(2,1)*element(i-1)%RT_j(1,2)
!        !    RL_j2(1,1)=element(1)%R_j(1,2)*element(1)%RT_j(2,1)
!        !    RL_j2(1,2)=element(1)%R_j(1,2)*element(1)%RT_j(2,2)
!        !    RL_j2(2,1)=element(1)%R_j(2,2)*element(1)%RT_j(2,1)
!        !    RL_j2(2,2)=element(1)%R_j(2,2)*element(1)%RT_j(2,2)
!        !    RL_j2p(1,1)=element(i-1)%R_j(1,2)*element(i-1)%RT_j(2,1)
!        !    RL_j2p(1,2)=element(i-1)%R_j(1,2)*element(i-1)%RT_j(2,2)
!        !    RL_j2p(2,1)=element(i-1)%R_j(2,2)*element(i-1)%RT_j(2,1)
!        !    RL_j2p(2,2)=element(i-1)%R_j(2,2)*element(i-1)%RT_j(2,2)
!        !else
!            !RL_j1(1,1)=element(i)%R_j(1,1)*element(i)%RT_j(1,1)
!            !RL_j1(1,2)=element(i)%R_j(1,1)*element(i)%RT_j(1,2)
!            !RL_j1(2,1)=element(i)%R_j(2,1)*element(i)%RT_j(1,1)
!            !RL_j1(2,2)=element(i)%R_j(2,1)*element(i)%RT_j(1,2)
!            !RL_j1p(1,1)=element(i-1)%R_j(1,1)*element(i-1)%RT_j(1,1)
!            !RL_j1p(1,2)=element(i-1)%R_j(1,1)*element(i-1)%RT_j(1,2)
!            !RL_j1p(2,1)=element(i-1)%R_j(2,1)*element(i-1)%RT_j(1,1)
!            !RL_j1p(2,2)=element(i-1)%R_j(2,1)*element(i-1)%RT_j(1,2)
!            !RL_j2(1,1)=element(i)%R_j(1,2)*element(i)%RT_j(2,1)
!            !RL_j2(1,2)=element(i)%R_j(1,2)*element(i)%RT_j(2,2)
!            !RL_j2(2,1)=element(i)%R_j(2,2)*element(i)%RT_j(2,1)
!            !RL_j2(2,2)=element(i)%R_j(2,2)*element(i)%RT_j(2,2)
!            !RL_j2p(1,1)=element(i-1)%R_j(1,2)*element(i-1)%RT_j(2,1)
!            !RL_j2p(1,2)=element(i-1)%R_j(1,2)*element(i-1)%RT_j(2,2)
!            !RL_j2p(2,1)=element(i-1)%R_j(2,2)*element(i-1)%RT_j(2,1)
!            !RL_j2p(2,2)=element(i-1)%R_j(2,2)*element(i-1)%RT_j(2,2)
!        !endif
!        
!        call compute_num_flux_LF( u11l,u12l,u11r,u12r,element_t1(i,io,iio)%xl,speed1(i),alpha,flux_LF1(:,i) )
!        call compute_num_flux_LF( u21l,u22l,u21r,u22r,element_t2(i,io,iio)%xl,speed2(i),alpha,flux_LF2(:,i) )
!    enddo
!
!    end subroutine rhs_int_flux