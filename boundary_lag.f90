    subroutine boundary_lag(io,iio)
    implicit none
    
    integer, intent(in) :: io,iio
    !write(*,*)1-1,io,iio,nk,nx,io,iio
    !pause
    element_t1(1-1,io,iio)%u1modal(1:nk+1) = element_t1(nx,io,iio)%u1modal(1:nk+1)
    element_t1(1-1,io,iio)%u2modal(1:nk+1) = element_t1(nx,io,iio)%u2modal(1:nk+1)
    !element_t1(1-1,io,iio)%v1modal(1:nk+1) = element_t1(nx,io,iio)%v1modal(1:nk+1)
    !element_t1(1-1,io,iio)%v2modal(1:nk+1) = element_t1(nx,io,iio)%v2modal(1:nk+1)
    element_t1(1-1,io,iio)%xr = element_t1(1,io,iio)%xl
    element_t1(1-1,io,iio)%xl = element_t1(1,io,iio)%xl-( element_t1(nx,io,iio)%xr-element_t1(nx,io,iio)%xl )
    element_t1(1-1,io,iio)%xc = ( element_t1(1-1,io,iio)%xl+element_t1(1-1,io,iio)%xr )/2.
    element_t1(1-1,io,iio)%dx = element_t1(1-1,io,iio)%xr-element_t1(1-1,io,iio)%xl
    !element_t1(1-1,io,iio)%weno_t(1)=element_t1(1,io,iio)%weno_t(1)-(element(nx)%xr-element(nx)%xl)
    !element_t1(1-1,io,iio)%weno_t(1)=element_t1(nx,io,iio)%weno_t(1)-(xright-xleft)
    !element_t1(1-1,io,iio)%weno_t(2)=element_t1(1,io,iio)%weno_t(2)-(element(nx)%xr-element(nx)%xl)
    !*******************
    element_t1(nx+1,io,iio)%u1modal(1:nk+1) = element_t1(1,io,iio)%u1modal(1:nk+1)
    element_t1(nx+1,io,iio)%u2modal(1:nk+1) = element_t1(1,io,iio)%u2modal(1:nk+1)
    !element_t1(nx+1,io,iio)%v1modal(1:nk+1) = element_t1(1,io,iio)%v1modal(1:nk+1)
    !element_t1(nx+1,io,iio)%v2modal(1:nk+1) = element_t1(1,io,iio)%v2modal(1:nk+1)
    element_t1(nx+1,io,iio)%xl = element_t1(nx,io,iio)%xr
    element_t1(nx+1,io,iio)%xr = element_t1(nx,io,iio)%xr+ ( element_t1(1,io,iio)%xr-element_t1(1,io,iio)%xl )
    element_t1(nx+1,io,iio)%xc = ( element_t1(nx+1,io,iio)%xl+element_t1(nx+1,io,iio)%xr )/2.
    element_t1(nx+1,io,iio)%dx = element_t1(nx+1,io,iio)%xr-element_t1(nx+1,io,iio)%xl
    !element_t1(nx+1,io,iio)%weno_t(1)=element_t1(nx,io,iio)%weno_t(1)+(element(nx)%xr-element(nx)%xl)
    !element_t1(nx+1,io,iio)%weno_t(1)=element_t1(1,io,iio)%weno_t(1)+(xright-xleft)
    !element_t1(nx+1,io,iio)%weno_t(2)=element_t1(nx,io,iio)%weno_t(2)+(element(nx)%xr-element(nx)%xl)
    !******
    
    element_t2(1-1,io,iio)%u1modal(1:nk+1) = element_t2(nx,io,iio)%u1modal(1:nk+1)
    element_t2(1-1,io,iio)%u2modal(1:nk+1) = element_t2(nx,io,iio)%u2modal(1:nk+1)
    !element_t2(1-1,io,iio)%v1modal(1:nk+1) = element_t2(nx,io,iio)%v1modal(1:nk+1)
    !element_t2(1-1,io,iio)%v2modal(1:nk+1) = element_t2(nx,io,iio)%v2modal(1:nk+1)
    element_t2(1-1,io,iio)%xr = element_t2(1,io,iio)%xl
    element_t2(1-1,io,iio)%xl = element_t2(1,io,iio)%xl-( element_t2(nx,io,iio)%xr-element_t2(nx,io,iio)%xl )
    element_t2(1-1,io,iio)%xc = ( element_t2(1-1,io,iio)%xl+element_t2(1-1,io,iio)%xr )/2.
    element_t2(1-1,io,iio)%dx = element_t2(1-1,io,iio)%xr-element_t2(1-1,io,iio)%xl
    !element_t2(1-1,io,iio)%weno_t(1)=element_t2(1,io,iio)%weno_t(1)-(element(nx)%xr-element(nx)%xl)
    !element_t2(1-1,io,iio)%weno_t(1)=element_t2(nx,io,iio)%weno_t(1)-(xright-xleft)
    !element_t2(1-1,io,iio)%weno_t(2)=element_t2(1,io,iio)%weno_t(2)-(element(nx)%xr-element(nx)%xl)
    !*******************
    element_t2(nx+1,io,iio)%u1modal(1:nk+1) = element_t2(1,io,iio)%u1modal(1:nk+1)
    element_t2(nx+1,io,iio)%u2modal(1:nk+1) = element_t2(1,io,iio)%u2modal(1:nk+1)
    !element_t2(nx+1,io,iio)%v1modal(1:nk+1) = element_t2(1,io,iio)%v1modal(1:nk+1)
    !element_t2(nx+1,io,iio)%v2modal(1:nk+1) = element_t2(1,io,iio)%v2modal(1:nk+1)
    element_t2(nx+1,io,iio)%xl = element_t2(nx,io,iio)%xr
    element_t2(nx+1,io,iio)%xr = element_t2(nx,io,iio)%xr+ ( element_t2(1,io,iio)%xr-element_t2(1,io,iio)%xl )
    element_t2(nx+1,io,iio)%xc = ( element_t2(nx+1,io,iio)%xl+element_t2(nx+1,io,iio)%xr )/2.
    element_t2(nx+1,io,iio)%dx = element_t2(nx+1,io,iio)%xr-element_t2(nx+1,io,iio)%xl
    !element_t2(nx+1,io,iio)%weno_t(1)=element_t2(nx,io,iio)%weno_t(1)+(element(nx)%xr-element(nx)%xl)
    !element_t2(nx+1,io,iio)%weno_t(1)=element_t2(1,io,iio)%weno_t(1)+(xright-xleft)
    !element_t2(nx+1,io,iio)%weno_t(2)=element_t2(nx,io,iio)%weno_t(2)+(element(nx)%xr-element(nx)%xl)
    !if(iexample==3)then
    !    element_t(1-1,io)%umodal(1:1) = 1.
    !    element_t(1-1,io)%umodal(2:nk+1) = 0.
    !    element_t(nx+1,io)%umodal(1:nk+1) = 0.
    !endif
    
    end subroutine boundary_lag