    !************************************************************************* 
    subroutine order_dg
    implicit none
    integer :: kk0,ig,igr
    real :: exact_ave
    real :: xrg
    real :: error11,error12,error13,error21,error22,error23,error31,error32,error33,error41,error42,error43,masspost1,masspost2
    real :: rr11,rr12,rr13,rr21,rr22,rr23,rr31,rr32,rr33,rr41,rr42,rr43
    real :: exe,exe2,exep,exep2
   

    do kk0 = 1 , nx
        !u1mean=u1mean+element(kk0)%u1modal(1)*dx
        !u2mean=u2mean+element(kk0)%u2modal(1)*dx
        write(1,*) x(kk0),element(kk0)%u1modal(1,0),element(kk0)%u2modal(1,0)
    enddo
    close(1)

    error11=0.0
    error12=0.0
    error13=0.0
    error21=0.0
    error22=0.0
    error23=0.0
    error31=0.0
    error32=0.0
    error33=0.0
    error41=0.0
    error42=0.0
    error43=0.0
    masspost1=0.0
    masspost2=0.0
    do kk0=1,nx
        !exact_ave = 0.
        !do ig = 1 , 6
        !    xrg = x(kk0) + dx * xg(ig)
        !    exact_ave = exact_ave + burgersriemann( xrg,time_final ) * wg(ig)
        !enddo

        do ig = 1 , 6
            igr=6*(kk0-1)+ig
            xrg =x(kk0) + dx * xg(ig)
            !write(2,*) xrg, exact(xrg,time_final), ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)
            !exe =exact(xrg,time_final)-ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)
            !write(2,*) xrg, exact1( xrg,time_final ),ortho_poly1d(element(kk0)%u1modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            write(2,*) xrg, exact2( xrg,time_final ),ortho_poly1d(element(kk0)%u2modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exe =exact1( xrg,time_final )-ortho_poly1d(element(kk0)%u1modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exep =exact1( xrg,time_final )-u1post(igr,2)
            masspost1=masspost1+u1post(igr,2)*wg(ig)
            masspost2=masspost2+u2post(igr,2)*wg(ig)
            exe2 =exact2( xrg,time_final )-ortho_poly1d(element(kk0)%u2modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exep2 =exact2( xrg,time_final )-u2post(igr,2)
            error11=error11+abs(exe) *wg(ig)
            error12=error12+(exe)**2 *wg(ig)
            error13=max(error13,abs(exe) )
            error21=error21+abs(exe2) *wg(ig)
            error22=error22+(exe2)**2 *wg(ig)
            error23=max(error23,abs(exe2) )
            error31=error31+abs(exep) *wg(ig)
            error32=error32+(exep)**2 *wg(ig)
            error33=max(error33,abs(exep) )
            error41=error41+abs(exep2) *wg(ig)
            error42=error42+(exep2)**2 *wg(ig)
            error43=max(error43,abs(exep2) )
        enddo
    enddo
    u1meanpost=masspost1*dx
    u2meanpost=masspost2*dx
    error11=error11/nx
    error12=sqrt(error12/nx)
    error21=error21/nx
    error22=sqrt(error22/nx)
    error31=error31/nx
    error32=sqrt(error32/nx)
    error41=error41/nx
    error42=sqrt(error42/nx)
    if(kkkk.eq.1) write(101,103) nx,error11,error12 ,error13,error31,error32 ,error33
    !write(*,*) 'error',error12,error11,error13,error32,error31,error33
    if(kkkk.gt.1) then
        rr11=log(er11/error11)/log(2.)
        rr12=log(er12/error12)/log(2.)
        rr13=log(er13/error13)/log(2.)
        !rr21=log(er21/error21)/log(2.)
        !rr22=log(er22/error22)/log(2.)
        !rr23=log(er23/error23)/log(2.)
        rr31=log(er31/error31)/log(2.)
        rr32=log(er32/error32)/log(2.)
        rr33=log(er33/error33)/log(2.)
        write(101,102) nx,error12,rr12,error11, rr11 ,error13,rr13
        !write(*,*) nx,rr12,rr11,rr13,rr22,rr21,rr23
        !write(*,*) nx,rr12,rr11,rr13,rr32,rr31,rr33
        !pause
    endif
    er11=error11
    er12=error12
    er13=error13
    !er21=error21
    !er22=error22
    !er23=error23
    er31=error31
    er32=error32
    er33=error33

102 format(i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'\\',1x,'\hline')
103 format(i6,1x,3('&',1x,es12.2E2,1x,'&',1x),'\\',1x,'\hline')
123 format( (1x,f16.6),3(1x,es12.5E2)  )

    !write(124,123)  cfl,error11,error12,error13
    !write(124,123)  cfl,error11,error31,error13
    !write(124,123)  nx,error11,rr11,error12,rr12,error13,rr13
    write(124,123)  nx,error31,rr31,error32,rr32,error33,rr33
    !write(3,*) cfl, error13, error23
    !%post-processing
    write(3,*) cfl, error33, error43
    end subroutine order_dg
    
   !************************************************************   
    subroutine Lerror_dg
    implicit none
    integer :: kk0,ig
    real :: exact_ave
    real :: xrg
    real :: error11,error12,error13,error21,error22,error23
    real :: rr11,rr12,rr13,rr21,rr22,rr23
    real :: exe,exe2
   

    error12=0.0
    do kk0=1,nx
        do ig = 1 , 6
            xrg =x(kk0) + dx * xg(ig)
            
            !write(2,*) xrg, exact2( xrg,time_final ),ortho_poly1d(element(kk0)%u2modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exe =exact1( xrg,time )-ortho_poly1d(element(kk0)%u1modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            !exe2 =exact2( xrg,time )-ortho_poly1d(element(kk0)%u2modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            !error11=error11+abs(exe) *wg(ig)
            error12=error12+(exe)**2 *wg(ig)
            !error13=max(error13,abs(exe) )
            !error21=error21+abs(exe2) *wg(ig)
            !error22=error22+(exe2)**2 *wg(ig)
            !error23=max(error23,abs(exe2) )
        enddo
    enddo
    !error11=error11/nx
    error12=sqrt(error12/nx)
    !error21=error21/nx
    !error22=sqrt(error22/nx)
    !if(kkkk.eq.1) write(101,103) nx,error11,error12 ,error13
   
    write(43,*) time, error12
    end subroutine Lerror_dg    
    
    subroutine error_dg(time)
    implicit none
    
    real,intent(in) :: time
    
    integer :: kk0,ig
    real :: exact_ave
    real :: xrg
    real :: error13,error23
    !real :: rr11,rr12,rr13,rr21,rr22,rr23
    real :: exe,exe2

    error13=0.0
    error23=0.0
    do kk0=1,nx
        !exact_ave = 0.
        !do ig = 1 , 6
        !    xrg = x(kk0) + dx * xg(ig)
        !    exact_ave = exact_ave + burgersriemann( xrg,time_final ) * wg(ig)
        !enddo

        do ig = 1 , 6
            xrg =x(kk0) + dx * xg(ig)
            !write(2,*) xrg, exact(xrg,time_final), ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)
            !exe =exact(xrg,time_final)-ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)
            !write(2,*) xrg, exact1( xrg,time_final ), ortho_poly1d(element(kk0)%u1modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exe =exact1( xrg,time )-ortho_poly1d(element(kk0)%u1modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            exe2 =exact2( xrg,time )-ortho_poly1d(element(kk0)%u2modal(1:nk+1,0),xrg,x(kk0),dx,nk)
            error13=max(error13,abs(exe) )
            error23=max(error23,abs(exe2) )
        enddo
    enddo
    !error11=error11/nx
    !error12=sqrt(error12/nx)
    !error21=error21/nx
    !error22=sqrt(error22/nx)
!    if(kkkk.eq.1) write(101,103) nx,error11,error12 ,error13
!    !write(*,*) 'error',error12,error11,error13,error22,error21,error23
!    if(kkkk.gt.1) then
!        rr11=log(er11/error11)/log(2.)
!        rr12=log(er12/error12)/log(2.)
!        rr13=log(er13/error13)/log(2.)
!        rr21=log(er21/error21)/log(2.)
!        rr22=log(er22/error22)/log(2.)
!        rr23=log(er23/error23)/log(2.)
!        write(101,102) nx,error12,rr12,error11, rr11 ,error13,rr13
!        !write(*,*) nx,rr12,rr11,rr13,rr22,rr21,rr23
!        !pause
!    endif
!    er11=error11
!    er12=error12
!    er13=error13
!    er21=error21
!    er22=error22
!    er23=error23
!
!102 format(i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'\\',1x,'\hline')
!103 format(i6,1x,3('&',1x,es12.2E2,1x,'&',1x),'\\',1x,'\hline')
!123 format( (1x,f16.6),3(1x,es12.5E2)  )

    !write(124,123)  cfl,error11,error12,error13
    write(*,*) time,error13,error23
    !write(31,*) time,error13,error23
    end subroutine error_dg