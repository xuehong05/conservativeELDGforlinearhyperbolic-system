    subroutine order
    implicit none
    integer :: kk0,ig
    real :: exact_ave
    real :: xrg
    real :: error1,error2,error3
    real :: rr1,rr2,rr3

    do kk0 = 1 , nx
        write(1,*) x(kk0),element(kk0)%umodal(1)
    enddo
    close(1)

    error1=0.0
    error2=0.0
    error3=0.0
    do kk0=1,nx
        exact_ave = 0.
        do ig = 1 , 6
            xrg = x(kk0) + dx * xg(ig)

            exact_ave = exact_ave + exact( xrg,time_final ) * wg(ig)

        enddo

        error1=error1+abs( exact_ave-element(kk0)%umodal(1)   )
        error3=max(error3,abs(exact_ave-element(kk0)%umodal(1)  ))
        error2=error2+( exact_ave-element(kk0)%umodal(1)  )**2
    enddo
    error1=error1/nx
    error2=sqrt(error2/nx)
    if(kkkk.eq.1) write(101,103) nx,error1,error2 ,error3
    write(*,*) 'error',error1,error2
    if(kkkk.gt.1) then
        rr1=log(er1/error1)/log(2.)
        rr2=log(er2/error2)/log(2.)
        rr3=log(er3/error3)/log(2.)
        write(101,102) nx,error1,rr1,error2, rr2 ,error3,rr3
        write(*,*) nx,rr1,rr2,rr3
    endif
    er1=error1
    er2=error2
    er3=error3

102 format(i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'\\',1x,'\hline')
103 format(i6,1x,3('&',1x,es12.2E2,1x,'&',1x),'\\',1x,'\hline')
123 format(4(1x,f16.6))


    end subroutine order