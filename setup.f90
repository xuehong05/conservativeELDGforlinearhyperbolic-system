    !*************************************************************************
    ! setup the coefficients of the scheme
    !   nx         ---> # of mesh
    !   nghost     ---> # of the ghost cells
    !   cfl        ---> CFL number
    !   time_final ---> the final time
    !   iexample   ---> for testing an example
    !                   e.g. iexample == 1: for u_t + u_x = 0.
    !                        iexample == 2: for u_t + (sin(x)u)_x = 0.
    !*************************************************************************
    subroutine setup
    implicit none

    nx = 10*2**kkkk
    !cfl=0.18
    
    !nghost = 11
    !cfl = 0.1*kkkk
    !nghost = int(3.*cfl)+2
    nghost = int(cfl)+10
    !time_final =  1.
    !time_final =  5!pi
    time_final =  2.85
    !time_final =  100.
    iexample = 1 !p2 cfl=0.15 for example BL
    nk = 1
    
    irk = 4
    irkdg = 0
    !irkdg = 1
    ! if we let speed =0, that is irkdg == 1
    ! and set the solution as initial condition at background elements,
    ! the code reduces to the normal RKDG method.
    ! setup computational domain [xleft,xright]!
    if(iexample == 1)then
        xleft = 0.
        xright = 2.*pi
        !xleft = 0.
        !xright = 40.*pi
    elseif(iexample ==2  )then
        !xleft = 0.
        !xright = 2.*pi
        xleft = -2.
        xright = 2.
    elseif(iexample ==3  )then
        !xleft = 0.
        !xright = 2.
        xleft = -10.
        xright = 10.
    elseif(iexample ==4  )then
        xleft = 0.
        xright = 2.*pi
    endif

    if(irk==1)then
        rk_d(0) = 0.
        rk_d(1) = 1.
    elseif(irk==2)then
        !!SSPRK 
        rk_d(0) = 0.
        rk_d(1) = 1.
        rk_d(2) = 1.
        !!RK
        !rk_d(0) = 0.
        !rk_d(1) = 1./2.
        !rk_d(2) = 1.
    elseif(irk==3)then
        !!SSPRK
        rk_d(0) = 0.
        rk_d(1) = 1.
        rk_d(2) = 0.5
        rk_d(3) = 1.
        ! !!RK
        !rk_d(0) = 0.
        !rk_d(1) = 1./3.
        !rk_d(2) = 2./3.
        !rk_d(3) = 1.
    elseif(irk==4)then
        rk_d(0) = 0.
        rk_d(1) = 0.5
        rk_d(2) = 0.5
        rk_d(3) = 1.
        rk_d(4) = 1.
    endif

    !*******************************************
    if(nk .eq. 1) then
        !  the points of 2th order Gauss-Lobatto quadrature
        gau_lob(1,1)=-0.5
        gau_lob(2,1)=0.5

        gau_lob(1,2)=0.5
        gau_lob(2,2)=0.5
    endif
    if(nk .eq. 2) then
        !  the points of 4th order Gauss-Lobatto quadrature
        gau_lob(1,1)=-0.5
        gau_lob(3,1)=0.5
        gau_lob(2,1)=0.0
        !   coefficients of 4h order Gauss-Lobatto quadrature
        gau_lob(1,2)=1.0/6.0
        gau_lob(3,2)=gau_lob(1,2)
        gau_lob(2,2)=2.0/3.0

    endif
    if( nk .eq. 3) then
        ! 6阶高斯积分点
        gau_lob(1,1)=-0.5
        gau_lob(2,1)=-sqrt(5.)/10.0
        gau_lob(3,1)= sqrt(5.0)/10.0
        gau_lob(4,1)=0.5
        !  6阶高斯积分系数
        gau_lob(1,2)=1.0/12.0
        gau_lob(2,2)=5.0/12.0
        gau_lob(3,2)=gau_lob(2,2)
        gau_lob(4,2)=gau_lob(1,2)
    endif
    if( nk .eq. 4) then
        !  8阶高斯积分点
        gau_lob(1,1)=-0.5

        gau_lob(2,1)=-sqrt(21.)/14.0
        gau_lob(3,1)=0.0
        gau_lob(4,1)= sqrt(21.0)/14.0

        gau_lob(5,1)=0.5
        !  8阶高斯积分系数
        gau_lob(1,2)=1.0/20.0

        gau_lob(2,2)=49.0/180.0
        gau_lob(3,2)=64.0/180.0
        gau_lob(4,2)=gau_lob(2,2)

        gau_lob(5,2)=gau_lob(1,2)
    endif
    if(nk.eq.5) then
        !  10阶高斯积分点
        gau_lob(1,1)=-0.5
        gau_lob(2,1)=-sqrt(147.+42.*sqrt(7.))/42.0
        gau_lob(3,1)=-sqrt(147.-42.*sqrt(7.))/42.0
        gau_lob(4,1)= sqrt(147.-42.*sqrt(7.))/42.0

        gau_lob(5,1)= sqrt(147.+42.*sqrt(7.))/42.0

        gau_lob(6,1)=0.5

        !  10阶高斯积分系数
        gau_lob(1,2)=1.0/30.0
        gau_lob(6,2)=gau_lob(1,2)
        gau_lob(2,2)=(-7.+5.*sqrt(7.))*sqrt(7.)*(7.+sqrt(7.))/840.0
        gau_lob(5,2)=gau_lob(2,2)
        gau_lob(3,2)=(7.+5.*sqrt(7.))*sqrt(7.)/(7.+sqrt(7.))/20.0
        gau_lob(4,2)=gau_lob(3,2)
    endif
    !******************************************************************

    !*******************************************
    if(nk .eq. 0)then
        gau(1,1) = 0.
        gau(1,2) = 1.
    endif
    if(nk .eq. 1) then
        !  the points of 4th order Gauss quadrature
        gau(1,1)=-sqrt(1./3.)*0.5
        gau(2,1)=sqrt(1./3.)*0.5

        gau(1,2)=0.5
        gau(2,2)=0.5
    endif
    if(nk .eq. 2) then
        !  the points of 6th order Gauss quadrature
        gau(1,1)=-sqrt(0.6)*0.5
        gau(2,1)=0.
        gau(3,1)=sqrt(0.6)*0.5
        !   coefficients of 6th order Gauss quadrature
        gau(1,2)=5./18.
        gau(2,2)=4./9.
        gau(3,2)=5./18.

    endif
    if( nk .eq. 3) then
        ! 8阶高斯积分点
        gau(1,1)=-sqrt( 3./7.+2./7.*sqrt(1.2) )*0.5
        gau(2,1)=-sqrt( 3./7.-2./7.*sqrt(1.2) )*0.5
        gau(3,1)= sqrt( 3./7.-2./7.*sqrt(1.2) )*0.5
        gau(4,1)=sqrt( 3./7.+2./7.*sqrt(1.2) )*0.5
        !  8阶高斯积分系数
        gau(1,2)=( 18.-sqrt(30.) )/36.*0.5
        gau(2,2)=( 18.+sqrt(30.) )/36.*0.5
        gau(3,2)=( 18.+sqrt(30.) )/36.*0.5
        gau(4,2)=( 18.-sqrt(30.) )/36.*0.5
    endif
    if( nk .eq. 4) then
        !  10阶高斯积分点
        gau(1,1)= -1./3.*sqrt( 5.+2.*sqrt(10./7.) )*0.5

        gau(2,1)= -1./3.*sqrt( 5.-2.*sqrt(10./7.) )*0.5
        gau(3,1)= 0.
        gau(4,1)= 1./3.*sqrt( 5.-2.*sqrt(10./7.) )*0.5

        gau(5,1)= 1./3.*sqrt( 5.+2.*sqrt(10./7.) )*0.5
        ! 10阶高斯积分系数
        gau(1,2)= ( 322.-13.*sqrt(70.) )/900.*0.5

        gau(2,2)= ( 322.+13.*sqrt(70.) )/900.*0.5
        gau(3,2)= 128./225.*0.5
        gau(4,2)= ( 322.+13.*sqrt(70.) )/900.*0.5

        gau(5,2)= ( 322.-13.*sqrt(70.) )/900.*0.5
    endif

    !************************************************************

    if(nk==1 .or. nk==0)then
        gau_f_phi_x(1,1) = -0.5
        gau_f_phi_x(2,1) = 0.0
        gau_f_phi_x(3,1) = 0.5
        !
        gau_f_phi_x(1,2) = 1.0/6.0
        gau_f_phi_x(2,2) = 2./3.
        gau_f_phi_x(3,2) = 1.0/6.0
    elseif(nk==2)then
        gau_f_phi_x(1,1) = -sqrt(0.6)*0.5
        gau_f_phi_x(2,1) = 0.
        gau_f_phi_x(3,1) = sqrt(0.6)*0.5
        !
        gau_f_phi_x(1,2) = 5./18.
        gau_f_phi_x(2,2) = 4./9.
        gau_f_phi_x(3,2) = 5./18.
    endif
    !************************************************************

    end subroutine setup