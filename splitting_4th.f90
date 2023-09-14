    subroutine splitting_4th
    implicit none
    !real :: umod( ng, 1 - nghost : nx_2d + nghost  )
    !real :: u( 1:ng, 1 - nghost : nx_2d + nghost  )
    real :: xrg
    real :: ftemp1,ftemp2,tnum
    integer :: ii
    
    tnum = time1 -dt1+cof_4thsplit(1)*dt1
    do je = 1, ny_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do ie = 1 ,nx_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi2( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !write(*,*)ftemp
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_x
    
            !*********************
    
            !do ii = 1-1,nel_x+1
            !    u( 1:nod+1,ii ) = elem( ii,je )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel_x,u(1:nod+1,1-1:nel_x+1) )
            !
            !call weno_limiter(nod,nel_x,u(1:nod+1,1-ighost:nel_x+ighost),ighost)
            !do ii = 1 , nel_x
            !    elem( ii,je )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_x
            !*********************
            do ii = 1,nx_2d
                element(ii)%u1modal(1:ng,0) = elem( ii,je )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ii,je )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 1,tnum, cof_4thsplit(1)*dt1, gy( 1  ,kk1,1,je) )
            !do ie = 1 ,nel_x
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 1,dx_2d,tnum, dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) , gy( 1  ,kk,ie,je),umod( 1:n_moment,ie ) )
            !
            !enddo !ie
    
            do ie = 1 , nx_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,ie)
                elem(ie,je)%split1_modal( 1:ng )=element(ie)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(ie)%u2modal(1:ng,0)
                !elem(ie,je)%split_modal( 1:ng )=element(ie)%umodal(1:ng)
                !************************compute nodal
                do ig = 1,ng
                    xrg = x_2d(ie) + dx_2d * x_g(ig)
                    elem(ie,je)%psi1(ig,kk1) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                    elem(ie,je)%psi2(ig,kk1) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! je
    
    
    tnum = time1 -dt1+cof_4thsplit(2)*dt1
    
    
    do ie = 1, nx_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do je = 1 ,ny_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi3( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !ftemp2 = ftemp2 + elem(ie,je)%psi2( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_y
    
            !*********************
            !do ii = 1-1,nel_y+1
            !    u( 1:nod+1,ii ) = elem( ie,ii )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_y
            !    elem( ie,ii )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_y
            !*********************
            do ii = 1,ny_2d
                element(ii)%u1modal(1:ng,0) = elem( ie,ii )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ie,ii )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 2,tnum, cof_4thsplit(2)*dt1, gx( kk1 ,1,ie,1) )
    
            !do je = 1 ,nel_y
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 2,dy,tnum,dt, gly( 1 ,1:n_gl,ie,je),gy( kk,1:n_g ,ie,je),gx( kk,1  ,ie,je) ,umod( 1:n_moment,je ) )
            !
            !enddo !je
    
            do je = 1 , ny_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,je)
                elem(ie,je)%split1_modal( 1:ng )=element(je)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(je)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = y_2d(je) + dy_2d * x_g(ig)
                    elem(ie,je)%psi1(kk1,ig) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    elem(ie,je)%psi3(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    !elem(ie,je)%psi2(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! ie
    
    !call output
    !pause
    !print *,1111111111
    !***************************************************************************************************
    !***************************************************************************************************
    tnum = time1 -dt1+(cof_4thsplit(1)+cof_4thsplit(3))*dt1
    do je = 1, ny_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do ie = 1 ,nx_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi2( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_x
    
            !!*********************
            !do ii = 1-1,nel_x+1
            !    u( 1:nod+1,ii ) = elem( ii,je )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_x
            !    elem( ii,je )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_x
            !!*********************
            do ii = 1,nx_2d
                element(ii)%u1modal(1:ng,0) = elem( ii,je )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ii,je )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 1,tnum, cof_4thsplit(3)*dt1, gy( 1  ,kk1,1,je) )
            !do ie = 1 ,nel_x
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 1,dx_2d,tnum,dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) ,gy( 1  ,kk,ie,je) ,umod( 1:n_moment,ie ) )
            !
            !enddo !ie
    
            do ie = 1 , nx_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,ie)
                elem(ie,je)%split1_modal( 1:ng )=element(ie)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(ie)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = x_2d(ie) + dx_2d * x_g(ig)
                    elem(ie,je)%psi1(ig,kk1) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                    elem(ie,je)%psi2(ig,kk1) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! je
     tnum = time1 -dt1+(cof_4thsplit(2)+cof_4thsplit(4))*dt1
    
    
    do ie = 1, nx_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do je = 1 ,ny_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi3( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !ftemp2 = ftemp2 + elem(ie,je)%psi2( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_y
    
            !*********************
            !do ii = 1-1,nel_y+1
            !    u( 1:nod+1,ii ) = elem( ie,ii )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_y
            !    elem( ie,ii )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_y
            !*********************
            do ii = 1,ny_2d
                element(ii)%u1modal(1:ng,0) = elem( ie,ii )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ie,ii )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 2,tnum, cof_4thsplit(4)*dt1, gx( kk1 ,1,ie,1) )
    
            !do je = 1 ,nel_y
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 2,dy,tnum,dt, gly( 1 ,1:n_gl,ie,je),gy( kk,1:n_g ,ie,je),gx( kk,1  ,ie,je) ,umod( 1:n_moment,je ) )
            !
            !enddo !je
    
            do je = 1 , ny_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,je)
                elem(ie,je)%split1_modal( 1:ng )=element(je)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(je)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = y_2d(je) + dy_2d * x_g(ig)
                    elem(ie,je)%psi1(kk1,ig) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    elem(ie,je)%psi3(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    !elem(ie,je)%psi2(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! ie
    
    !call output
    !pause
    !print *,1111111111
    !***************************************************************************************************
    !***************************************************************************************************
    tnum = time1 -dt1+(cof_4thsplit(1)+cof_4thsplit(3)+cof_4thsplit(5))*dt1
    do je = 1, ny_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do ie = 1 ,nx_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi2( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !write(*,*)ftemp
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_x
    
            !!*********************
            !do ii = 1-1,nel_x+1
            !    u( 1:nod+1,ii ) = elem( ii,je )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_x
            !    elem( ii,je )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_x
            !!*********************
            do ii = 1,nx_2d
                element(ii)%u1modal(1:ng,0) = elem( ii,je )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ii,je )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 1,tnum, cof_4thsplit(5)*dt1, gy( 1  ,kk1,1,je) )
            !do ie = 1 ,nel_x
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 1,dx_2d,tnum,dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) ,gy( 1  ,kk,ie,je) ,umod( 1:n_moment,ie ) )
            !
            !enddo !ie
    
            do ie = 1 , nx_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,ie)
                elem(ie,je)%split1_modal( 1:ng )=element(ie)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(ie)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = x_2d(ie) + dx_2d * x_g(ig)
                    elem(ie,je)%psi1(ig,kk1) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                    elem(ie,je)%psi2(ig,kk1) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! je
     tnum = time1 -dt1+(cof_4thsplit(2)+cof_4thsplit(4)+cof_4thsplit(6))*dt1
    
    
    do ie = 1, nx_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do je = 1 ,ny_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi3( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !ftemp2 = ftemp2 + elem(ie,je)%psi2( kk1,ig )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_y
    
            !*********************
            !do ii = 1-1,nel_y+1
            !    u( 1:nod+1,ii ) = elem( ie,ii )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_y
            !    elem( ie,ii )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_y
            !*********************
            do ii = 1,ny_2d
                element(ii)%u1modal(1:ng,0) = elem( ie,ii )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ie,ii )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 2,tnum, cof_4thsplit(6)*dt1, gx( kk1 ,1,ie,1) )
    
            !do je = 1 ,nel_y
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 2,dy,tnum,dt, gly( 1 ,1:n_gl,ie,je),gy( kk,1:n_g ,ie,je),gx( kk,1  ,ie,je) ,umod( 1:n_moment,je ) )
            !
            !enddo !je
    
            do je = 1 , ny_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,je)
                elem(ie,je)%split1_modal( 1:ng )=element(je)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(je)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = y_2d(je) + dy_2d * x_g(ig)
                    elem(ie,je)%psi1(kk1,ig) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    elem(ie,je)%psi3(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                    !elem(ie,je)%psi2(kk1,ig) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, y_2d(je) ,dy_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! ie
    
    !call output
    !pause
    !print *,1111111111
    !***************************************************************************************************
    !***************************************************************************************************
    tnum = time1 -dt1+(cof_4thsplit(1)+cof_4thsplit(3)+cof_4thsplit(5)+cof_4thsplit(7))*dt1
    do je = 1, ny_2d
    
        !call trans_to_split_modal
        do kk1= 1 ,ng
            do ie = 1 ,nx_2d
                !************************compute split modal
                do k1 = 1 , ng
                    ftemp1 = 0.
                    ftemp2 = 0.
                    do ig = 1,ng
                        ftemp1 = ftemp1 + elem(ie,je)%psi1( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        ftemp2 = ftemp2 + elem(ie,je)%psi2( ig,kk1 )*fle( k1-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split1_modal( k1 ) =  ftemp1 * ai(k1)
                    elem(ie,je)%split2_modal( k1 ) =  ftemp2 * ai(k1)
                enddo!k
                !************************end of computing split modal
            enddo
            !call bc_x
    
            !!*********************
            !do ii = 1-1,nel_x+1
            !    u( 1:nod+1,ii ) = elem( ii,je )%split_modal( 1:nod+1 )
            !enddo
            !call trouble_tvb(nod,nel,u(1:nod+1,1-1:nel+1) )
            !!
            !do ii = 1 , nel_x
            !    elem( ii,je )%split_modal( 1:nod+1 ) = u(1:nod+1,ii)
            !enddo
            !call bc_x
            !!*********************
            do ii = 1,nx_2d
                element(ii)%u1modal(1:ng,0) = elem( ii,je )%split1_modal( 1:ng )
                element(ii)%u2modal(1:ng,0) = elem( ii,je )%split2_modal( 1:ng )
            enddo
            call geldg_xy( 1,tnum, cof_4thsplit(7)*dt1, gy( 1  ,kk1,1,je) )
            !do ie = 1 ,nel_x
            !    ! call 1D SLDG subroutine
            !    call geldg_xy( 1,dx_2d,tnum,dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) ,gy( 1  ,kk,ie,je) ,umod( 1:n_moment,ie ) )
            !
            !enddo !ie
    
            do ie = 1 , nx_2d
                !elem(ie,je)%split_modal( 1:ng ) = umod(1:ng,ie)
                elem(ie,je)%split1_modal( 1:ng )=element(ie)%u1modal(1:ng,0)
                elem(ie,je)%split2_modal( 1:ng )=element(ie)%u2modal(1:ng,0)
                !************************compute nodal
                do ig = 1,ng
                    xrg = x_2d(ie) + dx_2d * x_g(ig)
                    elem(ie,je)%psi1(ig,kk1) = ortho_poly1d( elem(ie,je)%split1_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                    elem(ie,je)%psi2(ig,kk1) = ortho_poly1d( elem(ie,je)%split2_modal( 1:ng ), xrg, x_2d(ie) ,dx_2d,nk )
                enddo
                !************************end of computing nodal
            enddo
    
    
        enddo!kk
    
    enddo! je
    end subroutine splitting_4th