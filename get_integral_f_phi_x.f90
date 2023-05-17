    !subroutine get_integral_f_phi_x(pet,pet1,RT_j)
    !implicit none
    !type(element1d_dynamic), pointer :: pet,pet1
    !real,intent(in) :: RT_j(1:2,1:2)
    !real :: gau_tt(1:3),gau_tt1(1:3),flux,flux1,st,dx1,st1,stt,stt1
    !integer :: ig
    !integer :: nm,nm1,nm2
    !
    !pet%aa(1,1) = 1.;  pet%aa(2,1) = 0.;  pet%aa(3,1) = 0.;
    !pet%aa(1,2) = 0.;  pet%aa(2,2) = 1.;  pet%aa(3,2) = 0.;
    !pet%aa(1,3) = 0.;  pet%aa(2,3) = 0.;  pet%aa(3,3) = 1.;
    !
    !pet1%aa(1,1) = 1.;  pet1%aa(2,1) = 0.;  pet1%aa(3,1) = 0.;
    !pet1%aa(1,2) = 0.;  pet1%aa(2,2) = 1.;  pet1%aa(3,2) = 0.;
    !pet1%aa(1,3) = 0.;  pet1%aa(2,3) = 0.;  pet1%aa(3,3) = 1.;
    !
    !!dx1=xr-xl
    !gau_tt(1:3) = pet%xc + gau_f_phi_x(1:3,1)*pet%dx
    !gau_tt1(1:3) = pet1%xc + gau_f_phi_x(1:3,1)*pet1%dx
    !
    !do nm = 1,nk+1
    !    pet%sum_int(nm) = 0.
    !    pet1%sum_int(nm) = 0.
    !    stt=0.
    !    stt1=0.
    !    do ig = 1,3
    !        stt=stt + f_source(gau_tt(ig),pet%midtime)&
    !            *ortho_poly1d(pet%aa(1:nk+1,nm),gau_tt(ig),pet%xc,pet%dx,nk)*gau_f_phi_x(ig,2)
    !        stt1=stt1 + f_source(gau_tt1(ig),pet1%midtime)&
    !            *ortho_poly1d(pet1%aa(1:nk+1,nm),gau_tt1(ig),pet1%xc,pet1%dx,nk)*gau_f_phi_x(ig,2)
    !    enddo
    !    pet%source_int(nm) = RT_j(1,1)*stt*pet%dx
    !    pet1%source_int(nm) = RT_j(2,1)*stt1*pet1%dx
    !    if(nm>1)then
    !        st =0.
    !        st1 =0.
    !        do ig = 1,3
    !            !compute the valve of F=f(u,x)-a(i)u in Gauss points
    !            !call compute_flux_in_upstream( nk,pet%umodal(1:nk+1),gau_tt(ig),pet%xc,pet%dx,speed(i),speed(i+1),alpha,pet%xl,pet%xr,flux )
    !            call compute_flux_in_upstream(nk,pet%u1modal(1:nk+1),pet%u2modal(1:nk+1),RT_j(1,1),RT_j(1,2),gau_tt(ig),pet%xc,pet%dx,speed1(i),speed1(i+1),pet%xl,pet%xr,flux)
    !            call compute_flux_in_upstream(nk,pet1%u1modal(1:nk+1),pet1%u2modal(1:nk+1),RT_j(2,1),RT_j(2,2),gau_tt1(ig),pet1%xc,pet1%dx,speed2(i),speed2(i+1),pet1%xl,pet1%xr,flux1)
    !            !call compute_flux_in_upstream( nk,pet%umodal(1:nk+1),gau_tt(ig),pet%xc,pet%dx,speed(i),speed(i+1),alpha,xl,xr,flux )
    !            st = st + flux&
    !                !*diff_ortho_poly1d( pet%aa(1:nk+1,nm),gau_tt(ig),pet%xc,pet%dx,nk)*gau_f_phi_x(ig,2)
    !                *diff_ortho_poly1d( pet%aa(1:nk+1,nm),gau_tt(ig),pet%xc,pet%dx,nk)*gau_f_phi_x(ig,2)
    !            st1 = st1 + flux1&
    !                !*diff_ortho_poly1d( pet%aa(1:nk+1,nm),gau_tt(ig),pet%xc,pet%dx,nk)*gau_f_phi_x(ig,2)
    !                *diff_ortho_poly1d( pet1%aa(1:nk+1,nm),gau_tt1(ig),pet1%xc,pet1%dx,nk)*gau_f_phi_x(ig,2)
    !        enddo
    !        
    !        pet%sum_int(nm) =  st*pet%dx
    !        pet1%sum_int(nm) =  st1*pet1%dx
    !    endif
    !enddo
    !
    !
    ! !do nm1 = 1,nk+1
    ! !  do nm2 = 1,nk+1
    ! !      
    ! !       st1 =0.
    ! !       do ig = 1,3
    ! !          
    ! !           st1=st1+ortho_poly1d( pet%aa(1:nk+1,nm2),gau_tt(ig),pet%xc,pet%dx,nk)&
    ! !               *ortho_poly1d( pet%aa(1:nk+1,nm1),gau_tt(ig),pet%weno_t(1),dx1,nk)*gau_f_phi_x(ig,2)
    ! !       enddo
    ! !       pet%massmetrix(nm1,nm2) =  st1*pet%dx 
    ! !  enddo
    ! !enddo
    ! 
    !
    !end subroutine get_integral_f_phi_x