    !subroutine get_element_tmass(pet1,xl,xr)
    !implicit none
    !type(element1d_dynamic), pointer :: pet1
    !real,intent(in) :: xl,xr
    !real :: gau_tt(1:3), flux,st,dx1,st1
    !integer :: ig
    !integer :: nm,nm1,nm2
    !
    !pet1%aa(1,1) = 1.;  pet1%aa(2,1) = 0.;  pet1%aa(3,1) = 0.;
    !pet1%aa(1,2) = 0.;  pet1%aa(2,2) = 1.;  pet1%aa(3,2) = 0.;
    !pet1%aa(1,3) = 0.;  pet1%aa(2,3) = 0.;  pet1%aa(3,3) = 1.;
    !
    !dx1=xr-xl
    !gau_tt(1:3) = pet1%xc + gau_f_phi_x(1:3,1)*pet1%dx
    !
    !do nm1 = 1,nk+1
    !   do nm2 = 1,nk+1
    !       
    !        st1 =0.
    !        do ig = 1,3
    !           
    !            st1=st1+ortho_poly1d( pet1%aa(1:nk+1,nm2),gau_tt(ig),pet1%xc,pet1%dx,nk)&
    !                 *ortho_poly1d( pet1%aa(1:nk+1,nm1),gau_tt(ig),pet1%xc,pet1%dx,nk)*gau_f_phi_x(ig,2)
    !                !*ortho_poly1d( pet1%aa(1:nk+1,nm1),gau_tt(ig),pet1%weno_t(1),dx1,nk)*gau_f_phi_x(ig,2)
    !        enddo
    !        pet1%massmetrix(nm1,nm2) =  st1*pet1%dx 
    !   enddo
    ! enddo
    ! 
    !
    !end subroutine get_element_tmass