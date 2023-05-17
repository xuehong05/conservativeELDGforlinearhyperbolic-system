    subroutine ssp_runge_kutta
    implicit none
    type(element1d_dynamic), pointer :: pet,pet_new,pet1,pet1_new,pet0,pet10,pet01,pet101,pet02,pet102

    integer :: Lda,Ldb,Ipivot(1,nk+1),INFO,kk1,kk2
    real ::  M_hold(nk+1,nk+1),b(nk+1),b1(nk+1),tt1,tt0,tt2,tt3,tt4,tt5,tt01,tt02!,A0(nk+1,nx),A1(nk+1,nx)
    real :: dx1
    
   
  
    if(irk==1)then
        tt0=0.
        tt1=0.
        tt2=0.
        tt3=0.
        tt4=0.
        tt5=0.
        do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            !pet_new => element_t1(i,iio,iio)
            pet1 => element_t2(i,iio-1,iio)
            !pet1_new => element_t2(i,iio,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               tt0=tt0+pet%phi_int(kk1,2)+pet1%phi_int(kk1,2)
               tt1=tt1+pet%source_int(kk1,2)+pet1%source_int(kk1,2)
               tt2=tt2+pet%vector_int(kk1,2)+pet1%vector_int(kk1,2)
               tt3=tt3-pet%flux1(kk1,2)+pet%flux0(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2)
               !write(*,*)-pet1%flux1(kk1,1), pet1%flux0(kk1,1),tt3
               !-pet1%flux1(kk1,1)+pet1%flux0(kk1,1) pet%flux1(kk1,1)+pet%flux0(kk1,1)
               tt4=tt4+pet%sum_int(kk1,2)+pet1%sum_int(kk1,2)
               !tt5=tt5+pet%sum_int(kk1,2)+pet1%sum_int(kk1,2)
               b(kk1)= pet%sum_int(kk1,1)+dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet1%sum_int(kk1,1)+dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet%sum_int(kk1,2)+dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet1%sum_int(kk1,2)+dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            !tt4=tt4+element(i)%u1modal(1,iio)
            !tt5=tt5+element(i)%u2modal(1,iio)
        enddo
        !write(*,*)tt4,tt5
        !write(*,*)'time',time,'source',tt2,'flux',tt3,'u_int',tt4
        !write(*,*)'t0',tt0,'t1',tt1,'t2',tt2,'t3',tt3,'t4',tt4
        !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
     !SSPRk2
    elseif(irk==2)then
        if(iio==1)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            !pet_new => element_t1(i,iio,iio)
            pet1 => element_t2(i,iio-1,iio)
            !pet1_new => element_t2(i,iio,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               !tt0=tt0+pet%phi_int(kk1,1)+pet1%phi_int(kk1,1)
               !tt1=tt1+pet%source_int(kk1,1)+pet1%source_int(kk1,1)
               !tt2=tt2+pet%vector_int(kk1,1)+pet1%vector_int(kk1,1)
               !tt3=tt3-pet1%flux1(kk1,1)+pet1%flux0(kk1,1)
               !write(*,*)-pet1%flux1(kk1,1), pet1%flux0(kk1,1),tt3
               !-pet1%flux1(kk1,1)+pet1%flux0(kk1,1) pet%flux1(kk1,1)+pet%flux0(kk1,1)
               !tt4=tt4+pet%sum_int(kk1,1)+pet1%sum_int(kk1,1)
               !tt5=tt5+pet%sum_int(kk1,2)+pet1%sum_int(kk1,2)
               b(kk1)= pet%sum_int(kk1,1)+dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet1%sum_int(kk1,1)+dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet%sum_int(kk1,2)+dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet1%sum_int(kk1,2)+dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            !tt4=tt4+element(i)%u1modal(1,iio)
            !tt5=tt5+element(i)%u2modal(1,iio)
            enddo
        elseif(iio==2)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet0%sum_int(kk1,1)+0.5*dt*(pet0%phi_int(kk1,1)+pet0%source_int(kk1,1)+pet0%vector_int(kk1,1)-pet0%flux1(kk1,1)+pet0%flux0(kk1,1)+pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet10%sum_int(kk1,1)+0.5*dt*(pet10%phi_int(kk1,1)+pet10%source_int(kk1,1)+pet10%vector_int(kk1,1)-pet10%flux1(kk1,1)+pet10%flux0(kk1,1)+pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet0%sum_int(kk1,2)+0.5*dt*(pet0%phi_int(kk1,2)+pet0%source_int(kk1,2)+pet0%vector_int(kk1,2)-pet0%flux1(kk1,2)+pet0%flux0(kk1,2)+pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet10%sum_int(kk1,2)+0.5*dt*(pet10%phi_int(kk1,2)+pet10%source_int(kk1,2)+pet10%vector_int(kk1,2)-pet10%flux1(kk1,2)+pet10%flux0(kk1,2)+pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            enddo
        endif
    elseif(irk==3)then
        if(iio==1)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            !pet_new => element_t1(i,iio,iio)
            pet1 => element_t2(i,iio-1,iio)
            !pet1_new => element_t2(i,iio,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               !tt0=tt0+pet%phi_int(kk1,1)+pet1%phi_int(kk1,1)
               !tt1=tt1+pet%source_int(kk1,1)+pet1%source_int(kk1,1)
               !tt2=tt2+pet%vector_int(kk1,1)+pet1%vector_int(kk1,1)
               !tt3=tt3-pet1%flux1(kk1,1)+pet1%flux0(kk1,1)
               !write(*,*)-pet1%flux1(kk1,1), pet1%flux0(kk1,1),tt3
               !-pet1%flux1(kk1,1)+pet1%flux0(kk1,1) pet%flux1(kk1,1)+pet%flux0(kk1,1)
               !tt4=tt4+pet%sum_int(kk1,1)+pet1%sum_int(kk1,1)
               !tt5=tt5+pet%sum_int(kk1,2)+pet1%sum_int(kk1,2)
               b(kk1)= pet%sum_int(kk1,1)+dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet1%sum_int(kk1,1)+dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet%sum_int(kk1,2)+dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet1%sum_int(kk1,2)+dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            !tt4=tt4+element(i)%u1modal(1,iio)
            !tt5=tt5+element(i)%u2modal(1,iio)
            enddo
        elseif(iio==2)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet0%sum_int(kk1,1)+(1./4.)*dt*(pet0%phi_int(kk1,1)+pet0%source_int(kk1,1)+pet0%vector_int(kk1,1)-pet0%flux1(kk1,1)+pet0%flux0(kk1,1)+pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet10%sum_int(kk1,1)+(1./4.)*dt*(pet10%phi_int(kk1,1)+pet10%source_int(kk1,1)+pet10%vector_int(kk1,1)-pet10%flux1(kk1,1)+pet10%flux0(kk1,1)+pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet0%sum_int(kk1,2)+(1./4.)*dt*(pet0%phi_int(kk1,2)+pet0%source_int(kk1,2)+pet0%vector_int(kk1,2)-pet0%flux1(kk1,2)+pet0%flux0(kk1,2)+pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet10%sum_int(kk1,2)+(1./4.)*dt*(pet10%phi_int(kk1,2)+pet10%source_int(kk1,2)+pet10%vector_int(kk1,2)-pet10%flux1(kk1,2)+pet10%flux0(kk1,2)+pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            enddo
        elseif(iio==3)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet01 => element_t1(i,iio-3,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            pet101 => element_t2(i,iio-3,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet01%sum_int(kk1,1)+(1./6.)*dt*(pet01%phi_int(kk1,1)+pet01%source_int(kk1,1)+pet01%vector_int(kk1,1)-pet01%flux1(kk1,1)+pet01%flux0(kk1,1))&
                   +(1./6.)*dt*(pet0%phi_int(kk1,1)+pet0%source_int(kk1,1)+pet0%vector_int(kk1,1)-pet0%flux1(kk1,1)+pet0%flux0(kk1,1))&
                   +(2./3.)*dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet101%sum_int(kk1,1)+(1./6.)*dt*(pet101%phi_int(kk1,1)+pet101%source_int(kk1,1)+pet101%vector_int(kk1,1)-pet101%flux1(kk1,1)+pet101%flux0(kk1,1))&
                       +(1./6.)*dt*(pet10%phi_int(kk1,1)+pet10%source_int(kk1,1)+pet10%vector_int(kk1,1)-pet10%flux1(kk1,1)+pet10%flux0(kk1,1))&
                       +(2./3.)*dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet01%sum_int(kk1,2)+(1./6.)*dt*(pet01%phi_int(kk1,2)+pet01%source_int(kk1,2)+pet01%vector_int(kk1,2)-pet01%flux1(kk1,2)+pet01%flux0(kk1,2))&
                   +(1./6.)*dt*(pet0%phi_int(kk1,2)+pet0%source_int(kk1,2)+pet0%vector_int(kk1,2)-pet0%flux1(kk1,2)+pet0%flux0(kk1,2))&
                   +(2./3.)*dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet101%sum_int(kk1,2)+(1./6.)*dt*(pet101%phi_int(kk1,2)+pet101%source_int(kk1,2)+pet101%vector_int(kk1,2)-pet101%flux1(kk1,2)+pet101%flux0(kk1,2))&
                       +(1./6.)*dt*(pet10%phi_int(kk1,2)+pet10%source_int(kk1,2)+pet10%vector_int(kk1,2)-pet10%flux1(kk1,2)+pet10%flux0(kk1,2))&
                       +(2./3.)*dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            enddo
        endif
        
    elseif(irk==4)then
        if(iio==1)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            !pet_new => element_t1(i,iio,iio)
            pet1 => element_t2(i,iio-1,iio)
            !pet1_new => element_t2(i,iio,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               !tt0=tt0+pet%phi_int(kk1,1)+pet1%phi_int(kk1,1)
               !tt1=tt1+pet%source_int(kk1,1)+pet1%source_int(kk1,1)
               !tt2=tt2+pet%vector_int(kk1,1)+pet1%vector_int(kk1,1)
               !tt3=tt3-pet1%flux1(kk1,1)+pet1%flux0(kk1,1)
               !write(*,*)-pet1%flux1(kk1,1), pet1%flux0(kk1,1),tt3
               !-pet1%flux1(kk1,1)+pet1%flux0(kk1,1) pet%flux1(kk1,1)+pet%flux0(kk1,1)
               !tt4=tt4+pet%sum_int(kk1,1)+pet1%sum_int(kk1,1)
               !tt5=tt5+pet%sum_int(kk1,2)+pet1%sum_int(kk1,2)
               b(kk1)= pet%sum_int(kk1,1)+0.5*dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet1%sum_int(kk1,1)+0.5*dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet%sum_int(kk1,2)+0.5*dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet1%sum_int(kk1,2)+0.5*dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            !tt4=tt4+element(i)%u1modal(1,iio)
            !tt5=tt5+element(i)%u2modal(1,iio)
            enddo
        elseif(iio==2)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet0%sum_int(kk1,1)+0.5*dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet10%sum_int(kk1,1)+0.5*dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet0%sum_int(kk1,2)+0.5*dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet10%sum_int(kk1,2)+0.5*dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            enddo
        elseif(iio==3)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet01 => element_t1(i,iio-3,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            pet101 => element_t2(i,iio-3,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet01%sum_int(kk1,1)+dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet101%sum_int(kk1,1)+dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet01%sum_int(kk1,2)+dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet101%sum_int(kk1,2)+dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            !write(*,*)element(i)%u1modal(1:nk+1,iio),element(i)%u2modal(1:nk+1,iio)
            enddo
        elseif(iio==4)then
            do i=1,nx
            ! consider a dynamic element
            pet => element_t1(i,iio-1,iio)
            pet0 => element_t1(i,iio-2,iio)
            pet01 => element_t1(i,iio-3,iio)
            pet02 => element_t1(i,iio-4,iio)
            pet1 => element_t2(i,iio-1,iio)
            pet10 => element_t2(i,iio-2,iio)
            pet101 => element_t2(i,iio-3,iio)
            pet102 => element_t2(i,iio-4,iio)
            dx1=element(i)%xr-element(i)%xl
            do kk1 = 1,nk+1
               b(kk1)= pet02%sum_int(kk1,1)+(1./6.)*dt*(pet02%phi_int(kk1,1)+pet02%source_int(kk1,1)+pet02%vector_int(kk1,1)-pet02%flux1(kk1,1)+pet02%flux0(kk1,1))&
                   +(1./3.)*dt*(pet01%phi_int(kk1,1)+pet01%source_int(kk1,1)+pet01%vector_int(kk1,1)-pet01%flux1(kk1,1)+pet01%flux0(kk1,1))&
                   +(1./3.)*dt*(pet0%phi_int(kk1,1)+pet0%source_int(kk1,1)+pet0%vector_int(kk1,1)-pet0%flux1(kk1,1)+pet0%flux0(kk1,1))&
                   +(1./6.)*dt*(pet%phi_int(kk1,1)+pet%source_int(kk1,1)+pet%vector_int(kk1,1)-pet%flux1(kk1,1)+pet%flux0(kk1,1))&
                   +pet102%sum_int(kk1,1)+(1./6.)*dt*(pet102%phi_int(kk1,1)+pet102%source_int(kk1,1)+pet102%vector_int(kk1,1)-pet102%flux1(kk1,1)+pet102%flux0(kk1,1))&
                       +(1./3.)*dt*(pet101%phi_int(kk1,1)+pet101%source_int(kk1,1)+pet101%vector_int(kk1,1)-pet101%flux1(kk1,1)+pet101%flux0(kk1,1))&
                       +(1./3.)*dt*(pet10%phi_int(kk1,1)+pet10%source_int(kk1,1)+pet10%vector_int(kk1,1)-pet10%flux1(kk1,1)+pet10%flux0(kk1,1))&
                       +(1./6.)*dt*(pet1%phi_int(kk1,1)+pet1%source_int(kk1,1)+pet1%vector_int(kk1,1)-pet1%flux1(kk1,1)+pet1%flux0(kk1,1))
               b1(kk1)= pet02%sum_int(kk1,2)+(1./6.)*dt*(pet02%phi_int(kk1,2)+pet02%source_int(kk1,2)+pet02%vector_int(kk1,2)-pet02%flux1(kk1,2)+pet02%flux0(kk1,2))&
                   +(1./3.)*dt*(pet01%phi_int(kk1,2)+pet01%source_int(kk1,2)+pet01%vector_int(kk1,2)-pet01%flux1(kk1,2)+pet01%flux0(kk1,2))&
                   +(1./3.)*dt*(pet0%phi_int(kk1,2)+pet0%source_int(kk1,2)+pet0%vector_int(kk1,2)-pet0%flux1(kk1,2)+pet0%flux0(kk1,2))&
                   +(1./6.)*dt*(pet%phi_int(kk1,2)+pet%source_int(kk1,2)+pet%vector_int(kk1,2)-pet%flux1(kk1,2)+pet%flux0(kk1,2))&
                   +pet102%sum_int(kk1,2)+(1./6.)*dt*(pet102%phi_int(kk1,2)+pet102%source_int(kk1,2)+pet102%vector_int(kk1,2)-pet102%flux1(kk1,2)+pet102%flux0(kk1,2))&
                       +(1./3.)*dt*(pet101%phi_int(kk1,2)+pet101%source_int(kk1,2)+pet101%vector_int(kk1,2)-pet101%flux1(kk1,2)+pet101%flux0(kk1,2))&
                       +(1./3.)*dt*(pet10%phi_int(kk1,2)+pet10%source_int(kk1,2)+pet10%vector_int(kk1,2)-pet10%flux1(kk1,2)+pet10%flux0(kk1,2))&
                       +(1./6.)*dt*(pet1%phi_int(kk1,2)+pet1%source_int(kk1,2)+pet1%vector_int(kk1,2)-pet1%flux1(kk1,2)+pet1%flux0(kk1,2))
               b(kk1)=b(kk1)*ai(kk1)/dx1
               b1(kk1)=b1(kk1)*ai(kk1)/dx1
               !write(*,*)pet%sum_int(kk1,1),pet1%sum_int(kk1,1),pet%sum_int(kk1,2),pet1%sum_int(kk1,2)
            enddo
            element(i)%u1modal(1:nk+1,iio) = b(:)
            element(i)%u2modal(1:nk+1,iio) = b1(:)
            enddo
            endif !io  
            
            
    !!!RK2
    !!elseif(irk==2)then
    !!    if(iio==1)then
    !!        do i=1,nx
    !!        ! consider a dynamic element
    !!        pet => element_t1(i,iio-1,iio)
    !!        !pet0 => element_t1(i,iio-2,iio)
    !!        pet_new => element_t1(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           !tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !!              !tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !!           enddo
    !!           b(kk1)= tt1+0.5*dt*(pet%sum_int(kk1)+pet%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !!        
    !!        pet1 => element_t2(i,iio-1,iio)
    !!        !pet10 => element_t2(i,iio-2,iio)
    !!        pet1_new => element_t2(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        !dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !!           enddo
    !!           b(kk1)= tt1+0.5*dt*(pet1%sum_int(kk1)+pet1%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !!        element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        enddo
    !!        
    !!    elseif(iio==2)then
    !!        do i=1,nx
    !!        ! consider a dynamic element
    !!        pet => element_t1(i,iio-1,iio)
    !!        pet0 => element_t1(i,iio-2,iio)
    !!        pet_new => element_t1(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !!              tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !!           enddo
    !!           b(kk1)= tt0+dt*(pet%sum_int(kk1)+pet%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !!        
    !!        pet1 => element_t2(i,iio-1,iio)
    !!        pet10 => element_t2(i,iio-2,iio)
    !!        pet1_new => element_t2(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        !dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !!              tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !!           enddo
    !!           b(kk1)= tt0+dt*(pet1%sum_int(kk1)+pet1%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !!        element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        enddo
    !!    endif !io
    !!
    !  !elseif(irk==3)then
    !  !  if(iio==1)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      !pet0 => element_t1(i,iio-2,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         !tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            !tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= tt1+dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      !pet10 => element_t2(i,iio-2,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= tt1+dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !  !      
    !  !  elseif(iio==2)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      pet0 => element_t1(i,iio-2,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 0.75*tt0+0.25*tt1+0.25*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      pet10 => element_t2(i,iio-2,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !            tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 0.75*tt0+0.25*tt1+0.25*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !  !      
    !  !  elseif(iio==3)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      pet0 => element_t1(i,iio-2,iio)
    !  !      pet01 => element_t1(i,iio-3,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         tt01=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !            tt01=tt01+pet01%massmetrix(kk1,kk2)*pet01%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= (1./3.)*tt01+(2./3.)*tt1+(2./3.)*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      pet10 => element_t2(i,iio-2,iio)
    !  !      pet101 => element_t2(i,iio-3,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         tt01=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !            tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !  !            tt01=tt01+pet101%massmetrix(kk1,kk2)*pet101%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= (1./3.)*tt01+(2./3.)*tt1+(2./3.)*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !  !  endif !io
    !    !SSPRk3
    
    !!!    !Rk3
    !!elseif(irk==3)then
    !!    if(iio==1)then
    !!        do i=1,nx
    !!        ! consider a dynamic element
    !!        pet => element_t1(i,iio-1,iio)
    !!        !pet0 => element_t1(i,iio-2,iio)
    !!        pet_new => element_t1(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           !tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !!              !tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !!           enddo
    !!           b(kk1)= tt1+(1./3.)*dt*(pet%sum_int(kk1)+pet%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !!        
    !!        pet1 => element_t2(i,iio-1,iio)
    !!        !pet10 => element_t2(i,iio-2,iio)
    !!        pet1_new => element_t2(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        !dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !!           enddo
    !!           b(kk1)= tt1+(1./3.)*dt*(pet1%sum_int(kk1)+pet1%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !!        element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        enddo
    !!        
    !!    elseif(iio==2)then
    !!        do i=1,nx
    !!        ! consider a dynamic element
    !!        pet => element_t1(i,iio-1,iio)
    !!        pet0 => element_t1(i,iio-2,iio)
    !!        pet_new => element_t1(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !!              tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !!           enddo
    !!           b(kk1)= tt0+(2./3.)*dt*(pet%sum_int(kk1)+pet%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !!        
    !!        pet1 => element_t2(i,iio-1,iio)
    !!        pet10 => element_t2(i,iio-2,iio)
    !!        pet1_new => element_t2(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        !dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !!              tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !!           enddo
    !!           b(kk1)= tt0+(2./3.)*dt*(pet1%sum_int(kk1)+pet1%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !!        element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        enddo
    !!        
    !!    elseif(iio==3)then
    !!        do i=1,nx
    !!        ! consider a dynamic element
    !!        pet => element_t1(i,iio-1,iio)
    !!        pet0 => element_t1(i,iio-2,iio)
    !!        pet01 => element_t1(i,iio-3,iio)
    !!        pet_new => element_t1(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           tt01=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !!              tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !!              tt01=tt01+pet01%massmetrix(kk1,kk2)*pet01%v1modal(kk2)
    !!           enddo
    !!           b(kk1)= tt01+(1./4.)*dt*(pet01%sum_int(kk1)+pet01%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF01(1,i+1)+element(i)%RT_j(1,2)*flux_LF01(2,i+1))*ortho_poly1d( pet01%aa(1:nk+1,kk1),pet01%xr,pet01%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF01(1,i)+element(i)%RT_j(1,2)*flux_LF01(2,i))*ortho_poly1d( pet01%aa(1:nk+1,kk1),pet01%xl,pet01%weno_t(1),dx1,nk))+(3./4.)*dt*(pet%sum_int(kk1)+pet%source_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !!        enddo
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !!        
    !!        pet1 => element_t2(i,iio-1,iio)
    !!        pet10 => element_t2(i,iio-2,iio)
    !!        pet101 => element_t2(i,iio-3,iio)
    !!        pet1_new => element_t2(i,iio,iio)
    !!        Lda=nk+1
    !!        Ldb=nk+1
    !!        !dx1=element(i)%xr-element(i)%xl
    !!        do kk1 = 1,nk+1
    !!           tt1=0.
    !!           tt0=0.
    !!           tt01=0.
    !!           do kk2 = 1,nk+1
    !!              M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !!              tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !!              tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !!              tt01=tt01+pet101%massmetrix(kk1,kk2)*pet101%v2modal(kk2)
    !!           enddo
    !!           b(kk1)= tt01+(1./4.)*dt*(pet101%sum_int(kk1)+pet101%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF02(1,i+1)+element(i)%RT_j(2,2)*flux_LF02(2,i+1))*ortho_poly1d( pet101%aa(1:nk+1,kk1),pet101%xr,pet101%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF02(1,i)+element(i)%RT_j(2,2)*flux_LF02(2,i))*ortho_poly1d( pet101%aa(1:nk+1,kk1),pet101%xl,pet101%weno_t(1),dx1,nk))+(3./4.)*dt*(pet1%sum_int(kk1)+pet1%source_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !!               +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !!        enddo
    !!        
    !!        call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !!        element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !!        element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !!        enddo
    !!    endif !io
    !  !elseif(irk==4)then
    !  !  if(iio==1)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      !pet0 => element_t1(i,iio-2,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         !tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            !tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= tt1+0.5*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      !pet10 => element_t2(i,iio-2,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= tt1+0.5*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !  !      
    !  !  elseif(iio==2)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      pet0 => element_t1(i,iio-2,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 1.*tt0+0.*tt1+0.5*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      pet10 => element_t2(i,iio-2,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !            tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 1.*tt0+0.*tt1+0.5*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !  !      
    !  !  elseif(iio==3)then
    !  !      do i=1,nx
    !  !      ! consider a dynamic element
    !  !      pet => element_t1(i,iio-1,iio)
    !  !      pet0 => element_t1(i,iio-2,iio)
    !  !      pet01 => element_t1(i,iio-3,iio)
    !  !      pet_new => element_t1(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         tt01=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !  !            tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !  !            tt01=tt01+pet01%massmetrix(kk1,kk2)*pet01%v1modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 1.*tt01+0.*tt0+0.*tt1+1.*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !  !      
    !  !      pet1 => element_t2(i,iio-1,iio)
    !  !      pet10 => element_t2(i,iio-2,iio)
    !  !      pet101 => element_t2(i,iio-3,iio)
    !  !      pet1_new => element_t2(i,iio,iio)
    !  !      Lda=nk+1
    !  !      Ldb=nk+1
    !  !      !dx1=element(i)%xr-element(i)%xl
    !  !      do kk1 = 1,nk+1
    !  !         tt1=0.
    !  !         tt0=0.
    !  !         tt01=0.
    !  !         do kk2 = 1,nk+1
    !  !            M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !  !            tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !  !            tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !  !            tt01=tt01+pet101%massmetrix(kk1,kk2)*pet101%v2modal(kk2)
    !  !         enddo
    !  !         b(kk1)= 1.*tt01+0.*tt0+0.*tt1+1.*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !  !             +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !  !      enddo
    !  !      call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !  !      element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !  !      element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !  !      enddo
    !    !elseif(iio==4)then
    !    !    do i=1,nx
    !    !    ! consider a dynamic element
    !    !    pet => element_t1(i,iio-1,iio)
    !    !    pet0 => element_t1(i,iio-2,iio)
    !    !    pet01 => element_t1(i,iio-3,iio)
    !    !    pet02 => element_t1(i,iio-4,iio)
    !    !    pet_new => element_t1(i,iio,iio)
    !    !    Lda=nk+1
    !    !    Ldb=nk+1
    !    !    dx1=element(i)%xr-element(i)%xl
    !    !    do kk1 = 1,nk+1
    !    !       tt1=0.
    !    !       tt0=0.
    !    !       tt01=0.
    !    !       tt02=0.
    !    !       do kk2 = 1,nk+1
    !    !          M_hold(kk1,kk2)=pet_new%massmetrix(kk1,kk2)
    !    !          tt1=tt1+pet%massmetrix(kk1,kk2)*pet%v1modal(kk2)
    !    !          tt0=tt0+pet0%massmetrix(kk1,kk2)*pet0%v1modal(kk2)
    !    !          tt01=tt01+pet01%massmetrix(kk1,kk2)*pet01%v1modal(kk2)
    !    !          tt02=tt02+pet02%massmetrix(kk1,kk2)*pet02%v1modal(kk2)
    !    !       enddo
    !    !       b(kk1)= -(1./3.)*tt02+(1./3.)*tt01+(2./3.)*tt0+(1./3.)*tt1+(1./6.)*dt*(pet%sum_int(kk1)-(element(i)%RT_j(1,1)*flux_LF1(1,i+1)+element(i)%RT_j(1,2)*flux_LF1(2,i+1))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xr,pet%weno_t(1),dx1,nk)&
    !    !           +(element(i)%RT_j(1,1)*flux_LF1(1,i)+element(i)%RT_j(1,2)*flux_LF1(2,i))*ortho_poly1d( pet%aa(1:nk+1,kk1),pet%xl,pet%weno_t(1),dx1,nk))
    !    !    enddo
    !    !    call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !    !    element_t1(i,iio,iio)%v1modal(1:nk+1)=b
    !    !    
    !    !    pet1 => element_t2(i,iio-1,iio)
    !    !    pet10 => element_t2(i,iio-2,iio)
    !    !    pet101 => element_t2(i,iio-3,iio)
    !    !    pet102 => element_t2(i,iio-4,iio)
    !    !    pet1_new => element_t2(i,iio,iio)
    !    !    Lda=nk+1
    !    !    Ldb=nk+1
    !    !    !dx1=element(i)%xr-element(i)%xl
    !    !    do kk1 = 1,nk+1
    !    !       tt1=0.
    !    !       tt0=0.
    !    !       tt01=0.
    !    !       tt02=0.
    !    !       do kk2 = 1,nk+1
    !    !          M_hold(kk1,kk2)=pet1_new%massmetrix(kk1,kk2)
    !    !          tt1=tt1+pet1%massmetrix(kk1,kk2)*pet1%v2modal(kk2)
    !    !          tt0=tt0+pet10%massmetrix(kk1,kk2)*pet10%v2modal(kk2)
    !    !          tt01=tt01+pet101%massmetrix(kk1,kk2)*pet101%v2modal(kk2)
    !    !          tt02=tt02+pet102%massmetrix(kk1,kk2)*pet102%v2modal(kk2)
    !    !       enddo
    !    !       b(kk1)= -(1./3.)*tt02+(1./3.)*tt01+(2./3.)*tt0+(1./3.)*tt1+(1./6.)*dt*(pet1%sum_int(kk1)-(element(i)%RT_j(2,1)*flux_LF2(1,i+1)+element(i)%RT_j(2,2)*flux_LF2(2,i+1))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xr,pet1%weno_t(1),dx1,nk)&
    !    !           +(element(i)%RT_j(2,1)*flux_LF2(1,i)+element(i)%RT_j(2,2)*flux_LF2(2,i))*ortho_poly1d( pet1%aa(1:nk+1,kk1),pet1%xl,pet1%weno_t(1),dx1,nk))
    !    !    enddo
    !    !    call DGESV(nk+1,1,M_hold,Lda,Ipivot,b,Ldb,INFO)
    !    !    element_t2(i,iio,iio)%v2modal(1:nk+1)=b 
    !    !    element(i)%u1modal(1:nk+1,iio) = element(i)%R_j(1,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(1,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !    !    element(i)%u2modal(1:nk+1,iio) = element(i)%R_j(2,1)*element_t1(i,iio,iio)%v1modal(1:nk+1)+element(i)%R_j(2,2)*element_t2(i,iio,iio)%v2modal(1:nk+1)
    !    !    enddo
    !    !endif !io
   
      endif!irk
      
    
    
    end subroutine ssp_runge_kutta
    
    
    
    