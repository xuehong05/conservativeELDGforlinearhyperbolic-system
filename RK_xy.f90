    subroutine RK_xy
    implicit none
    
     call get_upstream_tn
     do iio = 1,irk
        if (iio >1 ) then
            call boundary_eulerian(iio-1) 
            do io=0,iio-1
                call projection(io,iio )
            enddo
        endif
        call ssp_runge_kutta
     enddo
     !!do io = 0,irk-1
     !!   call rhs_int_flux
     !!       
     !!   !if(io==0)then
     !!   !    call get_max_sl_wave
     !!   !endif
     !!   call boundary_lag
     !!       
     !!           !do i=0,nx+1
     !!           !    write(*,*) element_t(i,io)%umodal(1:nk+1)
     !!           !enddo
     !!           !pause
     !!           
     !!   call ssp_runge_kutta
     !!enddo !io
     !! update solution on background elements
     do i = 1,nx
        !u1mean=u1mean+element(i)%u1modal(1,irk)*dx
        !u2mean=u2mean+element(i)%u2modal(1,irk)*dx
        element(i)%u1modal(1:nk+1,0)=element(i)%u1modal(1:nk+1,irk)
        element(i)%u2modal(1:nk+1,0)=element(i)%u2modal(1:nk+1,irk)
                !write(*,*) element(i)%u1modal(1:nk+1)
     enddo
     
     !       do i = 1,nx
     !           element(i)%umodal(1:nk+1) = element_t(i,irk)%umodal(1:nk+1)
     !           !write(*,*) element(i)%umodal(1:nk+1)
     !       enddo
    end subroutine RK_xy