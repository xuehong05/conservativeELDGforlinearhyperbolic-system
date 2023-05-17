    subroutine get_max_sl_wave
    implicit none
    type(element1d_dynamic), pointer :: pet
    real :: speed_sl
    real :: alpha1
    real :: ul,ur
    ! compute alpha for global Lax-Friedrichs flux
    !speed_sl = 0.
    !do i = 1,nx+1
    !    pet=>element_t1(i-1,io)
    !    ul = ortho_poly1d(pet%umodal(1:nk+1),pet%xr,pet%xc,pet%dx,nk)
    !
    !    pet=>element_t1(i,io)
    !    ur = ortho_poly1d(pet%umodal(1:nk+1),pet%xl,pet%xc,pet%dx,nk)
    !    alpha1 = 1.0*max( abs(fp( ul,element_t(i-1,io)%xr  )-speed(i) ), &
    !        abs(fp( ur,element_t(i,io)%xl )-speed(i) ) )
    !
    !    
    !    speed_sl = max( alpha1,speed_sl ) 
    !enddo
    !
    !open(19,file='speed_sl_history.txt')
    !write(19,*)  time,speed_sl !/dx
    
    end subroutine get_max_sl_wave