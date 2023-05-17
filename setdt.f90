    subroutine setdt
    implicit none
    
    !call get_max_wave
   
    speed_max = 1. !for wave equation system 1
    
    ! set dt and get tn+1
    !dt = cfl*(dx**2.)/speed_max
    dt = cfl*(dx)/speed_max
    
    if(time+dt>time_final) dt = time_final- time

    time = time + dt

    !print *,time


    end subroutine setdt