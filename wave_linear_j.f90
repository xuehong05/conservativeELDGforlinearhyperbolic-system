
    !real function wave_linear_j(x,sl,sr,xl,xr,dxj,dt_b)
    !implicit none
    !real,intent(in) :: x,sl,sr,xl,xr,dxj,dt_b
    !real :: t1,t2
    !
    !t1 = ( sl*xr-sr*xl )/( xl-xr )
    !t2 = (sl-sr)/(xl-xr)
    ! 
    !wave_linear_j = -t1 + t2*x + t2*dt_b*( sl*xr-sr*xl )/( dt_b*sl-dt_b*sr-xl+xr )&
    !    - t2* x*dt_b *(sl-sr)/(dt_b*sl-dt_b*sr-xl+xr)
    !
    !
    !end function wave_linear_j
    
    real function wave_linear_j(x,sl,sr,xl,xr,dxj)
    implicit none
    real,intent(in) :: x,sl,sr,xl,xr,dxj
    !real :: t1,t2
    
    !t1 = ( sl*xr-sr*xl )/( xl-xr )
    !t2 = (sl-sr)/(xl-xr)
     
    wave_linear_j = sr*(x-xl)/dxj+sl*(xr-x)/dxj
    
    
    end function wave_linear_j