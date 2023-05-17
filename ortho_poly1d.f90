    real function ortho_poly1d(a,x ,xc ,dx ,k)
    implicit none
    real,intent(in) :: x ,xc ,dx 
    integer,intent(in) :: k
    real,intent(in) :: a(k+1)


    if(k==0)then
        ortho_poly1d = a(1)
    elseif(k==1)then
        ortho_poly1d = a(1) + a(2)*(x  - xc  )/dx 
    elseif(k==2)then
        ortho_poly1d = a(1) + a(2)*(x  - xc  )/dx  +  a(3)*(  ((x - xc ) /dx )**2 - 1./12. )
    elseif(k==3)then
        ortho_poly1d = a(1) + a(2)*(x - xc  )/dx  +  a(3)*(  ((x - xc ) /dx )**2 - 1./12. ) &
            +a(4) *( ((x - xc ) /dx )**3 - 3./20.*((x - xc ) /dx ) )
    endif

    end function ortho_poly1d
    !****************************************************************************************
    real function diff_ortho_poly1d(a,x ,xc ,dx ,k)
    implicit none
    real,intent(in) :: x ,xc ,dx 
    integer,intent(in) :: k
    real,intent(in) :: a(k+1)

    if(k==0)then
        diff_ortho_poly1d = 0.
    elseif(k==1)then
        diff_ortho_poly1d =  a(2)/dx 
    elseif(k==2)then
        diff_ortho_poly1d = a(2)/dx  +  a(3)*2.*( (x-xc)/dx )/dx
    endif

    end function diff_ortho_poly1d