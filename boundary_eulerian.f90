    !*************************************************************************
    ! periodic boundary conditions
    !*************************************************************************
    subroutine boundary_eulerian(m)
    implicit none
    
    integer :: m

    do i = 1,nghost
        !element(1-i)%u1modal(1:nk+1) = element(nx+1-i)%u1modal(1:nk+1)
        !element(1-i)%u2modal(1:nk+1) = element(nx+1-i)%u2modal(1:nk+1)
        !element(nx+i)%u1modal(1:nk+1) = element(0+i)%u1modal(1:nk+1)
        !element(nx+i)%u2modal(1:nk+1) = element(0+i)%u2modal(1:nk+1)
        element(1-i)%u1modal(1:nk+1,m) = element(nx+1-i)%u1modal(1:nk+1,m)
        element(1-i)%u2modal(1:nk+1,m) = element(nx+1-i)%u2modal(1:nk+1,m)
        element(nx+i)%u1modal(1:nk+1,m) = element(0+i)%u1modal(1:nk+1,m)
        element(nx+i)%u2modal(1:nk+1,m) = element(0+i)%u2modal(1:nk+1,m)
        !element(1-i)%weno_t(1) = element(nx+1-i)%weno_t(1)
        !element(nx+i)%weno_t(1) = element(0+i)%weno_t(1)
        !element(1-i)%weno_t(2) = element(nx+1-i)%weno_t(2)
        !element(nx+i)%weno_t(2) = element(0+i)%weno_t(2)
        !element(1-i)%R_j(1:2,1:2) = element(nx+1-i)%R_j(1:2,1:2)
        !element(nx+i)%R_j(1:2,1:2) = element(0+i)%R_j(1:2,1:2)
        !element(1-i)%RT_j(1:2,1:2) = element(nx+1-i)%RT_j(1:2,1:2)
        !element(nx+i)%RT_j(1:2,1:2) = element(0+i)%RT_j(1:2,1:2)
       
        !!Dirichlet boundary condition
        !if(iexample==3)then
        !    element(1-i)%umodal(1:1) = 1.
        !    element(1-i)%umodal(2:nk+1) = 0.
        !    element(nx+i)%umodal(1:nk+1) = 0.
        !    element(1-i)%weno_t(1) = element(nx+1-i)%weno_t(1)
        !    element(nx+i)%weno_t(1) = element(0+i)%weno_t(1)
        !endif
    enddo
    !do io=0,irk
    !    do i = 1,nghost
    !    element(1-i)%u1modal(1:nk+1,io) = element(nx+1-i)%u1modal(1:nk+1,io)
    !    element(1-i)%u2modal(1:nk+1,io) = element(nx+1-i)%u2modal(1:nk+1,io)
    !    element(nx+i)%u1modal(1:nk+1,io) = element(0+i)%u1modal(1:nk+1,io)
    !    element(nx+i)%u2modal(1:nk+1,io) = element(0+i)%u2modal(1:nk+1,io)
    !    enddo
    !enddo
    
            
    end subroutine boundary_eulerian