    subroutine allocate_var
    implicit none

    !*******************************
    !allocate_variables
    allocate( vertex(1:nx+1) )
    !allocate( vertex_star(1:nx+1) )
    allocate( vertex_star1(1:nx+1) )
    allocate( vertex_star2(1:nx+1) )
    !allocate( speed(1:nx+1) )
    allocate( speed1(1:nx+1) )
    allocate( speed2(1:nx+1) )
    allocate( u1post(nx*6,2) )
    allocate( u2post(nx*6,2) )
    allocate( x(1-nghost:nx+nghost) )
    allocate( xgrid(1-nghost:nx+1+nghost) )

    !allocate( element_star(1:nx) )
    allocate( element_star1(1:nx,0:4,1:4) )
    allocate( element_star2(1:nx,0:4,1:4) )

    allocate( element(1-nghost:nx+nghost) )
   
    !allocate( element_t(1-1:nx+1,0:4) )
    allocate( element_t1(1-1:nx+1,0:4,1:4) )
    allocate( element_t2(1-1:nx+1,0:4,1:4) )

    allocate( flux_LF(1:nx+1) )
    allocate( flux_LFp(1:nx+1) )
    allocate( flux_LF11(1:nx+1) )
    allocate( flux_LF21(1:nx+1) )
    allocate( flux_LF12(1:nx+1) )
    allocate( flux_LF22(1:nx+1) )
    allocate( flux_LF1(1:2,1:nx+1) )
    allocate( flux_LF2(1:2,1:nx+1) )
    allocate( flux_LF01(1:2,1:nx+1) )
    allocate( flux_LF02(1:2,1:nx+1) )
    allocate( flux_LF001(1:2,1:nx+1) )
    allocate( flux_LF002(1:2,1:nx+1) )
    allocate( flux_LF0001(1:2,1:nx+1) )
    allocate( flux_LF0002(1:2,1:nx+1) )
    
    allocate( itrouble(1:nx) )
    !end allocate_variables
    !*******************************
   allocate( A0(nk+1,nx),A1(nk+1,nx),A2(nk+1,nx) )
    end subroutine allocate_var
    !*********************************
    subroutine deallocate_var
    implicit none

    !*******************************
    !deallocate_variables
    deallocate( vertex )
    !deallocate( vertex_star )
    deallocate( vertex_star1 )
    deallocate( vertex_star2 )
    !deallocate( speed )
    deallocate( speed1 )
    deallocate( speed2 )
    deallocate( u1post )
    deallocate( u2post )
    deallocate( x )
    deallocate( xgrid )

    !deallocate( element_star )
    deallocate( element_star1 )
    deallocate( element_star2 )

    deallocate( element )

    !deallocate( element_t )
    deallocate( element_t1 )
    deallocate( element_t2 )

    deallocate( flux_LF )
    deallocate( flux_LFp )
    deallocate( flux_LF11 )
    deallocate( flux_LF21 )
    deallocate( flux_LF12 )
    deallocate( flux_LF22 )
    deallocate( flux_LF1 )
    deallocate( flux_LF2 )
    deallocate( flux_LF01 )
    deallocate( flux_LF02 )
    deallocate( flux_LF001 )
    deallocate( flux_LF002 )
    deallocate( flux_LF0001 )
    deallocate( flux_LF0002 )
    
    deallocate( itrouble )
    !end deallocate_variables
    !*******************************
    deallocate( A0 ,A1, A2  )
    end subroutine deallocate_var