    subroutine bc_x
    implicit none
    integer :: ibc
    
    do ibc = 1 , nghost
        elem( 1-ibc,je)%split1_modal( 1:ng ) = elem( nx_2d+1 -ibc,je)%split1_modal( 1:ng )
        elem( nx_2d+ibc,je)%split1_modal( 1:ng ) = elem(  ibc,je)%split1_modal( 1:ng )
        elem( 1-ibc,je)%split2_modal( 1:ng ) = elem( nx_2d+1 -ibc,je)%split2_modal( 1:ng )
        elem( nx_2d+ibc,je)%split2_modal( 1:ng ) = elem(  ibc,je)%split2_modal( 1:ng )
    enddo
    
    end subroutine bc_x
    !*********************************************************************
    subroutine bc_y
    implicit none
    integer :: ibc
    
    do ibc = 1 , nghost
        elem( ie, 1-ibc )%split1_modal( 1:ng ) = elem( ie,ny_2d+1 -ibc )%split1_modal( 1:ng )
        elem( ie,ny_2d+ibc )%split1_modal( 1:ng ) = elem(  ie,ibc )%split1_modal( 1:ng )
        elem( ie, 1-ibc )%split2_modal( 1:ng ) = elem( ie,ny_2d+1 -ibc )%split2_modal( 1:ng )
        elem( ie,ny_2d+ibc )%split2_modal( 1:ng ) = elem(  ie,ibc )%split2_modal( 1:ng )
    enddo
    
    end subroutine bc_y