    subroutine output_dg
    implicit none
    integer :: kk0
    integer :: ig
    real :: xrg

    open( 1,file='output.txt' )
    open( 2,file='exact.txt' )
    do kk0 = 1 , nx
        do ig = 1 , nk+1
            xrg =x(kk0) + dx * gau(ig,1)
            write(1,*) xrg,ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)
             
        enddo
    enddo
    close(1)
    close(2)

    end subroutine output_dg