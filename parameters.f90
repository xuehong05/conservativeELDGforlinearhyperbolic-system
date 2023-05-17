    !*************************************************************************
    !*************************************************************************
    subroutine parameters
    implicit none

    pi = 4.*atan(1.)
    eps = epsilon(0.)*10
    
    !***********************
    ! below is the six points Gaussian quadrature on [-1/2,1/2].
    ! where xg are nodes,
    !       wg are the cooresponding weights.
    xg(1)=-0.466234757101576013906150777246997304567e0
    xg(2)=-0.330604693233132256830699797509952673503e0
    xg(3)=-0.119309593041598454315250860840355967709e0
    xg(4)=-xg(3)
    xg(5)=-xg(2)
    xg(6)=-xg(1)
    wg(1)=1.71324492379170345040296142172733e-1/2e0
    wg(2)=3.60761573048138607569833513837716e-1/2e0
    wg(3)=4.67913934572691047389870343989551e-1/2e0
    wg(4)=wg(3)
    wg(5)=wg(2)
    wg(6)=wg(1)

    !****************************
    ai(1) = 1.
    ai(2) = 12.
    ai(3) = 180.
    ai(4) = 2800.
    ai(5) = 44100.

    end subroutine parameters