    !*************************************************************************
    !*************************************************************************
    module module_globals1d

    integer :: kkkk
    integer :: i,j
    integer :: nx,nghost
    integer :: nk ! nk is the number of degree of the polynomial of DG space!
    real :: dx,u10,u20,u1mean,u2mean,u210,uet,u1meanpost,u2meanpost
    real :: pi,eps
    real :: xg(6),wg(6)
    real :: ai(1:5)
    real :: xleft,xright
    real :: dt,time,time_final,cfl,time_begin,time_end
    real,allocatable :: x(:)
    real,allocatable :: xgrid(:)
    real,allocatable :: flux_LF(:),flux_LFp(:),flux_LF11(:),flux_LF21(:),flux_LF12(:),flux_LF22(:),flux_LF1(:,:),flux_LF2(:,:),flux_LF01(:,:),flux_LF02(:,:),flux_LF001(:,:),flux_LF002(:,:),flux_LF0001(:,:),flux_LF0002(:,:)
    real,allocatable :: A0(:,:),A1(:,:),A2(:,:)


    integer :: io,iio
    real :: er11,er12,er13,er21,er22,er23,er31,er32,er33,er41,er42,er43
 
    real,allocatable :: speed(:),speed1(:),speed2(:),u1post(:,:),u2post(:,:)
 
    real :: alpha !maximum wave speed!


    integer :: iexample
 
    real :: gau_lob(1:6,1:2),gau(1:6,1:2),gau_f_phi_x(1:3,1:2)
 
    integer :: kk

    !Runge-Kutta variables
    integer :: irk
    real :: rk_d(0:4)
    !*********************
    real :: dtt
    
    integer,allocatable :: itrouble(:)
    
    integer :: irkdg 
    
    real :: speed_max

    end module module_globals1d