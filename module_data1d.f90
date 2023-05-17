    !*************************************************************************
    !*************************************************************************
    module module_data1d
    ! data structure
    !*******************************
    type, public :: point1d
        sequence
        real :: coor
    end type

    type(point1d),allocatable,target,public :: vertex(:)
    !*******************************

    !*******************************
    type, public :: point_upstream
        sequence
        real :: coor
        integer :: id
    end type

    type(point_upstream),allocatable,target,public :: vertex_star(:)
    type(point_upstream),allocatable,target,public :: vertex_star1(:)
    type(point_upstream),allocatable,target,public :: vertex_star2(:)
    !*******************************

    type, public :: type_segment
        sequence
        type(point_upstream) :: porigin,pend
        integer :: id
    end type

    !*******************************
    type, public :: element1d_upstream
        sequence
        type(point_upstream) :: point_origin,point_end
        type(point_upstream) :: point_inter(0:11)
        type(type_segment) :: segment(10)
        integer :: nsub
        
        real :: xgl_star(1:5)
        real :: aa(1:5,1:3)
 
    end type

    type(element1d_upstream),allocatable,target,public :: element_star(:)
    type(element1d_upstream),allocatable,target,public :: element_star1(:,:,:)
    type(element1d_upstream),allocatable,target,public :: element_star2(:,:,:)
    !*******************************
    type,public :: element1d_Eulerian
        sequence
        real :: u1modal(1:6,0:4)
        real :: u2modal(1:6,0:4)
        real :: xl,xr
        
        real :: xgl(1:5)
        
    end type
    type(element1d_Eulerian),allocatable,target,public :: element(:)
    !*******************************
    type,public :: element1d_dynamic
        sequence
        
        real :: u1modal(1:6)
        real :: u2modal(1:6)
        real :: sum_int(6,1:2)
        real :: flux1(6,1:2)
        real :: flux0(6,1:2)
        real :: phi_int(6,1:2)
        real :: vector_int(6,1:2)
        real :: source_int(6,1:2)
        real :: aa(1:5,1:3)
        real :: xl,xr,xc
        real :: dx
        real :: midtime
        
        real :: massmetrix(1:5,1:5)
        
    end type
    
    type(element1d_dynamic),allocatable,target,public :: element_t(:,:)
    type(element1d_dynamic),allocatable,target,public :: element_t1(:,:,:)
    type(element1d_dynamic),allocatable,target,public :: element_t2(:,:,:)

    end module module_data1d