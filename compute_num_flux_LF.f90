    !subroutine compute_num_flux_LF(u1l,u2l,u1r,u2r,xs,speed_r,alpha,flux_LF11)
    !implicit none
    !real,intent(in) :: u1l,u2l,u1r,u2r,xs,alpha,speed_r
    !real,intent(out) :: flux_LF11(2)
    !
    !real :: u1bar,u2bar,u1jump,u2jump
    !
    !u1bar=0.5*(u1l+u1r)
    !u2bar=0.5*(u2l+u2r)
    !u1jump=u1r-u1l
    !u2jump=u2r-u2l
    !
    !!flux_LF11(1)= A11(xs)*u1bar+A12(xs)*u2bar-0.5*speed_r*(RL_j(1,1)*u1r+RL_j(1,2)*u2r+RL_jp(1,1)*u1l+RL_jp(1,2)*u2l)-0.5*alpha*u1jump
    !!flux_LF11(2)= A21(xs)*u1bar+A22(xs)*u2bar-0.5*speed_r*(RL_j(2,1)*u1r+RL_j(2,2)*u2r+RL_jp(2,1)*u1l+RL_jp(2,2)*u2l)-0.5*alpha*u2jump
    !
    !flux_LF11(1)= A11(xs)*u1bar+A12(xs)*u2bar-speed_r*u1bar-0.5*alpha*u1jump
    !flux_LF11(2)= A21(xs)*u1bar+A22(xs)*u2bar-speed_r*u2bar-0.5*alpha*u2jump
    !!
    !!flux_LF11(1)= 0.
    !!flux_LF11(2)= 0.
    !
    !!flux_LF11 = ((RT_j(1,1)*A11(xs)+RT_j(1,2)*A21(xs))*R_j(1,1)+(RT_j(1,1)*A12(xs)+RT_j(1,2)*A22(xs))*R_j(2,1)-speed)*v1bar+((RT_j(1,1)*A11(xs)+RT_j(1,2)*A21(xs))*R_j(1,2)+(RT_j(1,1)*A12(xs)+RT_j(1,2)*A22(xs))*R_j(2,2))*v2bar-alpha/2.*v1jump
    !!
    !!v1bar2=0.5*(v1l2+v1r2)
    !!v2bar2=0.5*(v2l2+v2r2)
    !!v2jump2=v2r2-v2l2
    !!
    !!flux_LF12 = ((RT_j(2,1)*A11(xs1)+RT_j(2,2)*A21(xs1))*R_j(1,1)+(RT_j(2,1)*A12(xs1)+RT_j(2,2)*A22(xs1))*R_j(2,1))*v1bar2+((RT_j(2,1)*A11(xs1)+RT_j(2,2)*A21(xs1))*R_j(1,2)+(RT_j(2,1)*A12(xs1)+RT_j(2,2)*A22(xs1))*R_j(2,2)-speed1)*v2bar2-alpha/2.*v2jump2
    !
    !end subroutine compute_num_flux_LF