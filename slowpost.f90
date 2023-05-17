    subroutine slowpost1
    
    !integer,intent(in) :: q,n,k
    !real,intent(in) :: u(n,k+1),y(n+1),x(n*q)
    !real,intent(out) :: upost(nx*6,2)
    
    integer :: sampleCount,kernelWidth,xIndex,i1,j1,num,allloc(2*nx),m,p,ig
    real :: bb,kernelBreaks(3*nk+2),compactSupport(2),halfWidth,width,localCompactSupport(2),offset,x1(nx*6)
    real :: localKernelBreaks(3*nk+2),intVal_1,intVal_2,a,b,quadZ(6),u1hatVal(6),u2hatVal(6),kernelVal(6),breaks(2*nx),brk(nk+2),coe(nk+1,nk+1),coeff(2*nk+1)
    
    !取H=max（h_j）
    !call gauss(xg,wg)
    sampleCount =nx*6
    do i1=1,nx
        do j1=1,6
            x1(6*(i1-1)+j1)=x(i1) + dx * xg(j1)
        end do
    end do
    u1post(:,1)=x1
    u2post(:,1)=x1
    !convolution
    call getbsbrk( nk+1,brk )
    call getbscoe( nk+1,coe )
    call getKernelCoeffs( nk,coeff )
    kernelWidth = (3 *nk+1)    
    halfWidth = kernelWidth * 0.5
    !call dxx(n,y,dx)
    !bb=0
    !do i=1,n
    !    if (dx(i)>bb) then 
    !        bb=dx(i)
    !    end if
    !end do
    !width=bb
    width=dx
    compactSupport = (/-halfWidth, halfWidth/) * width
    do i1=1,3*nk+2
        kernelBreaks(i1)=(i1-halfWidth-1)*width
    end do
    
    do xIndex=1,sampleCount
    
    ! find all the breakpoints for integration
        offset = u1post(xIndex,1)
        localCompactSupport = compactSupport + offset
        localKernelBreaks = kernelBreaks + offset
        call findbreaks(localKernelBreaks,breaks,allloc,num)
    ! integrate product of kernel and function
        intVal_1 = 0
        intVal_2 = 0
        do i1 = 1,(num - 1)
            a = breaks( i1 )
            b = breaks( i1 + 1 )
            !quadZ=(xg+1)*0.5*(b-a)+a
            quadZ=0.5*(b+a)+xg*(b-a)
            call UhEval(allloc(i1),quadZ,u1hatVal,u2hatVal)
            call evalKernelPP(nk,brk,coe,width,offset,coeff,quadZ,kernelVal)
            do j1=1,6
                !intVal_1=intVal_1+(b-a)*0.5*wg(j1)*kernelVal(j1)*u1hatVal(j1)
                !intVal_2=intVal_2+(b-a)*0.5*wg(j1)*kernelVal(j1)*u2hatVal(j1)
                intVal_1=intVal_1+(b-a)*wg(j1)*kernelVal(j1)*u1hatVal(j1)
                intVal_2=intVal_2+(b-a)*wg(j1)*kernelVal(j1)*u2hatVal(j1)
            end do
        end do
        u1post(xIndex,2)=intVal_1/width
        u2post(xIndex,2)=intVal_2/width
    end do
    !
    
    !call gauss(xg,wg)
    !upost(:,1)=x
    !!convolution
    !call getbsbrk( k+1,brk )
    !call getbscoe( k+1,coe )
    !call getKernelCoeffs( k,coeff )
    !kernelWidth = (3 *k+1)    
    !halfWidth = kernelWidth * 0.5
    !call dxx(n,y,dx)
    !do m=1,n
    !    width=dx(m)
    !    compactSupport = (/-halfWidth, halfWidth/) * width
    !    do i=1,3*k+2
    !        kernelBreaks(i)=(i-halfWidth-1)*width
    !    end do
    !    do p=1,q
    !        xIndex=q*(m-1)+p
    !        ! find all the breakpoints for integration
    !        offset = upost(xIndex,1)
    !        localCompactSupport = compactSupport + offset
    !        localKernelBreaks = kernelBreaks + offset
    !        call findbreaks(n,k,localKernelBreaks,y,breaks,allloc,num)
    !    ! integrate product of kernel and function
    !        intVal = 0
    !        do i = 1,(num - 1)
    !            a = breaks( i )
    !            b = breaks( i + 1 )
    !            quadZ=(xg+1)*0.5*(b-a)+a
    !            call UhEval(n,k,allloc(i),u,y,quadZ,uhatVal)
    !            call evalKernelPP(k,brk,coe,width,offset,coeff,quadZ,kernelVal)
    !            do j=1,6
    !                intVal=intVal+(b-a)*0.5*wg(j)*kernelVal(j)*uhatVal(j)
    !            end do
    !        end do
    !        upost(xIndex,2)=intVal/width
    !    end do
    !end do
    !
    end subroutine slowpost1
    !subroutine slowpost2(q,n,k,u,y,x,upost)
    !
    !integer,intent(in) :: q,n,k
    !real,intent(in) :: u(n,k+1),y(n+1),x(n*q)
    !real,intent(out) :: upost(n*q,2)
    !
    !integer :: sampleCount,kernelWidth,xIndex,i,j,num,allloc(2*n),m,p
    !real :: xg(6),wg(6),bb,kernelBreaks(3*k+2),compactSupport(2),halfWidth,dx(n),width,localCompactSupport(2),offset
    !real :: localKernelBreaks(3*k+2),intVal,a,b,quadZ(6),uhatVal(6),kernelVal(6),breaks(2*n),brk(k+2),coe(k+1,k+1),coeff(2*k+1)
    !
    !!!取H=max（h_j）
    !!call gauss(xg,wg)
    !!sampleCount =n*q
    !!upost(:,1)=x
    !!!convolution
    !!call getbsbrk( k+1,brk )
    !!call getbscoe( k+1,coe )
    !!call getKernelCoeffs( k,coeff )
    !!kernelWidth = (3 *k+1)    
    !!halfWidth = kernelWidth * 0.5
    !!call dxx(n,y,dx)
    !!bb=0
    !!do i=1,n
    !!    if (dx(i)>bb) then 
    !!        bb=dx(i)
    !!    end if
    !!end do
    !!width=bb
    !!compactSupport = (/-halfWidth, halfWidth/) * width
    !!do i=1,3*k+2
    !!    kernelBreaks(i)=(i-halfWidth-1)*width
    !!end do
    !!do xIndex=1,sampleCount
    !!
    !!! find all the breakpoints for integration
    !!    offset = upost(xIndex,1)
    !!    localCompactSupport = compactSupport + offset
    !!    localKernelBreaks = kernelBreaks + offset
    !!    call findbreaks(n,k,localKernelBreaks,y,breaks,allloc,num)
    !!! integrate product of kernel and function
    !!    intVal = 0
    !!    do i = 1,(num - 1)
    !!        a = breaks( i )
    !!        b = breaks( i + 1 )
    !!        quadZ=(xg+1)*0.5*(b-a)+a
    !!        call UhEval(n,k,allloc(i),u,y,quadZ,uhatVal)
    !!        call evalKernelPP(k,brk,coe,width,offset,coeff,quadZ,kernelVal)
    !!        do j=1,6
    !!            intVal=intVal+(b-a)*0.5*wg(j)*kernelVal(j)*uhatVal(j)
    !!        end do
    !!    end do
    !!    upost(xIndex,2)=intVal/width
    !!end do
    !!
    !
    !call gauss(xg,wg)
    !upost(:,1)=x
    !!convolution
    !call getbsbrk( k+1,brk )
    !call getbscoe( k+1,coe )
    !call getKernelCoeffs( k,coeff )
    !kernelWidth = (3 *k+1)    
    !halfWidth = kernelWidth * 0.5
    !call dxx(n,y,dx)
    !do m=1,n
    !    width=dx(m)
    !    compactSupport = (/-halfWidth, halfWidth/) * width
    !    do i=1,3*k+2
    !        kernelBreaks(i)=(i-halfWidth-1)*width
    !    end do
    !    do p=1,q
    !        xIndex=q*(m-1)+p
    !        ! find all the breakpoints for integration
    !        offset = upost(xIndex,1)
    !        localCompactSupport = compactSupport + offset
    !        localKernelBreaks = kernelBreaks + offset
    !        call findbreaks(n,k,localKernelBreaks,y,breaks,allloc,num)
    !    ! integrate product of kernel and function
    !        intVal = 0
    !        do i = 1,(num - 1)
    !            a = breaks( i )
    !            b = breaks( i + 1 )
    !            quadZ=(xg+1)*0.5*(b-a)+a
    !            call UhEval(n,k,allloc(i),u,y,quadZ,uhatVal)
    !            call evalKernelPP(k,brk,coe,width,offset,coeff,quadZ,kernelVal)
    !            do j=1,6
    !                intVal=intVal+(b-a)*0.5*wg(j)*kernelVal(j)*uhatVal(j)
    !            end do
    !        end do
    !        upost(xIndex,2)=intVal/width
    !    end do
    !end do
    !
    !
    !end subroutine slowpost2
    
    
    
    subroutine findbreaks(localKernelBreaks,allBreaks,allloc,l)
    
    !integer,intent(in) :: n,k
    real,intent(in) :: localKernelBreaks(3*nk+2)
    integer,intent(out) :: l,allloc(2*nx)
    real,intent(out) :: allBreaks(2*nx)
    
    integer :: i1,j1
    real :: a,b,yy(-nx+1:2*nx+1),pip
    
    !pi=4.*atan(1.)!domian perid
    !pip=8.*atan(1.)
    pip=xright-xleft
    allBreaks=0.
    allloc=0.
    do i1=1,nx
        yy(i1-nx)=xgrid(i1)-pip
        yy(i1)=xgrid(i1)
        yy(i1+nx)=xgrid(i1)+pip
    end do
    yy(2*nx+1)=xgrid(nx+1)+pip
    i1=1
    j1=1-nx
    l=1
    allBreaks(l)=localKernelBreaks(i1)
    do while (i1<3*nk+2)
        if (localKernelBreaks(i1)<yy(j1) .and. localKernelBreaks(i1+1)>=yy(j1)) then
            allloc(l)=j1-1
            l=l+1
            allBreaks(l)=yy(j1)
            j1=j1+1
        elseif (localKernelBreaks(i1+1)<yy(j1)) then
            allloc(l)=j1-1
            i1=i1+1
            l=l+1
            allBreaks(l)=localKernelBreaks(i1)
        else
            j1=j1+1
        end if
    end do          
    
    end subroutine findbreaks
    
    subroutine UhEval( n0, quadZ ,u1hatVal,u2hatVal)
    
    integer,intent(in) :: n0
    real,intent(in) :: quadZ(6)
    real,intent(out) :: u1hatVal(6),u2hatVal(6)
    
    integer :: i1,j1,n1
    real :: a,b,pip
    
    !pi=4.*atan(1.)为周期
    !pip=8.*atan(1.)
    pip=xright-xleft
    do i1=1,6
        !a=0
        if (n0<1) then
            n1=n0+nx
            !b=(quadZ(i1)+pip-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5!4.*atan(1.)为周期
            a=ortho_poly1d(element(n1)%u1modal(1:nk+1,0),quadZ(i1)+pip,x(n1),dx,nk)
            b=ortho_poly1d(element(n1)%u2modal(1:nk+1,0),quadZ(i1)+pip,x(n1),dx,nk)
        else if (n0>nx) then
            n1=n0-nx
            !b=(quadZ(i1)-pip-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5
            a=ortho_poly1d(element(n1)%u1modal(1:nk+1,0),quadZ(i1)-pip,x(n1),dx,nk)
            b=ortho_poly1d(element(n1)%u2modal(1:nk+1,0),quadZ(i1)-pip,x(n1),dx,nk)
        else
            n1=n0
            !b=(quadZ(i1)-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5
            a=ortho_poly1d(element(n1)%u1modal(1:nk+1,0),quadZ(i1),x(n1),dx,nk)
            b=ortho_poly1d(element(n1)%u2modal(1:nk+1,0),quadZ(i1),x(n1),dx,nk)
        end if
        u1hatVal(i1)=a
        u2hatVal(i1)=b
        !if (n0<1) then
        !    n1=n0+nx
        !    b=(quadZ(i1)+pip-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5!4.*atan(1.)为周期
        !else if (n0>nx) then
        !    n1=n0-nx
        !    b=(quadZ(i1)-pip-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5
        !else
        !    n1=n0
        !    b=(quadZ(i1)-xgrid(n1))/(xgrid(n1+1)-xgrid(n1))-0.5
        !end if
        !do j1=1,nk+1
        !    a=a+Uh(n1,j1)*(b**(j1-1))
        !end do
        !uhatVal(i1)=a
    end do
    
    end subroutine UhEval
    
    subroutine evalKernelPP(k,brk,coe,width,offset,coeff,quadZ,kernelVal)
    
    integer,intent(in) :: k
    real, intent(in) :: brk(k+2),coe(k+1,k+1),width,offset,coeff(2*k+1),quadZ(6)
    real,intent(out) :: kernelVal(6)
    
    integer :: i1,j1,cindex
    real :: val_p,off1
    
    kernelVal=0.
    do i1=1,6
        do j1=-k,k
            cindex=j1+k+1
            off1=offset+j1*width
            call evalMappedPP(k,brk,coe,width,off1,quadZ(i1),val_p)
            kernelVal(i1)=kernelVal(i1)+coeff(cindex)*val_p
        end do
    end do
    
    end subroutine evalKernelPP
    
    subroutine evalMappedPP(k,brk,coe,width,off1,x1,val)
    
    integer,intent(in) :: k
    real,intent(in) :: brk(k+2),coe(k+1,k+1),width,off1,x1
    real,intent(out) :: val
    
    real minVal,maxVal,xx
    
    minVal = brk(1)
    maxVal = brk(k+2)

    xx = x1-off1
    xx = xx / width
    if (xx >= minVal .and. xx <= maxVal) then
        call ppval( k,brk,coe,xx,val )
    else
        val=0
    end if
    
    end subroutine evalMappedPP
    
    subroutine ppval(k,brk,coe,x11,val)
    
    integer,intent(in) :: k
    real,intent(in) :: brk(k+2),coe(k+1,k+1),x11
    real,intent(out) :: val
    
    integer :: i1,j1
    real :: x1
    
    do i1=1,k+1
        if ((x11<brk(i1+1)).and.(x11>=brk(i1))) then
            val=0
            x1=x11-brk(i1)
            do j1=1,k+1
                val=val+coe(j1,i1)*(x1**(k+1-j1))
            end do
        end if
    end do
    if (x11==brk(k+2)) then
        val=0
    end if
    
    end subroutine ppval
    
    subroutine getbsbrk(k,brk)
    
    integer,intent(in) :: k
    real,intent(out) :: brk(k+1)
    
    if (k==1) then
        brk=(/-0.5, 0.5/)
    elseif(k==2) then
        brk=(/-1, 0, 1/)
    elseif(k==3) then
        brk=(/-1.5,-0.5,0.5,1.5/)
    elseif(k==4) then
        brk=(/-2,-1,0,1,2/)
    elseif(k==5) then
        brk=(/-2.5,-1.5,-0.5,0.5,1.5,2.5/)
    elseif(k==6) then
        brk=(/-3,-2,-1,0,1,2,3/)
    elseif(k==7) then
        brk=(/-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5/)
    elseif(k==8) then
        brk=(/-4,-3,-2,-1,0,1,2,3,4/)
    elseif(k==9) then
        brk=(/-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5/)
    elseif(k==10) then
        brk=(/-5,-4,-3,-2,-1,0,1,2,3,4,5/)
    end if
       
    end subroutine getbsbrk
      
    subroutine getbscoe(k,coe)
    
    integer,intent(in) :: k
    real,intent(out) :: coe(k,k)
    
    if (k==1) then
        coe=1
    elseif(k==2) then
        coe(:,1)=(/1.,0./)
        coe(:,2)=(/-1.,1./)
    elseif(k==3) then
        coe(:,1)=(/0.5, 0., 0./)
        coe(:,2)=(/-1., 1., 0.5/)
        coe(:,3)=(/0.5, -1., 0.5/)
    elseif(k==4) then
        coe(:,1)=(/1./6.,0.,0.,0./)
        coe(:,2)=(/-0.5,0.5,0.5,1./6./)
        coe(:,3)=(/0.5,-1.,0.,2./3./)
        coe(:,4)=(/-1./6.,0.5,-0.5,1./6./)
    elseif(k==5) then
        coe(:,1)=(/1./24.,0.,0.,0.,0./)
        coe(:,2)=(/-1./6.,1./6.,0.25,1./6.,1./24./)
        coe(:,3)=(/0.25,-0.5,-0.25,0.5,11./24./)
        coe(:,4)=(/-1./6.,0.5,-0.25,-0.5,11./ 24./)
        coe(:,5)=(/1./24.,-1./6.,0.25,-1./6.,1./24./)
    elseif(k==6) then
        coe(:,1)=(/1./120.,0.,0.,0.,0.,0./)
        coe(:,2)=(/-1./24., 1./24., 1./12., 1./12., 1./24., 1./120./)
        coe(:,3)=(/1./12., -1./6., -1./6., 1./6., 5./12., 13./60./)
        coe(:,4)=(/-1./12., 0.25, 0., -1./2., 0., 11./20./)
        coe(:,5)=(/1./24., -1./6., 1./6., 1./6., -5./12., 13./60./)
        coe(:,6)=(/-1./120., 1./24., -1./12., 1./12., -1./24., 1./120./)
    elseif(k==7) then
        coe(:,1)=(/1./720., 0., 0., 0., 0., 0., 0./)
        coe(:,2)=(/-1./120., 1./120., 1./48., 1./36., 1./48., 1./120., 1./720./)
        coe(:,3)=(/1./48., -1./24., -1./16., 1./36., 3./16., 5./24., 19./240./)
        coe(:,4)=(/-1./36., 1./12., 1./24., -2./9., -5./24., 1./3., 151./360./)
        coe(:,5)=(/1./48., -1./12., 1./24., 2./9., -5./24., -1./3., 151./360./)
        coe(:,6)=(/-1./120., 1./24., -1./16., -1./36., 3./16., -5./24., 19./240./)
        coe(:,7)=(/1./720., -1./120., 1./48., -1./36., 1./48., -1./120., 1./720./)
    elseif(k==8) then
        coe(:,1)=(/1./5040, 0., 0., 0., 0., 0., 0., 0./)
        coe(:,2)=(/-1./720, 1./720, 1./240, 1./144, 1./144, 1./240, 1./720, 1./5040/)
        coe(:,3)=(/1./240, -1./120, -1./60, 0., 1./18, 1./10, 7./90, 1./42/)
        coe(:,4)=(/-1./144, 1./48, 1./48, -1./16, -19./144, 1./16, 49./144, 397./1680/)
        coe(:,5)=(/1./144, -1./36, 0., 1./9, 0., -1./3, 0., 151./315/)
        coe(:,6)=(/-1./240, 1./48, -1./48, -1./16, 19./144, 1./16, -49./144, 397./1680/)
        coe(:,7)=(/1./720, -1./120, 1./60, 0., -1./18, 1./10, -7./90, 1./42/)
        coe(:,8)=(/-1./5040, 1./720, -1./240, 1./144, -1./144, 1./240, -1./720, 1./5040/)
    elseif(k==9) then
        coe(:,1)=(/1./40320, 0., 0., 0., 0., 0., 0., 0., 0./)
        coe(:,2)=(/-1./5040, -1./5040, 1./1440, 1./720, 1./576, 1./720, 1./1440, 1./5040, 1./40320/)
        coe(:,3)=(/1./1440, -1./720, -1./288, -1./720, 7./576, 23./720, 11./288, 17./720, 247./40320/)
        coe(:,4)=(/-1./720, 1./240, 1./160, -1./80, -3./64, -1./80, 21./160, 17./80, 477./4480/)
        coe(:,5)=(/1./576, -1./144, -1./288, 5./144, 19./576, -19./144, -49./288, 35./144, 15619./40320/)
        coe(:,6)=(/-1./720, 1./144, -1./288, -5./144, 19./576, 19./144, -49./288, -35./144, 15619./40320/)
        coe(:,7)=(/1./1440, -1./240, 1./160, 1./80, -3./64, 1./80, 21./160, -17./80, 477./4480/)
        coe(:,8)=(/-1./5040, 1./720, -1./288, 1./720, 7./576, -23./720, 11./288, -17./720, 247./40320/)
        coe(:,9)=(/1./40320, -1./5040, 1./1440, -1./720, 1./576, -1./720, 1./1440, -1./5040, 1./40320/)
    elseif(k==10) then
        coe(:,1)=(/1./362880, 0., 0., 0., 0., 0., 0., 0., 0.,0./)
        coe(:,2)=(/-1./40320, 1./40320, 1./10080, 1./4320, 1./2880, 1./2880, 1./4320, 1./10080, 1./40320, 1./362880/)
        coe(:,3)=(/1./10080, -1./5040, -1./1680, -1./2160, 1./480, 11./1440, 1./80, 59./5040, 41./6720, 251./181440/)
        coe(:,4)=(/-1./4320, 1./1440, 1./720, -1./540, -17./1440, -1./90, 67./2160, 17./180, 289./2880, 913./22680/)
        coe(:,5)=(/1./2880, -1./720, -1./720, 17./2160, 23./1440, -43./1440, -217./2160, 11./720, 809./2880, 44117./181440/)
        coe(:,6)=(/-1./2880, 1./576, 0., -5./432, 0., 19./288, 0., -35./144, 0., 15619./36288/)
        coe(:,7)=(/1./4320, -1./720, 1./720, 17./2160, -23./1440, -43./1440, 217./2160, 11./720, -809./2880, 44117./181440/)
        coe(:,8)=(/-1./10080, 1./1440, -1./720, -1./540, 17./1440, -1./90, -67./2160, 17./180, -289./2880, 913./2268/)
        coe(:,9)=(/1./40320, -1./5040, 1./1680, -1./2160, -1./480, 11./1440, -1./80, 59./5040, -41./6720, 251./181440/)
        coe(:,10)=(/-1./362880, 1./40320, -1./10080, 1./4320, -1./2880, 1./2880, -1./4320, 1./10080, -1./40320, 1./362880/)
    end if
     
    end subroutine getbscoe
    
    subroutine getKernelCoeffs( k,coeff )
    
    integer, intent(in) :: k
    real, intent(out) :: coeff(2*k+1)
    
    if (k==1) then
        coeff = (/-1./12, 7./6, -1./12/)
    elseif (k==2) then
            coeff = (/37./1920, -97./480, 437./320, -97./480, 37./1920/)
    elseif (k==3) then
            coeff = (/-41./7560, 311./5040, -919./2520, 12223./7560, -919./2520, 311./5040, -41./7560/)
    elseif (k==4) then
            coeff = (/153617./92897280, -35411./1658880, 3153959./23224320, -6803459./11612160, 18017975./9289728, -6803459./11612160, 3153959./23224320, -35411./1658880, 153617./92897280/)
    elseif (k==5) then
            coeff = (/-4201./7983360, 30773./3991680, -20813./380160, 2825./11088, -1179649./1330560, 1569217./665280, -1179649./1330560, 2825./11088, -20813./380160, 30773./3991680, -4201./7983360/)
    elseif (k==6) then
            coeff = (/13154671847./76517631590400, -18073154507./6376469299200, 287360344573./12752938598400, -2217732343517./19129407897600, 1240941746699./2833986355200, -275386671493./212548976640, 2648644782397./910924185600, -275386671493./212548976640, 1240941746699./2833986355200, -2217732343517./19129407897600, 287360344573./12752938598400, -18073154507./6376469299200, 13154671847./76517631590400/)
    elseif (k==7) then
            coeff = (/-800993./14010796800, 73587167./70053984000, -651305719./70053984000, 3714581677./70053984000, -3085236289./14010796800, 1426328231./2001542400, -43268401973./23351328000, 42401344373./11675664000, -43268401973./23351328000, 1426328231./2001542400, -3085236289./14010796800, 3714581677./70053984000, -651305719./70053984000, 73587167./70053984000, -800993./14010796800/)
    end if
    
    end subroutine getKernelCoeffs