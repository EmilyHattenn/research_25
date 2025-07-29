module usefulsubs
    use types
    use fftw
    use constants
    use helpers

    real(dp), private :: small = 1d-200
    real(dp), private :: kB    = 1.38064852d-23
    
    integer,  private, allocatable :: kmkp(:,:)	   	     ! ky - kyp index matrix
    integer,  private, allocatable :: kpkp(:,:)	   	     ! ky + kyp index matrix
    integer,  private, allocatable :: kmq(:,:)	   	     ! ky + kyp index matrix
    integer,  private, allocatable :: kpq(:,:)	   	     ! ky + kyp index matrix
    

    ! interface: CalcN2
    ! Provides an interface for calculating the nonlinear index of refraction.
    interface ReadIt
        module procedure ReadIt1D, ReadIt2D
    end interface ReadIt
    
    interface WriteIT
        module procedure WriteIt1D, WriteIt2D
    end interface WriteIt

    interface dfdy
        module procedure dfdy1D, dfdy2D
    end interface dfdy

    interface dfdx
        module procedure dfdx1D, dfdx2D
    end interface dfdx  

    interface dfdy_q
        module procedure dfdy1D_q, dfdy2D_q
    end interface dfdy_q

    interface dfdx_q
        module procedure dfdx1D_q, dfdx2D_q
    end interface dfdx_q

    interface GFFT
        module procedure GFFT_1D, GFFT_2D
    end interface GFFT

    interface GIFFT
        module procedure GIFFT_1D, GIFFT_2D
    end interface GIFFT   

    interface fflip
        module procedure fflip_dp, fflip_dpc
    end interface fflip 
contains


    function fflip_dp(f)
        real(dp) :: f(:)
        real(dp) :: fflip_dp(size(f))
        integer  :: N, i
        
        N = size(f)
        
        do i=1, N
            fflip_dp(i) = f(N+1-i)
        end do
    end function
    
    
    function fflip_dpc(f)
        complex(dp) :: f(:)
        complex(dp) :: fflip_dpc(size(f))
        integer     :: N, i
        
        N = size(f)
        
        do i=1, N
            fflip_dpc(i) = f(N+1-i)
        end do
    end function    


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate space of input array f(y) by FFT:
    ! df(z)/dz = IFFT( - ii q FFT(f(z)) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dfdy1D(f, qy)
        complex(dp) :: f(:)
        real(dp)    :: qy(:)
        complex*16 :: dfdy1D(size(f))

        if(size(qy)==1) then
            dfdy1D = 0d0
            return
        endif        
    
        dfdy1D = f
        call FFT(dfdy1D)
        dfdy1D = dfdy1D * (ii*qy)
        call IFFT(dfdy1D)

    end function


    function dfdy2D(f, qy)
        complex(dp) :: f(:,:)
        real(dp)    :: qy(:)
        complex(dp) :: dfdy2D(size(f,1), size(f,2))
        integer     :: i, Ny
        Nx = size(f,1)

        if(size(qy)==1) then
            dfdy2D = 0d0
            return
        endif
        
        !!$omp parallel do private(i)  
        do i=1, size(f,1)
            dfdy2D(i,:) = dfdy1D(f(i,:), qy)
        end do
        !!$omp end parallel do


        !dfdy2D = f
        !call FFT(dfdy2D)
        !forall(i=1:Nx)  dfdy2D(i,:) = dfdy2D(i,:) * (ii*qy(:))
        !call IFFT(dfdy2D)
        
    end function


    function dfdx1D(f, qx)
        complex(dp) :: f(:)
        real(dp)    :: qx(:)
        complex(dp) :: dfdx1D(size(f))
        integer     :: i, Nx

        if(size(qx)==1) then
            dfdx1D = 0d0
            return
        endif
        
        dfdx1D = f
        call FFT(dfdx1D)
        dfdx1D(:) = dfdx1D(:) * (ii*qx(:))
        call IFFT(dfdx1D)
        
    end function
    

    function dfdx2D(f, qx)
        complex(dp) :: f(:,:)
        real(dp)    :: qx(:)
        complex(dp) :: dfdx2D(size(f,1), size(f,2))
        integer     :: i, Nx
        Ny = size(f,2)

        if(size(qx)==1) then
            dfdx2D = 0d0
            return
        endif

        !!$omp parallel do private(i)  
        do i=1, size(f,2)
            dfdx2D(:,i) = dfdy1D(f(:,i), qx)
        end do
        !!$omp end parallel do
        
        !dfdx2D = f
        !call FFT(dfdx2D)
        !forall(i=1:Ny)  dfdx2D(:,i) = dfdx2D(:,i) * (ii*qx(:))
        !call IFFT(dfdx2D)
        
    end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate space of input array f(y) by FFT:
    ! df(z)/dz = IFFT( - ii q FFT(f(z)) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dfdy1D_q(f, qy)
        complex(dp) :: f(:)
        real(dp)    :: qy(:)
        complex*16 :: dfdy1D_q(size(f))

        if(size(qy)==1) then
            dfdy1D_q = 0d0
            return
        endif        
    
        dfdy1D_q = f * (ii*qy)
    end function


    function dfdy2D_q(f, qy)
        complex(dp) :: f(:,:)
        real(dp)    :: qy(:)
        complex(dp) :: dfdy2D_q(size(f,1), size(f,2))
        integer     :: i, Ny
        Nx = size(f,1)

        if(size(qy)==1) then
            dfdy2D_q = 0d0
            return
        endif
        
        !$omp parallel do private(i,j)  
        do j=1, size(f,2)
            do i=1, size(f,1)
                dfdy2D_q(i,j) = f(i,j) * (ii*qy(j))
            end do
        end do
        !$omp end parallel do        
    end function


    function dfdx1D_q(f, qx)
        complex(dp) :: f(:)
        real(dp)    :: qx(:)
        complex(dp) :: dfdx1D_q(size(f))
        integer     :: i, Nx

        if(size(qx)==1) then
            dfdx1D_q = 0d0
            return
        endif

        dfdx1D_q(:) = f(:) * (ii*qx(:))
    end function
    

    function dfdx2D_q(f, qx)
        complex(dp) :: f(:,:)
        real(dp)    :: qx(:)
        complex(dp) :: dfdx2D_q(size(f,1), size(f,2))
        integer     :: i, Nx
        Ny = size(f,2)

        if(size(qx)==1) then
            dfdx2D_q = 0d0
            return
        endif

        !$omp parallel do private(i,j)  
        do j=1, size(f,2)
            do i=1, size(f,1)
                dfdx2D_q(i,j) = f(i,j) * (ii*qx(i))
            end do
        end do
        !$omp end parallel do  
        
    end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine GFFT_1D(f, dx)
        complex(dp), intent(inout) :: f(:)
        real(dp)   , intent(in   ) :: dx

        call fft(f)
        f = f * dx / sqrt(twopi)
    end subroutine

    subroutine GIFFT_1D(f, dq)
        complex(dp), intent(inout) :: f(:)
        real(dp)   , intent(in   ) :: dq

        call ifft(f)
        f = f * dq / sqrt(twopi) * size(f)
    end subroutine

    subroutine GFFT_2D(f, dx, dy)
        complex(dp), intent(inout) :: f(:,:)
        real(dp)   , intent(in   ) :: dx, dy

        call fft(f)
        f = f * dx * dy / twopi
    end subroutine

    subroutine GIFFT_2D(f, dqx, dqy)
        complex(dp), intent(inout) :: f(:,:)
        real(dp)   , intent(in   ) :: dqx, dqy

        call ifft(f)
        f = f * dqx * dqy / twopi * size(f,1) * size(f,2)
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine test2Dfrom1D(x,y)
        real(dp)    :: x(:), y(:)
        integer     :: i, j
        complex(dp) :: gx(size(x)), gy(size(y)), gxy(size(x),size(y))
        real(dp)    :: dx, dy
        logical     :: xyonly = .true. 

        dx  = x(2) - x(1)
        dy  = y(2) - y(1)
        gx  = 0d0
        gy  = 0d0
        gxy = 0d0

        gy  = exp(-((y+8d-6)/2d-6)**12)
        gx  = exp(-((x+1d-6)/4d-7)**2) * exp(-((x+1d-6)/8d-7)**12)

        if(xyonly) then
            do j=1, size(y)
                do i=1, size(x)
                    gxy(i,j) =  gx(i) * gy(j)
                end do
            end do
        else
            call fft(gy)
            call fft(gx)
            do j=1, size(y)
                do i=1, size(x)
                    gxy(i,j) =  gx(i) * gy(j)
                end do
            end do
            call ifft(gxy)
        endif

        print*, sum(abs(gxy)**2) * dx * dy
        stop
    end subroutine


    subroutine testconv(y)
        real(dp) :: y(:)
        integer  :: i, j
        complex(dp) :: gate(size(y)), func(size(y)), conv(size(y))

        gate = 0d0
        func = 0d0

        gate = exp(-((y+8d-6)/2d-6)**12)
        func = exp(-((y+8d-6)/2.33d-6)**2)
        conv = 0d0


        call fft(func)
        call fft(gate)
        call nyquist_1D(func)
        call nyquist_1D(gate)
        func = cshift(func, Ny/2)
        gate = cshift(gate, Ny/2)
        conv = convolve(func, gate) * dyy
        conv = cshift(conv, -Ny/2)
        call nyquist_1D(conv)
        call ifft(conv)
 


        do j=1, size(y)
            print*, y(j), abs(func(j)), aimag(func(j))
            !print*, y(j), abs(gate(j)), aimag(gate(j))
            !print*, y(j), real(conv(j)), aimag(conv(j))
        end do

        stop
    end subroutine


    subroutine ApplyABC(Field, abc)
        complex(dp), intent(inout) :: Field(:,:)
        real(dp),    intent(in   ) :: abc(:,:)

        call ifft(Field)
        Field = Field * abc
        call fft(Field)
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    
!     function dgdx(g, dx)
!         complex(dp) :: g(:,:)
!         real(dp)    :: dx
!         complex(dp) :: dgdx(size(f,1), size(f,2))
!     
!         dgdx = (cshift(f,1,1) - cshift(f,-1,1)) / 2d0 / dx
!     end function
! 
! 
!     function dgdy(g, dy)
!         complex(dp) :: g(:,:)
!         real(dp)    :: dy
!         complex(dp) :: dgdx(size(f,1), size(f,2))
!     
!         dgdx = (cshift(f,1,2) - cshift(f,-1,2)) / 2d0 / dy
!     end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dEdx(E, dx)
        complex(dp) :: E(:,:)
        real(dp)    :: dx
        complex(dp) :: dEdx(size(E,1), size(E,2))
        integer     :: i,j
        integer     :: iu(size(E,1))
        
        dEdx = 0d0
        if(size(E,1)==1) return
        
        forall(i=1:size(E,1)) iu(i) = i
        iu = cshift(iu,+1)

        !$omp parallel do private(i,j) 
        do j=1, size(E,2)
            do i=1, size(E,1)
                dEdx(i,j) = (E(iu(i),j) - E(i,j)) / dx
            end do
        end do
        !$omp end parallel do
    end function


    function dEdy(E, dy)
        complex(dp) :: E(:,:)
        real(dp)    :: dy
        complex(dp) :: dEdy(size(E,1), size(E,2))
        integer     :: i,j
        integer     :: ju(size(E,2))
        
        dEdy = 0d0
        if(size(E,2)==1) return
        
        
        forall(j=1:size(E,2)) ju(j) = j
        ju = cshift(ju,+1)

        !$omp parallel do private(i,j) 
        do j=1, size(E,2)
            do i=1, size(E,1)
                dEdy(i,j) = (E(i,ju(j)) - E(i,j)) / dy
            end do
        end do
        !$omp end parallel do
    end function


    function dHdx(H, dx)
        complex(dp) :: H(:,:)
        real(dp)    :: dx
        complex(dp) :: dHdx(size(H,1), size(H,2))
        integer     :: i,j
        integer     :: id(size(H,1))
        
        dHdx = 0d0
        if(size(H,1)==1) return
        
        forall(i=1:size(H,1)) id(i) = i
        id = cshift(id,-1)

        !$omp parallel do private(i,j) 
        do j=1, size(H,2)
            do i=1, size(H,1)
                dHdx(i,j) = (H(i,j) - H(id(i),j)) / dx
            end do
        end do
        !$omp end parallel do
    end function



    function dHdy(H, dy)
        complex(dp) :: H(:,:)
        real(dp)    :: dy
        complex(dp) :: dHdy(size(H,1), size(H,2))
        integer     :: i,j
        integer     :: jd(size(H,2))
        
        dHdy = 0d0
        if(size(H,2)==1) return
        
        forall(j=1:size(H,2)) jd(j) = j
        jd = cshift(jd,-1)
        
        !$omp parallel do private(i,j) 
        do j=1, size(H,2)
            do i=1, size(H,1)
                dHdy(i,j) = (H(i,j) - H(i,jd(j))) / dy
            end do
        end do
        !$omp end parallel do
    end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    subroutine testing_conv(Ex, y, q)
    	complex(dp), intent(inout) :: Ex(:)
    	real(dp),    intent(in   ) :: y(:), q(:)
    	real(dp)				   :: Iy, Iq, Iy2
    	integer					   :: i
    	
    	Ex = exp(-(y/yw)**2)
    	
    	Iy = sum(abs(Ex)**2 * dy)
    	    	    	    	
		call FFTc(Ex)
		Ex = Ex * dy
		Iy = sum(abs(Ex)**2 * dq)
		Ex = Ex * Ex
		call IFFTc(Ex)
		Ex = Ex / dy / Iy
    	
    	do i=1, Ny
    		print*, real(Ex(i))
    	end do
		!print*, dy, dely/twopi, dq, delq

    	
    	Iy = sum(real(Ex * conjg(Ex))) * dy
    	    	    	
    	call IFFTc(Ex)
		Ex = Ex * dy * Ny
    	    	
    	Iq = sum(real(Ex * conjg(Ex))) * dq
    	    	
    	call FFTc(Ex)
    	Ex = Ex / dy / Ny
    	
    	Iy2 = sum(real(Ex * conjg(Ex))) * dy
    	
    	print*, Iy, Iq, Iy2
    	    	
    end subroutine


    subroutine testing_fftc(Ex, y, q)
    	complex(dp), intent(inout) :: Ex(:)
    	real(dp),    intent(in   ) :: y(:), q(:)
    	real(dp)				   :: Iy, Iq, Iy2
    	integer					   :: i
    	
    	Ex = exp(-(y/10d-6)**2)
    	
    	Iy = sum(abs(Ex)**2 * dy)
    	    	    	    	
		call FFTc(Ex)
		Ex = Ex / size(Ex)
		
		!call IFFTc(Ex)
		!Ex = Ex * size(Ex)


    	do i=1, size(ex)
    		print*, real(Ex(i))
    	end do

    end subroutine
    
    
    subroutine print2file(x,y,u,filename)
        real(dp), intent(in) :: x(:), y(:)
        integer,  intent(in) :: u
        character(len=*), intent(in) :: filename
        
        open(unit=u,file=filename)
        
        do i=1, size(x)
            write(u,*) x(i) , y(i), i
        end do
        
        close(u)
    end subroutine


    subroutine FT(y,x,q)
        complex(dp), intent(inout) :: y(:)
        real(dp),    intent(in   ) :: x(:), q(:)
        complex(dp)                :: ytmp(size(x))
        integer                    :: N, j

        N = size(x)
        
        forall(j=1:N) ytmp(j) = sum(y * exp(ii*x*q(j)))
        y = ytmp / N
        
!        forall(j=1:N) y(j) = y(j) * exp(-(q(j)/q(5*N/8))**8)
    end subroutine


    subroutine IFT(y,x,q)
        complex(dp), intent(inout) :: y(:)
        real(dp),    intent(in   ) :: x(:), q(:)
        complex(dp)                :: ytmp(size(x))
        integer                    :: N, j

        N = size(x)
                
        forall(j=1:N) ytmp(j) = sum(y * exp(-ii*x(j)*q))  
        
        !forall(j=1:N) y(j) = y(j) * (-1)**(j+1)
        
        y = ytmp
    end subroutine



    function Flip(x)
        complex(dp) :: x(:)
        complex(dp) :: Flip(size(x))
        integer     :: i, N

        Flip = 0d0
        N    = size(x)
        
        do i=0, N-1
            Flip(1+i) = x(N-i)
        end do
    end function Flip
    
    

    function GetArray0Index(x)
        real(dp) :: x(:)
        integer  :: GetArray0Index
        integer  :: xi(size(x))
        integer  :: i
        real(dp) :: dx

        dx = x(3)-x(2)
        xi = nint((x+0.1*dx)/dx)

        GetArray0Index = 0
        
        do i=1, size(x)
            GetArray0Index = i
            if(xi(i) == 0) return
        end do

        print*, "Error in GetArray0Index in usefulsubs.f90: N0 = ", GetArray0Index
        stop
    end function GetArray0Index
    
    elemental Function GaussDelta(a,b)
        real(dp), intent(in) :: a, b
        real(dp)             :: GaussDelta
        GaussDelta = 0d0
        GaussDelta = 1d0 / sqrt(pi) / b * Exp(-(a/b)**2)
    end Function


    ! The Dirac Delta function
    function delta(x)
        real(dp) :: x, delta
        delta = 0d0
        if(nint(x/dky) == 0) delta = 1
    end function

    ! The Kronecker Delta function
    function kdel(x)
        integer :: x, kdel
        kdel = 0
        if(x == 0) kdel = 1
    end function


    function delt(x)
        integer  :: x
        real(dp) :: delt
        delt = 0d0
        delt = 1d0 - abs(x) / (abs(x)+1d-100)
    end function

    integer function sgn(x)
        double precision :: x
        if(x < 0d0) then
            sgn = -1
        else
            sgn = +1
        endif
    end function

    function sgn2(x)
        double precision :: x(:)
        integer			 :: sgn2(size(x))
        integer          :: i
        do i=1, size(x)
            sgn2(i) = sgn(x(i))
        end do
    end function
    
    
    function TotalEnergy(n, E) result(Eng)
        complex*16	 :: n(:)
        double precision :: E(:)
        double precision :: Eng
        Eng = real(sum(n*E))
    end function


    function AvgEnergy(n, E) result(Eng)
        complex*16       :: n(:)
        double precision :: E(:)
        double precision :: Eng
        Eng = real(sum(n*E) / (sum(n)+small))
    end function


    function Temperature(n, E) result(temp)
        complex*16       :: n(:)
        double precision :: E(:)
        double precision :: temp
        temp = 2 * AvgEnergy(n, E) / kB
    end function Temperature


    ! The Lorentzian function
    elemental function Lrtz(a,b)
        real(dp), intent(in) :: a, b
        real(dp)             :: Lrtz

        Lrtz = (b/pi) / (a**2 + b**2)
    end function Lrtz



    real(dp) elemental function theta(x)
        real(dp), intent(in) :: x
        
        theta = (abs(x) + x) / 2d0 / (abs(x)+small)
    end function

    
    
    real(dp) elemental function softtheta(x, g)
        real(dp), intent(in) :: x, g
        
        softtheta = 0.5 * (1.0 + 2.0 / pi * atan(x/g))
    end function
    
    
    elemental function rad(degrees)
        double precision, intent(in) :: degrees
        double precision             :: rad
        rad = degrees * pi / 180
    end function
    
    
    subroutine RotateField(theta, Ex, Ey)
        real(dp),    intent(in   ) :: theta
        complex(dp), intent(inout) :: Ex(:,:), Ey(:,:)
        complex(dp)                :: Ex0(size(Ex,1),size(Ex,2))
        complex(dp)                :: Ey0(size(Ey,1),size(Ey,2))
        real(dp)                   :: R11, R22, R33, R44
        integer                    :: i, j
        
        Ex0 = Ex
        Ey0 = Ey
        
        R11 = cos(theta)
        R12 = - sin(theta)
        R21 = + sin(theta)
        R22 = cos(theta)
        
        !$omp parallel do private(j, i)
        do j=1, size(Ex,2)
            do i =1, size(Ex,1)
                Ex(i,j) = Ex0(i,j)*R11 + Ey0(i,j)*R12
                Ey(i,j) = Ex0(i,j)*R21 + Ey0(i,j)*R22
            end do
        end do
        !$omp end parallel do
    end subroutine
    
    
    subroutine ShiftField(Lx, Ly, dx, dy, Ex, Ey)
        real(dp),    intent(in   ) :: Lx, Ly
        real(dp),    intent(in   ) :: dx, dy
        complex(dp), intent(inout) :: Ex(:,:), Ey(:,:)
        integer                    :: nx, ny
 
        nx = nint(Lx / dx)
        ny = nint(Ly / dy)
        
        Ex = Cshift(Ex,nx,1)
        Ex = Cshift(Ex,ny,2)
        
        Ey = Cshift(Ey,nx,1)
        Ey = Cshift(Ey,ny,2)
    end subroutine
    
    
    subroutine RotateShiftEField(theta, qx, qy, Ex, Ey)
        real(dp),    intent(in   ) :: theta, qx(:), qy(:)
        complex(dp), intent(inout) :: Ex(:,:), Ey(:,:)
        real(dp)                   :: dqx, dqy
        integer                    :: i, j, nqx, nqy
        
        dqx = qx(2) - qx(1)
        dqy = qy(2) - qy(1)

        
        call FFTC(Ex)
        call FFTC(Ey)
        
        dumb = sum(abs(Ex)**2 + abs(Ey)**2)
        
        do j=1, size(Ex,2)
            do i=1, size(Ex,1)
                qx0 = qx0 + qx(i) * (abs(Ex(i,j))**2 + abs(Ey(i,j))**2) / dumb
                qy0 = qy0 + qy(i) * (abs(Ex(i,j))**2 + abs(Ey(i,j))**2) / dumb
            end do
        end do
        
        nqx = nint(q0x * sin(theta) / dqx)
        nqy = nint(q0y * cos(theta) / dqy)
        
        Ex = Cshift(Ex,nqx,1)
        Ex = Cshift(Ex,nqy,2)
        
        Ey = Cshift(Ey,nqx,1)
        Ey = Cshift(Ey,nqy,2)

        call IFFTC(Ex)
        call IFFTC(Ey)

    end subroutine
    

subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

!*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z), 
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, 
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), 
!    and K1'(z).
!
  implicit none

  real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
    0.125D+00,           7.03125D-02,&
    7.32421875D-02,      1.1215209960938D-01,&
    2.2710800170898D-01, 5.7250142097473D-01,&
    1.7277275025845D+00, 6.0740420012735D+00,&
    2.4380529699556D+01, 1.1001714026925D+02,&
    5.5133589612202D+02, 3.0380905109224D+03 /)
  real ( kind = 8 ) a0
  real ( kind = 8 ), save, dimension ( 10 ) :: a1 = (/ &
    0.125D+00,            0.2109375D+00, &
    1.0986328125D+00,     1.1775970458984D+01, &
    2.1461706161499D+002, 5.9511522710323D+03, &
    2.3347645606175D+05,  1.2312234987631D+07, &
    8.401390346421D+08,   7.2031420482627D+10 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
   -0.375D+00,           -1.171875D-01, &
   -1.025390625D-01,     -1.4419555664063D-01, &
   -2.7757644653320D-01, -6.7659258842468D-01, &
   -1.9935317337513D+00, -6.8839142681099D+00, &
   -2.7248827311269D+01, -1.2159789187654D+02, &
   -6.0384407670507D+02, -3.3022722944809D+03 /)
  complex ( kind = 8 ) ca
  complex ( kind = 8 ) cb
  complex ( kind = 8 ) cbi0
  complex ( kind = 8 ) cbi1
  complex ( kind = 8 ) cbk0
  complex ( kind = 8 ) cbk1
  complex ( kind = 8 ) cdi0
  complex ( kind = 8 ) cdi1
  complex ( kind = 8 ) cdk0
  complex ( kind = 8 ) cdk1
  complex ( kind = 8 ) ci
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) cs
  complex ( kind = 8 ) ct
  complex ( kind = 8 ) cw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  !real ( kind = 8 ) pi
  real ( kind = 8 ) w0
  complex ( kind = 8 ) z
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2
  complex ( kind = 8 ) zr
  complex ( kind = 8 ) zr2

  !pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0D+00 ) then
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cbi1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
    cbk0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cbk1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    z1 = -z
  end if

  if ( a0 <= 18.0D+00 ) then

    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * k )
      cbi0 = cbi0 + cr
      if ( abs ( cr / cbi0 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * ( k + 1 ) )
      cbi1 = cbi1 + cr
      if ( abs ( cr / cbi1 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = 0.5D+00 * z1 * cbi1

  else

    if ( a0 < 35.0D+00 ) then
      k0 = 12
    else if ( a0 < 50.0D+00 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    zr = 1.0D+00 / z1
    do k = 1, k0
      cbi0 = cbi0 + a(k) * zr ** k
    end do
    cbi0 = ca * cbi0
    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, k0
      cbi1 = cbi1 + b(k) * zr ** k
    end do
    cbi1 = ca * cbi1

  end if

  if ( a0 <= 9.0D+00 ) then

    cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
    w0 = 0.0D+00
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      w0 = w0 + 1.0D+00 / k
      cr = 0.25D+00 * cr / ( k * k ) * z2
      cs = cs + cr * ( w0 + ct )
      if ( abs ( ( cs - cw ) / cs ) < 1.0D-15 ) then
        exit
      end if
      cw = cs
    end do

    cbk0 = ct + cs

  else

    cb = 0.5D+00 / z1
    zr2 = 1.0D+00 / z2
    cbk0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 10
      cbk0 = cbk0 + a1(k) * zr2 ** k
    end do
    cbk0 = cb * cbk0 / cbi0

  end if

  cbk1 = ( 1.0D+00 / z1 - cbi1 * cbk0 ) / cbi0

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then

    if ( imag ( z ) < 0.0D+00 ) then
      cbk0 = cbk0 + ci * pi * cbi0
      cbk1 = - cbk1 + ci * pi * cbi1
    else
      cbk0 = cbk0 - ci * pi * cbi0
      cbk1 = - cbk1 - ci * pi * cbi1
    end if

    cbi1 = - cbi1

  end if

  cdi0 = cbi1
  cdi1 = cbi0 - 1.0D+00 / z * cbi1
  cdk0 = - cbk1
  cdk1 = - cbk0 - 1.0D+00 / z * cbk1

  return
end subroutine


function K03(x)
        real(dp)    :: x, K03
        complex(dp) :: z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1

        K03 = 0d0
        if(x > 1d2) return

        z = 0d0
        z = x

        call cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

        K03 = real(cbk0)
end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Write 2D field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine WriteIT2D(V, file)
        real(dp),    intent(in) :: V(:,:)
        character(len=*)        :: file 
        character(len=25)       :: filename, fmt
        integer                 :: i, j, u

        fmt = '(I5.5)'

        filename = 'dataQW/'//trim(file)//'.dat'

        u = 169

        open(unit=u, file=filename)

        do i=1, size(V,1)
            do j=1, size(V,2)
                write(u,*) V(i,j)
            end do
        end do

        close(u)
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read 2D field from file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ReadIT2D(V, file)
        real(dp), intent(inout) :: V(:,:)
        character(len=*)        :: file 
        character(len=25)       :: filename, fmt
        integer                 :: i, j, u

        fmt = '(I5.5)'

        filename = 'dataQW/'//trim(file)//'.dat'

        u = 169

        open(unit=u, file=filename)

        do i=1, size(V,1)
            do j=1, size(V,2)
                read(u,*) V(i,j)
            end do
        end do

        close(u)
    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Write 2D field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine WriteIT1D(V, file)
        real(dp),    intent(in) :: V(:)
        character(len=*)        :: file 
        character(len=25)       :: filename, fmt
        integer                 :: i, j, u

        fmt = '(I5.5)'

        filename = 'dataQW/'//trim(file)//'.dat'

        u = 169

        open(unit=u, file=filename)

        do i=1, size(V)
            write(u,*) V(i)
        end do

        close(u)
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read 2D field from file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ReadIT1D(V, file)
        real(dp), intent(inout) :: V(:)
        character(len=*)        :: file 
        character(len=25)       :: filename, fmt
        integer                 :: i, j, u

        fmt = '(I5.5)'

        filename = 'dataQW/'//trim(file)//'.dat'

        u = 169

        open(unit=u, file=filename)

        do i=1, size(V)
            read(u,*) V(i)
        end do

        close(u)
    end subroutine
    
    
    
    
    
    ! function: FieldBiInterp
    ! Computes a linear interpolation of a
    ! complex 1D array at a specified position.
    !
    ! Parameters:
    ! f  - 1D array to be interpolated
    ! x  - 1D position array corresponding to 'f'
    ! x0 - Position at which 'f' is to be interpolated
	function EAtX(f, x, x0) result(f0)
		complex(dp), intent(in) :: f(:,:)
		real(dp),    intent(in) :: x(:)
		real(dp),    intent(in) :: x0
        complex(dp)             :: f0(size(f,2))
		
		integer  :: i, j
		
        f0 = 0d0
        
        if(size(x)==1) then
            f0 = f(1,:)
            return
        endif
                
		i = locator(x, x0)
		
        if(i<1 .or. i>=size(x)) return
       
        do j=1, size(f,2)
    		f0(j) = (f(i,j) * (x(i+1) - x0) + f(i+1,j) * (x0 - x(i)) )  / (x(i+1) - x(i))
    	end do
	end function	


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printIT(Dx, z, n, file)
        complex(dp), intent(in) :: Dx(:)
        real(dp),    intent(in) :: z(:)
        integer                 :: n
        character(len=*)        :: file 
        character(len=50)       :: filename, fmt
        integer                 :: i, u



            fmt = '(I6.6)'

            write(filename,fmt) n
            filename = 'dataQW/'//trim(file)//trim(filename)//'.dat'

            u = n+20

            open(unit=u, file=trim(filename))

            do i=1, size(z)
                write(u,*) sngl(z(i)), sngl(real(Dx(i))), sngl(aimag(Dx(i)))
            end do

            close(u)


    end subroutine
    
    
    subroutine printITR(Dx, z, n, file)
        real(dp), intent(in)    :: Dx(:)
        real(dp),    intent(in) :: z(:)
        integer                 :: n
        character(len=*)        :: file 
        character(len=50)       :: filename, fmt
        integer                 :: i, u



            fmt = '(I6.6)'

            write(filename,fmt) n
            filename = 'dataQW/'//trim(file)//trim(filename)//'.dat'

            u = n+20

            open(unit=u, file=trim(filename))

            do i=1, size(z)
                write(u,*) sngl(z(i)), sngl(real(Dx(i)))
            end do

            close(u)


    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printIT2D(Dx, z, n, file)
        complex(dp), intent(in) :: Dx(:,:)
        real(dp),    intent(in) :: z(:)
        integer                 :: n
        character(len=*)        :: file 
        character(len=50)       :: filename, fmt
        integer                 :: i, u

            fmt = '(I7.7)'

            write(filename,fmt) n
            filename = 'dataQW/'//trim(file)//trim(filename)//'.dat'

            u = n+20

            open(unit=u, file=trim(filename))

            do j=1, size(Dx,2)
                do i=1, size(Dx,1)
                    write(u,*) sngl(abs(Dx(i,j)))!, sngl(aimag(Dx(i,j)))
                end do
            end do

            close(u)
    end subroutine


    elemental function gaussian(x,x0)
        double precision, intent(in) :: x, x0
        double precision :: gaussian

        gaussian = exp(-x**2 / x0**2)
    end function gaussian


  function convolve(x, h)
        implicit none
        
        !x is the signal array
        !h is the noise/impulse array
        complex(dp), dimension(:), allocatable :: convolve, y
        complex(dp), dimension(:) :: x, h
        integer :: kernelsize, datasize
        integer :: i,j,k
        
        datasize = size(x)
        kernelsize = size(h)
        
        allocate(y(datasize))
        allocate(convolve(datasize))
       
        !last part
        do i=kernelsize,datasize
            y(i) = 0.0
            j=i
            do k=1,kernelsize
                y(i) = y(i) + x(j)*h(k)
                j = j-1
            end do
        end do
        
        !first part
        do i=1,kernelsize
            y(i) = 0.0
            j=i
            k=1
            do while (j > 0)
                y(i) = y(i) + x(j)*h(k)
                j = j-1
                k = k+1
            end do
        end do
        
        convolve = y
               
    end function convolve


    subroutine FFTG(F)
        complex(dp), intent(inout) :: F(:)
        integer                 :: Nf
        Nf = size(F)
        call FFTC(F)
        F = - F / (Nf*1d0)
    end subroutine
    
    subroutine iFFTG(F)
        complex(dp), intent(inout) :: F(:)
        integer                 :: Nf
        Nf = size(F)
        call iFFTC(F)
        F = - F * Nf
    end subroutine

end module
