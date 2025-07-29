module phost
    use types
    use constants
    use helpers
    use fftw

    implicit none

    integer,     private              :: osc = 2
    real(dp),    private              :: q = -e0
    real(dp),    private, allocatable :: w(:), gam(:)
    real(dp),    private, allocatable :: B(:), C(:), chi1(:), Nf(:)
    real(dp),    private              :: A0
    complex(dp), private              :: chi3 = 0d-3
    real(dp),    private              :: epsr_0 = 10.0d0
    real(dp),    private              :: epsr_infty = 8.2
    real(dp),    private              :: w1 = 0d0
    real(dp),    private              :: w2 = 1d20
    
    real(dp),    private              :: w0, lambda0
    integer,     private              :: N1, N2
    character(len=4),private          :: material = 'AlAs'
    
    complex(dp), private, allocatable :: Px_before(:,:,:), Py_before(:,:,:) 
    complex(dp), private, allocatable :: Px_now(:,:,:)   , Py_now(:,:,:)
    complex(dp), private, allocatable :: Px_after(:,:,:) , Py_after(:,:,:)
    complex(dp), private, allocatable :: omega_q(:,:)
    complex(dp), private, allocatable :: EpsrWq(:,:)



    
contains
	
    subroutine CalcPHost(Ex, Ey, dt, m, epsb, Px, Py)
        complex*16, intent(in   ) :: Ex(:,:), Ey(:,:)
	    real(dp),   intent(in   ) :: dt
	    integer,    intent(in   ) :: m
	    real(dp),   intent(inout) :: epsb
        complex*16, intent(inout) :: Px(:,:), Py(:,:)


        Px_before = Px_now
        Py_before = Py_now

		Px_now    = Px_after
		Py_now    = Py_after
            
		Px_after  = CalcNextP(Px_before, Px_now, Ex, dt)
		Py_after  = CalcNextP(Py_before, Py_now, Ey, dt)
		    
        epsb      = A0
        Px(:,:)   = sum(Px_after,3)
        Py(:,:)   = sum(Py_after,3)

   !Px = 0d0
   !Py = 0d0
   !epsb = nw2_no_gam(w0)
    end subroutine CalcPHost
	
	
	
    subroutine CalcPHostOld(Ex, Ey, dt, m, epsb, Px, Py)
        complex*16, intent(in   ) :: Ex(:,:), Ey(:,:)
	    real(dp),   intent(in   ) :: dt
	    integer,    intent(in   ) :: m
	    real(dp),   intent(inout) :: epsb
        complex*16, intent(inout) :: Px(:,:), Py(:,:)


        Px_before = Px_now
        Py_before = Py_now
		
        if(m > 2) then
		    Px_now    = Px_after
		    Px_after  = CalcNextP(Px_before, Px_now, Ex, dt)
            
		    Py_now    = Py_after
		    Py_after  = CalcNextP(Py_before, Py_now, Ey, dt)
		    
            epsb      = A0
        elseif(m>=2) then
		    Px_now    = CalcMonoP(Ex)
		    Px_after  = CalcNextP(Px_before, Px_now, Ex, dt)
            
		    Py_now    = CalcMonoP(Ey)
		    Py_after  = CalcNextP(Py_before, Py_now, Ey, dt)
        
            epsb     = A0
        else
            Px_now    = CalcMonoP(Ex)
            Px_after  = 0d0
            
            Py_now    = CalcMonoP(Ey)
            Py_after  = 0d0
            
            epsb      = nw2_no_gam(w0)
        end if
        
        Px(:,:)   = sum(Px_after,3)
        Py(:,:)   = sum(Py_after,3)

    end subroutine CalcPHostOld
	
	
	
	
	function CalcNextP(P1, P2, E, dt)
	    complex(dp), intent(in   ) :: P1(:,:,:), P2(:,:,:), E(:,:)
	    real(dp)   , intent(in   ) :: dt
	    complex(dp)                :: CalcNextP(size(E,1), size(E,2), osc)
	    real(dp)                   :: f1(osc), f2(osc), f3(osc)
	    integer                    :: i,j,n
	
		f1(:) = - (1d0  - gam(:)  * dt   ) / (gam(:) * dt + 1d0)
	    f2(:) =   (2d0  - w(:)**2 * dt**2) / (gam(:) * dt + 1d0)
		f3(:) =   (B(:) * w(:)**2 * dt**2) / (gam(:) * dt + 1d0) * eps0
		
		
        !$omp parallel do private(n, j, i)   
        do n = 1, osc
            do j=1, size(E,2)
                do i=1, size(E,1)
                    CalcNextP(i,j,n) = f1(n) * P1(i,j,n) + f2(n) * P2(i,j,n) + f3(n) * E(i,j)
                end do
            end do
        end do
        !$omp end parallel do
	end function
	
	
	
		
	function CalcMonoP(E)
	    complex(dp), intent(in) :: E(:,:)
	    complex(dp)             :: CalcMonoP(size(E,1), size(E,2),osc)
	    integer                 :: i,j,n
	    
        !$omp parallel do private(n, j, i)   
        do n = 1, osc
            do j=1, size(E,2)
                do i=1, size(E,1)
                    CalcMonoP(i,j,n) = eps0 * E(i,j) * real(chi1(n))
                end do
            end do
        end do
        !$omp end parallel do
	end function
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!! All Remaining Code is for Initializing this Module !!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	subroutine SetHostMaterial(host, mat, lam, epsr, n0)
	    logical,          intent(in   ) :: host
        character(len=*), intent(in   ) :: mat
        real(dp),         intent(in   ) :: lam
        real(dp),         intent(inout) :: epsr, n0

        lambda0 = lam

        if    (trim(mat)=='AlAs') then
            call SetParamsAlAs()
        elseif(trim(mat)=='fsil') then
            call SetParamsSilica()
        elseif(trim(mat)=='GaAs') then
            call SetParamsGaAs()
        elseif(trim(mat)=='none') then
            call SetParamsNone()
        else
            print*, "ERROR: Host Material = ", mat, "Is Not Included In phost.f90 Code"
            stop
        end if


        material = trim(mat)

        allocate(chi1(osc))
        
        chi1 = B * lam**2 / (lam**2 - C)

        w0 = twopi*c0/lam
        
        if(host) then
            epsr = real(nw2_no_gam(w0))
            n0   = real(sqrt(epsr))
            call WriteHostDispersion()
        endif
  
        epsr = n0**2
    end subroutine SetHostMaterial
	
	
	
    subroutine InitializeHost(Nx, Ny, n0, qsq, host)
        integer,     intent(in   ) :: Nx, Ny
        real(dp),    intent(in   ) :: n0
        complex(dp), intent(in   ) :: qsq(:,:)
        logical,     intent(in   ) :: host

        allocate(omega_q(Nx,Ny))
        allocate(EpsrWq(Nx,Ny))
        
        if(host) then
            allocate(Px_before(Nx, Ny, osc), Px_now(Nx, Ny, osc), Px_after(Nx, Ny, osc))
            allocate(Py_before(Nx, Ny, osc), Py_now(Nx, Ny, osc), Py_after(Nx, Ny, osc))
        
            Px_before = 0d0
            Px_now    = 0d0
            Px_after  = 0d0
            Py_before = 0d0
            Py_now    = 0d0
            Py_after  = 0d0
            omega_q   = 0d0
            call CalcWq(sqrt(qsq))
            call CalcEpsrWq(sqrt(qsq))
        else
            omega_q = sqrt(real(qsq)) * c0 / n0
            EpsrWq  = n0**2
        end if
    end subroutine InitializeHost
	
	
	
    subroutine CalcWq(q)
        complex(dp), intent(in   ) :: q(:,:)
        complex(dp)                :: n0, np0, x
        integer                    :: i,j
	
	    n0  = sqrt(nw2_no_gam(w0))
	    np0 = nwp_no_gam(w0)
	    x   = n0 - np0 * w0
	
	    !w   = (c0 * q + np0 * w0**2) / x
	    
	    omega_q   = (- x + sqrt(x**2 + 4*np0*q*c0)) / (2*np0)

	    open(unit=750, file="fields/host/w.q.dat")
	    do j=1, max(size(q,2)/2, 1)
	        do i=1, max(size(q,1)/2, 1)
	            write(750,*) real(q(i,j))*1d-7, real(omega_q(i,j))*1d-15, aimag(omega_q(i,j))*1d-15
	        end do
	    end do
	    close(750)
    end subroutine CalcWq



    subroutine CalcEpsrWq(q)
        complex(dp), intent(in   ) :: q(:,:)
        real(dp)                   :: aw(2), bw(2)
        integer                    :: i,j

        call DetermineCoeffs(aw, bw)

        do j=1, size(q,2)
            do i=1, size(q,1)

                call CalcEpsrWq_ij(abs(omega_q(i,j)), aw, bw, EpsrWq(i,j))

                !EpsrWq(i,j) = 1d0 - nw2_no_gam(real(omega_q(i,j)))
                !ChiWq(i,j) = min(abs(ChiWq(i,j)), 15d0)
       	       	!ChiWq(i,j) = max(abs(ChiWq(i,j)), 1d0)
            end do
        end do
    end subroutine


    subroutine CalcEpsrWq_ij(w_ij,  aw, bw, epsr_ij)
        real(dp),    intent(in   ) :: w_ij, aw(2), bw(2)
        complex(dp), intent(inout) :: epsr_ij

        if    (w_ij < w1) then

            epsr_ij = epsr_0 + aw(1) * w_ij**2 + aw(2) * w_ij**3

        elseif(w_ij > w2) then

            epsr_ij = epsr_infty + bw(1) / w_ij**2 + bw(2) / w_ij**3

        else
            epsr_ij = nw2_no_gam(w_ij)
        end if
        
    end subroutine



    subroutine DetermineCoeffs(aw, bw)
        real(dp), intent(inout) :: aw(2), bw(2)
        real(dp)                :: e1, e2, ep1, ep2
        
        e1 = nw2_no_gam(w1)
        e2 = nw2_no_gam(w2)

        ep1 = epsrwp_no_gam(w1)
        ep2 = epsrwp_no_gam(w2)

        aw(1) = + 3/w1**2 * (e1 - epsr_0) - ep1 / w1
        aw(2) = - 2/w1**3 * (e1 - epsr_0) + ep1 / w1**2

        bw(1) = + 3*w2**2 * (e2 - epsr_infty) + ep2 * w2**3
        bw(2) = - 2*w2**3 * (e2 - epsr_infty) - ep2 * w2**4

    end subroutine

    function Epsr_q(q)
        complex(dp)    :: q(:,:)
        complex(dp) :: Epsr_q(size(q,1),size(q,2))
 
        Epsr_q(:,:) = EpsrWq(:,:)
    end function


    function Epsr_qij(i,j)
        integer     :: i,j
        complex(dp) :: Epsr_qij
 
        Epsr_qij = EpsrWq(i,j)
    end function


    subroutine FDTD_Dispersion(qx, qy, dx, dy, dt, n0)
        real(dp), intent(in) :: qx(:), qy(:), dx, dy, dt, n0
        integer              :: i, j

        omega_q = 0d0
        
        do j=1, size(qy)
            do i=1, size(qx)
                omega_q(i,j) =  sqrt( Sin(qx(i) * dx / 2_dp)**2 / dx**2 + Sin(qy(j) * dy / 2_dp)**2 / dy**2 )
                omega_q(i,j) =  2_dp / dt * ASin( (c0/n0) * dt * real(omega_q(i,j)) ) 
            end do
        end do
    end subroutine FDTD_Dispersion


    function wq(i,j)
        integer     :: i,j
        complex(dp) :: wq
        wq = omega_q(i,j)
    end function
	
	
	
    subroutine SetInitialP(Ex, Ey, qx, qy, qsq, dt, Px, Py, epsb)
        complex(dp), intent(inout) :: Ex(:,:), Ey(:,:)
        real(dp),    intent(in)    :: qx(:), qy(:), qsq(:,:), dt
        complex(dp), intent(inout) :: Px(:,:), Py(:,:)
        real(dp),    intent(inout) :: epsb
        integer                    :: n

        forall(n=1:osc) Px_after(:,:,n) = eps0 * Ex(:,:) * B(n)*w(n)**2 / (w(n)**2 - omega_q(:,:)**2)
        forall(n=1:osc) Py_after(:,:,n) = eps0 * Ey(:,:) * B(n)*w(n)**2 / (w(n)**2 - omega_q(:,:)**2)

        forall(n=1:osc) Px_now(:,:,n) = Px_after(:,:,n) * exp(-ii * omega_q(:,:) * (-dt))
        forall(n=1:osc) Py_now(:,:,n) = Py_after(:,:,n) * exp(-ii * omega_q(:,:) * (-dt))

        call MakeTransverse(Ex      , Ey      , qx, qy, qsq)
        do n=1, osc
            call MakeTransverse(Px_now(:,:,n)  , Py_now(:,:,n)  , qx, qy, qsq)
            call MakeTransverse(Px_after(:,:,n), Py_after(:,:,n), qx, qy, qsq)
        end do

        do n=1, osc
            call IFFT(Px_now(:,:,n))
            call IFFT(Py_now(:,:,n))
            call IFFT(Px_after(:,:,n))
            call IFFT(Py_after(:,:,n))
        end do

        Px_now    = real(Px_now)
        Py_now    = real(Py_now)
        Px_after  = real(Px_after)
        Py_after  = real(Py_after)

        epsb      = A0
        Px(:,:)   = sum(Px_after,3)
        Py(:,:)   = sum(Py_after,3)
    end subroutine SetInitialP
    
	
    subroutine MakeTransverse(Ex, Ey, qx, qy, qsq)
        complex(dp), intent(Inout) :: Ex(:,:), Ey(:,:)
        real(dp),    intent(in   ) :: qx(:), qy(:), qsq(:,:)
        integer                    :: j, Ny

        Ny = size(Ex,2)

        forall(j=1:Ny) Ex(:,j) = Ex(:,j) - qx(:) * (qx(:) * Ex(:,j) + qy(j) * Ey(:,j)) / qsq(:,j)
        forall(j=1:Ny) Ey(:,j) = Ey(:,j) - qy(j) * (qx(:) * Ex(:,j) + qy(j) * Ey(:,j)) / qsq(:,j)
    end subroutine MakeTransverse
	
	
    subroutine SetParamsSilica()
        osc   = 3
        allocate(B(osc), C(osc), w(osc), gam(osc), Nf(osc))
    
        A0    = 1d0
        B(1)  = 0.696166300
        B(2)  = 0.407942600
        B(3)  = 0.897479400
        C(1)  = 4.67914826d-3 * 1d-12
        C(2)  = 1.35120631d-2 * 1d-12
        C(3)  = 97.9340025    * 1d-12


        w   = twopi * c0 / sqrt(C)
        gam = w / 25_dp 


        Nf = B * w**2 * eps0 * me0 / q**2
        
        epsr_0     = A0 + sum(B*lambda0**2 / (lambda0**2 - C))
        epsr_infty = A0 + sum(B*lambda0**2 / (lambda0**2 - C))
    end subroutine SetParamsSilica



    subroutine SetParamsGaAs()
        osc   = 3
        allocate(B(osc), C(osc), w(osc), gam(osc), Nf(osc))
    
        A0    = 4.37251400
        B(1)  = 5.46674200
        B(2)  = 0.02429960
        B(3)  = 1.95752200
        C(1)  = 0.4431307**2 * 1d-12
        C(2)  = 0.8746453**2 * 1d-12
        C(3)  = 36.916600**2 * 1d-12

        w    = twopi * c0 / sqrt(C)
        gam  = w / 10_dp 

        Nf = B * w**2 * eps0 * me0 / q**2
        
        epsr_0     = A0 + sum(B*lambda0**2 / (lambda0**2 - C))
        epsr_infty = A0 + sum(B*lambda0**2 / (lambda0**2 - C))
    end subroutine SetParamsGaAs
		


    subroutine SetParamsAlAs()
        osc   = 2
        allocate(B(osc), C(osc), w(osc), gam(osc), Nf(osc))
    
        A0    = 2.0792           
        B(1)  = 6.0840           
        B(2)  = 1.9000            
        C(1)  = 0.2822**2 * 1d-12
        C(2)  = 27.620**2 * 1d-12

        w     = twopi * c0 / sqrt(C)
        !gam = w / 1000_dp 
        gam = 0d0

        Nf = B * w**2 * eps0 * me0 / q**2

        w1 = twopi * c0 / (2.2d-6)
        w2 = twopi * c0 / (0.56d-6)
        
        epsr_0     = 10.0d0
        epsr_infty = 8.2d0
    end subroutine  

    subroutine SetParamsNone()
        osc   = 1
        allocate(B(osc), C(osc), w(osc), gam(osc), Nf(osc))
    
        A0    = 1d0          
        B(1)  = 0d0           
        C(1)  = 1d0

        w     = twopi * c0 / sqrt(C)
        gam   = 0d0 
        Nf = 0d0
        
        epsr_0     = 1d0
        epsr_infty = 1d0

    end subroutine  


!   function nl2(lam)
!       real(dp) :: lam, nl2
!       
!       if(lam < lammin) then
!           nl2 = epsr_infty
!       elseif(lam > lammax) then
!           nl2 = epsr_0
!       else 
!           nl2 = A0 + sum(B*lam**2 / (lam**2 - Ci))
!       endif
!   end function


   function nw2_no_gam(wL)
       real(dp)    :: wL
       complex(dp) :: nw2_no_gam
       nw2_no_gam = A0 + sum(B*w**2 / (w**2 - wL**2))
   end function
   
   
   function nw2(wL)
       real(dp)    :: wL
       complex(dp) :: nw2
       nw2 = A0 + sum(B(:)*w(:)**2 / (w(:)**2 - ii*2*gam*wL- wL**2))
   end function
   
   
   function nwp_no_gam(wL)
       real(dp)    :: wL
       complex(dp) :: nwp_no_gam, nw
       nw  = sqrt(nw2_no_gam(wL))
       
       nwp_no_gam = sum(B(:)*w(:)**2 * (wL) / (w(:)**2 - wL**2)**2) / nw
   end function   
   

   function epsrwp_no_gam(wL)
       real(dp)    :: wL
       complex(dp) :: epsrwp_no_gam
       
       epsrwp_no_gam = sum(B(:)*w(:)**2 * (2*wL) / (w(:)**2 - wL**2)**2)
   end function  

   
   function nwp(wL)
       real(dp)    :: wL
       complex(dp) :: nwp, nw
       nw  = sqrt(nw2(wL))
       
       nwp = sum(B(:)*w(:)**2 * (ii*gam + wL) / (w(:)**2 - ii*2*gam*wL- wL**2)**2) / nw
   end function
   
   
   function nl2_no_gam(lam)
       real(dp)    :: lam
       complex(dp) :: nl2_no_gam
       nl2_no_gam = A0 + sum(B*lam**2 / (lam**2 - C))
   end function
    
    
   function nl2(lam)
       real(dp)    :: lam
       complex(dp) :: nl2
       complex(dp) :: wL
       wL = twopi*c0 / (lam+1d-100)
       nl2 = A0 + sum(B(:)*w(:)**2 / (w(:)**2 - ii*2*gam*wL- wL**2))
   end function
    
    subroutine WriteHostDispersion()
        real(dp)  			 :: x, dx, x0, xf
        complex(dp)          :: n2, n2X, n, nX
        integer              :: l, Nxx			


        Nxx = 10000
        
        x0 = 0d0
        xf = maxval(w)*3d0
        dx = (xf-x0) / Nxx

        open(unit=745,file="fields/host/n.w.real.dat")
        open(unit=746,file="fields/host/n.w.imag.dat")
        open(unit=747,file="fields/host/epsr.w.real.dat")
        open(unit=748,file="fields/host/epsr.w.imag.dat")
        open(unit=749,file="fields/host/nogam/n.w.real.dat")
        open(unit=750,file="fields/host/nogam/n.w.imag.dat")
        open(unit=751,file="fields/host/nogam/epsr.w.real.dat")
        open(unit=752,file="fields/host/nogam/epsr.w.imag.dat")    
        open(unit=753,file="fields/host/q2.w.real.dat")
        open(unit=754,file="fields/host/q2.w.imag.dat")
        do l=1, Nxx
            x = x0 + (l)*dx
            !n2  = nw2(twopi*c0/lam)
            n2   = nw2(x)
            n2X  = nw2_no_gam(x)
            n    = sqrt(n2)
            nX   = sqrt(n2X)

            write(745,*) x, real(n)
            write(746,*) x, aimag(n)
            write(747,*) x, real(n2)
            write(748,*) x, aimag(n2)
            write(749,*) x, real(nX)   * exp(-((abs(n2X)-9)/8)**12)
            write(750,*) x, aimag(nX)  * exp(-((abs(n2X)-9)/8)**12)
            write(751,*) x, real(n2X)  * exp(-((abs(n2X)-9)/8)**12)
            write(752,*) x, aimag(n2X) * exp(-((abs(n2X)-9)/8)**12)
        end do
        close(745)
        close(746)
        close(747)
        close(748)
        close(749)
        close(750)
        close(751)
        close(752)
        close(753)
        close(754)


        x0 = 0d0
        xf = twopi * c0 / (minval(w)+1d-100) * 3d0
        dx = (xf-x0) / Nxx

        open(unit=745,file="fields/host/n.l.real.dat")
        open(unit=746,file="fields/host/n.l.imag.dat")
        open(unit=747,file="fields/host/epsr.l.real.dat")
        open(unit=748,file="fields/host/epsr.l.imag.dat")
        open(unit=749,file="fields/host/nogam/n.l.real.dat")
        open(unit=750,file="fields/host/nogam/n.l.imag.dat")
        open(unit=751,file="fields/host/nogam/epsr.l.real.dat")
        open(unit=752,file="fields/host/nogam/epsr.l.imag.dat")    
        do l=1, Nxx
            x = x0 + (l-1)*dx
            n2   = nl2(x)
            n2X  = nl2_no_gam(x)
            n    = sqrt(n2)
            nX   = sqrt(n2X)
            
            ! Convert to micrometers
            x    = x*1d6

            write(745,*) x, real(n)
            write(746,*) x, aimag(n)
            write(747,*) x, real(n2)
            write(748,*) x, aimag(n2)
            write(749,*) x, real(nX)   * exp(-((abs(n2X)-9)/8)**12)
            write(750,*) x, aimag(nX)  * exp(-((abs(n2X)-9)/8)**12)
            write(751,*) x, real(n2X)  * exp(-((abs(n2X)-9)/8)**12)
            write(752,*) x, aimag(n2X) * exp(-((abs(n2X)-9)/8)**12)
        end do
        close(745)
        close(746)
        close(747)
        close(748)
        close(749)
        close(750)
        close(751)
        close(752)



    end subroutine
	
end module
