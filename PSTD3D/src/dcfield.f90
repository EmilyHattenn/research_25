! This module calculates the dc field carrier transport
! contributions to the Semiconductor Bloch equations in
! support of propagation simulations for a quantum wire.
module dcfield

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers
    use usefulsubs
    use spliner
    implicit none



    ! Pre-calculated QW arrays required for solving SBEs.
    real(dp), private              :: small = 1d-200  ! Smallest # worthy of consideration
    real(dp), private, allocatable :: Y(:), QY(:)     ! fft array in y-space for derivative evals
    real(dp), private, allocatable :: xe(:), xh(:)    ! k-dependent delta-func. coefficients
    real(dp), private, allocatable :: qinv(:)         !
    real(dp), private              :: ERate, HRate    ! The temp damping rates for elecs & holes
    real(dp), private              :: VEDrift, VHDrift! Stored drift velocities
    real(dp), private              :: dY, dQy, dkk    ! Y and QY       
    real(dp), private              :: kmin, kmax      ! 
    integer,  private              :: NQ, Nk          ! ky and QY point numbers
    logical,  private              :: WithPhns=.true. ! Couple the damping rate with phonons or no
   


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    subroutine InitializeDC(ky, me, mh)
        real(dp), intent(in   ) :: ky(:)
        real(dp), intent(in   ) :: me, mh
        real(dp)                :: Lk
        real(dp)                :: dky, dffty
        integer                 :: Nk, k, k1
        
        Nk  = size(ky)
        dky = ky(3)-ky(2)
        dkk = dky

        ERate = 0d0
        HRate = 0d0

        allocate(Y(Nk))

        Y    = GetKArray(Nk, (Nk-1) *dky)

        allocate(xe(Nk), xh(Nk))

        xe =  me / hbar**2 * abs(ky) / (abs(ky)+dky*1d-5)**2 / dky
        xh =  mh / hbar**2 * abs(ky) / (abs(ky)+dky*1d-5)**2 / dky

        allocate(qinv(0:Nk+1))
        qinv = 0d0
        qinv(1:Nk) = ky / (abs(ky)+dky*1d-5)**2
        
        kmin = ky(1 ) - 2*dky
        kmax = ky(Nk) + 2*dky

        open(file='dataQW/Fe.dat',unit=931)
        open(file='dataQW/Fh.dat',unit=932)
    end subroutine InitializeDC
    

    subroutine CalcDCE2(DCTrans, ky, Cq2, Edc, me, ge, Ephn, N0, ne, Ee, Vee, n,j, DC)
        logical,     intent(in   ) :: DCTrans
        real(dp),    intent(in   ) :: ky(:)
        real(dp),    intent(in   ) :: Cq2(:)
        real(dp),    intent(in   ) :: Edc, me, ge
        real(dp),    intent(in   ) :: Ephn, N0
        complex(dp), intent(in   ) :: ne(:)
        real(dp),    intent(in   ) :: Ee(:), Vee(:,:) 
integer         ,    intent(in   ) :: n,j       
        real(dp),    intent(inout) :: DC(:)
        complex(dp)                :: DC0(size(ky))
        real(dp)                   :: Eec(size(ky))
        real(dp)                   :: Fd(size(ky))
        real(dp)                   :: gate(size(ky))
        real(dp)                   :: v, dk
        integer                    :: k, Nk
        
        DC   = 0d0
        gate = exp(-(ky / 0.6d9)**6) 
gate=1d0
        dk   = ky(2) - ky(1) 
       
        Eec = EkReNorm(real(ne(:)),  Ee(:), Vee(:,:))

        !if(.not.DCTrans) return

        v   = DriftVt( real(ne(:)), Eec(:))        

        DC0 = 0d0
        
        if(WithPhns) then
           !Fd = FDrift( Ephn, me, ky, real(dndk), Cq2, v, N0, xe)
           !Fd = FDrift2(Ephn, me, ky, real(ne), Cq2, v, N0, xe)
            Fd = FDrift2(Ephn, me, ge, ky, real(ne), Cq2, v, N0, xe)
        else
           !Fd = CalcPD(ky, me, ne) * 1d13
           !Fd = hbar * ky * 1d12
            Fd = 0d0
        endif


        if(j==25) call printIT(Fd*(1d0,0d0), ky, n, "Wire/Fe/Fe.k.")
 
        ERate   = sum(Fd(:) / hbar / (abs(ky) + 1d-5)**2 * ky * ne) / sum(ne)

        Fd = sum(Fd) / (sum(abs(ne))+1d-20) * 2

        write(931,*) n, Fd(1)

        if(.not.DCTrans) return

        DC0 = -(- e0 * Edc - Fd) * gate / hbar * ne


        VEDrift = v

         !call FFT(DC0)
         !   DC0 = DC0 * (ii*Y)
         !call IFFT(DC0)

        DC0 = (cshift(DC0,1) - cshift(DC0,-1)) / 2d0 / dk

        DC = real(DC0)

        if(j==25) call printIT(DC0, ky, n, "Wire/Fe/De.k.")

    end subroutine CalcDCE2
    
    
    subroutine CalcDCH2(DCTrans, ky, Cq2, Edc, mh, gh, Ephn, N0, nh, Eh, Vhh, n,j,DC)
        logical,     intent(in   ) :: DCTrans
        real(dp),    intent(in   ) :: ky(:)
        real(dp),    intent(in   ) :: Cq2(:)
        real(dp),    intent(in   ) :: Edc, mh, gh
        real(dp),    intent(in   ) :: Ephn, N0
        complex(dp), intent(in   ) :: nh(:)
        real(dp),    intent(in   ) :: Eh(:), Vhh(:,:)
        integer,     intent(in   ) :: n,j       
        real(dp),    intent(inout) :: DC(:)
        complex(dp)                :: DC0(size(ky))
        real(dp)                   :: Ehc(size(ky))
        real(dp)                   :: Fd(size(ky))
        real(dp)                   :: gate(size(ky))
        real(dp)                   :: v, dk
        integer                    :: k, Nk
        
        DC   = 0d0

!gate = exp(-(ky / 1d9)**8)
gate = 1d0
        dk   = ky(2) - ky(1)

        Ehc = EkReNorm(real(nh(:)),  Eh(:), Vhh(:,:))

        !if(.not.DCTrans) return

        v   = DriftVt( real(nh(:)), Ehc(:))
        
        DC0 = 0d0

        if(WithPhns) then
           !Fd = FDrift( Ephn, mh, ky, real(dndk), Cq2, v, N0, xh)
           !Fd = FDrift2(Ephn, mh, ky, real(nh), Cq2, v, N0, xh)
            Fd = FDrift2(Ephn, mh, gh, ky, real(nh), Cq2, v, N0, xh)
        else
            !Fd = CalcPD(ky, mh, nh) * 1d13
            !Fd = hbar * ky * 1d12
             Fd = 0d0
        endif

 HRate   = sum(Fd(:) / hbar / (abs(ky) + 1d-5)**2 * ky * nh) / sum(nh)

 if(j==25)call printIT(Fd*(1d0,0d0), ky, n, "Wire/Fh/Fh.k.")

 Fd = sum(Fd) / (sum(abs(nh))+1d-20) * 2

        write(932,*) n, Fd(1)

        if(.not.DCTrans) return

        DC0 = -(- e0 * Edc - Fd) * gate / hbar * nh


        VHDrift = v

         !call FFT(DC0)
         !    DC0 = DC0 * (ii*Y)
         !call IFFT(DC0)

DC0 = (cshift(DC0,1) - cshift(DC0,-1)) / 2d0 / dk

        DC = real(DC0)

if(j==25) call printIT(DC0, ky, n, "Wire/Fh/Dh.k.")
    end subroutine CalcDCH2


    
    subroutine CalcDCE(DCTrans, ky, Cq2, Edc, me, ge, Ephn, N0, ne, Ee, Vee, DC)
        logical,     intent(in   ) :: DCTrans
        real(dp),    intent(in   ) :: ky(:)
        real(dp),    intent(in   ) :: Cq2(:)
        real(dp),    intent(in   ) :: Edc, me, ge
        real(dp),    intent(in   ) :: Ephn, N0
        complex(dp), intent(in   ) :: ne(:)
        real(dp),    intent(in   ) :: Ee(:), Vee(:,:)        
        real(dp),    intent(inout) :: DC(:)
        complex(dp)                :: dndk(size(ky))
        real(dp)                   :: Eec(size(ky))
        real(dp)                   :: v, Fd(size(ky))
        integer                    :: k, Nk
        
        DC = 0d0
        
        Eec = EkReNorm(real(ne(:)),  Ee(:), Vee(:,:))

        if(.not.DCTrans) return

        v   = DriftVt( real(ne(:)), Eec(:))        

        dndk = 0d0
        dndk = real(ne(:))
        
        call FFT(dndk)
            dndk = dndk * (ii*Y)
        call IFFT(dndk)

        if(WithPhns) then
           !Fd = FDrift( Ephn, me, ky, real(dndk), Cq2, v, N0, xe)
            Fd = FDrift2(Ephn, me, ge, ky, real(ne), Cq2, v, N0, xe)
        else
           !Fd = CalcPD(ky, me, ne) * 1d13
            Fd = 0d0
        endif

        DC = -(- e0 * Edc - Fd) / hbar * real(dndk)

       !ERate   = Fd / me / (abs(v) + 1d0)**2 * v
        ERate   = sum(Fd(:) / hbar / (abs(ky) + 1d-5)**2 * ky * ne) / sum(ne)

        VEDrift = v
    end subroutine CalcDCE
    
    
    subroutine CalcDCH(DCTrans, ky, Cq2, Edc, mh, gh, Ephn, N0, nh, Eh, Vhh, DC)
        logical,     intent(in   ) :: DCTrans
        real(dp),    intent(in   ) :: ky(:)
        real(dp),    intent(in   ) :: Cq2(:)
        real(dp),    intent(in   ) :: Edc, mh, gh
        real(dp),    intent(in   ) :: Ephn, N0
        complex(dp), intent(in   ) :: nh(:)
        real(dp),    intent(in   ) :: Eh(:), Vhh(:,:)
        real(dp),    intent(inout) :: DC(:)
        complex(dp)                :: dndk(size(ky))
        real(dp)                   :: Ehc(size(ky))
        real(dp)                   :: v, Fd(size(ky))
        integer                    :: k, Nk
        
        DC = 0d0

        Ehc = EkReNorm(real(nh(:)),  Eh(:), Vhh(:,:))

        if(.not.DCTrans) return

        v   = DriftVt( real(nh(:)), Ehc(:))
        
        dndk = 0d0
        dndk = real(nh)
        
        call FFT(dndk)
            dndk = dndk * (ii*Y)
        call IFFT(dndk)

        if(WithPhns) then
           !Fd = FDrift( Ephn, mh, ky, real(dndk), Cq2, v, N0, xh)
            Fd = FDrift2(Ephn, mh, gh, ky, real(nh), Cq2, v, N0, xh)
        else        
           !Fd = CalcPD(ky, mh, nh) * 1d13
            Fd = 0d0
        endif

        DC = -(+ e0 * Edc - Fd) / hbar * real(dndk)

       !HRate   = Fd / mh / (abs(v) + 1d0)**2 * v
        HRate   = sum(Fd(:) / hbar / (abs(ky) + 1d-5)**2 * ky * nh) / sum(nh)
        
        VHDrift = v
    end subroutine CalcDCH

    function CalcI0n(ne, me, ky)
        complex(dp), intent(in   ) :: ne(:)
        real(dp),    intent(in   ) :: me, ky(:)
        real(dp)                   :: CalcI0n
        real(dp)                   :: dk
        real(dp)                   :: v, ve, Ie
        real(dp)                   :: Ec(size(ne))

        dk = ky(2)-ky(1)

        Ie    = -e0 * sum(ne*ky*hbar/me)* 2 * dk

        CalcI0n = Ie

        !CalcI0n = sum(ne*ky)/dk
    end function


    subroutine CalcI0(ne, nh, Ee, Eh, VC, dk, ky, I0)
        complex(dp), intent(in   ) :: ne(:), nh(:)
        real(dp),    intent(in   ) :: Ee(:), Eh(:)
        real(dp),    intent(in   ) :: VC(:,:,:), dk, ky(:)
        real(dp),    intent(inout) :: I0
        real(dp)                   :: v, ve, vh	
        real(dp)                   :: Ec(size(ne))

dkk  = dk

        Ec(:) = EkReNorm(real(ne(:)), Ee(:), VC(:,:,2))
        ve    = DriftVt( real(ne(:)), Ec(:))

!print*, "Ee", maxval(Ee), maxval(Ec)

        Ec(:) = EkReNorm(real(nh(:)), Eh(:), VC(:,:,3))
        vh    = DriftVt( real(nh(:)), Ec(:))

!print*, "Eh", maxval(Eh), maxval(Ec)

       !v     = abs(ve) + abs(vh)
        v     = ve + vh

!print*, "ve", ve, real(sum(ne*ky) / (sum(abs(ne))+1d-100) / 0.07d0/ me0) * hbar
!print*, "vh", vh, real(sum(nh*ky) / (sum(abs(nh))+1d-100) / 0.45d0/ me0) * hbar

        I0    = - e0 * v * sum(ne) * dk * 2
    end subroutine
    

    function EkReNorm(n, En, V) result(Ec)
        real(dp), intent(in) :: n(:), En(:), V(:,:)
        real(dp)             :: Ec(size(n))
        integer              :: k

        do k=1, size(n)
            Ec(k) = En(k) + sum( n(:) * (V(k,k) - V(k,:)) ) / 2d0
        end do
    end function

    

    real(dp) function DriftVt(n, Ec) result(v)
        real(dp), intent(in) :: n(:), Ec(:)
        complex(dp)          :: dEdk(size(n))
        integer :: i, Nk

        dEdk = 0d0

        Nk=size(n)

        !dEcdk(:) = Ec(:)
        !call FFT(dEdk)
        !    dEdk = dEdk * (ii*Y)
        !call IFFT(dEdk)

        do i=2,Nk-1
            dEdk(i) = real(Ec(i+1)-Ec(i-1)) / 2d0 / dkk
        end do
        dEdk(1)    = 2*dEdk(2 ) - dEdk(3)
        dEdk(Nk)   = 2*dEdk(Nk-1) - dEdk(Nk-2)

        v = sum(dEdk(:) * n(:)) / (1d-100 + sum(n(:))) / hbar
    end function DriftVt

    

    function FDrift2(Ephn, m, g, ky, n, Cq2, v, N0, x)
        real(dp), intent(in) :: Ephn, m, g             ! Average phonon energy, carrier mass
        real(dp), intent(in) :: ky(:), n(:), Cq2(:)    ! K-space array, carrier occupation #, carrier constant
        real(dp), intent(in) :: v, N0                 ! 
        real(dp), intent(in) :: x(:)
        real(dp)             :: FDrift2(size(ky))
        real(dp)             :: EM(  size(ky), size(ky))
        real(dp)             :: ABSB(size(ky), size(ky))
        integer              :: k, Nk, q

        Nk = size(ky)

        !$omp parallel do private(q, k)
        do k = 1, Nk
            do q = 1, Nk
                EM(q,k)   =  ThetaEM( Ephn, m, g, ky, n, Cq2, v, N0, q, k)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(q, k)
        do k = 1, Nk
            do q = 1, Nk        
                ABSB(q,k) =  ThetaABS(Ephn, m, g, ky, n, Cq2, v, N0, q, k)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(k)
        do k = 1, Nk
            FDrift2(k) = sum( hbar*ky(:) * (EM(:,k) - ABSB(:,k)) )
        end do
        !$omp end parallel do
        
!        FDrift2 = sum(FDrift2) / (sum(n)+1d-10) * 2d0

       !FDrift2 = sum( hbar*q(:) * (EM(:)*(N0+1) - ABSB(:)*N0) * x(:))
    end function FDrift2

    

    function FDrift(Ephn, m, q, dndk, Cq2, v, N0, x)
        real(dp), intent(in) :: Ephn, m
        real(dp), intent(in) :: q(:), dndk(:), Cq2(:)
        real(dp), intent(in) :: v, N0
        real(dp), intent(in) :: x(:)
        real(dp)             :: FDrift
        real(dp)             :: dq

        dq = q(2)-q(1)

        FDrift = sum( hbar*q(:) * ThetaEMABS(Ephn, m, q, dndk, Cq2, v) * (2*N0 + 1) * x(:)) / 2d0
    end function FDrift


    
    function ThetaEM(Ephn, m, g, ky, n, Cq2, v, N0, q, k)
        real(dp), intent(in) :: Ephn, m, g
        real(dp), intent(in) :: ky(:), n(:), Cq2(:)
        real(dp), intent(in) :: v, N0
        integer,  intent(in) :: q, k
        real(dp)             :: ThetaEM
        real(dp)             :: xq, dk, Ek, Ekmq
        integer              :: kmq, Nk0, Nk

        ThetaEM = 0d0

        ! Find the central index for ky-array, where ky(Nk0) = 0.0
        Nk0 = ceiling(size(ky)/2d0)

        ! Size of ky-array
        Nk  = size(ky)

        ! Get ky-array differential
        dk  = ky(2)-ky(1)

        !  Find the k-q array index within the ky-array
        kmq = Nk0 + nint(ky(k)/dk) - nint(ky(q)/dk)

        !  We don't need to do this calculation if k-q is out of bounds.
        !  So if it is, srew it!  Just return zero.
        if(kmq < 1 .or. kmq > Nk) return

        !  Short-hand variables for a common arguments
        xq    = Ephn  -  hbar * ky(q) * v
        Ek    = hbar**2 * ky(k  )**2 / 2d0 / m
        Ekmq  = hbar**2 * ky(kmq)**2 / 2d0 / m

        ThetaEM = 4*pi/hbar * Cq2(q) * n(k) * (1d0 - n(kmq)) * (N0+1d0) * &
                  Lrtz(Ekmq - Ek + xq, hbar * g) * theta(xq)
    end function ThetaEM



    function ThetaABS(Ephn, m, g, ky, n, Cq2, v, N0, q, k)
        real(dp), intent(in) :: Ephn, m, g
        real(dp), intent(in) :: ky(:), n(:), Cq2(:)
        real(dp), intent(in) :: v, N0
        integer,  intent(in) :: q, k
        real(dp)             :: ThetaABS
        real(dp)             :: xq, dk, Ek, Ekmq
        integer              :: kmq, Nk0, Nk

        ThetaABS = 0d0

        ! Find the central index for ky-array, where ky(Nk0) = 0.0
        Nk0 = ceiling(size(ky)/2d0)

        ! Size of ky-array
        Nk  = size(ky)

        ! Get ky-array differential
        dk  = ky(2)-ky(1)

        !  Find the k-q array index within the ky-array
        kmq = Nk0 + nint(ky(k)/dk) - nint(ky(q)/dk)

        !  We don't need to do this calculation if k-q is out of bounds.
        !  So if it is, srew it!  Just return zero.
        if(kmq < 1 .or. kmq > Nk) return

        !  Short-hand variables for a common arguments
        xq    = Ephn  -  hbar * ky(q) * v
        Ek    = hbar**2 * ky(k  )**2 / 2d0 / m
        Ekmq  = hbar**2 * ky(kmq)**2 / 2d0 / m
        
        ThetaABS = 4*pi/hbar * Cq2(q) * n(kmq) * (1d0 - n(k)) * N0 * &
                   Lrtz(Ek - Ekmq - xq, hbar * g) * theta(xq)
    end function ThetaABS



    subroutine CalcAvgCoeff(ky, dk, k1, k2, i1, i2, x1, x2, x3, x4)
        real(dp), intent(in   ) :: ky(:), dk, k1, k2
        integer,  intent(in   ) :: i1, i2
        real(dp), intent(inout) :: x1, x2, x3, x4
        real(dp)                :: k(-1:size(ky)+2)
        integer                 :: Nk

        Nk = size(ky)
        
        k = 0d0
        k(1:size(ky)) = ky(:)
        
        k(0 ) = ky(1) - dk
        k(-1) = ky(1) - 2*dk

        k(Nk+1) = ky(Nk) + dk
        k(Nk+2) = ky(Nk) + 2*dk

        x1 = ( k(i1+1) - k1    ) * ( k(i2+1) - k2    ) / dk**2
        x2 = ( k1      - k(i1) ) * ( k(i2+1) - k2    ) / dk**2
        x3 = ( k(i1+1) - k1    ) * ( k2      - k(i2) ) / dk**2        
        x4 = ( k1      - k(i1) ) * ( k2      - k(i2) ) / dk**2



    end subroutine
    

    function ThetaEMABS(Ephn, m, q, dndk, Cq2, v)
        real(dp), intent(in) :: Ephn, m
        real(dp), intent(in) :: q(:), dndk(:), Cq2(:)
        real(dp), intent(in) :: v
        real(dp)             :: ThetaEMABS(size(q))

        ThetaEMABS(:) = 4*pi/hbar * Cq2(:) * hbar*q(:) * v * (-dndEk(Ephn, m, q(:), dndk(:))) 
    end function ThetaEMABS
    
    
    function dndEk(Ephn, m, q, dndq)
        real(dp), intent(in) :: Ephn, m
        real(dp), intent(in) :: q(:), dndq(:)
        real(dp)             :: dndEk(size(q))
        real(dp)             :: x0

        x0 = Ephn * 2 * m / hbar**2
        
        dndEk(:) = dndq(:) * 4 * q(:)**3 / (q(:)**4 - x0**2) * m / hbar**2
    end function dndEk
    
    
    function CalcVD(ky, m, n)
        real(dp)    :: ky(:), m
        complex(dp) :: n(:)
        real(dp)    :: CalcVd
    
        CalcVd = sum( real(n) * hbar * ky / m) / sum(real(n) + small)
    end function
    
    
    function CalcPD(ky, m, n)
        real(dp)    :: ky(:), m
        complex(dp) :: n(:)
        real(dp)    :: CalcPd
    
        CalcPd = sum( abs(n) * hbar * ky) / sum(abs(n) + small)
    end function



    function GetEDrift()
        real(dp) :: GetEDrift
        GetEDrift = ERate
    end function


    function GetHDrift()
        real(dp) :: GetHDrift
        GetHDrift = HRate
    end function


    function GetVEDrift()
        real(dp) :: GetVEDrift
        GetVEDrift = VEDrift
    end function


    function GetVHDrift()
        real(dp) :: GetVHDrift
        GetVHDrift = VHDrift
    end function
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!    Old Legacy Code   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






   subroutine DC_Step_Scale(ne, nh, ky, Edc, dt)
        complex(dp), intent(inout) :: ne(:), nh(:)
        real(dp),    intent(in   ) :: ky(:), Edc, dt
        real(dp)                   :: dky
        complex(dp)                :: ned(size(ne)) , nhd(size(nh))
        complex(dp)                :: ned2(size(ne)), nhd2(size(nh))
        integer                    :: s, Nk, k

        Nk  = size(ky)
        dky = ky(3)-ky(2)

        ned = ne
        nhd = nh
        call rescale_1D(ky - e0*Edc/hbar*dt, ned, ky, ne)
        call rescale_1D(ky + e0*Edc/hbar*dt, ned, ky, nh)
    end subroutine
    
    
    subroutine DC_Step_FD(ne, nh, nemid, nhmid, ky, Edc, dt, me, mh)
        complex(dp), intent(inout) :: ne(:), nh(:)
        complex(dp), intent(in   ) :: nemid(:), nhmid(:)
        real(dp),    intent(in   ) :: ky(:), Edc, dt, me, mh
        real(dp)                   :: dky
        complex(dp)                :: neu(size(ne)) , nhu(size(nh))
        complex(dp)                :: ned(size(ne)),  nhd(size(nh))
        integer                    :: s, Nk, k

        Nk  = size(ky)
        dky = ky(3)-ky(2)
        s   = nint(Edc / (Edc+small))
        
        ned = cshift(nemid,  1)
        nhd = cshift(nhmid,  1)
        
        neu = cshift(nemid, -1)
        nhu = cshift(nhmid, -1)
        

        ne = ne + (- e0 * Edc - CalcPD(ky, me, ne) * 1d13) / hbar * (neu-ned) / dky * dt
        nh = nh + (+ e0 * Edc - CalcPD(ky, mh, nh) * 1d13) / hbar * (nhu-nhd) / dky * dt
    end subroutine


    
    subroutine ShiftN1D(ne, dk)
        complex(dp), intent(inout) :: ne(:)
        real(dp),    intent(in   ) :: dk
        integer                    :: k, Nk
        Nk = size(ne)
        
        call FFT(ne)
            ne = ne * exp(-ii * Y * dk)
        call IFFT(ne)
    end subroutine ShiftN1D



    subroutine ShiftN2D(C, dk)
        complex(dp), intent(inout) :: C(:,:)
        real(dp),    intent(in   ) :: dk
        integer                    :: k1, k2, Nk
        real(dp)                   :: Y2(size(Y))

        Nk = size(C,1)
        Y2 = Y
        
        call FFT(C)
        !$omp parallel do private(k1, k2)
        do k1=1, Nk
            do k2=1, Nk
                C(k1,k2) = C(k1,k2) * exp(-ii * (Y(k1)+Y2(k2)) * dk)
            end do
        end do
        !$omp end parallel do
        call IFFT(C) 
    end subroutine ShiftN2D


    subroutine Transport(C , Edc, Eac, dt, DCTrans, k1nek2)
        complex(dp), intent(inout) :: C(:,:)
        real(dp),    intent(in   ) :: Edc, Eac, dt
        logical,     intent(in   ) :: DCTrans, k1nek2
        real(dp)                   :: dk = 0d0
        complex(dp)                :: nk(size(C,1))
        integer                    :: k

        if(.not.DCTrans) return

        dk = - e0 * (Edc + Eac) / hbar * dt

        if(k1nek2) then
            call ShiftN2D(C,dk)
        else
            forall(k=1:size(C,1)) nk(k) = C(k,k)
            call ShiftN1D(nk,dk)
            forall(k=1:size(C,1)) C(k,k) = nk(k)
        end if
    end subroutine


    
end module
