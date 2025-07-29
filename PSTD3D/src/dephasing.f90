! This module calculates dephasing rates for the SBEs
! in support of propagation simulations for a quantum wire.
module dephasing

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use helpers
    use usefulsubs
    implicit none


    ! Arrays needed for dephasing calculations
    integer, private, allocatable :: k_p_q(:,:)        ! k+q array
    integer, private, allocatable :: k_m_q(:,:)        ! k-q array
    integer, private, allocatable :: k1_m_q(:,:)       ! k1+q array for h-dephasing
    integer, private, allocatable :: k1p_m_q(:,:)      ! k1p+q array for e-dephasing
    integer, private, allocatable :: k1(:,:)           ! k1 array for eh-dephasing of holes
    integer, private, allocatable :: k1p(:,:)          ! k1p array for eh-dephasing of electrons
    real(dp),private, allocatable :: xe(:), xh(:)      ! coeffients for delta functions


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!    Calculations for SBEs.f90     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    


    subroutine InitializeDephasing(ky, me, mh)
        real(dp),    intent(in) :: ky(:), me, mh
        real(dp)                :: dk, kmax, kmin, k1_0, k1p_0
        integer                 :: Nk, Nk0, k, kp, q

        Nk = size(ky)
        dk = ky(2)-ky(1)

        kmax = ky(Nk) + dk
        kmin = ky(1)  - dk

       !Nk0  = ceiling(Nk/2d0)
        Nk0  = (Nk-1)/2 + 1

        allocate(k_p_q(Nk,Nk))
        allocate(k_m_q(Nk,Nk))
        allocate(k1_m_q(Nk,Nk))
        allocate(k1p_m_q(Nk,Nk))
        allocate(k1(Nk,Nk))
        allocate(k1p(Nk,Nk))

        allocate(xe(Nk), xh(Nk))

        xe =  me / hbar**2 * abs(ky) / (abs(ky)+dk*1d-5)**2 / dk
        xh =  mh / hbar**2 * abs(ky) / (abs(ky)+dk*1d-5)**2 / dk
        
        k_p_q   = 0d0
        k_m_q   = 0d0
        k1_m_q  = 0d0
        k1p_m_q = 0d0
        k1      = 0d0
        k1p     = 0d0
        
        do q=1, Nk
            do k=1, Nk
                k_p_q(k,q) = nint( max( min(ky(k) + ky(q), kmax),  kmin ) / dk ) + Nk0
                k_m_q(k,q) = nint( max( min(ky(k) - ky(q), kmax),  kmin ) / dk ) + NK0
            end do
        end do

        do q=1, Nk
            do k=1, Nk
                k1p_0  = ( (me + mh) * ky(q) - 2*mh * ky(k) ) / 2d0 / me

                k1p(k,q)     = nint( max( min(k1p_0        , kmax),  kmin ) / dk ) + Nk0
                k1p_m_q(k,q) = nint( max( min(k1p_0 - ky(q), kmax),  kmin ) / dk ) + Nk0
            end do
        end do

        do q=1, Nk
            do kp=1, Nk
                k1_0  = ( (me + mh) * ky(q) - 2*me * ky(kp) ) / 2d0 / mh

                k1(kp,q)     = nint( max( min(k1_0        , kmax),  kmin ) / dk ) + Nk0
                k1_m_q(kp,q) = nint( max( min(k1_0 - ky(q), kmax),  kmin ) / dk ) + Nk0
            end do
        end do
        
        open(unit=383,file="dataQW/Wire/info/MaxOffDiag.e.dat")
        open(unit=384,file="dataQW/Wire/info/MaxOffDiag.h.dat")

    end subroutine 
    

    subroutine CalcGammaE(ky, ne0, nh0, VC, GammaE)
        real(dp),    intent(in)    :: ky(:)
        complex(dp), intent(in)    :: ne0(:), nh0(:)
        real(dp),    intent(in)    :: VC(:,:,:)
        real(dp),    intent(inout) :: GammaE(:)
        real(dp)                   :: Veh2(size(ky))
        real(dp)                   :: Vee2(size(ky))
        real(dp)                   :: ne(0:size(ky)+1), nh(0:size(ky)+1)
        real(dp)                   :: se(0:size(ky)+1), sh(0:size(ky)+1)
        real(dp)                   :: dk
        integer                    :: Nk, k, q, kp
        
        Veh2(:) = Vxx2(ky, VC(:,:,1))
        Vee2(:) = Vxx2(ky, VC(:,:,2))
        
        Nk = size(ky)
        dk = ky(2)-ky(1)
                
        ne     = 0d0
        nh     = 0d0
        GammaE = 0d0

        ne(1:Nk) = real(ne0(:))
        nh(1:Nk) = real(nh0(:))

        se(1:Nk)    = 1d0 - ne(1:Nk)
        sh(1:Nk)    = 1d0 - nh(1:Nk)

        ! Electron-Electron dephasing of electons
        ! (fully simplified by delta functions and parabolic bands)
        do q=1, Nk
            do k=1, Nk
                GammaE(k) = GammaE(k) + pi/hbar * Vee2(q) *  ne(k_p_q(k,q))  * se(k_p_q(k,q)) * abs(xe(q))
            end do
        end do

        ! Electron-Hole dephasing of electrons
        ! (simplified by delat functions and parabolic bands)
        do q=1, Nk
            do k=1, Nk
                GammaE(k) = GammaE(k) + pi/hbar * Veh2(q) * &
                            ( + nh(k1p_m_q(k,q)) * sh(k1p(k,q)) * ne(k_m_q(k,q)) &
                              + sh(k1p_m_q(k,q)) * nh(k1p(k,q)) * se(k_m_q(k,q)) &
                            ) * abs(xh(q))
            end do
        end do        
    end subroutine






    subroutine CalcGammaH(ky, ne0, nh0, VC, GammaH)
        real(dp),    intent(in)    :: ky(:)
        complex(dp), intent(in)    :: ne0(:), nh0(:)
        real(dp),    intent(in)    :: VC(:,:,:)
        real(dp),    intent(inout) :: GammaH(:)
        real(dp)                   :: Veh2(size(ky))
        real(dp)                   :: Vhh2(size(ky))
        real(dp)                   :: ne(0:size(ne0)+1), nh(0:size(ne0)+1)
        real(dp)                   :: se(0:size(ne0)+1), sh(0:size(ne0)+1)
        real(dp)                   :: dk
        integer                    :: Nk, k, q, kp
        
        Veh2(:) = Vxx2(ky, VC(:,:,1))
        Vhh2(:) = Vxx2(ky, VC(:,:,3))
        
        Nk      = size(ky)
        dk      = ky(2)-ky(1)
        
        ne      = 0d0
        nh      = 0d0
        GammaH  = 0d0

        ne(1:Nk) = real(ne0(:))
        nh(1:Nk) = real(nh0(:))

        se(1:Nk)    = 1d0 - ne(1:Nk)
        sh(1:Nk)    = 1d0 - nh(1:Nk)
        
        ! Electron-Electron dephasing of electons
        ! (fully simplified by delta functions and parabolic bands)
        do q=1, Nk
            do kp=1, Nk
                GammaH(kp) = GammaH(kp) + pi/hbar * Vhh2(q) *  nh(k_p_q(kp,q))  * sh(k_p_q(kp,q)) * abs(xh(q))
            end do
        end do

        ! Electron-Hole dephasing of electrons
        ! (simplified by delat functions and parabolic bands)
        do q=1, Nk
            do kp=1, Nk
                GammaH(kp) = GammaH(kp) + pi/hbar * Veh2(q) * &
                             ( + ne(k1_m_q(kp,q)) * se(k1(kp,q)) * nh(k_m_q(kp,q)) &
                               + se(k1_m_q(kp,q)) * ne(k1(kp,q)) * sh(k_m_q(kp,q)) &
                             ) * abs(xe(q))
            end do
        end do
    end subroutine 



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!! BEGIN OFF DIAGONAL DEPHASING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine OffDiagDephasing(ne, nh, p, ky, Ee, Eh, g, VC, x)
        complex(dp), intent(in   ) :: ne(:), nh(:)
        complex(dp), intent(in   ) :: p(:,:)
        real(dp)   , intent(in   ) :: ky(:), Ee(:), Eh(:), g(:)
        real(dp),    intent(in   ) :: VC(:,:,:)
        complex(dp), intent(inout) ::  x(size(ky),size(ky))
        complex(dp)                :: De(size(ky),size(ky))
        complex(dp)                :: Dh(size(ky),size(ky))
        complex(dp)                :: pp(0:size(ky)+1,0:size(ky)+1)
        complex(dp)                :: pt(0:size(ky)+1,0:size(ky)+1)
        integer                    :: k, kp, q, qp, Nk
        real(dp)                   :: undel(size(ky))

        x  = 0d0

        Nk = size(ky)
        undel = abs(ky) / (abs(ky)+1d-10)

        De = 0d0
        Dh = 0d0
        pp = 0d0
        pt = 0d0

        De = CalcOffDiagDeph_E(ne, nh, ky, Ee, Eh, g(1), g(3), VC)
        Dh = CalcOffDiagDeph_H(ne, nh, ky, Ee, Eh, g(2), g(3), VC)
        De = transpose(De)
        Dh = transpose(Dh)

        ! Note that the "new" code transposed p from p(ke,kh) to p(kh,ke)
        ! So the old code commended below was replace with an unintuitive version (Feb 08 2023)
        !pp(1:Nk,1:Nk) = p
        !pt(1:Nk,1:Nk) = transpose(p)
         pt(1:Nk,1:Nk) = p
         pp(1:Nk,1:Nk) = transpose(p)






        ! NOTE IN EVERYTHING BELOW THAT THE MATRIX k_p_q IS THE TRANSPOSE OF ITSELF!!!!!
       

        !$omp parallel do private(k, kp, qp)
        do k=1, Nk
            do kp=1, Nk
                do qp=1, Nk
                    x(kp,k) = x(kp,k) + Dh(qp,kp) * pt(k_p_q(qp,kp), k) * undel(qp)
                end do
            end do
        end do
        !$omp end parallel do
        
        x = transpose(x)

        !$omp parallel do private(kp, k, q)
        do kp=1, Nk
            do k=1, Nk
                do q=1, Nk
                    x(k,kp) = x(k,kp) + De(q,k) * pp(k_p_q(q,k),kp) * undel(q)
                end do
            end do
        end do
        !$omp end parallel do
        
        x =  ii * hbar * x
    end subroutine
    
    
    
    Function CalcOffDiagDeph_E(ne0, nh0, ky, Ee0, Eh0, gee, geh, VC) result(D)
        complex(dp), intent(in   ) :: ne0(:), nh0(:)
        real(dp),    intent(in   ) :: ky(:), Ee0(:), Eh0(:), gee, geh
        real(dp),    intent(in   ) :: VC(:,:,:)
        real(dp)                   :: D(size(ky),size(ky))
        real(dp)                   :: ne(0:size(ne0)+1)
        real(dp)                   :: nh(0:size(nh0)+1)
        real(dp)                   :: Veh2(0:size(ne)+1, 0:size(nh)+1)
        real(dp)                   :: Vee2(0:size(ne)+1, 0:size(nh)+1)
        real(dp)                   :: Ee(0:size(ky)+1), Eh(0:size(ky)+1)
        integer                    :: k, q, k1, kmq, Nk
        integer                    :: kpq, k1mq, k1pq
    
        D    = 0d0
        Veh2 = 0d0
        Vee2 = 0d0
        ne   = 0d0
        nh   = 0d0
        Ee   = 0d0
        Eh   = 0d0
        
        Nk   = size(ne0)
        
        Veh2(1:Nk,1:Nk) = VC(:,:,1)**2
        Vee2(1:Nk,1:Nk) = VC(:,:,2)**2
        
        ne(1:Nk) = abs(ne0)
        nh(1:Nk) = abs(nh0)
        Ee(1:Nk) = Ee0
        Eh(1:Nk) = Eh0
        
        !$omp parallel do private(q, k, kpq, k1, k1pq)
        do q=1, Nk
            do k=1, Nk
                kpq = k_p_q(k,q)
                do k1=1, Nk
                    k1pq = k_p_q(k1,q)
            
                    D(k,q) = D(k,q) + Vee2(k1,k1pq) * Lrtz(Ee(k1pq)+Ee(k)-Ee(k1)-Ee(kpq), hbar*gee) &
                                    * ( ne(k1pq) * ne(k) * (1d0-ne(k1)) + (1d0-ne(k1pq)) * (1d0-ne(k)) * ne(k1) )
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(q, k, kpq, k1, k1mq)
        do q=1, Nk
            do k=1, Nk
                kpq = k_p_q(k,q)
                do k1=1, Nk
                    k1mq = k_m_q(k1,q)
           
                    D(k,q) = D(k,q) + Veh2(k,kpq) * Lrtz(Eh(k1mq)+Ee(k)-Eh(k1)-Ee(kpq), hbar*geh) &
                                    * ( nh(k1mq) * (1d0-nh(k1)) * ne(k) + (1d0-nh(k1mq)) * nh(k1) * (1d0-ne(k)) )
                end do
            end do
        end do
        !$omp end parallel do

        D = D * pi / hbar
    end Function



    Function CalcOffDiagDeph_H(ne0, nh0, ky, Ee0, Eh0, ghh, geh, VC) result(D)
        complex(dp), intent(in   ) :: ne0(:), nh0(:)
        real(dp),    intent(in   ) :: ky(:), Ee0(:), Eh0(:), ghh, geh
        real(dp),    intent(in   ) :: VC(:,:,:)
        real(dp)                   :: D(size(ky),size(ky))
        real(dp)                   :: ne(0:size(ne0)+1)
        real(dp)                   :: nh(0:size(nh0)+1)
        real(dp)                   :: Veh2(0:size(ne)+1, 0:size(nh)+1)
        real(dp)                   :: Vhh2(0:size(ne)+1, 0:size(nh)+1)
        real(dp)                   :: Ee(0:size(ky)+1), Eh(0:size(ky)+1)
        integer                    :: k, q, k1, kmq, Nk
        integer                    :: kpq, k1mq, k1pq
    
        D    = 0d0
        Veh2 = 0d0
        Vhh2 = 0d0
        ne   = 0d0
        nh   = 0d0
        Ee   = 0d0
        Eh   = 0d0
        
        Nk   = size(ne0)
        
        Veh2(1:Nk,1:Nk) = VC(:,:,1)**2
        Vhh2(1:Nk,1:Nk) = VC(:,:,3)**2
        
        ne(1:Nk) = abs(ne0)
        nh(1:Nk) = abs(nh0)
        Ee(1:Nk) = Ee0
        Eh(1:Nk) = Eh0


        !$omp parallel do private(q, k, kpq, k1, k1pq)        
        do q=1, Nk
            do k=1, Nk
                kpq = k_p_q(k,q)
                do k1=1, Nk

                    D(k,q) = D(k,q) + Vhh2(k1,k1pq) * Lrtz(Eh(k1pq)+Eh(k)-Eh(k1)-Eh(kpq), hbar*ghh) &
                                    * ( nh(k1pq) * nh(k) * (1d0-nh(k1)) + (1d0-nh(k1pq)) * (1d0-nh(k)) * nh(k1) )
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(q, k, kpq, k1, k1mq)        
        do q=1, Nk
            do k=1, Nk
                kpq = k_p_q(k,q)
                do k1=1, Nk
                    k1mq = k_m_q(k1,q)
           
                    D(k,q) = D(k,q) + Veh2(k1,k1mq) * Lrtz(Ee(k1mq)+Eh(k)-Ee(k1)-Eh(kpq), hbar*geh) &
                                    * ( ne(k1mq) * (1d0-ne(k1)) * nh(k) + (1d0-ne(k1mq)) * ne(k1) * (1d0-nh(k)) )
                end do
            end do
        end do
        !$omp end parallel do
    
        D = D * pi / hbar 
    end Function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     ATTEMPT #2 FOR OFF-DIAGONAL DEPHASING
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    subroutine OffDiagDephasing2(ne, nh, p, ky, Ee, Eh, g, VC, t, x)
        complex(dp), intent(in   ) :: ne(:), nh(:)
        complex(dp), intent(in   ) :: p(:,:)
        real(dp)   , intent(in   ) :: ky(:), Ee(:), Eh(:), g(:)
        real(dp),    intent(in   ) :: VC(:,:,:), t
        complex(dp), intent(inout) ::  x(size(ky),size(ky))
        real(dp)                   :: De(-size(ky):size(ky),size(ky))
        real(dp)                   :: Dh(-size(ky):size(ky),size(ky))
        complex(dp)                :: pt(size(ky),size(ky))
        integer                    :: ke, kh, q, Nk
        real(dp)                   :: undel(-size(ky):size(ky))

        x  = 0d0
        Nk = size(ky)
       
        undel(-Nk:Nk)    = 1d0
        undel(0) = 0

        Nk = size(ky)
        !undel = abs(ky) / (abs(ky)+1d-10)

        De = 0d0
        Dh = 0d0


        De = CalcOffDiagDeph_E2(ne, nh, ky, Ee, Eh, g(1), g(3), VC, Nk)
        Dh = CalcOffDiagDeph_H2(ne, nh, ky, Ee, Eh, g(2), g(3), VC, Nk)

        pt = transpose(p)

        ! NOTE IN EVERYTHING BELOW THAT THE MATRIX k_p_q IS THE TRANSPOSE OF ITSELF!!!!!

        !$omp parallel do private(ke, kh, q)
        do ke=1, Nk
            do kh=1, Nk
                do q=1-kh, Nk-kh
                    x(kh,ke) = x(kh,ke) + Dh(q,kh) * p(kh+q, ke) * undel(q)
                end do
            end do
        end do
        !$omp end parallel do


        !$omp parallel do private(ke, kh, q)
        do ke=1, Nk
            do kh=1, Nk
                do q=1-ke, Nk-ke
                    x(kh,ke) = x(kh,ke) + De(q,ke) * pt(ke+q,kh) * undel(q)
                end do
            end do
        end do
        !$omp end parallel do

        write(383,*) real(t), real(MaxVal(De)), real(MinVal(De))
        write(384,*) real(t), real(MaxVal(Dh)), real(MinVal(Dh))
        
        x =  ii * hbar * x
    end subroutine




    Function CalcOffDiagDeph_E2(ne, nh, ky, Ee, Eh, gee, geh, VC, Nk) result(D)
        complex(dp), intent(in   ) :: ne(:), nh(:)
        real(dp),    intent(in   ) :: ky(:), Ee(:), Eh(:), gee, geh
        real(dp),    intent(in   ) :: VC(:,:,:)
        integer,     intent(in   ) :: Nk
        real(dp)                   :: D(-Nk:Nk, Nk)
        integer                    :: k, q, p
        real(dp)                   :: Vsq(-Nk:Nk,3)

        Vsq  = 0d0


        ! Set up 1D Vsq calculation.
        ! Note that Vsq = 0 for q = 0
        forall(q=1:Nk-1) Vsq(+q,:) = VC(1+q,1,:)**2
        forall(q=1:Nk-1) Vsq(-q,:) = VC(1+q,1,:)**2
    
        D    = 0d0

        !$omp parallel do private(k, p, q)
        do k=1, Nk
            do p=1, Nk
                do q=max(p-Nk,1-k), min(p-1,Nk-k)
                    D(q,k) = D(q,k) + Vsq(q,2) * Lrtz(Ee(k+q)+Ee(p-q)-Ee(p)-Ee(k), hbar*gee) &
                                    * real( ne(p-q) * (1-ne(p)) * (1-ne(k))  +  (1-ne(p-q)) * ne(p) * ne(k) )
                end do
            end do
        end do
        !$omp end parallel do


        !$omp parallel do private(k, p, q)
        do k=1, Nk
            do p=1, Nk
                do q=max(1-p,1-k), min(Nk-p,Nk-k)
                    D(q,k) = D(q,k) + Vsq(q,1) * Lrtz(Ee(k+q)+Eh(p+q)-Eh(p)-Ee(k), hbar*geh) &
                                    * real( nh(p+q) * (1-nh(p)) * (1-ne(k))  +  (1-nh(p+q)) * nh(p) * ne(k) )
                end do
            end do
        end do
        !$omp end parallel do
    
        D = D * pi / hbar
    end Function



    Function CalcOffDiagDeph_H2(ne, nh, ky, Ee, Eh, ghh, geh, VC, Nk) result(D)
        complex(dp), intent(in   ) :: ne(:), nh(:)
        real(dp),    intent(in   ) :: ky(:), Ee(:), Eh(:), ghh, geh
        real(dp),    intent(in   ) :: VC(:,:,:)
        integer,     intent(in   ) :: Nk
        real(dp)                   :: D(-Nk:Nk, Nk)
        integer                    :: k, q, p
        real(dp)                   :: Vsq(-Nk:Nk,3)


        Vsq  = 0d0

        ! Set up 1D Vsq calculation.
        ! Note that Vsq = 0 for q = 0
        forall(q=1:Nk-1) Vsq(+q,:) = VC(1+q,1,:)**2
        forall(q=1:Nk-1) Vsq(-q,:) = VC(1+q,1,:)**2
    
        D    = 0d0

        !$omp parallel do private(k, p, q)
        do k=1, Nk
            do p=1, Nk
                do q=max(p-Nk,1-k), min(p-1,Nk-k)
                    D(q,k) = D(q,k) + Vsq(q,3) * Lrtz(Eh(k+q)+Eh(p-q)-Eh(p)-Eh(k), hbar*ghh) &
                                    * real( nh(p-q) * (1-nh(p)) * (1-nh(k))  +  (1-nh(p-q)) * nh(p) * nh(k) )
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(k, p, q)
        do k=1, Nk
            do p=1, Nk
                do q=max(1-p,1-k), min(Nk-p,Nk-k)
                    D(q,k) = D(q,k) + Vsq(q,1) * Lrtz(Eh(k+q)+Ee(p+q)-Ee(p)-Eh(k), hbar*geh) &
                                    * real( ne(p+q) * (1-ne(p)) * (1-nh(k))  +  (1-ne(p+q)) * ne(p) * nh(k) )
                end do
            end do
        end do
        !$omp end parallel do
    
        D = D * pi / hbar
    end Function






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!! END OFF DIAGONAL DEPHASING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    
    ! Calculates Cq for use in the DC Field module
    function Vxx2(q, V)
        real(dp), intent(in) :: q(:), V(:,:)
        real(dp)             :: Vxx2(size(q))
        real(dp)             :: dq
        integer              :: i, iq(size(q))

        dq = q(2) - q(1)
        
        iq = nint( abs(q / dq) )
        
        do i=1, size(q)
            Vxx2(i) = V(1 + iq(i), 1)**2
        end do
    end function Vxx2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   Printing statements for SBEs.f90     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine WriteDephasing(ky, gamE, gamH, w, xxx)
        real(dp),         intent(in) :: ky(:)      ! QW momemtum
        real(dp),         intent(in) :: gamE(:)    ! Elec diagonal dephasing rate
        real(dp),         intent(in) :: gamH(:)    ! Elec diagonal dephasing rate
        integer,          intent(in) :: w          ! Wire index
        integer,          intent(in) :: xxx        ! Time index
        character(len=50)            :: fmt, wire
        
        fmt = '(I2.2)'
        write(wire,fmt) w

        call printGam( GamE, ky, xxx, 'Wire/Ge/Ge.'//trim(wire)//'.k.')
        call printGam( GamH, ky, xxx, 'Wire/Gh/Gh.'//trim(wire)//'.k.')
    end subroutine
    
    

    subroutine printGam(Dx, z, n, file)
        real(dp), intent(in)    :: Dx(:)
        real(dp),    intent(in) :: z(:)
        integer                 :: n
        character(len=*)        :: file 
        character(len=50)       :: filename, fmt
        integer                 :: i, u



            fmt = '(I5.5)'

            write(filename,fmt) n
            filename = 'dataQW/'//trim(file)//trim(filename)//'.dat'

            u = n+20

            open(unit=u, file=trim(filename))

            do i=1, size(z)
                write(u,*) sngl(z(i)), sngl(Dx(i))
            end do

            close(u)


    end subroutine
    
end module
