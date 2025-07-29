! This module calculates the electron-hole, electron-electron,
! and hole-hole collision integrals and other carrier-optical
! related calculations required for the Semiconductor Bloch
! Equations module (SBEs.f90) in support of simulations of pulse
! propagation (pstd.f90) through a quantum wire.
module coulomb 

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers
    use spliner
    use usefulsubs
    implicit none

    ! Pre-calculated QW arrays required for solving SBEs.
    real(dp),    private, allocatable :: Vee0(:,:)  ! Unscreened Electron-Electron collision arrays (J)
    real(dp),    private, allocatable :: Vhh0(:,:)  ! Unscreened Hole-hole collision arrays (J)
    real(dp),    private, allocatable :: Veh0(:,:)  ! Unscreened Electron-hole collision arrays (J)
    real(dp),    private, allocatable :: UnDel(:,:) ! Reverse Delta function array (1 - delta_ij)
    real(dp),    private, allocatable :: Chi1De(:,:) ! The screening 1D-Chi for QW electrons
    real(dp),    private, allocatable :: Chi1Dh(:,:) ! The screening 1D-Chi for QW holes

    
    integer,     private, allocatable :: k3(:,:,:)   ! k3 = k1+k2-k4 indexing array
    real(dp),    private, allocatable :: Ceh(:,:,:)  ! 
    real(dp),    private, allocatable :: Ceh2(:,:,:) ! 
    real(dp),    private, allocatable :: Cee(:,:,:)  ! 
    real(dp),    private, allocatable :: Chh(:,:,:)  !
    real(dp),    private, allocatable :: qe(:,:)
    real(dp),    private, allocatable :: qh(:,:)

    ! Initial Numerical parameters
    real(dp),    private              :: small = 1d-200  ! Smallest # worthy of consideration

    ! Should we read in precalculated potential arrays (Vee, Vhh, and Veh)?
    logical,     private              :: ReadArrays = .false.
    logical,     private              :: ScrewThis = .false.
    logical,     private              :: LorentzDelta=.false.

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! Construction of Unscreened Coulomb Arrays !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    ! Initialize this module and all of its needed quantities
    subroutine InitializeCoulomb(y, ky, L, Delta0, me, mh, Ee, Eh, ge, gh, alphae, alphah, er, Qy, kkp, screened)
        real(dp), intent(in) :: y(:), ky(:)    ! Length and momentum coordinates of QW (m, 1/m)
        real(dp), intent(in) :: L, Delta0      ! The length and thickness of the QW (m)
        real(dp), intent(in) :: me, mh         ! Effective electron/hole masses (kg)
        real(dp), intent(in) :: Ee(:), Eh(:)   ! Electron and Hole energies (J)
        real(dp), intent(in) :: ge, gh         ! Inverse electron and hole lifetives (Hz)
        real(dp), intent(in) :: alphae, alphah ! Level separation between ground and 1st excited state (1/m)
        real(dp), intent(in) :: er             ! Background dielectric constant
        real(dp), intent(in) :: Qy(:)          !
        integer,  intent(in) :: kkp(:,:)       !
        logical,  intent(in) :: screened      !

        ! Make the inverse delta function
        if(.not.allocated(UnDel)) call MakeUnDel(ky)

        ! Make the 3D index array for k3 = k1 + k2 - k4
        if(.not.allocated(k3))    call MakeK3(ky)

        ! Make the q, qe, and qh arrays
        if(.not.allocated(qe))    call MakeQs(ky, alphae, alphah)
        
        ! Calculate the 3D many-body interaction arrays
        if(.not.allocated(Ceh))   call CalcMBArrays(ky, Ee, Eh, ge, gh)

        ! Calculate the unscreened Coulomb collision arrays
        if(.not.allocated(Veh0))  call CalcCoulombArrays(y, ky, er, alphae, alphah, L, Delta0, Qy, kkp)

        ! Calculate the Chi1D array for Coulomb screening
        if(.not.allocated(Chi1De)) call CalcChi1D(ky, alphae, alphah, Delta0, L, er, me, mh) 
    end subroutine InitializeCoulomb

    
    ! Construct the unscreened Coulomb collision arrays
    subroutine CalcCoulombArrays(y, ky, er, alphae, alphah, L, Delta0, Qy, kkp)
        real(dp), intent(in) :: y(:), ky(:)    ! Length and momentum coordinates of QW (m, 1/m)
        real(dp), intent(in) :: er             ! Background dielectric constant 
        real(dp), intent(in) :: alphae, alphah ! Level separation between ground and 1st excited state (1/m)
        real(dp), intent(in) :: L, Delta0      ! The length and thickness of the QW (m)
        real(dp), intent(in) :: Qy(:)          !
        integer,  intent(in) :: kkp(:,:)       !
        integer              :: k, q, N, N0, NQ
        real(dp)             :: eh(size(Qy))
        real(dp)             :: ee(size(Qy))
        real(dp)             :: hh(size(Qy))



        N  = size(ky)
        N0 = (N-1)/2 + 1 !GetArray0Index(ky)
        NQ = size(Qy)
        
        allocate(Veh0(N,N), Vee0(N,N), Vhh0(N,N))
        
        ! Initialized to zero for now
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Veh0 = 0d0
        Vee0 = 0d0
        Vhh0 = 0d0
        eh   = 0d0
        ee   = 0d0
        hh   = 0d0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        if(ReadArrays) then
            call ReadIt(Veh0,'Veh')
            call ReadIt(Vee0,'Vee')
            call ReadIt(Vhh0,'Vhh')
            
        elseif(ScrewThis) then
            ! Do nothing
        else
            print*, "Calculating Coulomb Arrays"

            do k=1, NQ
                eh(k) = e0**2 / (twopi*eps0*er*L) * Vint(Qy(k), y, alphae, alphah, Delta0)
                ee(k) = e0**2 / (twopi*eps0*er*L) * Vint(Qy(k), y, alphae, alphae, Delta0)
                hh(k) = e0**2 / (twopi*eps0*er*L) * Vint(Qy(k), y, alphah, alphah, Delta0)
            end do
            
            
            do k=1, N
                do q=1, N
                    if(kkp(k,q) >= 0) then
                        Veh0(k,q) = eh(kkp(k,q))
                        Vee0(k,q) = ee(kkp(k,q))                   
                        Vhh0(k,q) = hh(kkp(k,q))
                    end if
                end do
            end do

            
            call WriteIt(Veh0,'Veh')
            call WriteIt(Vee0,'Vee')
            call WriteIt(Vhh0,'Vhh')

            print*, "Finished Calculating Unscreened Coulomb Arrays"
        end if
    end subroutine CalcCoulombArrays


    
    function Vehint(k, q, y, ky, alphae, alphah, Delta0)
        integer  :: k, q
        real(dp) :: y(:), ky(:)
        real(dp) :: alphae, alphah, Delta0
        real(dp) :: Vehint
        integer  :: i, j, Ny, N1, N2
        real(dp) :: aey2(size(y)), ahy2(size(y))
        real(dp) :: multconst, kmin, dk

        Vehint    = 0d0
        aey2      = (alphae * y)**2
        ahy2      = (alphah * y)**2
        kmin      = (alphae + alphah) / 4d0
        dk        = max(abs(ky(k) - ky(q)), kmin)
        
        multconst = alphae * alphah / pi * (y(2)-y(1))**2
        Ny        = size(y)
        N1        = Ny/4
        N2        = 3*Ny/4
        
        do i=N1, N2
            do j=N1, N2
                Vehint = Vehint + exp(- aey2(i) - ahy2(j)) * multconst * &
                         K03(dk * sqrt( (y(i)-y(j))**2 + Delta0**2 ) )
            end do
        end do
    end function Vehint
    

    function Vint(Qyk, y, alphae, alphah, Delta0)
        real(dp) :: Qyk                           ! Momentum Difference
        real(dp) :: y(:)
        real(dp) :: alphae, alphah, Delta0
        real(dp) :: Vint
        integer  :: i, j, Ny, N1, N2
        real(dp) :: aey2(size(y)), ahy2(size(y))
        real(dp) :: multconst, kmin, dk

        Vint      = 0d0
        aey2      = (alphae * y)**2
        ahy2      = (alphah * y)**2
        kmin      = (alphae + alphah) / 4d0
        dk        = max(abs(Qyk), kmin)
        
        multconst = alphae * alphah / pi * (y(2)-y(1))**2
        Ny        = size(y)
        N1        = Ny/4
        N2        = 3*Ny/4
        
        do i=N1, N2
            do j=N1, N2
                Vint = Vint + exp(- aey2(i) - ahy2(j)) * multconst * &
                       K03(dk * sqrt( (y(i)-y(j))**2 + Delta0**2 ) )
            end do
        end do
    end function Vint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! Setup Calculations for this module !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
    subroutine MakeK3(ky)
        real(dp), intent(in) :: ky(:)
        integer              :: N, N0, k1, k2, k3i, k4

        N   = size(ky)
        allocate(k3(N,N,N))
        
        do k4=1, N
            do k2=1, N
                do k1=1, N
                    k3i = k1+k2-k4
                    
                    if(k3i < 1 .or. N < k3i) k3i = 0

                    k3(k1, k2, k4) = k3i
                end do
            end do
        end do      
    end subroutine MakeK3
    

    subroutine MakeQs(ky, ae, ah)
        real(dp), intent(in) :: ky(:), ae, ah
        integer              :: k1, k2, Nk

        Nk = size(ky)

        allocate(qe(Nk,Nk))
        allocate(qh(Nk,Nk))

        qe = 0d0
        qh = 0d0

        do k2=1, Nk
            do k1=1, Nk
                qe(k1, k2) = max(abs( ky(k2) - ky(k1) ), ae/2d0)
            end do
        end do

        do k2=1, Nk
            do k1=1, Nk
                qh(k1, k2) = max(abs( ky(k2) - ky(k1) ), ah/2d0)
            end do
        end do
    end subroutine MakeQs
    

        
    
    subroutine MakeUnDel(ky)
        real(dp), intent(in) :: ky(:)
        integer              :: N, i, j

        N = size(ky)

        allocate(unDel(0:N,0:N))

        UnDel = 1d0

        Undel(0,:) = 0d0
        Undel(:,0) = 0d0

        forall(i=1:N) unDel(i,i) = 0d0
        
    end subroutine MakeUnDel


    subroutine CalcMBArrays(ky, Ee, Eh, ge, gh)
        real(dp), intent(in) :: ky(:), Ee(:), Eh(:), ge, gh
        real(dp)             :: a, dky, const
        real(dp)             :: geh, hge2, hgh2, hgeh2
        integer              :: N, k1, k2, k4, k30, k40

        N   = size(ky)
        dky = ky(2) - ky(1)

        geh   = (ge+gh) / 2d0
        hge2  = (hbar * ge)**2
        hgh2  = (hbar * gh)**2
        hgeh2 = (hbar * geh)**2

        allocate(Ceh(0:N,0:N,0:N))
        allocate(Cee(0:N,0:N,0:N))
        allocate(Chh(0:N,0:N,0:N))
        !allocate(Deh(N,N,N))

        Ceh = 0d0
        Cee = 0d0
        Chh = 0d0
        !Deh = 0d0
        
        if(LorentzDelta) then
        do k1=1, N
            do k2=1, N
                do k4=1, N
                    ! Veh(k1,k2,k30,k4)
                    k30 = k3(k4,k2,k1)
                    Ceh(k1,k2,k4) = 2d0 * geh * UnDel(k1,k4)  * Undel(k2,k30) / &
                                    ( (Ee(k1) + Eh(k2) - Eh(k30) - Ee(k4))**2 + hgeh2)

                    !! Veh(k1,k2,k4,k40)
                    !k40 = k3(k1,k4,k2)
                    !Deh(k1,k2,k4) = 2d0 * geh * kdel(k40) * UnDel(k1,k40) * UnDel(k2,k4)  / &
                    !                ( (Ee(k1) + Eh(k2) - Eh(k4) - Ee(k40))**2 + hgeh2)

                    ! Vee(k1,k2,k30,k4)
                    k30 = k3(k1,k2,k4)
                    Cee(k1,k2,k4) = 2d0 * ge  * UnDel(k1,k4)  * Undel(k2,k30) / &
                                    ( (Ee(k1) + Ee(k2) - Ee(k30) - Ee(k4))**2 + hge2)

                    ! Vhh(k1,k2,k30,k4)
                    k30 = k3(k1,k2,k4)
                    Chh(k1,k2,k4) = 2d0 * gh  * UnDel(k1,k4)  * Undel(k2,k30) / &
                                    ( (Eh(k1) + Eh(k2) - Eh(k30) - Eh(k4))**2 + hgh2)
                end do
            end do
        end do   
        else
        do k1=1, N
            do k2=1, N
                do k4=1, N
                    ! Veh(k1,k2,k30,k4)
                    k30 = k3(k4,k2,k1)
                    Ceh(k1,k2,k4) = twopi/hbar * UnDel(k1,k4)  * Undel(k2,k30) * &
                                    GaussDelta(Ee(k1) + Eh(k2) - Eh(k30) - Ee(k4), hbar*geh)

                    !! Veh(k1,k2,k4,k40)
                    !k40 = k3(k1,k4,k2)
                    !Deh(k1,k2,k4) = 2d0 * geh * kdel(k40) * UnDel(k1,k40) * UnDel(k2,k4)  / &
                    !                ( (Ee(k1) + Eh(k2) - Eh(k4) - Ee(k40))**2 + hgeh2)

                    ! Vee(k1,k2,k30,k4)
                    k30 = k3(k1,k2,k4)
                    Cee(k1,k2,k4) = twopi/hbar * UnDel(k1,k4)  * Undel(k2,k30) * &
                                    GaussDelta(Ee(k1) + Ee(k2) - Ee(k30) - Ee(k4), hbar*ge )

                    ! Vhh(k1,k2,k30,k4)
                    k30 = k3(k1,k2,k4)
                    Chh(k1,k2,k4) = twopi/hbar * UnDel(k1,k4)  * Undel(k2,k30) * &
                                    GaussDelta(Eh(k1) + Eh(k2) - Eh(k30) - Eh(k4), hbar*gh )
                end do
            end do
        end do 
        end if
        
    end subroutine CalcMBArrays


    !elemental Function GaussDelta(a,b)
    !    real(dp), intent(in) :: a, b
    !    real(dp)             :: GaussDelta
    !    GaussDelta = 0d0
    !    GaussDelta = 1d0 / sqrt(pi) / b * Exp(-(a/b)**2)
    !end Function
    
    subroutine SetLorentzDelta(boolean)
        logical, intent(in) :: boolean
        LorentzDelta = boolean
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! Coulomb Screening Calculations !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine GetChi1Dqw(alphae, alphah, Delta0, L, epsr, game, gamh, ky, Ee, Eh, ne, nh, qq, w, chir, chii)
        real(dp), intent(in)    :: alphae, alphah, Delta0, L, epsr, gamE(:), gamH(:)
        real(dp), intent(in)    :: ky(:), Ee(:), Eh(:), ne(:), nh(:), qq, w
        real(dp), intent(inout) :: chir, chii

        real(dp)             :: qmine, qminh
        real(dp)             :: Re, Rh, hge(size(ky)), hgh(size(ky))
        complex(dp)          :: chi
        real(dp)             :: ql, Ce, Ch, Ke, Kh, beta, dk
        real(dp)             :: OhmEp, OhmEm, OhmHp, OhmHm
        integer              :: k, q, Nk
        !real(dp)             :: hge, hgh

        Nk = size(ky)
        dk = ky(2) - ky(1)


        hge = hbar*game
        hgh = hbar*gamh

        forall(k=1:Nk) hge(k) = max(hge(k),1d-4*eV)
        forall(k=1:Nk) hgh(k) = max(hgh(k),1d-4*eV)

        qmine = alphae / 2d0
        qminh = alphah / 2d0
        
        Re    = sqrt((2d0/alphae)**2 + Delta0**2)
        Rh    = sqrt((2d0/alphah)**2 + Delta0**2)

        beta  = e0**2 / (4*pi*eps0*epsr)

        Ke = K03(max(ql,qminh) * Re) * beta*2
        Kh = K03(max(ql,qmine) * Rh) * beta*2

        chii = 0d0
        chir = 0d0
        chi  = 0d0
        q    = nint(qq / dk)
        
        do k=max(1,1-q), min(Nk,Nk-q)
            chi = chi - Ke/pi*(ne(k) - ne(k+q)) / (hbar*w - Ee(k+q) + Ee(k) + ii*(hge(k+q)+hge(k))) * dk &
                      - Kh/pi*(nh(k) - nh(k+q)) / (hbar*w - Eh(k+q) + Eh(k) + ii*(hgh(k+q)+hgh(k))) * dk
        end do

        chir = real(chi)
        chii = aimag(chi)
    end subroutine



    subroutine GetEps1Dqw(alphae, alphah, Delta0, L, epsr, me, mh, n1D, q, w, epr, epi)
        real(dp), intent(in) :: alphae, alphah, Delta0, L, epsr, me, mh, n1D
        real(dp)             :: qmine, qminh
        real(dp)             :: Re, Rh
        real(dp)             :: q, w, epr, epi
        real(dp)             :: ql, Ce, Ch, Ke, Kh, beta
        real(dp)             :: OhmEp, OhmEm, OhmHp, OhmHm
        integer              :: k1,k2

        qmine = alphae / 2d0
        qminh = alphah / 2d0
        
        Re    = sqrt((2d0/alphae)**2 + Delta0**2)
        Rh    = sqrt((2d0/alphah)**2 + Delta0**2)

        beta  = e0**2 / (4*pi*eps0*epsr)

        if(abs(q)<1d0) q = 1d0
        ql = abs(q)
        
        Ce = 2*beta*me/pi/hbar**2/ql
        Ch = 2*beta*mh/pi/hbar**2/ql

        Ke = K03(max(ql,qminh) * Re)
        Kh = K03(max(ql,qmine) * Rh)

        OhmEp = hbar*ql/2d0/me * abs(ql + pi*n1D)
        OhmEm = hbar*ql/2d0/me * abs(ql - pi*n1D)
        OhmHp = hbar*ql/2d0/mh * abs(ql + pi*n1D)
        OhmHm = hbar*ql/2d0/mh * abs(ql - pi*n1D)

        epi = 0d0
        epr = 1d0
        if(abs(w) < OhmEm) epr = epr - Ce * Ke * log((w**2 - OhmEm**2)/(w**2 - OhmEp**2))
        if(abs(w) < OhmHm) epr = epr - Ch * Kh * log((w**2 - OhmHm**2)/(w**2 - OhmHp**2))
        if(abs(w) > OhmEp) epr = epr - Ce * Ke * log((w**2 - OhmEm**2)/(w**2 - OhmEp**2))
        if(abs(w) > OhmHp) epr = epr - Ch * Kh * log((w**2 - OhmHm**2)/(w**2 - OhmHp**2))

        epr = 1d0 - Ce * Ke * log(abs((w**2 - OhmEm**2)/(w**2 - OhmEp**2))) &
                  - Ch * Kh * log(abs((w**2 - OhmHm**2)/(w**2 - OhmHp**2)))

!print*, "OEp, OEm	", OhmEp, OhmEm
!print*, "OHp, OHm	", OhmHp, OhmHm
!print*, "Ke & Kh	", Ke, Kh
!print*, "Ce & Ch	", Ce, Ch
!print*, "ql & w		", ql, w
!print*, "Em, Ep		", (w**2 - OhmEm**2),(w**2 - OhmEp**2)
!print*, "Hm, Hp		", (w**2 - OhmHm**2),(w**2 - OhmHp**2)
!print*, "logE, logH	", log((w**2 - OhmEm**2)/(w**2 - OhmEp**2)), log((w**2 - OhmHm**2)/(w**2 - OhmHp**2))
!print*, "epr & epi	", epr, epi
!stop
        
        epi = 0d0
        if(min(OhmEm,OhmEp) < w .and. w < max(OhmEm,OhmEp)) then
            if(min(OhmHm,OhmHp) < w .and. w < max(OhmHm,OhmHp)) then
                epi = - Ce*Ke - Ch*Kh
            endif
        endif

    if (isnan(epr) .or. isnan(epi)) then
        print*, "NAN in EpsL(q,w) at (q,w) =	", q, w
        print*, "Re[EpsL(q,w)] =		", epr
        print*, "Im[EpsL(q,w)] =		", epi
        ! stop
    end if
    end subroutine



    subroutine CalcChi1D(ky, alphae, alphah, Delta0, L, epsr, me, mh)
        real(dp), intent(in) :: ky(:)
        real(dp), intent(in) :: alphae, alphah, Delta0, L, epsr, me, mh
        real(dp)             :: qmine, qminh
        real(dp)             :: Re, Rh
        integer              :: k1,k2

        allocate(Chi1De(size(ky),size(ky)))
        allocate(Chi1Dh(size(ky),size(ky)))

        Chi1De = 0d0
        Chi1Dh = 0d0

        
        qmine = alphae / 2d0
        qminh = alphah / 2d0
        
        Re    = sqrt((2d0/alphae)**2 + Delta0**2)
        Rh    = sqrt((2d0/alphah)**2 + Delta0**2)
        
        do k2=1, size(ky)
            do k1=1, size(ky)                
                Chi1De(k1, k2) = me * K03(qe(k1,k2) * Re) / qe(k1,k2)
            end do
        end do
        
        do k2=1, size(ky)
            do k1=1, size(ky)                
                Chi1Dh(k1, k2) = mh * K03(qh(k1,k2) * Rh) / qh(k1,k2)
            end do
        end do
        
        Chi1De = Chi1De * e0**2 / (twopi*eps0*epsr*hbar**2)
        Chi1Dh = Chi1Dh * e0**2 / (twopi*eps0*epsr*hbar**2)
    end subroutine


    function Eps1D(n1D, Nk)
        real(dp), intent(in   ) :: n1D
        integer,  intent(in   ) :: Nk
        real(dp)                :: Eps1D(Nk,Nk)
                
        Eps1D(:,:) = 1d0 - Chi1De(:,:) * 2*log( abs(qe(:,:) - pi*n1D) / (qe(:,:) + n1D) ) &
                         - Chi1Dh(:,:) * 2*log( abs(qh(:,:) - pi*n1D) / (qh(:,:) + n1D) )
    end function

    
    subroutine CalcScreenedArrays(screened, L, ne, nh,  VC, E1D)
        logical,     intent(in   ) :: screened
        real(dp),    intent(in   ) :: L
        complex(dp), intent(in   ) :: ne(:), nh(:)
        real(dp),    intent(inout) :: VC(:,:,:)
        real(dp)                   :: E1D(:,:)
        real(dp)                   :: density_1D
        real(dp)                   :: density_max

        VC(:,:,1) = Veh0
        VC(:,:,2) = Vee0
        VC(:,:,3) = Vhh0
        E1D(:,: ) = 1d0 

        density_max = min(qe(2,2), qh(2,2)) / pi * 0.99
        density_1D  = sum(real(ne) + real(nh)) / 2d0 / L
        density_1D  = min(density_1D, density_max)
        
        if(screened) then
            E1D(:,:)  = Eps1D(density_1D, size(ne))
            
            VC(:,:,1) = VC(:,:,1) / E1D(:,:)
            VC(:,:,2) = VC(:,:,2) / E1D(:,:)
            VC(:,:,3) = VC(:,:,3) / E1D(:,:)
        end if            
    end subroutine CalcScreenedArrays
    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! Coulomb Calculations for the SBEs !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
    subroutine CalcMVeh(p, VC, MVeh)
        complex(dp), intent(in   ) :: p(:,:,:)
        real(dp),    intent(in   ) :: VC(:,:,:)
        complex(dp), intent(inout) :: MVeh(:, :, :)
        real(dp)                   :: Veh(size(p,1), size(p,2))
        integer                    :: k, kp, N, f
        integer                    :: q, qp, qmin, qmax
        integer                    :: kkpq0, kkpq

        MVeh = 0d0
        N    = size(p,1)

        Veh = VC(:,:,1)
        
        !$omp parallel do private(f, kp, k, q, qp)
        do f=1, size(p,3)
            do kp=1, N
                do k=1, N
                    do q=1, N
                        qp = k3(kp,q,k)
                        MVeh(k,kp,f) = MVeh(k,kp,f) + p(q, qp,f) * Veh(k,q) * UnDel(k,q) * UnDel(kp,qp)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine


    
    function undell(k,q)
        integer :: k,q
        double precision :: undell
        undell = undel(k,q)
        !undell = 1d0
        !undell = 1d0*(k-q)/(1d0*k-1d0*q + small)
    end function undell


    
    subroutine BGRenorm(C, D, VC, BGR)
        complex(dp), intent(in   ) :: C(:,:), D(:,:)
        real(dp),    intent(in   ) :: VC(:,:,:)
        complex(dp), intent(inout) :: BGR(:,:)
        integer                    :: k, kp
        real(dp)                   :: Veh(size(C,1), size(D,2))
        real(dp)                   :: Vee(size(C,1), size(C,2))
        real(dp)                   :: Vhh(size(D,1), size(D,2))
        complex(dp)                :: ne(size(C,1)), nh(size(D,1))

        forall(k=1:size(ne)) ne(k) = C(k,k)
        forall(k=1:size(nh)) nh(k) = D(k,k)

        Veh = VC(:,:,1)
        Vee = VC(:,:,2)
        Vhh = VC(:,:,3)
        
        BGR = 0d0
        
        !$omp parallel do private(kp, k)
        do kp=1, size(nh)
            do k=1, size(ne)
               !BGR(k,kp) = + 2d0 * sum(nh(:) * Vhh(kp,kp)) - sum(nh(:) * Vhh(:,kp) * UnDel(:,kp)) &
               !            + 2d0 * sum(ne(:) * Vee(k ,k )) - sum(ne(:) * Vee(:,k ) * UnDel(:,k )) &
               !            - 2d0 * sum(ne(:) * Veh(kp,kp)) - sum(nh(:) * Veh(k,k))
                BGR(k,kp) = - sum(nh(:) * Vhh(:,kp) * UnDel(:,kp)) &
                            - sum(ne(:) * Vee(:,k ) * UnDel(:,k )) 
            end do
        end do
        !$omp end parallel do

    end subroutine BGRenorm
    

  
    subroutine EeRenorm(ne, VC, BGR)
        complex(dp), intent(in   ) :: ne(:)
        real(dp),    intent(in   ) :: VC(:,:,:)
        complex(dp), intent(inout) :: BGR(:,:)
        integer                    :: k, kp
        real(dp)                   :: Vee(size(ne), size(ne))

        Vee = VC(:,:,2)
        BGR = 0d0
        
        !$omp parallel do private(kp, k)
        do kp=1, size(ne)
            do k=1, size(ne)
                BGR(k,kp) = + 2d0 * sum(ne(:) * Vee(kp,kp)) - sum(ne(:) * Vee(:,kp) * UnDel(:,kp)) &
                            + 2d0 * sum(ne(:) * Vee(k ,k )) - sum(ne(:) * Vee(:,k ) * UnDel(:,k )) &
                            - 2d0 * sum(ne(:) * Vee(kp,kp)) - sum(ne(:) * Vee(k,k))
            end do
        end do
        !$omp end parallel do
    end subroutine EeRenorm


    subroutine EhRenorm(nh, VC, BGR)
        complex(dp), intent(in   ) :: nh(:)
        real(dp),    intent(in   ) :: VC(:,:,:)
        complex(dp), intent(inout) :: BGR(:,:)
        integer                    :: k, kp
        real(dp)                   :: Vhh(size(nh), size(nh))

        Vhh = VC(:,:,3)
        BGR = 0d0
        
        !$omp parallel do private(kp, k)
        do kp=1, size(nh)
            do k=1, size(nh)
                BGR(k,kp) = + 2d0 * sum(nh(:) * Vhh(kp,kp)) - sum(nh(:) * Vhh(:,kp) * UnDel(:,kp)) &
                            + 2d0 * sum(nh(:) * Vhh(k ,k )) - sum(nh(:) * Vhh(:,k ) * UnDel(:,k )) &
                            - 2d0 * sum(nh(:) * Vhh(kp,kp)) - sum(nh(:) * Vhh(k,k))
            end do
        end do
        !$omp end parallel do
    end subroutine EhRenorm
 

  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! MANY BODY RELAXATION EFFECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! i.e., Non-Hartree-Foch Terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Calculate the Many-body Coulomb In/Out rates for electrons
    subroutine MBCE2(ne0, nh0, ky, Ee, Eh, VC, geh, ge, Win, Wout)
        real(dp), intent(in   ) :: ne0(:), nh0(:)         ! Carrier populations
        real(dp), intent(in   ) :: ky(:), Ee(:), Eh(:)
        real(dp), intent(in   ) :: VC(:,:,:)
        real(dp), intent(in   ) :: geh, ge
        real(dp), intent(inout) :: Win(:), Wout(:)
        real(dp)                :: ne(0:size(ne0))
        real(dp)                :: nh(0:size(nh0))
        integer                 :: k1, k2, k30, k4
        integer                 :: kp, k1p
        integer                 :: Nk, k, q1, q2
        real(dp)                :: u(size(ne0))
        real(dp)                :: Veh2(size(ne), size(nh))
        real(dp)                :: Vee2(size(ne), size(nh))

        Nk   = size(ne0)
        
        Veh2 = VC(:,:,1)**2
        Vee2 = VC(:,:,2)**2

        ne = 0d0
        nh = 0d0
        ne(1:Nk) = abs(ne0)
        nh(1:Nk) = abs(nh0)

        !$omp parallel do private(k, q1, q2, kp, k1, k2, k30, k4, k1p)        
        do k=1, Nk
            do q1=1, Nk
                do q2=1, Nk

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Electron-Hole In/Out rates, Veh(k,kp,k1p,k1), Veh(k1,kp,kp1,k) !!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! k + k1p = kp + k1, k1 + k1p = kp + k !!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    kp      = q1
                    k1      = q2
                    
                    k1p     = k3(kp,k1,k)
                    Win(k)  = Win(k)  + Veh2(k,k1) * (1d0 - nh(kp)) * nh(k1p)        * ne(k1)  * Ceh(k,kp,k1)

                    k1p     = k3(kp,k,k1)
                    Wout(k) = Wout(k) + Veh2(k1,k) * (1d0 - ne(k1)) * (1d0 - nh(kp)) * nh(k1p) * Ceh(k1,kp,k)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Elec-Elec In/Out rates, Vee(k,k2,k30,k4), Vee(k4,k2,k30,k) !!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!! k + k2 = k30 + k4, k4 + k2 = k30 + k !!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
                    k2      = q1
                    k4      = q2
                    
                    k30     = k3(k,k2,k4)
                    Win(k)  = Win(k)  + Vee2(k,k4) * (1d0 - ne(k2)) * ne(k30)        * ne(k4)  * Cee(k,k2,k4)

                    k30     = k3(k4,k2,k)
                    Wout(k) = Wout(k) + Vee2(k4,k) * (1d0 - ne(k4)) * (1d0 - ne(k2)) * ne(k30) * Cee(k4,k2,k)
                end do
            end do
        end do
        !$omp end parallel do
        

    end subroutine MBCE2



    ! Calculate the Many-body Coulomb In/Out rates for electrons
    subroutine MBCE(ne0, nh0, ky, Ee, Eh, VC, geh, ge, Win, Wout)
        real(dp), intent(in   ) :: ne0(:), nh0(:)         ! Carrier populations
        real(dp), intent(in   ) :: ky(:), Ee(:), Eh(:)
        real(dp), intent(in   ) :: VC(:,:,:)
        real(dp), intent(in   ) :: geh, ge
        real(dp), intent(inout) :: Win(:), Wout(:)
        real(dp)                :: ne(0:size(ne0))
        real(dp)                :: nh(0:size(nh0))
        integer                 :: k1, k2, k30, k4
        integer                 :: kp, k1p
        integer                 :: Nk, k, q1, q2
        real(dp)                :: u(size(ne0))
        real(dp)                :: Veh2(size(ne), size(nh))
        real(dp)                :: Vee2(size(ne), size(nh))

        Nk   = size(ne0)
        
        Veh2 = VC(:,:,1)**2
        Vee2 = VC(:,:,2)**2

        ne = 0d0
        nh = 0d0
        ne(1:Nk) = abs(ne0)
        nh(1:Nk) = abs(nh0)

        !$omp parallel do private(k, q1, q2, kp, k1, k2, k30, k4, k1p)        
        do k=1, Nk
            do q1=1, Nk
                do q2=1, Nk

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Electron-Hole In/Out rates, Veh(k,kp,k1p,k1), Veh(k1,kp,kp1,k) !!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! k + k1p = kp + k1, k1 + k1p = kp + k !!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    kp      = q1
                    k1      = q2
                    
                    k1p     = k3(kp,k1,k)
                    Win(k)  = Win(k)  + Veh2(k,k1) * (1d0 - nh(kp)) * nh(k1p)        * ne(k1)  * Ceh(k,kp,k1)

                    k1p     = k3(kp,k,k1)
                    Wout(k) = Wout(k) + Veh2(k1,k) * (1d0 - ne(k1)) * (1d0 - nh(kp)) * nh(k1p) * Ceh(k1,kp,k)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Elec-Elec In/Out rates, Vee(k,k2,k30,k4), Vee(k4,k2,k30,k) !!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!! k + k2 = k30 + k4, k4 + k2 = k30 + k !!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
                    k2      = q1
                    k4      = q2
                    
                    k30     = k3(k,k2,k4)
                    Win(k)  = Win(k)  + Vee2(k,k4) * (1d0 - ne(k2)) * ne(k30)        * ne(k4)  * Cee(k,k2,k4)

                    k30     = k3(k4,k2,k)
                    Wout(k) = Wout(k) + Vee2(k4,k) * (1d0 - ne(k4)) * (1d0 - ne(k2)) * ne(k30) * Cee(k4,k2,k)
                end do
            end do
        end do
        !$omp end parallel do
        

    end subroutine MBCE


    ! Calculate the Many-body Coulomb In/Out rates for electrons
    subroutine MBCH(ne0, nh0, ky, Ee, Eh, VC, geh, gh, Win, Wout)
        real(dp), intent(in   ) :: ne0(:), nh0(:)         ! Carrier populations
        real(dp), intent(in   ) :: ky(:), Ee(:), Eh(:)
        real(dp), intent(in   ) :: VC(:,:,:)
        real(dp), intent(in   ) :: geh, gh
        real(dp), intent(inout) :: Win(:), Wout(:)
        real(dp)                :: ne(0:size(ne0))
        real(dp)                :: nh(0:size(nh0))
        integer                 :: k1p, k2p, k3p, k4p
        integer                 :: k, k1
        integer                 :: Nk, kp, q1, q2
        real(dp)                :: u
        real(dp)                :: Veh2(size(ne), size(nh))
        real(dp)                :: Vhh2(size(ne), size(nh))

        Nk   = size(ne0)
        
        Veh2 = VC(:,:,1)**2
        Vhh2 = VC(:,:,2)**2

        ne = 0d0
        nh = 0d0
        ne(1:Nk) = abs(ne0)
        nh(1:Nk) = abs(nh0)

        !$omp parallel do private(k, q1, q2, kp, k1p, k2p, k3p, k4p, k1)        
        do kp=1, Nk
            do q1=1, Nk
                do q2=1, Nk

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Electron-Hole In/Out rates, Veh(k,kp,k1p,k1), Veh(k,k1p,kp,k1) !!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! k + k1p = kp + k1, k + kp = k1p + k1 !!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    k       = q1
                    k1      = q2
                    
                    k1p     = k3(kp,k1,k)
                    Win(kp)  = Win(kp)  + Veh2(k,k1) * (1d0 - ne(k)) * nh(k1p)         * ne(k1)  * Ceh(k,kp,k1)

                    k1p     = k3(kp,k,k1)
                    Wout(kp) = Wout(kp) + Veh2(k1,k) * (1d0 - ne(k)) * (1d0 - nh(k1p)) * ne(k1) * Ceh(k,k1p,k1)


                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!! Hole-Hole In/Out rates, Vhh(kp,k2p,k3p,k4p), Vhh(k4p,k2p,k3p,kp) !!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!! kp + k2p = k3p + k4p, k4p + k2p = k3p + kp !!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
                    k2p     = q1
                    k4p     = q2
                    
                    k3p      = k3(kp,k2p,k4p)
                    Win(kp)  = Win(kp)  + Vhh2(kp,k4p) * (1d0 - nh(k2p)) * nh(k3p)         * nh(k4p) * Chh(kp,k2p,k4p)

                    k3p      = k3(k4p,k2p,kp)
                    Wout(kp) = Wout(kp) + Vhh2(k4p,kp) * (1d0 - nh(k4p)) * (1d0 - nh(k2p)) * nh(k3p) * Chh(k4p,k2p,kp)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine MBCH    
               
    
end module coulomb
