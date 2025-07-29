! This module calculates the many-body electron-hole phonon
! contribution to the Semiconductor Bloch equations in
! support of propagation simulations for a quantum wire.
module emission

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers
    use spliner
    use usefulsubs
    implicit none



    ! Pre-calculated QW arrays required for solving SBEs.
    real(dp), private              :: RScale
    real(dp), private              :: kB    = 1.3806504d-23    ! Boltzmann Constant (J/K)
    real(dp), private              :: Temp  = 77d0            ! Temperature of QW solid (K)
    real(dp), private, allocatable :: HOmega(:)
    real(dp), private, allocatable :: square(:)
    integer,  private, allocatable :: idel(:,:)

    integer, private               :: kkk


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine InitializeEmission(ky, Ee, Eh, dcv, epsr, geh, ehint)
        real(dp), intent(in   ) :: ky(:), Ee(:), Eh(:) ! Momentum, Electron and Hole energies
        real(dp), intent(in   ) :: dcv, epsr             !
        real(dp), intent(in   ) :: geh, ehint            !
        integer                 :: Nk, k, k1

        Nk = size(ky)
        
        allocate(idel(Nk, Nk))
        idel = 0d0
        forall(k=1:Nk) idel(k,k) = 1
        idel = 1 - idel

        RScale = 3d0 * dcv**2 / eps0 / sqrt(epsr) * ehint**2

        call CalcHOmega(kB*Temp, hbar*geh)

        allocate(Square(size(HOmega)))

        square(:) = 3d0 * dcv**2 * ehint / eps0 / sqrt(epsr) * ehint / hbar * &
                    Lrtz(HOmega(:), hbar*geh) * Exp(- HOmega(:) / kB / Temp)

    end subroutine InitializeEmission
    
        

    subroutine SpontEmission(ne, nh, Ee, Eh, gap, geh, VC, Rsp)
        complex(dp), intent(in   ) :: ne(:), nh(:)    ! Electron and Hole occupation numbers
        real(dp), intent(in   )    :: Ee(:), Eh(:)    ! Electron and Hole energies
        real(dp), intent(in   )    :: gap, geh        !
        real(dp), intent(in   )    :: VC(:,:,:)       !
        real(dp), intent(inout)    :: RSP(size(Ee))   !
        real(dp)                   :: Ek(size(ne))    !
        integer                    :: k               !
        
        RSP = 0d0

        Ek = gap + Ee + Eh + Ec(real(ne), real(nh), VC)
      
        do k=1, size(ne)
            RSP(k) = SpontIntegral(Ek(k))
        end do        
    end subroutine SpontEmission

    

    function Ec(ne, nh, VC)
        real(dp), intent(in) :: ne(:), nh(:)
        real(dp), intent(in) :: VC(:,:,:)
        real(dp)             :: Ec(size(ne))
        real(dp)             :: Veh(size(ne), size(ne))
        real(dp)             :: Vee(size(ne), size(ne))
        real(dp)             :: Vhh(size(ne), size(ne))
        integer              :: k

        Ec = 0d0
        
        Veh = VC(:,:,1)
        Vee = VC(:,:,2)
        Vhh = VC(:,:,3)
        
        do k=1, size(ne)
            Ec(k) = sum(ne(:) * (Vee(k,k) - Vee(:,k)) + nh(:) * (Vhh(k,k) - Vhh(:,k))) + &
                    sum(ne(:) *  Veh(k,k) * idel(:,k) - nh(:) *  Veh(k,k) * idel(:,k)) - Veh(k,k)
        end do

    end function Ec




    ! Numerical solution to the integral
    function SpontIntegral(Ek)
        real(dp), intent(in) :: Ek
        real(dp)             :: SpontIntegral     
        real(dp)             :: Integrand(size(HOmega))
        real(dp)             :: dhw

        dhw   = HOmega(2) - HOmega(1)

        Integrand(:) = (HOmega(:) + Ek) * rho0(HOmega(:) + Ek) * square(:)
    
        SpontIntegral = sum(Integrand) * dhw
        
    end function SpontIntegral


    
    ! Photon density of states as a function of hbar*w
    elemental function rho0(hw)
        real(dp), intent(in) :: hw
        real(dp)             :: rho0     

        rho0 = hw**2 / (c0**3 * pi**2 * hbar**3)
        
    end function rho0



    
        
!    ! Analytic solution to the integral
!    ! int_{a}^{3a} \frac{b x^3}{(x-a)^2 + b^2}
!    ! Here a = Ek, b = hbar*geh, and x = hbar*omega
!    function SpontIntegral(a, b0)
!        real(dp), intent(in) :: a, b0
!        real(dp)             :: SpontIntegral     
!        real(dp)             :: b, b2, a2
!
!        
!        b  = b0 + hbar*(1d6)
!        a2 = a**2
!        b2 = b**2
!        
!        SpontIntegral = b * ((3d0*a2 - b2) * log((4*a2 + b2) / b2) - 16d0*a2) / 2d0  + &
!                        a * (a2 - 3d0*b2) * ATAN(2*a/b)
!    end function SpontIntegral



    subroutine CalcHOmega(kBT, hg)
        real(dp)              :: kBT, hg
        real(dp)              :: dhw, hwmax
        integer               :: Nw, n

        hwmax = (kBT + hg) * 4d0
        dhw   = min(kBT, hg) / 20d0
        Nw    = ceiling(hwmax / dhw)

        if(Nw < 10) then
            print*, "Error: temperature is too low in emission.f90"
            stop
        end if
        
        allocate(HOmega(Nw))
        
        forall(n=1:Nw) HOmega(n) = 0d0 + (n-0.5)*dhw
    end subroutine CalcHOmega


    subroutine Calchw(hw, PLS, Estart, Emax)
        real(dp) :: hw(:), PLS(:), Estart, Emax
        real(dp) :: dhw, hw_start
        integer  :: w, Nw

        Nw  = size(hw)
        hw  = 0d0
        PLS = 0d0

        dhw = (Emax-Estart) / (1.0*Nw)
        
        do w=1, Nw
            hw(w) = Estart + (w-1) * dhw
        end do
    end subroutine Calchw
    

    subroutine PLSpectrum(ne, nh, Ee, Eh, gap, geh, VC, hw, t, PLS)
        complex(dp), intent(in   ) :: ne(:), nh(:)  ! Electron and Hole occupation numbers
        real(dp), intent(in   )    :: Ee(:), Eh(:)         ! Electron and Hole energies
        real(dp), intent(in   )    :: gap, geh             !
        real(dp), intent(in   )    :: VC(:,:,:)            !
        real(dp), intent(in   )    :: hw(:), t             !
        real(dp), intent(inout)    :: PLS(:) !
        real(dp)                   :: Ek(size(ne))         !
        real(dp), allocatable      :: E(:), nenh(:)        !
        real(dp), allocatable      :: ky(:), qy(:)         !
        real(dp)                   :: tempavg, Te, Th
        integer                    :: w, k, Nk             !
        integer                    :: X = 100
        
        Nk = size(ne)

        Ek = gap + Ee + Eh + Ec(abs(ne), abs(nh), VC)

	Te = Temperature(ne, Ee)
	Th = Temperature(nh, Eh)
        
        tempavg = (Te + Th)/2d0

        allocate(E(X*Nk), nenh(X*Nk), qy(X*Nk), ky(Nk))

        forall(k=1:Nk  ) ky(k) = (k-1)*1d0
        forall(k=1:Nk*X) qy(k) = (k-1)*1d0 / (X*1d0)

        do k=1, X*Nk
            E(k)     = LinearInterp_dp( Ek(:), ky(:), qy(k))
            nenh(k)  = abs(LinearInterp_dpc(ne*nh, ky(:), qy(k)))
        end do
        
        
        !$omp parallel do private(w)        
        do w=1, size(hw)
            PLS(w) = PLS(w) + Rscale * Sum( hw(w) * rho0(hw(w)) * nenh(:)      * &
                                          Exp(- abs(hw(w) - E(:)) / kB / tempavg ) * &
                                          Lrtz(     hw(w) - E(:), hbar*geh) * &
                                          softtheta(hw(w) - E(X*Nk/2), hbar*geh)  &
                                          )
	end do
        !$omp end parallel do

        PLS = PLS * softtheta(t,geh)
    end subroutine PLSpectrum

    
end module
