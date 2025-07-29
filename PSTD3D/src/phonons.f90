! This module calculates the many-body electron-hole phonon
! contribution to the Semiconductor Bloch equations in
! support of propagation simulations for a quantum wire.
module phonons

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers
    use spliner
    use usefulsubs
    implicit none



    ! Pre-calculated QW arrays required for solving SBEs.
    real(dp), private              :: small = 1d-200           ! Smallest # worthy of consideration
    real(dp), private              :: Temp  = 77d0             ! Temperature of QW solid (K)
    real(dp), private              :: kB    = 1.3806504d-23    ! Boltzmann Constant (J/K)
    real(dp), private              :: epsr0 = 10.0             ! Dielectric constant in host AlAs at w=0
    real(dp), private              :: epsrINF = 8.2            ! Dielectric constant in host AlAs at w=INFINITY
    real(dp), private              :: Vscale                   ! Scaling constant e-e and h-h Coulomb arrays
    real(dp), private              :: NO                       ! Boze function for thermal equilibrium
                                                               ! longitudinal-optical phonons in host
    
    real(dp), private, allocatable :: EP(:,:), EPT(:,:)        !
    real(dp), private, allocatable :: HP(:,:), HPT(:,:)        !
    integer,  private, allocatable :: idel(:,:)


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine InitializePhonons(ky, Ee, Eh, L, epsr, Gph, Oph)
        real(dp), intent(in   ) :: ky(:), Ee(:), Eh(:) ! Momentum, Electron and Hole energies
        real(dp), intent(in   ) :: L, epsr             !
        real(dp), intent(in   ) :: Gph, Oph            !
        integer                 :: Nk, k, k1
        
        NO = 1d0 / (exp(hbar*Oph / kB / Temp) - 1d0)
        Nk = size(ky)

        allocate(idel(Nk, Nk))
        idel = 0d0
        forall(k=1:Nk) idel(k,k) = 1
        idel = 1 - idel
        
        allocate(EP(Nk, Nk), EPT(Nk, Nk))
        allocate(HP(Nk, Nk), HPT(Nk, Nk))

        EP  = 0d0
        HP  = 0d0
        EPT = 0d0
        HPT = 0d0
        
        do k1=1, Nk            
            do k=1, Nk
                EP(k,k1) = (NO    ) / ( (Ee(k) - Ee(k1) - hbar*Oph)**2 + (hbar * Gph)**2) + &
                           (NO+1d0) / ( (Ee(k) - Ee(k1) + hbar*Oph)**2 + (hbar * Gph)**2) 

                HP(k,k1) = (NO    ) / ( (Eh(k) - Eh(k1) - hbar*Oph)**2 + (hbar * Gph)**2) + &
                           (NO+1d0) / ( (Eh(k) - Eh(k1) + hbar*Oph)**2 + (hbar * Gph)**2) 
                
            end do
        end do
        
        ! Multiply by scaling factor of (2*pi/hbar)*(hbar*Gph/pi)
        ! Multiplication by the idel array ensures that k =/= kp
        EP = EP * 2d0 * Gph * idel
        HP = HP * 2d0 * Gph * idel

        EPT = Transpose(EP)
        HPT = Transpose(HP)

        Vscale = hbar * Oph * epsr * (1d0/epsrINF - 1d0/epsr0 )

    end subroutine InitializePhonons
    
        
        

    
    subroutine MBPE(ne, VC, E1D, Win, Wout)
        real(dp), intent(in   ) :: ne(:)                   ! Carrier populations
        real(dp), intent(in   ) :: VC(:,:,:), E1D(:,:)
        real(dp), intent(inout) :: Win(:), Wout(:)
        real(dp)                :: Vep(size(ne), size(ne))
        integer                 :: k, Nk
        
        Nk  = size(ne)
        Vep = Vc(:,:,2) / E1D(:,:) * Vscale

        ! Note that since Vep = Transpose(Vep) I
        ! use Vep(:,k) instead of Vep(k,:) below
        !
        !$omp parallel do private(k)        
        do k=1, Nk
            Win(k)  = Win(k)  + sum(Vep(k,:) * (   ne(:)   ) * EPT(:,k))
            Wout(k) = Wout(k) + sum(Vep(k,:) * (1d0 - ne(:)) *  EP(:,k))
        end do
        !$omp end parallel do
    end subroutine MBPE


    subroutine MBPH(nh, VC, E1D, Win, Wout)
        real(dp), intent(in   ) :: nh(:)                 ! Carrier populations
        real(dp), intent(in   ) :: VC(:,:,:), E1D(:,:)
        real(dp), intent(inout) :: Win(:), Wout(:)
        real(dp)                :: Vhp(size(nh), size(nh))
        integer                 :: kp, Nk

        Nk  = size(nh)    
        Vhp = Vc(:,:,3) / E1D(:,:) * Vscale

        ! Note that since Vep = Transpose(Vep) I
        ! use Vep(:,k) instead of Vhp(k,:) below
        !
        !$omp parallel do private(kp)        
        do kp=1, Nk
            Win(kp)  = Win(kp)  + sum(Vhp(:,kp) * (      nh(:)) * HPT(:,kp))
            Wout(kp) = Wout(kp) + sum(Vhp(:,kp) * (1d0 - nh(:)) *  HP(:,kp))
        end do
        !$omp end parallel do
    end subroutine
    
    
    !function id(k,kp)
    !    integer :: k,kp,id
    !    id = idel(k,kp)
    !end function


    ! Calculates Cq for use in the DC Field module
    function Cq2(q, V, E1D)
        real(dp), intent(in) :: q(:), V(:,:), E1D(:,:)
        real(dp)             :: Cq2(size(q))
        real(dp)             :: dq
        integer              :: i, iq(size(q))

        dq = q(2) - q(1)
        
        iq = nint( abs(q / dq) )
        
        do i=1, size(q)
            Cq2(i) = V(1 + iq(i), 1) / E1D(1 + iq(i), 1) * Vscale
        end do
    end function Cq2
    

    ! Calculate Fermi-Dirac Distribution assuming
    ! host temperature and the Fermi Energy = 0
    elemental function FermiDistr(En)
        real(dp),    intent(in   ) :: En
        complex(dp)                :: FermiDistr

        FermiDistr = 1d0 / (Exp(En / kB / Temp) + 1d0)
    end function FermiDistr


    elemental function BoseDistr(En)
        real(dp), intent(in) :: En
        real(dp)             :: BoseDistr

        BoseDistr =  1d0 / (Exp(En / kB / Temp) - 1d0)
    end function

    
   function N00()
        real(dp) :: N00

        N00 = NO
    end function N00
    

end module phonons
