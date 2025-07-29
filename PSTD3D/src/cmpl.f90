!===========================================================
! Module: cpml
! Implements 3D Convolutional PML (CPML) for a pseudo-spectral
! time-domain (PSTD) 
!
! This module:
!   - Allocates and initializes CPML coefficients
!   - Updates auxiliary psi fields for absorbing boundary layers
!   - Applies CPML corrections in E- and H-field updates
!
! Variables:
!   npml_x/y/z : PML thickness (# grid points)
!   sigma_max  : maximum conductivity at PML outer edge
!   kappa_max  : scaling for stretched coordinates
!   alpha_max  : low-frequency term (avoids DC divergence)
!   be,ce,bh,ch: CPML coefficient arrays
!   psi_***    : Auxiliary fields for each field component
!
! Author: Emily Hatten
!===========================================================
module cpml
    use constants
    use fftw
    implicit none

#include "config.f"

    ! === CPML Parameters === !
    integer                     :: NPML_x, NPML_y, NPML_z
    real(dp)                    :: sigma_max, alpha_max, kappa_max
    real(dp)                    :: dx, dy, dz, dt

    ! === Arrays for CPML Coefficents === !
    real(dp), allocatable       :: bEx(:), cEx(:), bHx(:), cHx(:)
    real(dp), allocatable       :: bEy(:), cEy(:), bHy(:), cHy(:)
    real(dp), allocatable       :: bEz(:), cEz(:), bHz(:), cHz(:)

    ! === Arrays for Auxiliary CPML  === !
    real(dp), allocatable       :: psi_Exy(:, :, :), psi_Exz(:, :, :)
    real(dp), allocatable       :: psi_Eyx(:, :, :), psi_Eyz(:, :, :)
    real(dp), allocatable       :: psi_Ezx(:, :, :), psi_Ezy(:, :, :)

    real(dp), allocatable       :: psi_Hxy(:, :, :), psi_Hxz(:, :, :)
    real(dp), allocatable       :: psi_Hyx(:, :, :), psi_Hyz(:, :, :)
    real(dp), allocatable       :: psi_Hzx(:, :, :), psi_Hzy(:, :, :)
   



contains

! ==================================================
!   subroutine: InitCPML
!   Inialaizes CPML Arrays and Coefficents
!   Input:
!       Nx, Ny, Nz      :: Grid Dimensions
!       Dx, Dy, Dz, Dt  :: Spatial/Temporal Steps
!       NPML            :: PML Thickness
! =================================================
    subroutine InitCPML(Nx, Ny, Nz, NPML, Dx, Dy, Dz, Dt)
        integer, intent(in)     :: Nx, Ny, Nz, NPML
        real(dp), intent(in)    :: Dx, Dy, Dz, Dt
        real(dp)                :: sigma_max, kappa_max, alpha_max


        
        

        ! CPML Params (Recommend from Literature)
        sigma_max = 
        kappa_max =
        alpha_max = 

        ! Allocate Coeffient Arrays
        allocate(bEx(npml), bHx(npml), cEx(npml), cHx(npml))
        allocate(bEy(npml), bHy(npml), cEy(npml), cHy(npml))
        allocate(bEz(npml), bHz(npml), cEz(npml), cHz(npml))

        ! Compute CPML Coefficents
        call Calc_CoefficentsCPML()
        call Calc_CoefficentsCPML()
        call Calc_CoefficentsCPML()
        call Calc_CoefficentsCPML()

        ! Allocate Auxiliary Fields
        allocate(psi_Exy(Nx, Ny, Nz), psi_Exz(Nx, Ny, Nz))
        allocate(psi_Eyx(Nx, Ny, Nz), psi_Eyz(Nx, Ny, Nz))
        allocate(psi_Ezx(Nx, Ny, Nz), psi_Ezy(Nx, Ny, Nz))

        allocate(psi_Hxy(Nx, Ny, Nz), psi_Hxz(Nx, Ny, Nz))
        allocate(psi_Hyx(Nx, Ny, Nz), psi_Hyz(Nx, Ny, Nz))
        allocate(psi_Hzx(Nx, Ny, Nz), psi_Hzy(Nx, Ny, Nz))

        psi_Exy = 0d0
        psi_Exz = 0d0
        psi_Eyx = 0d0
        psi_Eyz = 0d0
        psi_Ezx = 0d0
        psi_Ezy = 0d0
        psi_Hxy = 0d0
        psi_Hxz = 0d0
        psi_Hyx = 0d0
        psi_Hyz = 0d0
        psi_Hzx = 0d0
        psi_Hzy = 0d0

    end subroutine

! ==================================================
!   subroutine: UpdateCPML_E
!   
!   Input:
!       Nx, Ny, Nz      :: Grid Dimensions
!       Dx, Dy, Dz, Dt  :: Spatial/Temporal Steps
!       NPML            :: PML Thickness
! =================================================
    subroutine UpdateCPML_E()

        do k = 1, Nz -1
            do j = 1, Ny -1
                do i = 1, Nx -1
                    Ex(i,j,k) = 
                end do
            end do
        end do

    end subroutine


! ==================================================
!   subroutine: UpdateCPML_H
!  
!   Input:
!       Nx, Ny, Nz      :: Grid Dimensions
!       Dx, Dy, Dz, Dt  :: Spatial/Temporal Steps
!       NPML            :: PML Thickness
! =================================================
    subroutine UpdateCPML_H()



        ! Hx Component 
        do k = 1, Nz -1
            do j = 1, Ny -1
                do i = 1, Nx -1
                    psi_Hxy(i, j, k) = bHy(j) * psi_Hxy(i, j, k) + (cHy(j) * ((Ez(i, j+1, k) - Ez(i, j, k)) / dy))
                    psi_Hxz(i, j, k) = bHz(k) * psi_Hyz(i, j, k) + (cHz(k) * ((Ey(i, j, k+1) - Ey(i, j, k)) / dz))

                    Hx(i,j,k) = DaX(i, j, k) * Hx(i, j, k) - DbX(i, j, k) * (((Ez(i,j,k) - Ez(i,j-1,k)) / dy  ) - ((Ey(i,j,k)-Ey(i,j,k-1)) / dz) * (dt / epsilon0 )) - psi_Hxy(i, j, k) + psi_Hxz(i, j, k)
                end do
            end do
        end do

        ! Hy Component 
        do k = 1, Nz -1
            do j = 1, Ny -1
                do i = 1, Nx -1
                    psi_Hyx(i, j, k) = bHx(i) * psi_Hyx(i, j, k) + cHx(i) * (Ez(i+1, j, k) - Ez(i, j, k)) / dx
                    psi_Hyz(i, j, k) = bHz(k) * psi_Hyz(i, j, k) + cHz(k) * (Ex(i, j, k+1) - Ex(i, j, k)) / dz
                    Hy(i,j,k) = Hy(i, j, k) -              - psi_Hyx(i, j, k) + psi_Hyz(i, j, k)
                end do
            end do
        end do


        ! Hz Component 
        do k = 1, Nz -1
            do j = 1, Ny -1
                do i = 1, Nx -1
                    psi_Hzx(i, j, k) = bHx(i) * psi_Hzx(i, j, k) + cHx(i) * (Ex(i+1, j, k) - Ex(i, j, k)) / dy
                    psi_Hzy(i, j, k) = bHy(j) * psi_Hzy(i, j, k) + cHy(j) * (Ey(i, j+1, k) - Ey(i, j, k)) / dx
                    Hz(i,j,k) = Hz(i, j, k) +             - psi_Hzx(i, j, k) + psi_Hzy(i, j, k)
                end do
            end do
        end do

    end subroutine


! ==================================================
!   subroutine: Calc_CoefficentsCPML
!   Computes: 
!       Graded Conductivity (sigma)
!       Stretching Factor   (kappa)
!       Low Frequency Loss  (alpha)
!  
!   Uses to compute:
!       bE/cE: E-Field Auxiliary Terms
!       bH/cH: H-Field Auxiliary Terms 
!               
!   Input:
!       N               :: 
!       NPML            :: PML Thickness
!       d               ::
!       bE, cE          ::
!       bH, cH          ::
! =================================================
    subroutine Calc_CoefficentsCPML(N, NPML, d, bE, cE, bH, cH)

        real(dp), allocatable       :: 
        real(dp), allocatable       :: bE(:), cE(:), bH(:), cH(:)
        real(dp)                    :: epsilon0, c0
        integer                     :: m                                ! 
        integer                     :: d                                ! Spaital Step in 
        integer                     :: dPML                             ! Thickness of PML
        real(dp)                    :: pos                              ! Position in the PML


        m = 3
        sigma_max = real(0.8d0 * (m + 1) * epsilon0 * c0, dp) / 2 * dPML
        kappa_max = 5d0
        alpha_max = 0.1d0

        ! Loop through Grid Points within the PML
        do d = 1, nPML

            sigma = sigma_max * (d / dPML)**m
            kappa = 1 + (kappa_max - 1) * (d/dPML)**m
            alpha = alpha_max * (1 - (d/dPML))

            ! Coefficents for E-Field Updates 

            ! Coefficents for H-Field Updtes 

            ! Mirror at other end
        end do

    end subroutine



end module 