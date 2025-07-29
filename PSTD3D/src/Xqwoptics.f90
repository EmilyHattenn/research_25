! This module serves the SBEs.f90 module.  It contains functions
! and subroutines to interpolate the macroscopic electric fields from
! the propagation program to the field within the quantum wire (QW).
! After the SBEs have been solved, this module also uses SBE
! solutions to calculate the QW polariztion, charge densities, and
! current densities in the QW k- and y-spaces.  It then interpolates
! these source terms back into the Y-space of the propagation simulation.
! This module also writes the QW data to files as requested by the user.
module qwoptics

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers     ! 
    use spliner     !
    use usefulsubs  ! 
    implicit none

    real(dp),    private              :: small = 1d-100    ! Smallest # worthy of consideration
    real(dp),    private              :: epsr  = 9.1 
    real(dp),    private, allocatable :: QWWindow(:)
    complex(dp), private, allocatable :: Expikr(:,:), Expikrc(:,:)
    complex(dp), private              :: dcv0
    complex(dp), private, allocatable :: Xcv0(:,:), Ycv0(:,:), Zcv0(:,:)
    complex(dp), private, allocatable :: Xvc0(:,:), Yvc0(:,:), Zvc0(:,:)
    logical,     private              :: firsttime = .true.
    real(dp),    private              :: Vol
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Subroutines called from the SBEs.f90 module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
    ! Converts the Maxwell electric fields (Eii[:] and Vrr[:]) in
    ! the propagation space (RR[:]) into the the QW electric 
    ! fields (Ei[:] and Vr[:]) in the QW electric field space (r[:] and Qr[:]).
    ! contributions to be converted.  After interpolating from the YY-
    ! space to the y-space, the FFT from the y- to Qy- space is taken.
    ! The average QW field value is also calculated and returned.
    subroutine Prop2QW(RR, Exx, Eyy, Ezz, Vrr, Edc, R, Ex, Ey, Ez, Vr, t, xxx) 
        real(dp),   intent(in   ) :: RR(:)                 ! Maxwell RR spatial array
        complex*16, intent(in   ) :: Exx(:), Eyy(:)        ! Maxwell X & Y Electric field
        complex*16, intent(in   ) :: Ezz(:), Vrr(:)        ! Maxwell Z Electric Field and Free Charge Potential
        real(dp),   intent(inout) :: Edc                   ! QW spatial array
        real(dp),   intent(in   ) :: R(:)                  ! QW spatial array
        complex*16, intent(inout) :: Ex(:), Ey(:)          ! QW X & Y Electric fields
        complex*16, intent(inout) :: Ez(:), Vr(:)          ! QW Z Electric Field and Free Charge Potential
        real(dp),   intent(in   ) :: t                     ! Current time
        integer,    intent(in   ) :: xxx                   ! Time index
        integer                   :: q, Ny

        Ny = size(EY)

        ! Initialize quantum wire fields
        Ex = 0d0
        Ey = 0d0
        Ez = 0d0
        Vr = 0d0
        
        ! Take the fields from the propagation (Exx & Eyy)
        ! and produce the fields in the QW only (Ex & Ey)
        call rescale_1D(RR, Exx(:), R, Ex(:))
        call rescale_1D(RR, Eyy(:), R, Ey(:))
        call rescale_1D(RR, Ezz(:), R, Ez(:))
        call rescale_1D(RR, Vrr(:), R, Vr(:))

        
        !! Calculate Vr(r) from Ey 1D field:  Obselete
        !forall(q=2:size(Ex)-1) Vr(q) = - sum(Ey(1:q-1)+Ey(1:q)) * (YY(2)-YY(1)) / 2d0
        
      
        
        ! Calculate the fields only between -L/2 to L/2
        Ex(:)   = Ex(:)   * QWWindow
        Ey(:)   = Ey(:)   * QWWindow
        Ez(:)   = Ez(:)   * QWWindow
        Vr(:)   = Vr(:)   * QWWindow

        Edc = sum(Ey) / (size(Ey)*0.5)

        ! Make sure the QW y-space fields are real
        Ex = real(Ex)
        Ey = real(Ey)
        Ez = real(Ez)
        Vr = real(Vr)
        
        
        
        ! Perform a Centered Fast Fourier Transform on the fields... centered Gulley-style!
        call FFTG(Ex(:))
        call FFTG(Ey(:))
        call FFTG(Ez(:))
        call FFTG(Vr(:))
    end subroutine Prop2QW



    subroutine QW2Prop(r, Qr, Ex, Ey, Ez, Vr, Px, Py, Pz, re, rh, RR, Pxx, Pyy, Pzz, RhoE, RhoH, w, xxx, WriteFields, Plasmonics)
        real(dp),   intent(in)    :: r(:)                   ! QW Y-spaces
        real(dp),   intent(in   ) :: Qr(:)                  ! QW momentum space
        complex*16, intent(inout) :: Ex(:) , Ey(:), Ez(:)   ! X & Y components QW Electric fields
        complex*16, intent(inout) :: Vr(:)                  ! Free charge potential
        complex*16, intent(inout) :: Px(:) , Py(:), Pz(:)   ! X & Y components QW Polarization fields
        complex*16, intent(inout) :: re(:) , rh(:)          ! X & Y components QW Polarization fields
        real(dp),   intent(in)    :: RR(:)                  ! QW and Propagation Y-spaces
        complex*16, intent(inout) :: Pxx(:), Pyy(:), Pzz(:) ! X & Y propagation Polarizations
        complex*16, intent(inout) :: RhoE(:), RhoH(:)       ! Free charge density
        integer,    intent(in   ) :: w, xxx                 ! Wire and time indeces
        logical,    intent(in   ) :: WriteFields            ! Record fields?
        logical,    intent(in   ) :: Plasmonics             ! Calculate Charge Densities?
        real(dp)                  :: total, dr, dRR
        integer                   :: i

        ! Constants needed for integration below
        dr  = r(2)  - r(1)
        dRR = RR(2) - RR(1)
    
    
        !! Record the Field arrays in the Fourier Space
        !if(WriteFields) call WriteQWFields(Qr, Ex, Ey, Ez, Px, Py, Pz, Re, Rh, 'k', w, xxx)     
        
        ! Inverse Fourier Transform the fields back into the QW y-space
        call iFFTG(Ex(:))
        call iFFTG(Ey(:))
        call iFFTG(Ez(:))
        call iFFTG(Vr(:))
        call iFFTG(Px(:))
        call iFFTG(Py(:))
        call iFFTG(Pz(:))
        call iFFTG(re(:))
        call iFFTG(rh(:))


        ! Make sure the charge densities are well behaved (real, balanced, etc.)
        if(Plasmonics) then
            re = abs(re)
            rh = abs(rh)
            total = (sum(re)*dr + sum(rh)*dr + small) / 2d0
            re = re * total / (sum(re)*dr + small)
            rh = rh * total / (sum(rh)*dr + small)
        end if


        ! Record the QW Field arrays in the real Space
        if(WriteFields) call WriteQWFields(r, Ex, Ey, Ez, Vr, Px, Py, Ez, Re, Rh, 'r', w, xxx)  

        
        ! Rescale the fields to the propagation YY-space scale
        call rescale_1D(r, Px(:), RR, Pxx(:))
        call rescale_1D(r, Py(:), RR, Pyy(:))
        call rescale_1D(r, Pz(:), RR, Pzz(:))
        call rescale_1D(r, re(:), RR, rhoe(:))
        call rescale_1D(r, rh(:), RR, rhoh(:))

        
        if(Plasmonics) then
            ! Make sure the charge densities are well behaved (real, balanced, etc.)
            rhoe = abs(rhoe)
            rhoh = abs(rhoh)
            rhoe = rhoe * total / (sum(rhoe)*dRR + small)
            rhoh = rhoh * total / (sum(rhoh)*dRR + small)
        end if         
    end subroutine QW2Prop
    
    

   subroutine QWPolarization3(y, ky, p, ehint, area, L, Px, Py, Pz, xxx)
        real(dp),    intent(in   ) :: y(:), ky(:)
        complex(dp), intent(in   ) :: p(:,:)
        real(dp),    intent(in   ) :: ehint, area, L
        complex(dp), intent(inout) :: Px(:), Py(:), Pz(:)
        integer                    :: xxx, r, ke, kh

        Px    = 0d0
        Py    = 0d0
        Pz    = 0d0

        !$omp parallel do private(r, ke, kh)
        do r=1, size(y)
          do ke=1, size(ky)
            do kh=1, size(ky)
              Px(r) =  Px(r) + real((p(kh,ke)) * (+Xvc0(kh,ke)) *  expikr(ke,r) * expikrc(kh,r) ) * QWWindow(r) !&
            end do
          end do
        end do
        !$omp end parallel do
        !$omp parallel do private(r, ke, kh)
        do r=1, size(y)
          do ke=1, size(ky)
            do kh=1, size(ky)
              Py(r) =  Py(r) + real((p(kh,ke)) * (Yvc0(kh,ke)) *  expikr(ke,r) * expikrc(kh,r) ) * QWWindow(r) !&
            end do
          end do
        end do
        !$omp end parallel do
        !$omp parallel do private(r, ke, kh)
        do r=1, size(y)
          do ke=1, size(ky)
            do kh=1, size(ky)
              Pz(r) =  Pz(r) + real((p(kh,ke)) * (Zvc0(kh,ke)) *  expikr(ke,r) * expikrc(kh,r) ) * QWWindow(r) !&
            end do
          end do
        end do
        !$omp end parallel do

        Px(:) =  Px * 2 * ehint / area / (L) 
        Py(:) =  Py * 2 * ehint / area / (L) 
        Pz(:) =  Pz * 2 * ehint / area / (L) 

        call FFTG(Px)
        call FFTG(Py)
        call FFTG(Pz)
    end subroutine QWPolarization3



    subroutine WriteSBESolns(ky, ne, nh, C, D, P, Ee, Eh, w, xxx)
        real(dp),         intent(in) :: ky(:)         ! QW momemtum
        complex*16,       intent(in) :: ne(:), nh(:)  ! QW electron/hole occupation #s
        complex*16,       intent(in) :: C(:,:)        ! QW electron/electron coherence
        complex*16,	  intent(in) :: D(:,:)        ! QW hole/hole coherence
        complex*16,	  intent(in) :: P(:,:)        ! QW electron/hole coherence
        complex*16,	  intent(in) :: Ee(:), Eh(:)  ! QW electron/hole energies
        integer,          intent(in) :: w             ! Wire index
        integer,          intent(in) :: xxx           ! Time index
        character(len=50)            :: fmt, wire
        integer                      :: k
        complex(dp)                  :: pkk(size(ky))        
        complex(dp)                  :: px(size(ky))        
        complex(dp)                  :: py(size(ky))        
        
        fmt = '(I2.2)'
        write(wire,fmt) w

        !forall(k=1:size(ky)) pkk(k) =  p(k,k)
        !call printIT(pkk, ky, xxx, 'Wire/p/p.'//trim(wire)//'.kk.')
        call printIT2D(C(:,:), ky, xxx, 'Wire/C/C.'//trim(wire)//'.k.kp.')
        call printIT2D(D(:,:), ky, xxx, 'Wire/D/D.'//trim(wire)//'.k.kp.')
        call printIT2D(P(:,:), ky, xxx, 'Wire/P/P.'//trim(wire)//'.k.kp.')
        
        call printIT( ne, ky, xxx, 'Wire/ne/ne.'//trim(wire)//'.k.')
        call printIT( nh, ky, xxx, 'Wire/nh/nh.'//trim(wire)//'.k.')

        call printIT( Ee, ky, xxx, 'Wire/Ee/Ee.'//trim(wire)//'.k.')
        call printIT( Eh, ky, xxx, 'Wire/Eh/Eh.'//trim(wire)//'.k.')
    end subroutine WriteSBESolns
    

    subroutine WritePLSpectrum(hw, PLS, w, xxx)
        real(dp),         intent(in) :: hw(:)         ! Photon energy
        real(dp),         intent(in) :: PLS(:)        ! Photoluminescence spectrum
        integer,          intent(in) :: w             ! Max wire index
        integer,          intent(in) :: xxx           ! Time index
        complex(dp)                  :: PLS0(size(hw))
        character(len=50)            :: fmt, wire

        PLS0 = 0d0

        PLS0 = PLS
        
        fmt = '(I2.2)'
        write(wire,fmt) w
        
        call printIT( PLS0, hw/e0, xxx, 'Wire/PL/pl.'//trim(wire)//'.hw.')
    end subroutine WritePLSpectrum

    
    subroutine WriteQWFields(QY, Ex, Ey, Ez, Vr, Px, Py, Pz, Re, Rh, sp, w, xxx)
        real(dp),   intent(in) :: QY(:)               ! QY momentum/y-space array
        complex*16, intent(in) :: Ex(:), Ey(:), Ez(:) ! X & Y components QW Electric fields
        complex*16, intent(in) :: Vr(:)               ! Free charge potential
        complex*16, intent(in) :: Px(:), Py(:), Pz(:) ! X & Y components QW Polarization fields
        complex*16, intent(in) :: Re(:), Rh(:)        ! X & Y components QW Polarization fields
        character(len=*), intent(in) :: sp            ! Domain (k or y) label for file name
        integer,    intent(in) :: w                   ! Wire index
        integer,    intent(in) :: xxx                 ! Time index
        character(len=50)      :: fmt, wire
        
        fmt = '(I2.2)'
        write(wire,fmt) w

        call printIT(Ex(:),  Qy, xxx, 'Wire/Ex/Ex.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Ey(:),  Qy, xxx, 'Wire/Ey/Ey.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Ez(:),  Qy, xxx, 'Wire/Ez/Ez.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Vr(:),  Qy, xxx, 'Wire/Vr/Vr.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Px(:),  Qy, xxx, 'Wire/Px/Px.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Py(:),  Qy, xxx, 'Wire/Py/Py.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Pz(:),  Qy, xxx, 'Wire/Pz/Pz.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Re(:),  Qy, xxx, 'Wire/Re/Re.'//trim(wire)//'.'//trim(sp)//'.')
        call printIT(Rh(:),  Qy, xxx, 'Wire/Rh/Rh.'//trim(wire)//'.'//trim(sp)//'.')
        
        call printIT(Rh(:)-Re(:),   Qy, xxx, 'Wire/Rho/Rho.'//trim(wire)//'.'//trim(sp)//'.')


    end subroutine WriteQWFields



    subroutine WritePropFields(y, Ex, Ey, Ez, Vr, Px, Py, Pz, Re, Rh, sp, w, xxx)
        real(dp),   intent(in) :: y(:)                ! QY momentum/y-space array
        complex*16, intent(in) :: Ex(:), Ey(:), Ez(:) ! X & Y components QW Electric fields
        complex*16, intent(in) :: Vr(:)               ! Free charge potential
        complex*16, intent(in) :: Px(:), Py(:), Pz(:) ! X & Y components QW Polarization fields
        complex*16, intent(in) :: Re(:), Rh(:)        ! X & Y components QW Polarization fields
        character(len=*), intent(in) :: sp            ! space label for file name
        integer,    intent(in) :: w                   ! Wire index
        integer,    intent(in) :: xxx                 ! Time index
        
        character(len=50)      :: fmt, wire
        
        fmt = '(I2.2)'
        write(wire,fmt) w

        call printITReal2(Ex(:) ,  y, xxx, 'Prop/Ex/Ex.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Ey(:) ,  y, xxx, 'Prop/Ey/Ey.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Ez(:) ,  y, xxx, 'Prop/Ez/Ez.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Vr(:) ,  y, xxx, 'Prop/Vr/Vr.'//trim(wire)//'.'//trim(sp)//'.')        
        call printITReal2(Px(:) ,  y, xxx, 'Prop/Px/Px.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Py(:) ,  y, xxx, 'Prop/Py/Py.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Pz(:) ,  y, xxx, 'Prop/Pz/Pz.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Re(:),   y, xxx, 'Prop/Re/Re.'//trim(wire)//'.'//trim(sp)//'.')
        call printITReal2(Rh(:),   y, xxx, 'Prop/Rh/Rh.'//trim(wire)//'.'//trim(sp)//'.')
        
        call printIT(Rh(:)-Re(:), y, xxx, 'Prop/Rho/Rho.'//trim(wire)//'.'//trim(sp)//'.')
    end subroutine WritePropFields
  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Subroutines and functions called from within this module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine QWRho5(Qr, kr, R, L, kkp, p, CC, DD, ne, nh, re, rh, xxx, jjj)
        real(dp),    intent(in)    :: Qr(:), kr(:), R(:), L
        integer,     intent(in)    :: kkp(:,:)
        complex(dp), intent(in)    :: ne(:), nh(:)
        complex(dp), intent(in)    :: p(:,:), CC(:,:), DD(:,:)
        complex(dp), intent(inout) :: re(size(Qr))
        complex(dp), intent(inout) :: rh(size(Qr))
        integer,     intent(in   ) :: xxx, jjj
        real(dp)                   :: dkr, NeTotal, NhTotal, dr
        integer			   :: ri, k1, k2, Nk, NQ0, Nk0, Nr

        Nr  = size(R)
        Nk  = size(kr)
        NQ0 = GetArray0Index(Qr)
        Nk0 = GetArray0Index(kr)
        dkr = kr(2)-kr(1)
        dr  =  R(2)-R(1)

        NeTotal  = sum(abs(ne)) + 1d-50
        NhTotal  = sum(abs(nh)) + 1d-50
        
        re = 0d0
        rh = 0d0

        !$omp parallel do private(ri, k1, k2)                
        do ri=1, Nr
            do k2=1, Nk
                do k1=1, Nk
                   re(ri) = re(ri) +       CC(k1,k2)  * expikrc(k1,ri) * expikr(k2,ri) * QWWindow(ri) / (2*L)
                   rh(ri) = rh(ri) + conjg(DD(k1,k2)) * expikrc(k1,ri) * expikr(k2,ri) * QWWindow(ri) / (2*L)
                end do
            end do
        end do
        !$omp end parallel do
        
        re = re - real( sum(re(1:10)) + sum(re(Nr-9:Nr)) ) / 20d0
        rh = rh - real( sum(rh(1:10)) + sum(rh(Nr-9:Nr)) ) / 20d0
       
        re = re * Netotal / (sum(abs(re))*dr + small)
        rh = rh * Nhtotal / (sum(abs(rh))*dr + small)
        
        call FFTG(re)
        Call FFTG(rh)
    end subroutine

   

    subroutine printIT3D(Dx, z, n, file)
        complex(dp), intent(in) :: Dx(:,:,:)
        real(dp),    intent(in) :: z(:)
        integer                 :: n
        character(len=*)        :: file 
        character(len=50)       :: filename, fmt
        integer                 :: i, j, k, u



            fmt = '(I6.6)'

            write(filename,fmt) n
            filename = 'dataQW/'//trim(file)//trim(filename)//'.dat'

            u = n+20

            open(unit=u, file=trim(filename))

            do k=1, size(Dx, 3)
                do j=1, size(Dx,2)
                    do i=1, size(Dx,1)
                        write(u,*) sngl(real(Dx(i,j,k))), sngl(aimag(Dx(i,j,k)))
                    end do
                end do
            end do

            close(u)


    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printITReal(Dx, z, n, file)
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
                write(u,*) sngl(z(i)), sngl(real(Dx(i)))
            end do

            close(u)

    end subroutine

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print field to file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printITReal2(Dx, z, n, file)
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

            if(firsttime) then
                do i=1, size(z)
                    write(u,*) sngl(real(Dx(i))), sngl(z(i))
                end do
                firsttime = .false.
            else
                do i=1, size(z)
                    write(u,*) sngl(real(Dx(i)))
                end do
            end if

            close(u)

    end subroutine


    !! This is just in case I want to print the QW's linear chi 
    function QWChi1(lam, dky, Ee, Eh, area, geh, dcv)
        real(dp)    :: lam, dky, Ee(:), Eh(:), area, geh
        complex(dp) :: dcv
        complex(dp) :: QWChi1
        real(dp)    :: ww
        
        ww = twopi * c0 / lam
        

        QWChi1 =  4 * dcv**2 / eps0 / area * dky / twopi &                      
                    * sum((Ee+Eh) / (Ee + Eh - ii*hbar*geh - hbar*ww)  &
                                  / (Ee + Eh + ii*hbar*geh + hbar*ww))
        
    end function QWChi1
    





    subroutine CalcQWWindow(YY,L)
        real(dp), intent(in) :: YY(:), L
        integer              :: k, Ny
        
        Ny = size(YY)
        
        allocate(QWWindow(Ny))
        
        QWWindow = 1d0
        
        do k=1, Ny
            if(abs(YY(k)) > L/2d0) QWWindow(k) = 0d0
        end do
        
        QWWindow = exp(-(YY/(L/2))**150)

        call printIT(QWWindow*(ii)**0d0,  yy, 0, 'Envl.y.')
    end subroutine
    
   
 
       
    subroutine InitializeQWOptics(RR, L, dcv, kr, Qr, Ee, Eh, ehint, area, gap)
        real(dp), intent(in)    :: RR(:), L
        complex(dp), intent(in) :: dcv
        real(dp), intent(in)    :: kr(:), Qr(:)
        real(dp), intent(in)    :: Ee(:), Eh(:), ehint, area, gap
        integer                 :: ke, kh, Nk, Nr
        real(dp)                :: R0
        
        Nr = size(RR)
        Nk = size(kr)
        
        call CalcQWWindow(RR,L)
        call CalcExpikr(RR,kr)
        !call Calcg(kr, Qr)
        
        dcv0 = dcv
        R0   = sqrt(area/twopi)
        Vol  = L * area / ehint

        allocate(Xcv0(Nk,Nk), Ycv0(Nk,Nk), Zcv0(Nk,Nk))
        allocate(Xvc0(Nk,Nk), Yvc0(Nk,Nk), Zvc0(Nk,Nk))

        do kh=1, Nk
            do ke=1, Nk
                Xcv0(ke,kh) = dcv * (-1)**kh ! gap / (Ee(k)+Eh(kp)+gap)
                Ycv0(ke,kh) = dcv * (-1)**kh ! gap / (Ee(k)+Eh(kp)+gap)
                Zcv0(ke,kh) = 0d0            ! gap / (Ee(k)+Eh(kp)+gap)
            end do
        end do
       
        Xvc0 = conjg(transpose(Xcv0))
        Yvc0 = conjg(transpose(Ycv0))
        Zvc0 = conjg(transpose(Zcv0))
    end subroutine
    

    function Xcv(k,kp)
        integer, intent(in) :: k,kp
        complex(dp)         :: Xcv
        Xcv = Xcv0(k,kp) 
    end function 

    function Ycv(k,kp)
        integer, intent(in) :: k,kp
        complex(dp)         :: Ycv
        Ycv = Ycv0(k,kp) 
    end function 
    
    function Zcv(k,kp)
        integer, intent(in) :: k,kp
        complex(dp)         :: Zcv
        Zcv = Zcv0(k,kp) 
    end function 
       


    subroutine CalcExpikr(y,ky)
        real(dp), intent(in) :: y(:), ky(:)
        integer :: k, r

        allocate(Expikr( size(ky),size(y)))
        allocate(Expikrc(size(ky),size(y)))

        do r=1, size(y)
            do k=1, size(ky)
                Expikr(k,r) = exp(ii * y(r) * ky(k))
            end do
        end do

        Expikrc = conjg(Expikr)
    end subroutine
    
    
   subroutine GetVn1n2(kr, rcv, Hcc, Hhh, Hcv, Vcc, Vvv, Vcv, Vvc)
        real(dp),   intent(in   ) :: kr(:)
        complex*16, intent(in   ) :: rcv(:)
        complex*16, intent(in   ) :: Hcc(:,:), Hhh(:,:), Hcv(:,:)
        complex*16, intent(inout) :: Vcc(:,:), Vvv(:,:), Vcv(:,:), Vvc(:,:)
        
        complex*16       :: Hvv(size(kr),size(kr))
        complex*16       :: Hvc(size(kr),size(kr))
        complex*16       :: rvc(size(rcv))
        integer          :: Nk, k1, k2

        Nk = size(kr)
        
        Vcc = 0d0
        Vvv = 0d0
        Vcv = 0d0
        Vvc = 0d0
        
        Hvv = - transpose(Hhh)
        Hvc = conjg(transpose(Hcv))
        rvc = conjg(rcv)
        
        !$omp parallel do private(k1, k2)                
        do k2=1, Nk
            do k1=1, Nk
                Vcv(k1, k2) = (-ii/hbar) * (rcv(k1) * Hvv(k1,k2) - Hcc(k1,k2) * rcv(k2))
            end do 
        end do
        !$omp end parallel do

        Vvc = conjg(transpose(Vcv))

        !$omp parallel do private(k1, k2)                
        do k2=1, Nk
            do k1=1, Nk
                Vcc(k1, k2) = (-ii/hbar) * ( rcv(k1) * Hvc(k1,k2) - Hcv(k1,k2) * rvc(k2))
            end do 
        end do        
        !$omp end parallel do

        !$omp parallel do private(k1, k2)                
        do k2=1, Nk
            do k1=1, Nk
                Vvv(k1, k2) = (-ii/hbar) * (rvc(k1) * Hcv(k1,k2) - Hvc(k1,k2) * rcv(k2))
            end do 
        end do         
        !$omp end parallel do

    end subroutine 



!    subroutine GetJ(Rr, kr, Volume, Ccc, Chh, Cvc, Hcc, Hhh, Hcv, Jx, Jy, Jz)
!        double precision :: Rr(:), kr(:), Volume
!        complex*16       :: Ccc(:,:), Chh(:,:), Cvc(:,:)
!        complex*16       :: Hcc(:,:), Hhh(:,:), Hcv(:,:), Jx(:), Jy(:), Jz(:) 
!        complex*16       :: Vcc(size(kr), size(kr)), Vvv(size(kr), size(kr))
!        complex*16       :: Vcv(size(kr), size(kr)), Vvc(size(kr), size(kr))
!        complex*16       :: Ccv(size(kr), size(kr)), Cvv(size(kr), size(kr))
!        complex*16       :: Identity(size(kr), size(kr))
!        integer          :: Nk, Nr, k1, k2, k, r
!        
!        Nk = size(kr)
!        Nr = size(Rr)     

!        Identity = 0d0
!        forall(k=1:Nk) Identity(k,k) = 1

!        Cvc = conjg(transpose(Ccv))
!        Cvv = Identity(:,:) - transpose(Chh(:,:))
!        
!        Jx = 0d0
!        Jy = 0d0
!        Jz = 0d0
!        
!        Vcc = 0d0
!        Vcv = 0d0
!        Vvc = 0d0
!        Vvv = 0d0
!       
!        call GetVn1n2(kr, Xcv0(1,:), Hcc, Hhh, Hcv, Vcc, Vvv, Vcv, Vvc)
!        !$omp parallel do private(r)                
!        do r=1, Nr
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            ! THE FOLLOWING NEEDS ITS OWN FUNCTION         
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            do k2=1, Nk
!                do k1=1, Nk
!                     Jx(r) = Jx(r) + sum( expikr(:, r) * expikrc(k1, r) / Volume * QWWindow(r) * &
!                                          ( + Vcc(:, k2)  * Ccc(k1, k2) + Vvv(:, k2)  * Cvv(k1, k2) &
!                                            + Vcv(:, k2)  * Ccv(k1, k2) + Vvc(:, k2)  * Cvc(k1, k2) &
!                                          ) &
!                                        )
!                end do
!            end do
!        end do        
!        !$omp end parallel do
!        Jx = real(Jx)
!        

!        !call GetVn1n2(kr, Ycv0(1,:), Hcc, Hhh, Hcv, Vcc, Vvv, Vcv, Vvc)
!        !$omp parallel do private(r)                
!        do r=1, Nr
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            !! THE FOLLOWING NEEDS ITS OWN FUNCTION         
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            do k2=1, Nk
!               do k1=1, Nk
!                     Jy(r) = Jy(r) + sum( expikr(:, r) * expikrc(k1, r) / Volume * QWWindow(r) * &
!                                          ( + Vcc(:, k2)  * Ccc(k1, k2) + Vvv(:, k2)  * Cvv(k1, k2) &
!                                            + Vcv(:, k2)  * Ccv(k1, k2) + Vvc(:, k2)  * Cvc(k1, k2) &
!                                          ) &
!                                        )
!                end do
!            end do
!        end do 
!        !$omp end parallel do       
!        Jy = real(Jy)
!        
!    end subroutine
!    
!    
!   subroutine QWPolarization4(Rr, kr, ehint, area, L, dt, Ccc, Dhh, phe, Hcc, Hhh, Heh, Px, Py, Pz, xxx)
!        real(dp),    intent(in   ) :: Rr(:), kr(:)
!        real(dp),    intent(in   ) :: ehint, area, L, dt
!        complex(dp), intent(in   ) :: Ccc(:,:), Dhh(:,:), phe(:,:)
!        complex(dp), intent(in   ) :: Hcc(:,:), Hhh(:,:), Heh(:,:)
!        complex(dp), intent(inout) :: Px(:), Py(:), Pz(:)
!        complex(dp)                :: Jx(size(Rr)), Jy(size(Rr)), Jz(size(Rr))
!        real(dp)                   :: V
!        integer                    :: xxx, r, ke, kh

!        V = area * 2*L * ehint
!        
!        call GetJ(Rr, kr, V, Ccc, Dhh, phe, Hcc, Hhh, Heh, Jx, Jy, Jz)        

!        Px = Px + Jx * dt
!        Py = Py + Jy * dt
!        Pz = Pz + Jz * dt

!        call FFTG(Px)
!        call FFTG(Py)
!        call FFTG(Pz)
!        
!    end subroutine QWPolarization4   
    
end module
