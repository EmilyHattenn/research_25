module rhoPJ
    use types
    use fftw
    use constants
    use helpers
    use phost
    use SBEs
    use usefulsubs
    use typespace
    use typetime
    use typepulse

    implicit none


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TAKE CARE OF INITIAL SETUP PROCEDURES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! reading qwarray.params
    integer , private                 :: Nw   = 4            ! Number of quantum wires
    real(dp), private                 :: y0   = 0.0d-6        ! Wire placement along x-axis [m]
    real(dp), private                 :: x0   = 0.0d-9        ! X-Center of QW array (m)
    real(dp), private                 :: z0   = 0.0d-6        ! Z-Center of QW array (m)
    real(dp), private                 :: dxqw = 0.0d-9      ! Wire spacing along x [m]
    real(dp), private                 :: ay   = 5.0d-9        ! SHO oscillator length in y [m]
    real(dp) , private                :: az   = 5.0d-9        ! SHO oscillator length in z [m]

    logical, private                  :: QW = .true.        ! Turn on the quantum wire(s)
    logical, private                  :: host = .false.       ! Turn on the host dispersion
    logical, private                  :: propagate = .true.   ! Is the propagation on?

    ! dummies
    integer, private         :: Ny, Nx, Nz, i, j, k  ! Number of space x, y, z points
    ! real(dp)                 :: dyy , dely           ! Y-Space step & window (m)
    ! real(dp)                 :: dxx , delx           ! X-Space step & window (m)
    ! real(dp)                 :: dzz, delz            ! Z-Space step & window (m)
    ! real(dp)                 :: tf, t                  ! Time step, final time, current time (s)

    ! Maxwell grid accumulator from Quantum Wire, intermediate vals
    complex(dp), private, allocatable :: PxxOld(:,:), PyyOld(:,:), PzzOld(:,:)

    ! Free charge density storage for continuity eq.
    complex(dp), private, allocatable :: RhoOld(:,:)

    contains



    subroutine InitializeSources(space, time, pulse)
        type(ss), intent(in) :: space
        type(ts), intent(in) :: time
        type(ps), intent(in) :: pulse
        integer              :: iostat
    
        ! Read parameters at the start of the subroutine
        call read_qw_parameters('params/qwarray.params', iostat)
        if (iostat /= 0) then
            print *, "Error reading quantum wire parameters, using defaults"
            QW = .false.  ! Disable QW if there was an error
            return
        end if

        if (.not. allocated(PxxOld)) then
            allocate( PxxOld(GetNx(space), Nw),  &
                      PyyOld(GetNx(space), Nw),  &
                      PzzOld(GetNx(space), Nw),  &
                      RhoOld(GetNx(space), Nw))
            PxxOld = 0d0
            PyyOld = 0d0
            PzzOld = 0d0
            RhoOld = 0d0
        endif        
        
        call InitializeSBE(GetKxArray(space), GetXArray(space), x0, GetAmp(pulse), GetLambda(pulse), Nw, QW)
    end subroutine


    ! Subroutine to read quantum wire parameters from file
    subroutine read_qw_parameters(filename, iostat)
        character(len=*), intent(in)  :: filename
        integer,          intent(out) :: iostat
        character(len=256)           :: errmsg

        open(unit=495, file=filename, status='old', action='read', &
             iostat=iostat, iomsg=errmsg)
        if (iostat /= 0) then
            print *, "Error opening quantum wire parameters file: ", trim(errmsg)
            return
        end if

        read(495, *, iostat=iostat) QW
        read(495, *, iostat=iostat) Nw
        read(495, *, iostat=iostat) x0
        read(495, *, iostat=iostat) dxqw
        read(495, *, iostat=iostat) ay
        read(495, *, iostat=iostat) az

        close(495)
    end subroutine read_qw_parameters



    ! (Porting to 3D)
    subroutine CalcJ(space, time, Ex, Ey, Ez, Jx, Jy, Jz)
        type(ss),    intent(in   ) :: space
        type(ts),    intent(in   ) :: time
        complex(dp), intent(in   ) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        complex(dp), intent(inout) :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
        
        ! Arrays for QWPlacement aiding E longitudinal computations with Rho
        complex(dp) :: Rho(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: ExlfromRho(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: EylfromRho(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: EzlfromRho(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: Px(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: Py(GetNx(space), GetNy(space), GetNz(space))
        complex(dp) :: Pz(GetNx(space), GetNy(space), GetNz(space))


        ! Local Work arrays (not dummy arguments)
        ! 1D arrays for the SBE code
        complex(dp)                :: Exx(GetNx(space)), Eyy(GetNx(space)), Ezz(GetNx(space))   ! No. of points E field points on each axis
        complex(dp)                :: Jxx(GetNx(space)), Jyy(GetNx(space)), Jzz(GetNx(space))   ! No. of points Jfield points on each axis
        complex(dp)                :: Jfx(GetNx(space))                                         ! Temp. array for free current
        ! complex(dp)              :: Vrr(GetNx(space), GetNy(space), GetNz(space))
        complex(dp)                :: Vrr(GetNx(space))

        ! Arrays to be turned 1D after being called through two wires
        complex(dp)                :: Pxx(GetNx(space),Nw), Pyy(GetNx(space),Nw), Pzz(GetNx(space),Nw)
        complex(dp)                :: RhoEH(GetNx(space),Nw)

        ! Real-space and Q-space arrays
        real(dp)                   :: rx(GetNx(space)), ry(GetNy(space)), rz(GetNz(space))        ! Maxwell domain
        real(dp)                   :: qx(GetNx(space)), qy(GetNy(space)), qz(GetNz(space))        ! SBE domain
        real(dp)                   :: x(GetNx(space)) , y(GetNy(space)) , z(GetNz(space))         ! X, Y, and Z Space Arrays (m)

        complex(dp)                :: gx(GetNx(space))
        complex(dp)                :: gy(GetNy(space))
        complex(dp)                :: gz(GetNz(space))
        real(dp)                   :: xqw(Nw)
        integer                    :: NQ0, k, w

        logical                    :: DoQWP  = .true.
        logical                    :: DoQWDl = .true.
        logical                    :: DoQWCurr = .true.

        real(dp) :: ps = 1d-12
        integer  :: iostat, Ny1, Ny2, Nz1, Nz2, Nyz, n
        real(dp) :: dt

        dt = GetDt(time)
        n  = GetN(time)
        
        xqw = 0d0

        ! Read parameters at the start of the subroutine
        call read_qw_parameters('params/qwarray.params', iostat)
        if (iostat /= 0) then
            print *, "Error reading quantum wire parameters, using defaults"
            QW = .false.  ! Disable QW if there was an error
            return
        end if

        if (.not. allocated(PxxOld)) then
            allocate( PxxOld(GetNx(space), Nw),  &
                      PyyOld(GetNx(space), Nw),  &
                      PzzOld(GetNx(space), Nw),  &
                      RhoOld(GetNx(space), Nw))
            PxxOld = 0d0
            PyyOld = 0d0
            PzzOld = 0d0
            RhoOld = 0d0
        endif

        Vrr = 0d0

        rx = GetXArray(space)
        ry = GetYArray(space)
        rz = GetZArray(space)
        qx = GetKXArray(space)
        qy = GetKYArray(space)
        qz = GetKZArray(space)

        Pxx   = 0d0
        Pyy   = 0d0
        Pzz   = 0d0
        Px    = 0d0
        Py    = 0d0
        Pz    = 0d0
        RhoEH = 0d0

        Ny1 = ceiling(0.5*GetNy(space))
        Ny2 =   floor(0.5*GetNy(space)) + 1
        Ny1 = ceiling(0.5*GetNz(space))
        Ny2 =   floor(0.5*GetNz(space)) + 1 
        Nyz = (Ny2-Ny1+1) * (Nz2-Nz1+1)

        ! the full 3D field at the wire center not along the whole wire (xqw(w), y0, z0)
        do w = 1, Nw
            Exx = sum(sum( Ex(:, Ny1:Ny2, Nz1:Nz2), 3), 2) / (1d0*Nyz)
            Eyy = sum(sum( Ey(:, Ny1:Ny2, Nz1:Nz2), 3), 2) / (1d0*Nyz)
            Ezz = sum(sum( Ez(:, Ny1:Ny2, Nz1:Nz2), 3), 2) / (1d0*Nyz)

            ! spatial solve per wire, and temporal solve for internal quantum states (solving in 1D)
            call QWCalculator(Exx, Eyy, Ezz, Vrr, rx, qx, dt, w, Pxx(:,w), Pyy(:,w), Pzz(:,w), RhoEH(:,w), DoQWP, DoQWDl)

        end do

        ! PxxOld handling time evolutions
		forall(i=1:GetNx(space)) Jxx(i) = sum(Pxx(i,:) - PxxOld(i,:)) / dt
		forall(i=1:GetNx(space)) Jyy(i) = sum(Pyy(i,:) - PyyOld(i,:)) / dt
		forall(i=1:GetNx(space)) Jzz(i) = sum(Pzz(i,:) - PzzOld(i,:)) / dt


		! To account for two wires, moved down and fixed all errors? I guess!
        Jxx = Jxx / (1d0*Nw) + CalcJfx(sum(RhoEH,2), sum(RhoOld,2), dt, GetDx(space)) / (1d0*Nw)
        Jyy = Jyy / (1d0*Nw)
        Jzz = Jzz / (1d0*Nw)

        ! gate profiles, where is real gx?  In QWOptics.f90 -> supergaussian in qwwindow
		gx  = 1d0
        gy  = exp( -((ry - y0) ** 2) / (2d0 * ay ** 2) )
        gz  = exp( -((rz - z0) ** 2) / (2d0 * az ** 2) )

        ! computing ElongfromRho using the given subroutine in 3D
        if(DoQWP .and. propagate) then
           call QWPlacement((sum(Pxx, 2) / (1d0*Nw))      , gx, gy, gz, Px)
           call QWPlacement((sum(Pyy, 2) / (1d0*Nw))      , gx, gy, gz, Py)
           call QWPlacement((sum(Pzz, 2) / (1d0*Nw))      , gx, gy, gz, Pz)

           call ElongfromRho(space, Rho, Px, Py, Pz, ExlfromRho, EylfromRho, EzlfromRho)

        end if

        if (DoQWCurr .and. propagate) then
            call QWPlacement( Jxx     , gx, gy, gz,  Jx)
            call QWPlacement( Jyy     , gx, gy, gz,  Jy)
            call QWPlacement( Jzz     , gx, gy, gz,  Jz)
        end if

        if(DoQWDl .and. propagate) then
            call QWPlacement((sum(rhoEH, NW) / (1d0*Nw))     , gx, gy, gz, Rho)
        end if

        PxxOld = Pxx
        PyyOld = Pyy
        PzzOld = Pzz
        RhoOld = RhoEH
    end subroutine CalcJ


    ! computes J_free = - ∫(∂ρ/∂t)​dx
	pure function CalcJfx(RhoNew, RhoPrev, dt, dx)
		complex(dp), intent(in) :: RhoNew(:), RhoPrev(:)
		real(dp),    intent(in) :: dt, dx
	    complex(dp)             :: CalcJfx(size(RhoNew))
	    complex(dp)             :: drhodt(size(RhoNew))
	    integer                 :: i

        ! finite differencing
	    drhodt = (RhoNew - RhoPrev) / dt

	    CalcJfx = 0d0

        ! Using Riemann sums approx. integrals
	    do i=2, size(RhoNew)
			CalcJfx(i) = CalcJfx(i) - drhodt(i) * dx
	    end do
	end function


    ! Fwire(i), Fxy (before)   : complex polarization along the x-axis (length Nx)
    ! gx(i),gy(j),gz(k) : real gate profiles in x,y,z
    ! Fgrid(i,j,k), Fy (before) : complex 3-D Maxwell polarization array to deposit into
    subroutine QWPlacement(Fwire, gx, gy, gz, Fgrid)
        complex(dp), intent(in   ) :: Fwire(:), gx(:), gy(:), gz(:)

        complex(dp), intent(inout) :: Fgrid(:,:,:)

        integer :: Nx, Ny, Nz
        Nx = size(Fgrid,1)
        Ny = size(Fgrid,2)
        Nz = size(Fgrid,3)

        Fgrid = 0d0

        !$omp parallel do private(k, j, i)
        do k=1, Nz
            do j=1, Ny
                do i=1, Nx
                    Fgrid(i,j,k) = Fgrid(i,j,k) + gx(i) * gy(j) * gz(k) * Fwire(i)                  ! with SHO gate in y and z
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine QWPlacement

    ! ElongSeparate

    subroutine ElongSeparate(space, Ex, Ey, Ez, Exl, Eyl, Ezl)
        use typespace
        implicit none

        type(ss),    intent(in)         :: space
        complex(dp), intent(inout)      :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        complex(dp), intent(inout)      :: Exl(:,:,:), Eyl(:,:,:), Ezl(:,:,:)

        integer                         :: i, j, k, Nx, Ny, Nz
        real(dp)                        :: qx(GetNx(space)), qy(GetNy(space)), qz(GetNz(space))
        real(dp)                        :: qsquare
        complex(dp)                     :: qdotE
        real(dp)                        :: small

        Exl = 0d0
        Eyl = 0d0
        Ezl = 0d0

        qx = GetKXArray(space)
        qy = GetKYArray(space)
        qz = GetKZArray(space)

        small = (qx(2)-qx(1)) / 1d4

        Nx = GetNx(space)
        Ny = GetNy(space)
        Nz = GetNz(space)

call fft(Ex)
call fft(Ey)
call fft(Ez)


        !$omp parallel do collapse(3) private(i,j,k,qsquare,qdotE)
        do k=1, Nz
            do j=1, Ny
                do i=1, Nx
                    qsquare    = qx(i) ** 2 + qy(j) ** 2 + qz(k) ** 2 + small
                    qdotE      = qx(i)*Ex(i,j,k) + qy(j)*Ey(i,j,k) + qz(k)*Ez(i,j,k)
                    Exl(i,j,k) = (qdotE * qx(i)) / qsquare
                    Eyl(i,j,k) = (qdotE * qy(j)) / qsquare
                    Ezl(i,j,k) = (qdotE * qz(k)) / qsquare
                end do
            end do
        end do
        !$omp end parallel do

call ifft(Exl)
call ifft(Eyl)
call ifft(Ezl)
    end subroutine ElongSeparate

    ! subroutine ElongfromRho
    subroutine ElongfromRho(space, Rho, Px, Py, Pz, ExlfromRho, EylfromRho, EzlfromRho)
        use typespace
        implicit none

        type(ss),    intent(in)         :: space
        complex(dp), intent(inout)      :: Px(:,:,:), Py(:,:,:), Pz(:,:,:)
        complex(dp), intent(inout)      :: Rho(:,:,:)

        complex(dp), intent(inout)      :: ExlfromRho(:,:,:), EylfromRho(:,:,:), EzlfromRho(:,:,:)

        integer                         :: i, j, k, Nx, Ny, Nz
        real(dp)                        :: qx(GetNx(space)), qy(GetNy(space)), qz(GetNz(space))
        complex(dp)                     :: qdotP, qsquare, rhoq, rhoPqsquare, rho_term, P_term
        real(dp)                        :: small

        qx = GetKXArray(space)
        qy = GetKYArray(space)
        qz = GetKZArray(space)

        Nx = GetNx(space)
        Ny = GetNy(space)
        Nz = GetNz(space)

        ExlfromRho = 0d0
        EylfromRho = 0d0
        EzlfromRho = 0d0

call fft(Rho)
call fft(Px)
call fft(Py)
call fft(Pz)

        small = (qx(2)-qx(1)) / 1d4

        !$omp parallel do collapse(3) private(i, j, k, qsquare, qdotP, rho_term, P_term)
        do k=1, Nz
            do j=1, Ny
                do i=1, Nx
                    qsquare    = qx(i) ** 2 + qy(j) ** 2 + qz(k) ** 2 + small
                    qdotP      = qx(i)*Px(i,j,k) + qy(j)*Py(i,j,k) + qz(k)*Pz(i,j,k)
                    rho_term   = - ii * rho(i,j,k) / (qsquare * eps0 * GetEpsr(space))
                    P_term     = -qdotP / (qsquare * eps0 * GetEpsr(space))

                    ExlfromRho = rho_term * qx(i) + P_term * qx(i)
                    EylfromRho = rho_term * qy(j) + P_term * qy(j)
                    EzlfromRho = rho_term * qz(k) + P_term * qz(k)
                end do
            end do
        end do
        !$omp end parallel do

call ifft(ExlfromRho)
call ifft(EylfromRho)
call ifft(EzlfromRho)

    end subroutine ElongfromRho

end module
