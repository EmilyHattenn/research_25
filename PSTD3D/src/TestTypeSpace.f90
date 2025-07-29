! program TestTypeSpace
program TestTypeSpace
    use types
    use fftw
    use constants
    use helpers
    use typespace
    implicit none
    

    ! === Space Variables === !
    type(ss)                        :: space
    real(dp),   allocatable         :: Kx(:), Ky(:), Kz(:)
    real(dp),   allocatable         :: x(:), y(:), z(:)
    real(dp)                        :: Dx, Dy, Dz, Vol, dqVol
    real(dp)                        :: Qx, Qy, Qz
    integer                         :: Nx, Ny, Nz   
    
    ! === Field Variables === !
    complex(dp), allocatable        :: e(:,:,:)
    integer                         :: i, j, k
    character(len=256)              :: filename, fnspace


!!! ================================================================ typesapce.f90 tests ================================================================ !!!

    ! === Read from params file === !
    fnspace = 'params/space.params'
    call ReadSpaceParams(fnspace, space)
    
    Nx = GetNx(space)
    Ny = GetNy(space)
    Nz = GetNz(space)
    Dx = GetDx(space)
    Dy = GetDy(space)
    Dz = GetDz(space)

    ! === Print Space Params === !
    print*, "Nx =", Nx
    print*, "Ny =", Ny
    print*, "Nz =", Nz
    print*, "Dx =", Dx
    print*, "Dy =", Dy
    print*, "Dz =", Dz
    
    ! === Width of Each Cartesian Coordinate (Computational Space) === !
    print*, "X Width =", GetXWidth(space)
    print*, "Y Width =", GetYWidth(space)
    print*, "Z Width =", GetZWidth(space)
    

    ! === Cartesian Space Arrays === !
    allocate(x(Nx))
    x = GetXArray(space)
    allocate(y(Ny))
    y = GetYArray(space)
    allocate(z(Nz))
    z = GetZArray(space)

    print*, "First 5 of X Array =", x(1:5)          ! Only prints out first five values in array
    print*, "Last 5 of X Array =", x(Nx-4:Nx)       ! Only prints out last five values in array
    print*, "First 5 of Y Array =", y(1:5)          ! Only prints out first five values in array
    print*, "Last 5 of Y Array =", y(Ny-4:Ny)       ! Only prints out last five values in array
    print*, "First 5 of Z Array =", z(1:5)          ! Only prints out first five values in array
    print*, "Last 5 of Z Array =", z(Nz-4:Nz)       ! Only prints out first five values in array

    
    ! === Fourier Space Arrays === !
    allocate(Kx(Nx))
    Kx = GetKxArray(space)
    allocate(Ky(Ny))
    Ky = GetKyArray(space)
    allocate(Kz(Nz))
    Kz = GetKzArray(space)


    Qx = GetdQx(space)
    Qy = GetdQy(space)
    Qz = GetdQz(space)
    print*, "Qx = ", Qx
    print*, "Qy = ", Qy
    print*, "Qz = ", Qz

    print*, "First 5 of Kx Array =", Kx(1:5)          ! Only prints out first five values in array
    print*, "Last 5 of Kx Array =",  Kx(Nx-4:Nx)      ! Only prints out last five values in array
    print*, "First 5 of Ky Array =", Ky(1:5)          ! Only prints out first 5 values in array
    print*, "Last 5 of Ky Array =",  Ky(Ny-5:Ny)      ! Only prints out last five values in array
    print*, "First 5 of Kz Array =", Kz(1:5)          ! Only prints out first 5 values in array
    print*, "Last 5 of Kz Array =",  Kz(Nz-2:Nz)      ! Only prints out last five values in array
    
    ! === Differential for Conjugate Coordinate System === !
    print*, "dQx =", GetDQx(space)
    print*, "dQy =", GetDQy(space)
    print*, "dQz =", GetDQz(space)
    
    ! === Volume of each Space === !
    print*, "D Vol =", GetDVol(space)
    print*, "dQ Vol =", GetDQVol(space)

    
    ! ===  Create 3DnField, Planes and Lines for gPulse3D === !
    call writeNread(e, 'gPulse3D', space)
    ! === Creates File for Point in 3D Pulse === !
    open(unit=77,file='output/gPulse3D_points.dat')
    do k=1,Nz
        do j=1,Ny
            do i=1,Nx
                write(77,'(4E20.10)') x(i), y(j), z(k), real(e(i,j,k))
            end do
        end do
    end do
    close(77)
    call printLine(e, x, 'gPulse3D', 'x')
    call printLine(e, y,'gPulse3D', 'y')
    call printLine(e, z,'gPulse3D', 'z')
    call printPlane(e, 'gPulse3D', 'xy')
    call printPlane(e,'gPulse3D', 'xz')
    call PrintPlane(e, 'gPulse3D', 'yz')


    ! ===  Create 3DnField, Planes and Lines for gPulseCreation === !
    call writeNread(e, 'gPulseCreation', space)
    call printLine(e, x, 'gPulseCreation', 'x')
    call printLine(e, y,'gPulseCreation', 'y')
    call printLine(e, z,'gPulseCreation', 'z')
    call printPlane(e, 'gPulseCreation', 'xy')
    call printPlane(e,'gPulseCreation', 'xz')
    call PrintPlane(e, 'gPulseCreation', 'yz')

    ! ===  Create 3DnField, Planes and Lines for sech2Pulse === !
    call writeNread(e, 'sech2Pulse', space)
    call printLine(e, x, 'sech2Pulse', 'x')
    call printLine(e, y,'sech2Pulse', 'y')
    call printLine(e, z,'sech2Pulse', 'z')
    call printPlane(e, 'sech2Pulse', 'xy')
    call printPlane(e,'sech2Pulse', 'xz')
    call PrintPlane(e, 'sech2Pulse', 'yz')

    print*, " All Pulse Fields, Planes & Lines Created!!!"


    ! ! ===  Create 3DnField, Planes and Lines for rectanglePulse === !
    ! call writeNread(e, 'rectanglePulse', space)
    ! call printLine(e, x, 'rectanglePulse', 'x')
    ! call printLine(e, y,'rectanglePulse', 'y')
    ! call printLine(e, z,'rectanglePulse', 'z')
    ! call printPlane(e, 'rectanglePulse', 'xy')
    ! call printPlane(e,'rectanglePulse', 'xz')
    ! call PrintPlane(e, 'rectanglePulse', 'yz')


    ! ! ===  Create 3DnField, Planes and Lines for trianglePulse === !
    ! call writeNread(e, 'trianglePulse', space)
    ! call printLine(e, x, 'trianglePulse', 'x')
    ! call printLine(e, y,'trianglePulse', 'y')
    ! call printLine(e, z,'trianglePulse', 'z')
    ! call printPlane(e, 'trianglePulse', 'xy')
    ! call printPlane(e,'trianglePulse', 'xz')
    ! call PrintPlane(e, 'trianglePulse', 'yz')

    ! === Clean Up === !
    deallocate(x, y, z, Kx, Ky, Kz, e)

    contains

!!! ================================================================ typesapce.f90 subroutines ================================================================ !!!

    ! === Read and Writes Space Structure and Fields === !
    subroutine writeNread(e, pulseType, space)
        character(len=*),              intent(in)       :: pulseType
        type(ss),                      intent(inout)    :: space
        complex(dp),      allocatable, intent(inout)    :: e(:,:,:)
        real(dp),         allocatable                   :: x(:), y(:), z(:)
        real(dp)                                        :: EE0, sigma, lam
        real(dp)                                        :: x0, y0, z0
        integer                                         :: Nx, Ny, Nz   
        character(len=256)                              :: filename


        Nx = GetNx(space)
        Ny = GetNy(space)
        Nz = GetNz(space)

        allocate(x(Nx))
        x = GetXArray(space)
        allocate(y(Ny))
        y = GetYArray(space)
        allocate(z(Nz))
        z = GetZArray(space)

        ! === Create Gaussian Field === !
        if (allocated(e)) then 
            deallocate(e)
        end if

        allocate(e(Nx,Ny,Nz))

        ! Read Pulse Params
        !call ReadPulseParms() function for attribute of pulse
        EE0 = 1.25d8     ! Copied from 'params/pulse.params"
        sigma = 1.0d-6   ! Copied from 'params/pulse.params"
        lam = 800d-9     ! Centeral Wavelength
        x0 = 0d0
        y0 = 0d0
        z0 = 0d0


        if (pulseType == 'gPulse3D') then 
            call gPulse3D(e, x, y, z, EE0, sigma, x0, y0, z0, lam)
        else if (pulseType == 'gPulseCreation') then
            call gPulseCreation(e, x, EE0, sigma, x0, lam)
        else if (pulseType == 'rectanglePulse') then
            call rectanglePulse(e, x, EE0, sigma, x0, lam) 
        else if (pulseType == 'trianglePulse') then
            call trianglePulse(e, x, EE0, sigma, x0, lam)
        else if (pulseType == 'sech2Pulse') then
            call sech2Pulse(e, x, EE0, sigma, x0, lam)
        else 
            print*, "Invalid Pulse: ", trim(pulseType)
        end if
        
        ! === Write Field and Space to File === !
        filename = 'output/fields/'//trim(pulseType)//'.dat'
        call writefield(filename, e, space=space, binmode=.false., single=.true., fnspace=fnspace)

        print*, "Field's written to files: ", trim(filename)

        ! === Read Space from Created File === !
        call readspace_only(filename, space=space, binmode=.false., single=.true.)
        print*, "Reading Space from File: ", trim(filename)
        print*, "Nx =", GetNx(space)
        print*, "Ny =", GetNy(space)
        print*, "Nz =", GetNz(space)
        print*, "Dx =", GetDx(space)
        print*, "Dy =", GetDy(space)
        print*, "Dz =", GetDz(space)
        print*, "X Width =", GetXWidth(space)
        print*, "Y Width =", GetYWidth(space)
        print*, "Z Width =", GetZWidth(space)


        ! === Read Field from Created File === !
        call readfield(filename, e, space=space, binmode=.false., single=.true.)
        print*, "Reading Field from File: ", trim(filename)
        print*, e(1,1,1)
        print*, e(Nx/2,Ny/2,Nz/2)
        print*, e(Nx,Ny,Nz)

    end subroutine


    ! === Creates 3D Guassian Pulse === !
    subroutine gPulse3D(e, x, y, z, EE0, sigma, x0, y0, z0, lam)
        complex(dp),    intent(out)         :: e(:,:,:)
        real(dp),       intent(in)          :: x(:), y(:), z(:)
        real(dp),       intent(in)          :: EE0               ! Wave Amplitude
        real(dp),       intent(in)          :: sigma             ! Spatial Width of Wave
        real(dp),       intent(in)          :: lam               ! Wavelength
        real(dp),       intent(in)          :: x0, y0, z0        ! Center of Wave Pulse 
        real(dp)                            :: k0, r2
        integer                             :: Nx, Ny, Nz, i, j, k
        complex(dp)                         :: phase

        Nx = size(e,1)
        Ny = size(e,2)
        Nz = size(e,3)

        k0 = 2* pi /lam                                         ! Centeral Wavevector

        do k=1, Nz
            do j=1, Ny
                do i=1, Nx
                    r2 = ((x(i)-x0)**2 +(y(j)-y0)**2 + (z(k)-z0)**2)
                    e(i,j,k) = EE0 * exp(-(r2 / (2.0d0 * sigma**2))) * cos(k0*x(i))
                end do
            end do
        end do
    end subroutine

    ! === Creates 1D Guassian Pulse === !
    subroutine gPulseCreation(e, x, EE0, sigma, x0, lam)
        complex(dp),    intent(out)         :: e(:,:,:)
        real(dp),       intent(in)          :: x(:)
        real(dp),       intent(in)          :: EE0               ! Wave Amplitude
        real(dp),       intent(in)          :: sigma             ! Spatial Width of Wave
        real(dp),       intent(in)          :: lam               ! Wavelength
        real(dp),       intent(in)          :: x0                ! Center of Wave Pulse 
        real(dp)                            :: k0, r2, wx, c
        integer                             :: Nx, i

        Nx = size(e,1)
        
        c = 3d8
        wx = sigma
        k0 = 2* pi /lam                                          ! Centeral Wavevector

        do i=1, Nx
            r2 = ((x(i)-x0)**2 / wx**2)
            e(i,:,:) = EE0 * exp(-(r2)) * cos(k0*x(i))
        end do
    end subroutine
    
    ! === Creates 1D Rectangular Pulse === !
    subroutine rectanglePulse(e, x, EE0, sigma, x0, lam)
        complex(dp),    intent(out)         :: e(:,:,:)
        real(dp),       intent(in)          :: x(:)
        real(dp),       intent(in)          :: EE0               ! Wave Amplitude
        real(dp),       intent(in)          :: sigma             ! Spatial Width of Wave
        real(dp),       intent(in)          :: lam               ! Wavelength
        real(dp),       intent(in)          :: x0                ! Center of Wave Pulse 
        real(dp)                            :: k0, r, envelope
        integer                             :: Nx, i


        Nx = size(e,1)
        k0 = 2* pi /lam   

        do i=1, Nx
            r = (x(i)- x0)
            if (r <= sigma) then 
                envelope = 1.0 
            else 
                envelope = 0.0 
            end if
            e(i,:,:) = EE0 * envelope * cos(k0*x(i))
        end do
    end subroutine

    ! ! === Creates 1D Triangular Pulse === !
    subroutine trianglePulse(e, x, EE0, sigma, x0, lam)
        complex(dp),    intent(out)         :: e(:,:,:)
        real(dp),       intent(in)          :: x(:)
        real(dp),       intent(in)          :: EE0               ! Wave Amplitude
        real(dp),       intent(in)          :: sigma             ! Spatial Width of Wave
        real(dp),       intent(in)          :: lam               ! Wavelength
        real(dp),       intent(in)          :: x0                ! Center of Wave Pulse 
        real(dp)                            :: k0, r, envelope
        integer                             :: Nx, i

        Nx = size(e,1)
        k0 = 2* pi /lam   

        do i=1, Nx
            r = (x(i)- x0)
            if (r <= sigma) then 
                envelope = 1.0d0 - r / sigma 
            else 
                envelope = 0.0d0
            end if
            e(i,:,:) = EE0 * envelope * cos(k0*x(i))
        end do
    end subroutine

    ! === Creates 1D sech2 Pulse === !
    subroutine sech2Pulse(e, x, EE0, sigma, x0, lam)
        complex(dp),    intent(out)         :: e(:,:,:)
        real(dp),       intent(in)          :: x(:)
        real(dp),       intent(in)          :: EE0               ! Wave Amplitude
        real(dp),       intent(in)          :: sigma             ! Spatial Width of Wave
        real(dp),       intent(in)          :: lam               ! Wavelength
        real(dp),       intent(in)          :: x0                ! Center of Wave Pulse 
        real(dp)                            :: k0, r
        integer                             :: Nx, i

        Nx = size(e,1)
        k0 = 2* pi /lam   

        do i=1, Nx
            r = (x(i)- x0)
            e(i,:,:) = EE0 * (1.0 / cosh(r / sigma))**2 * cos(k0*x(i))
        end do
    end subroutine


    ! === Prints out Lines === !
    subroutine printLine(e, r, pulseType, axis)
        complex(dp),        intent(inout)       :: e(:,:,:)
        real(dp),           intent(inout)       :: r(:)
        character(len=*),   intent(in)          :: pulseType
        character(len=*),   intent(in)          :: axis
        character(len=50)                       :: filename
        integer                                 :: Nx, Ny, Nz
        integer                                 :: i, j, k, u

        Nx = size(e,1)
        Ny = size(e,2)
        Nz = size(e,3)


        if (Nx==1) then 
            Nx = 2
        end if
        if (Ny==1) then 
            Ny = 2
        end if
        if (Nz==1) then 
            Nz = 2
        end if

        u = 101
        if (axis == 'x') then
            filename = 'output/lines/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/x.dat')
            do i=1, Nx 
                write(u,*) r(i), real(e(i, Ny/2, Nz/2))
            end do
        else if (axis == 'y') then 
            
            filename = 'output/lines/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/y.dat')
            do j=1, Ny
                write(u,*) r(j), real(e(Nx/2, j, Nz/2))
            end do
        else if (axis == 'z') then 
            filename = 'output/lines/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/z.dat')
            do k=1, Nz 
                write(u,*) r(k), real(e(Nx/2, Ny/2, k))
            end do
        else 
            print*, "Invalid Axis"
        end if
        close(u)
    end subroutine printLine


    ! === Prints out Planes === !
    subroutine printPlane(e, pulseType, plane)
        complex(dp),        intent(inout)       :: e(:,:,:)
        character(len=*),   intent(in)          :: pulseType
        character(len=*),   intent(in)          :: plane
        integer                                 :: Nx, Ny, Nz
        integer                                 :: i, j, k, u
        character(len=50)                       :: filename

        Nx = size(e,1)
        Ny = size(e,2)
        Nz = size(e,3)


        if (Nx==1) then 
            Nx = 2
        end if
        if (Ny==1) then 
            Ny = 2
        end if
        if (Nz==1) then 
            Nz = 2
        end if

        u = 99

        if (plane == 'xy') then
            filename = 'output/planes/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/xy.dat')
            do j=1, Ny
                do i=1, Nx
                    write(u,*) real(e(i, j, Nz/2))
                end do
            end do
        else if (plane == 'yz') then 
            filename = 'output/planes/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/yz.dat')
            do k=1, Nz
                do j=1, Ny
                    write(u,*) real(e(Nx/2, j, k))
                end do
            end do
        else if (plane == 'xz') then 
            filename = 'output/planes/'//trim(pulseType)
            open(unit=u, file=trim(filename)//'/xz.dat')
            do k=1, Nz
                do i=1, Nx
                    write(u, *) real(e(i, Ny/2, k))
                end do
            end do
        else 
            print*, "Invalid Plane"
        end if
        close(u)
    end subroutine


    
end program TestTypeSpace

! subroutine: usage
! Provides a short usage statement to remind the user of the options.
! Also used at the top of the help section.
subroutine usage
    use pscommandline
    character(len=256) :: nm
end subroutine usage

! subroutine: help
! Provides complete help for the program.
subroutine help
    use pscommandline
    use propagator
    use stdlog
    character(len=256) :: nm
    write(*,"(A)")"Report bugs to Jeremy Gulley <jgulley@furman.edu>"
end subroutine help

! subroutine: version
! Ouputs the version information.
! Uses CVS RCSfile and Revision tags to generate file name and version.
subroutine version
    use pscommandline
end subroutine
