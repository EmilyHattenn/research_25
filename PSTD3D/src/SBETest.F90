program SBETest
    use types
    use fftw
    use constants
    use helpers
    use phost
    use SBEs
    use usefulsubs

    implicit none


    integer                  :: Nr  = 100               ! Number pixels along the quantum wire direction
    real(dp)                 :: drr = 10d-9             ! Pixel size along the quantum wire direction (m)
    real(dp)                 :: n0  = 3.1               ! Background refractive index
    
    integer                  :: Nt  = 10000             ! Number pixels in time
    real(dp)                 :: dt  = 10d-18            ! Pixel size in time (s)
    real(dp)                 :: t   = 0d0               ! Time variable, starting at t=0 (s)
    
    real(dp)                 :: E0x = 1d7               ! Peak Ex-field value (V/m)
    real(dp)                 :: twx = 10d-15            ! Ex-field pulsewidth (s)
    real(dp)                 :: tpx = 50d-15            ! Ex-field time of pulse peak value at origin (s)
    real(dp)                 :: lamX= 800d-9            ! Field wavelength for Ex (m)

    real(dp)                 :: E0y = 0                 ! Peak Ey-field value (V/m)
    real(dp)                 :: twy = 1                 ! Ey-field pulsewidth (s)
    real(dp)                 :: tpy = 0                 ! Ey-field time of pulse peak value at origin (s)
    real(dp)                 :: lamY= 0                 ! Field wavelength for Ey (m)
    
    real(dp)                 :: E0z = 0                 ! Peak Ez-field value (V/m)
    real(dp)                 :: twz = 1                 ! Ez-field pulsewidth (s)
    real(dp)                 :: tpz = 0                 ! Ez-field time of pulse peak value at origin (s)
    real(dp)                 :: lamZ= 0                 ! Field wavelength for Ez (m)    

    complex(dp), allocatable :: Exx(:), Eyy(:), Ezz(:)  ! E-Field vector components along wire (V/m)
    complex(dp), allocatable :: Pxx(:), Pyy(:), Pzz(:)  ! Polarizations along wire (C/m^2)
    complex(dp), allocatable :: Rho(:), Vrr(:)          ! Charge Density (C/m^2) and Potential (V) along wire

    real(dp),    allocatable :: rr(:)                   ! Spatial Array along the wire (m)
    real(dp),    allocatable :: qrr(:)                  ! Momentum arrays (rad/m)

    real(dp)                 :: w0x, w0y, w0z           ! Angular time frequencies (rad/s), calculated later
    real(dp)                 :: k0x, k0y, k0z           ! Angular space frequencies (rad/m), calculated later
    real(dp)                 :: Tcx, Tcy, Tcz           ! Optical cycles (s), calculated later
    real(dp)                 :: Emax0 = 0d0             ! Peak Electric field value (V/m), calculated later

    integer                  :: n                       ! Time index
    logical                  :: boolF = .false.         ! Dummy boolean value for false
    logical                  :: boolT = .true.          ! Dummy boolean value for false
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate all fields and initialize to zero
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(Exx(Nr), Eyy(Nr), Ezz(Nr))
    allocate(Pxx(Nr), Pyy(Nr), Pzz(Nr))
    allocate(Rho(Nr),  rr(Nr), qrr(Nr), Vrr(Nr))
    call initializefields()


    ! Calculate angular frequencies, & optical cycle (for X-direction only for now...)
    w0x = twopi * c0 / lamX 
    k0x = twopi / lamX * n0
    Tcx = lamX / c0


    ! Calculate the maximum field possible during the simulation:
    Emax0 = sqrt(E0x**2 + E0y**2 + E0z**2)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the real-space & q-space array
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rr  = GetSpaceArray(Nr, (Nr-1) * drr)
    qrr = GetKArray(Nr, Nr * drr)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize the SBEs in SBEs.f90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    Call InitializeSBE(qrr, rr, 0d0, Emax0, lamX, 1,.true.)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Open a files to record data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    open(unit=444,file='fields/Ex.dat')
    open(unit=445,file='fields/Px.dat')


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin Time-Evolving the SBEs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n=1, Nt
    
        ! Update the user on the command line
        print*, n,  Nt

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate E-fields
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Exx = E0x * exp(-(t-tpx)**2 / (twx)**2) * cos(w0x*(t-tpx)) * exp(-(t-tpx)**20 / (2*twx)**20)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Time-Evolve the SBEs from t(n) to t(n+1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call QWCalculator(Exx, Eyy, Ezz, Vrr, rr, qrr, dt, 1, Pxx, Pyy, Pzz, Rho, boolT, boolF)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Print the electric field for the record
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(444,*) t, real(Exx(Nr/2)) ! Record to file in 'fields/Ex.dat'
        write(445,*) t, real(Pxx(Nr/2)) ! Record to file in 'fields/Px.dat'

        t = t + dt
    
    end do

    close(444)
    close(445)
    deallocate(Exx,Eyy,Ezz,Pxx,Pyy,Pzz,Rho,rr,qrr) 
    
contains


    subroutine initializefields()        
        Exx  = 0d0    
        Eyy  = 0d0    
        Ezz  = 0d0    
        Pxx  = 0d0    
        Pyy  = 0d0    
        Pzz  = 0d0    
        Rho  = 0d0
        Vrr  = 0d0    
        rr   = 0d0    
        qrr  = 0d0
    end subroutine initializefields
    
      
    
end program



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
    write(*,"(A)")"Report bugs to Jeremy Gulley <jgulley@kennesaw.edu>"
end subroutine help

! subroutine: version
! Ouputs the version information.
! Uses CVS RCSfile and Revision tags to generate file name and version.
subroutine version
    use pscommandline
end subroutine
