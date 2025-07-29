! module: constants
! Some very useful constants
module constants
    use types
    implicit none

    ! const: pi
    ! The value of pi (3.1415926535897932384626433832795).
    real(dp), parameter :: pi       = 3.141592653589793238462643383279502884197_dp
    ! const: const_e
    ! The value of e (2.7182818284590452353602874713527)
    real(dp), parameter :: const_e  = 2.7182818284590452353602874713527_dp
    ! const: pio2
    ! pi/2
    real(dp), parameter :: pio2     = 1.57079632679489661923132169163975144209858_dp
    ! const: twopi
    ! 2*pi
    real(dp), parameter :: twopi    = 6.283185307179586476925286766559005768394_dp
    ! const: sqrt2
    ! sqrt(2)
    real(dp), parameter :: sqrt2    = 1.41421356237309504880168872420969807856967_dp
    ! const: euler
    real(dp), parameter :: euler    = 0.5772156649015328606065120900824024310422_dp

    ! const: ii
    ! i (sqrt(-1))
    complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)     ! sqrt(-1)

    ! const: c0
    ! The speed of light (299792458.0 m/s)
    real(dp), parameter :: c0       = 299792458.0_dp
    ! const: eps0
    ! The value of eps0 (8.8541878176203898505365630317107e-12 F/m)
    real(dp), parameter :: eps0     = 8.8541878176203898505365630317107e-12_dp
    ! const: mu0
    ! The value of mu0 (1.2566370614359172953850573533118e-6 H/m)
    real(dp), parameter :: mu0      = 1.2566370614359172953850573533118e-6_dp
    ! const: hplank
    ! The Plank constant (6.62606876e-34 J*s)
    real(dp), parameter :: hplank   = 6.62606876e-34_dp
    ! const: hbar
    ! The value of h-bar (1.05457159e-34 J*s)
    real(dp), parameter :: hbar     = 1.05457159e-34_dp
    ! const: e0
    ! The electronic charge determined in 1991
    ! 1.60217733e-19 +- 0.00000049e-19 C
    real(dp), parameter :: e0       = 1.60217733e-19_dp
    ! const: eV
    ! The electron Volt determined in 1991
    ! 1.60217733e-19 +- 0.00000049e-19 J
    real(dp), parameter :: eV       = 1.60217733e-19_dp
    ! const: me0
    ! The rest mass of the electron (kg)
    real(dp), parameter :: me0      = 9.109534e-31_dp

    ! const: efrmt
    ! Format for default display of numerics.
    character(len=*), parameter :: efrmt   = '(SP,1PE23.15E3,1X)'
    ! const: pfrmt
    ! Format for default parameter display of numerics.
    character(len=*), parameter :: pfrmt   = '(ES25.14E3)'
    ! const: efrmt2x
    ! Format for two numeric values
    character(len=*), parameter :: efrmt2x = '(2'//efrmt//')'
    ! const: efrmtA
    ! Format for a numeric followed by a string.
    character(len=*), parameter :: efrmtA  = '('//efrmt//',A)'
    ! const: pfrmtA
    ! Format for a numeric parameter followed by a string.
    character(len=*), parameter :: pfrmtA  = '('//pfrmt//',A)'
    ! const: ifrmtA
    ! Format for a integer parameter followed by a string.
    character(len=*), parameter :: ifrmtA  = '((I25),A)'
    ! const: Aefrmt
    ! Format for a string followed by a numeric.
    character(len=*), parameter :: Aefrmt  = '(A,'//efrmt//')'

    ! const: stdin
    ! File for the standard input.
    character(len=*), parameter :: stdin  = "/dev/stdin"
    ! const: stdout
    ! File for standard output.
    character(len=*), parameter :: stdout = "/dev/stdout"
    ! const: stderr
    ! File for standard error.
    character(len=*), parameter :: stderr = "/dev/stderr"

end module constants
