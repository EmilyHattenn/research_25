! module: helpers
! Mostly useful functions that did not fit anywhere else.
module helpers
    use constants
    implicit none

    ! interface: sech
    ! Single and double precision sech function.
    !
    ! Implemented in:
    !   <sech_sp>, <sech_dp>
    interface sech
        module procedure sech_sp, sech_dp
    end interface

    ! interface: arg
    ! Returns the angle of a complex number.
    !
    ! Implemented in:
    !   <arg_sp>, <arg_dp>
    interface arg
        module procedure arg_sp, arg_dp
    end interface

    ! interface: gauss
    ! Returns exp(-x^2)
    !
    ! Implemented in:
    !   <gauss_sp>, <gauss_dp>
    interface gauss
        module procedure gauss_sp, gauss_dp
    end interface

    ! interface: magsq
    ! Provides a faster abs(Z)**2
    !
    ! Implemented in:
    !   <magsq_dp>, <magsq_sp>
    interface magsq
        module procedure magsq_dp, magsq_sp
    end interface

    ! interface: constrain
    ! Constrains values between limits.
    !
    ! Implemented in:
    !   <constrain_dp>, <constrain_int>
    interface constrain
        module procedure constrain_dp, constrain_int
    end interface

    ! interface: LinearInterp
    ! Returns the linearaly interplated value of
	! the array f(:) at the position x0.
    !
    ! Implemented in:
    !   <LinearInterp_dp>, <LinearInterp_dpc>
    interface LinearInterp
        module procedure LinearInterp_dp, LinearInterp_dpc
    end interface
	
    ! interface: BilinearInterp
    ! Returns the linearaly interplated value of
	! the array f(:,:) at the position (x0,y0).
    !
    ! Implemented in:
    !   <BilinearInterp_dp>, <BilinearInterp_dpc>
    interface BilinearInterp
        module procedure BilinearInterp_dp, BilinearInterp_dpc
    end interface

    ! interface: TrilinearInterp
    ! Returns the linearaly interplated value of
	! the array f(:,:,:) at the position (x0,y0,z0).
    !
    ! Implemented in:
    !   <TrilinearInterp_dp>, <TrilinearInterp_dpc>
    interface TrilinearInterp
        module procedure TrilinearInterp_dp, TrilinearInterp_dpc
    end interface

    ! interface: dfdt
    ! Returns the first derivative value of
    ! the array f(:) with respect to t at 
    ! the index k with a five-point stencil method.
    !
    ! Implemented in:
    !   <dfdt_dp>, <dfdt_dpc>
    interface dfdt
        module procedure dfdt_dp, dfdt_dpc, dfdt_1D_dp, dfdt_1D_dpc
    end interface


    ! interface: isnan
    ! Provides an isnan function for those compilers that lack it.
    !
    ! Implemented in:
    !   <isnan_sp>, <isnan_dp>
#ifdef NAG
    interface isnan
        module procedure isnan_sp, isnan_dp
    end interface
#endif

#ifdef XLF
    interface isnan
        module procedure isnan_sp, isnan_dp
    end interface
#endif

#ifdef __GNUC__
    interface isnan
        module procedure isnan_sp, isnan_dp
    end interface
#endif



    ! const: ec2
    ! PRIVATE: Stores the result of 2.0 * eps0 * c0 for later use.
    real(dp), private, parameter :: ec2 = 2.0_dp * eps0 * c0

contains
    ! function: AmpToInten
    ! Converts a real amplitute into the intensity
    ! in a medium with refractive index n0.
    ! Default n0 = 1.0
    !
    ! Parameters:
    !   e  - the amplitute
    !   n0 - the linear refractive index
    !
    ! Returns:
    !   The corresponding intensity.
    elemental real(dp) function AmpToInten(e, n0) result(inten)
        real(dp), intent(in) :: e
        real(dp), intent(in), optional :: n0

        if(present(n0)) then
            inten = n0 * ec2 * e**2
        else
            inten = ec2 * e**2
        end if
    end function

    ! function: FldToInten
    ! Converts a complex amplitute into the intensity
    ! in a medium with refractive index n0.
    ! Default n0 = 1.0
    !
    ! Parameters:
    !   e  - the complex amplitute
    !   n0 - the linear refractive index
    !
    ! Returns:
    !   The corresponding intensity.
    elemental real(dp) function FldToInten(e, n0) result(inten)
        complex(dp), intent(in) :: e
        real(dp), intent(in), optional :: n0

        if(present(n0)) then
            inten = n0 * ec2 * magsq_dp(e)
        else
            inten = ec2 * magsq_dp(e)
        end if
    end function


    ! function: IntenToAmp
    ! Converts the medium intensity to a real amplitude.
    ! Default n0 = 1.0
    !
    ! Parameters:
    !   inten - the intensity
    !   n0    - the linear refractive index
    !
    ! Returns:
    !   The corresponding amplitude.
    elemental real(dp) function IntenToAmp(inten, n0) result(amp)
        real(dp), intent(in) :: inten
        real(dp), intent(in), optional :: n0

        if(present(n0)) then
            amp = sqrt(inten / ec2 / n0)
        else
            amp = sqrt(inten / ec2)
        end if
    end function

    ! function: arg_dp
    ! Returns the angle of a complex number wrt the real axis.
    elemental real(dp) function arg_dp(Z) result(arg)
        complex(dp), intent(in) :: Z
        arg=atan2(aimag(Z), real(Z))
    end function

    ! function: arg_sp
    ! Returns the angle of a complex number wrt the real axis.
    elemental real(sp) function arg_sp(Z) result(arg)
        complex(sp), intent(in) :: Z
        arg=atan2(aimag(Z), real(Z))
    end function

    ! function: sech_sp
    ! Single precision sech
    elemental function sech_sp(t)
        implicit none
        real (sp) :: sech_sp
        real (sp), intent (in) :: t
        sech_sp = 1.0_sp / cosh (t)
    end function sech_sp

    ! function: sech_dp
    ! Double precision sech
    elemental function sech_dp(t)
        implicit none
        real (dp) :: sech_dp
        real (dp), intent (in) :: t
        sech_dp = 1.0_dp / cosh (t)
    end function sech_dp

    ! function: gauss_dp
    ! The Gaussian function.
    elemental real(dp) function gauss_dp(x)
        real(dp), intent(in) :: x
        gauss_dp = exp(-x**2)
    end function

    ! function: gauss_sp
    ! The Gaussian function.
    elemental real(sp) function gauss_sp(x)
        real(sp), intent(in) :: x
        gauss_sp = exp(-x**2)
    end function

    ! function: magsq_dp
    ! Computes the magnitude squared of a complex number.
    !
    ! real(Z)^2 + imag(Z)^2
    elemental real(dp) function magsq_dp(Z)
        implicit none
        complex(dp), intent(in) :: Z

        magsq_dp = real(Z)**2 + aimag(Z)**2
    end function magsq_dp

    ! function: magsq_sp
    ! Computes the magnitude squared of a complex number.
    !
    ! real(Z)^2 + imag(Z)^2
    elemental real(sp) function magsq_sp(Z)
        implicit none
        complex(sp), intent(in) :: Z

        magsq_sp = real(Z)**2 + aimag(Z)**2
    end function magsq_sp

    ! function: LAX
    ! Implements the Lax mathod in finite differencing.
    !
    ! Improves the stability of explicit methods.
    pure complex(dp) function LAX(u, i, j, k)
        implicit none
        complex(dp), intent(in) :: u(:,:,:)
        integer,     intent(in) :: i,j,k

        LAX = (u(i-1,j,k) + u(i+1,j,k) &
            + u(i,j-1,k) + u(i,j+1,k) &
            + u(i,j,k-1) + u(i,j,k+1)) / 6.0

    end function

    ! function: noLAX
    ! Function with the same signature as <LAX> for a direct replacement.
    !
    ! Does not implement the Lax method.
    pure complex(dp) function noLAX(u, i, j, k)
        implicit none
        complex(dp), intent(in) :: u(:,:,:)
        integer,     intent(in) :: i,j,k

        noLAX = u(i,j,k)

    end function

    ! function: constrain_dp
    ! Constrains a number between a high and a low value.
    elemental real(dp) function constrain_dp(x, H, L)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: H, L

        constrain_dp = max(min(L,H), min(max(L,H), x))
    end function constrain_dp

    ! function: constrain_int
    ! Constrains a number between a high and a low value.
    elemental integer function constrain_int(x, H, L)
        integer, intent(in) :: x
        integer, intent(in) :: H, L

        constrain_int = max(min(L,H), min(max(L,H), x))
    end function constrain_int

    ! function: l2f
    ! Converts a wavelength into its corresponding frequency
    elemental real(dp) function l2f(lam)
        real(dp), intent(in) :: lam

        l2f = c0 / lam
    end function l2f

    ! function: l2w
    ! Converts a wavelength into its corresponding angular frequency.
    elemental real(dp) function l2w(lam)
        real(dp), intent(in) :: lam

        l2w = 2.0_dp * pi * c0 / lam
    end function l2w

    ! function: w2l
    ! Converts an angular frequency into its corresponding wavelength.
    elemental real(dp) function w2l(w)
        real(dp), intent(in) :: w

        w2l = 2.0_dp * pi * c0 / w
    end function w2l

    ! function: GetSpaceArray
    ! Helper function to calculate the spatial array.
    ! Change this function to change the centering of the field.
    !
    ! TODO: Move this to somewhere more appropriate.
    function GetSpaceArray(N, length) result(X)
        integer,  intent(in) :: N
        real(dp), intent(in) :: length

        real(dp) :: X(N)
        real(dp) :: dl
        integer  :: i

        if (N == 1) then
            X = 0.0
            return
        end if

        dl = length / (N-1)     ! 0 is between the middle two points

        forall (i=1:N) X(i) = real(i-1) * dl - length / 2.0
    end function GetSpaceArray

    ! function: GetKArray
    ! Calculates the wavevector array for the FFT.
    ! TODO: Move this to somewhere more appropriate.
    function GetKArray(N, length)
        integer,  intent(in) :: N
        real(dp), intent(in) :: length

        real(dp) :: GetKArray(N)
        real(dp) :: dl
        integer  :: i

        dl = 2.0 * pi / length

        do i=1, N
            if (i .le. N/2+1) then
                GetKArray(i) = real(i-1) * dl
            else
                GetKArray(i) = real(i-N-1) * dl
            endif
        end do
    end function GetKArray

    ! function: unwrap
    ! Reconstructs a smoth phase from one that has been wrapped by modula 2 pi.
    function unwrap(phase) result(unwrapped)
        real(dp), intent(in) :: phase(:)

        real(dp) :: unwrapped(size(phase))
        real(dp) :: pm1, p0, thr, pi2, cp, dpp
        integer  :: N,i

        N = size(phase)

        unwrapped = 0.0_dp

        pm1   = phase(1)
        unwrapped(1) = pm1
        p0    = 0.0_dp
        thr   = pi - epsilon(pi)
        pi2   = 2.0_dp * pi
        do i=2, N
            cp  = phase(i) +p0
            dpp = cp - pm1
            pm1 = cp
            if (dpp>thr) then
                do while (dpp>thr)
                    p0  = p0 - pi2
                    dpp = dpp - pi2
                end do
            end if
            if (dpp < -thr) then
                do while (dpp < -thr)
                    p0  = p0 + pi2
                    dpp = dpp + pi2
                end do
            end if
            cp    = phase(i) + p0
            pm1   = cp
            unwrapped(i) = cp
        end do
    end function

    ! function: factorial
    ! Computes the factorial of an integer.
    integer elemental function factorial(p) result(l)
        integer, intent(in) :: p
        integer :: i
        l=1
        if (p <= 1) then
            return
        else
            do i = 2, p
                l = l * i
            end do
        end if
    end function
    

    ! function: LinearInterp_dpc
    ! Computes a linear interpolation of a
    ! complex 1D array at a specified position.
    !
    ! Parameters:
    ! f  - 1D array to be interpolated
    ! x  - 1D position array corresponding to 'f'
    ! x0 - Position at which 'f' is to be interpolated
	complex function LinearInterp_dpc(f, x, x0) result(f0)
		complex(dp), intent(in) :: f(:)
		real(dp),    intent(in) :: x(:)
		real(dp),    intent(in) :: x0
		
		integer  :: i
		
		i = locator(x, x0)
		
		f0 = (f(i) * (x(i+1) - x0) + f(i+1) * (x0 - x(i)) )  / (x(i+1) - x(i))
	end function	

    ! function: BilinearInterp_dpc
    ! Computes a linear interpolation of a
    ! complex 2D array at a specified position.
    !
    ! Parameters:
    ! f  - 2D array to be interpolated
    ! x  - 1D X position array corresponding to 'f'
    ! y  - 1D Y position array corresponding to 'f'
    ! x0 - X position at which 'f' is to be interpolated
    ! y0 - Y position at which 'f' is to be interpolated	
	complex function BilinearInterp_dpc(f, x, y, x0, y0) result(f0)
		complex(dp), intent(in) :: f(:,:)
		real(dp),    intent(in) :: x(:), y(:)
		real(dp),    intent(in) :: x0, y0
		
		integer  :: j
		real(dp) :: f1, f2
		
		j  = locator(y, y0)
		f1 = LinearInterp_dpc(f(:,j  ), x, x0)
		f2 = LinearInterp_dpc(f(:,j+1), x, x0)
		
		f0 = (f1 * (y(j+1) - y0) + f2 * (y0 - y(j)) ) / (y(j+1) - y(j))
	end function

    ! function: TrilinearInterp_dpc
    ! Computes a linear interpolation of a
    ! complex 2D array at a specified position.
    !
    ! Parameters:
    ! f  - 2D array to be interpolated
    ! x  - 1D X position array corresponding to 'f'
    ! y  - 1D Y position array corresponding to 'f'
    ! Z  - 1D Z position array corresponding to 'f'
    ! x0 - X position at which 'f' is to be interpolated
    ! y0 - Y position at which 'f' is to be interpolated
    ! z0 - Z position at which 'f' is to be interpolated
	complex function TrilinearInterp_dpc(f, x, y, z, x0, y0, z0) result(f0)
		complex(dp), intent(in) :: f(:,:,:)
		real(dp),    intent(in) :: x(:), y(:), z(:)
		real(dp),    intent(in) :: x0, y0, z0
		
		integer  :: k
		real(dp) :: f1, f2
		
		k  = locator(z, z0)
		f1 = BilinearInterp_dpc(f(:,:,k  ), x, y, x0, y0)
		f2 = BilinearInterp_dpc(f(:,:,k+1), x, y, x0, y0)
		
		f0 = (f1 * (z(k+1) - z0) + f2 * (z0 - z(k)) ) / (z(k+1) - z(k))
	end function

    ! function: LinearInterp_dpc
    ! Computes a linear interpolation of a
    ! complex 1D array at a specified position.
    !
    ! Parameters:
    ! f  - 1D array to be interpolated
    ! x  - 1D position array corresponding to 'f'
    ! x0 - Position at which 'f' is to be interpolated	
	real function LinearInterp_dp(f, x, x0) result(f0)
		real(dp), intent(in) :: f(:)
		real(dp), intent(in) :: x(:)
		real(dp), intent(in) :: x0
		
		integer  :: i
		
		i = locator(x, x0)
		
		f0 = (f(i) * (x(i+1) - x0) + f(i+1) * (x0 - x(i)) )  / (x(i+1) - x(i))
	end function	

    ! function: BilinearInterp_dpc
    ! Computes a linear interpolation of a
    ! real 2D array at a specified position.
    !
    ! Parameters:
    ! f  - 2D array to be interpolated
    ! x  - 1D X position array corresponding to 'f'
    ! y  - 1D Y position array corresponding to 'f'
    ! x0 - X position at which 'f' is to be interpolated
    ! y0 - Y position at which 'f' is to be interpolated	
	real function BilinearInterp_dp(f, x, y, x0, y0) result(f0)
		real(dp), intent(in) :: f(:,:)
		real(dp), intent(in) :: x(:), y(:)
		real(dp), intent(in) :: x0, y0
		
		integer  :: j
		real(dp) :: f1, f2
		
		j  = locator(y, y0)
		f1 = LinearInterp_dp(f(:,j  ), x, x0)
		f2 = LinearInterp_dp(f(:,j+1), x, x0)
		
		f0 = (f1 * (y(j+1) - y0) + f2 * (y0 - y(j)) ) / (y(j+1) - y(j))
	end function

    ! function: TrilinearInterp_dpc
    ! Computes a linear interpolation of a
    ! real 2D array at a specified position.
    !
    ! Parameters:
    ! f  - 2D array to be interpolated
    ! x  - 1D X position array corresponding to 'f'
    ! y  - 1D Y position array corresponding to 'f'
    ! Z  - 1D Z position array corresponding to 'f'
    ! x0 - X position at which 'f' is to be interpolated
    ! y0 - Y position at which 'f' is to be interpolated
    ! z0 - Z position at which 'f' is to be interpolated	
	real function TrilinearInterp_dp(f, x, y, z, x0, y0, z0) result(f0)
		real(dp), intent(in) :: f(:,:,:)
		real(dp), intent(in) :: x(:), y(:), z(:)
		real(dp), intent(in) :: x0, y0, z0
		
		integer  :: k
		real(dp) :: f1, f2
		
		k  = locator(z, z0)
		f1 = BilinearInterp_dp(f(:,:,k  ), x, y, x0, y0)
		f2 = BilinearInterp_dp(f(:,:,k+1), x, y, x0, y0)
		
		f0 = (f1 * (z(k+1) - z0) + f2 * (z0 - z(k)) ) / (z(k+1) - z(k))
	end function

    ! function: locator
    ! Locates the position in the array braketing a position.
    !
    ! Returns the position, i, in x where x(i) and x(i+1)
    ! braket u.  If a subsequent call is in the same interval
    ! nosearching is done.  Uses a binary search.
    !
    ! Parameters:
    !   x - the x array
    !   u - the position to braket
    integer function locator(x, u)
        implicit none
        real(dp), intent(in) :: x(:), u

        integer, save :: i = 1
        integer :: j, k, n

        n = size(x)

        if ( i >= n ) i = n-1

        if ( u < x(i) .or. u > x(i+1) ) then

            ! binary search
            i = 1
            j = n+1
            do
                k = (i+j)/2
                if ( u < x(k) ) j = k
                if ( u >= x(k) ) i = k
                if ( j-i == 1 ) exit
            end do
        endif

        locator = i
    end function locator


    function dfdt_dp(f, dt, k)
	real(dp)    :: f(:)     ! Function of t
	real(dp)    :: dt       ! t differential
	integer     :: k        ! t-index for dfdt
	complex(dp) :: dfdt_dp
	integer     :: N

	N = size(f)

	if(k>2 .and. k<N-1) then
	    dfdt_dp = (-f(k+2) + 8.0_dp*f(k+1) - 8.0_dp*f(k-1) + f(k-2)) / dt / 12.0_dp
	elseif(k==1) then
	    dfdt_dp = (f(2) - 0.0_dp) / 2.0_dp / dt
	elseif(k==2) then
	    dfdt_dp = (f(3) - f(1)  ) / 2.0_dp / dt
	elseif(k==N)then
	    dfdt_dp = (0.0_dp - f(N-1)) / 2.0_dp / dt
	elseif(k==N-1) then
	    dfdt_dp = (f(N) - f(N-2)) / 2.0_dp / dt
	else
	    dfdt_dp = 0.0_dp
	endif
    end function

    function dfdt_dpc(f, dt, k)
	complex(dp) :: f(:)     ! Function of t
	real(dp)    :: dt       ! t differential
	integer     :: k        ! t-index for dfdt
	complex(dp) :: dfdt_dpc
	integer     :: N

	N = size(f)

	if(k>2 .and. k<N-1) then
	    dfdt_dpc = (-f(k+2) + 8.0_dp*f(k+1) - 8.0_dp*f(k-1) + f(k-2)) / dt / 12.0_dp
	elseif(k==1) then
	    dfdt_dpc = (f(2) - 0.0_dp) / 2.0_dp / dt
	elseif(k==2) then
	    dfdt_dpc = (f(3) - f(1)  ) / 2.0_dp / dt
	elseif(k==N)then
	    dfdt_dpc = (0.0_dp - f(N-1)) / 2.0_dp / dt
	elseif(k==N-1) then
	    dfdt_dpc = (f(N) - f(N-2)) / 2.0_dp / dt
	else
	    dfdt_dpc = 0.0_dp
	endif
    end function


    function dfdt_1D_dp(f, dt)
	real(dp)    :: f(:)     ! Function of t
	real(dp)    :: dt       ! t differential
	integer     :: k        ! t-index for dfdt
	real(dp)    :: dfdt_1D_dp
	integer     :: N

	N = size(f)

        do k=1,N
            dfdt_1D_dp = dfdt_dp(f,dt,k)
        end do
    end function


    function dfdt_1D_dpc(f, dt)
	complex(dp) :: f(:)     ! Function of t
	real(dp)    :: dt       ! t differential
	integer     :: k        ! t-index for dfdt
	complex(dp) :: dfdt_1D_dpc
	integer     :: N

	N = size(f)

        do k=1,N
            dfdt_1D_dpc = dfdt_dpc(f,dt,k)
        end do
    end function


	
#ifdef NAG
    ! function: isnan_dp
    ! A double precision isnan.
    !
    ! Wraps the native isnan.
    elemental logical function isnan_dp(X) result(Y)
        use, intrinsic :: ieee_arithmetic
        real(dp), intent(in) :: X
        Y = ieee_is_nan(X)
    end function

    ! function: isnan_sp
    ! A single precision isnan.
    !
    ! Wraps the native isnan.
    elemental logical function isnan_sp(X) result(Y)
        use, intrinsic :: ieee_arithmetic
        real(sp), intent(in) :: X
        Y = ieee_is_nan(X)
    end function
#endif

#ifdef XLF
    elemental logical function isnan_dp(X) result(Y)
        use ieee_arithmetic
        real(dp), intent(in) :: X
        Y = ieee_is_nan(X)
    end function

    elemental logical function isnan_sp(X) result(Y)
        use ieee_arithmetic
        real(sp), intent(in) :: X
        Y = ieee_is_nan(X)
    end function
#endif

#ifdef __GNUC__
    elemental logical function isnan_dp(X) result(Y)
        real(dp), intent(in) :: X
        Y = .false.                 ! TODO: Fix this.
    end function

    elemental logical function isnan_sp(X) result(Y)
        real(sp), intent(in) :: X
        Y = .false.                 ! TODO: Fix this.
    end function
#endif
end module helpers
