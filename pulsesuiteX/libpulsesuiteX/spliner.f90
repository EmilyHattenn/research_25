! module: spliner
! Holds routines for interpolation of arrays.
module spliner
    use types
    implicit none

    ! interface: spline
    ! Overloads the spline functions.
    !
    ! See Also:
    !   <spline_dp>, <spline_dpc>, <spline2_dp>, <spline2_dpc>
    interface spline
        module procedure spline_dp, spline_dpc, spline2_dp, spline2_dpc
    end interface

    ! interface: seval
    ! Overloads the seval functions.
    !
    ! See Also:
    !   <seval_dp>, <seval_dpc>, <seval2_dp>, <seval2_dpc>
    interface seval
        module procedure seval_dp, seval_dpc, seval2_dp, seval2_dpc
    end interface

    interface rescale_1D
        module procedure rescale_1D_dp, rescale_1D_dpc
    end interface

    interface rescale_2D
        module procedure rescale_2D_dp, rescale_2D_dpc
    end interface

contains
    ! subroutine: rescale_2D_dp
    subroutine rescale_2D_dp(x0, y0, z0, x1, y1, z1)
        implicit none
        real(dp), intent(in   ) :: x0(:), y0(:)
        real(dp), intent(in   ) :: x1(:), y1(:)
        real(dp), intent(in   ) :: z0(:,:)
        real(dp), intent(  out) :: z1(:,:)

        real(dp) :: zt(size(x1), size(y0))
        integer :: i

        do i = 1, size(y0)
            call rescale_1D(x0, z0(:,i), x1, zt(:,i))
        enddo

        do i = 1, size(x1)
            call rescale_1D(y0, zt(i,:), y1, z1(i,:))
        enddo
    end subroutine rescale_2D_dp

    ! subroutine: rescale_2D_dpc
    subroutine rescale_2D_dpc(x0, y0, z0, x1, y1, z1)
        implicit none
        real(dp),    intent(in   ) :: x0(:), y0(:)
        real(dp),    intent(in   ) :: x1(:), y1(:)
        complex(dp), intent(in   ) :: z0(:,:)
        complex(dp), intent(  out) :: z1(:,:)

        complex(dp) :: zt(size(x1), size(y0))
        integer :: i

        do i = 1, size(y0)
            call rescale_1D(x0, z0(:,i), x1, zt(:,i))
        enddo

        do i = 1, size(x1)
            call rescale_1D(y0, zt(i,:), y1, z1(i,:))
        enddo

    end subroutine rescale_2D_dpc

    ! subroutine: rescale_1D_dpc
    subroutine rescale_1D_dpc(x0, z0, x1, z1)
        implicit none
        real(dp),     intent(in   ) :: x0(:), x1(:)
        complex(dp), intent(in   ) :: z0(:)
        complex(dp), intent(  out) :: z1(:)

        complex(dp) :: b(size(x0)), c(size(x0)), d(size(x0))
        integer :: i

        if (size(x0) /= size(z0)) stop "Bad sizes in rescale_1D_dp"

        call spline_dpc(x0, z0, b, c, d)

        do i = 1, size(x1)
            if (x1(i) >= minval(x0) .and. x1(i) <= maxval(x0)) then
                z1(i) = seval_dpc(x1(i), x0, z0, b, c, d)
            else
                z1(i) = 0
            end if
        enddo
    end subroutine rescale_1D_dpc

    subroutine rescale_1D_cyl_dpc(x0, z0, x1, z1)
        implicit none
        real(dp),     intent(in   ) :: x0(:), x1(:)
        complex(dp), intent(in   ) :: z0(:)
        complex(dp), intent(  out) :: z1(:)

        complex(dp) :: b(size(x0)), c(size(x0)), d(size(x0))
        integer :: i

        if (size(x0) /= size(z0)) stop "Bad sizes in rescale_1D_dp"

        call spline_dpc(x0, z0, b, c, d)

        z1(1) = z0(1)

        do i = 2, size(x1)
            if (x1(i) >= minval(x0) .and. x1(i) <= maxval(x0)) then
                z1(i) = seval_dpc(x1(i), x0, z0, b, c, d)
            else
                z1(i) = 0
            end if
        enddo
    end subroutine rescale_1D_cyl_dpc

    ! subroutine: rescale_1D_dp
    subroutine rescale_1D_dp(x0, y0, x1, y1)
        implicit none
        real(dp), intent(in   ) :: x0(:), y0(:), x1(:)
        real(dp), intent(  out) :: y1(:)

        real(dp) :: b(size(x0)), c(size(x0)), d(size(x0))
        integer :: i

        if (size(x0) /= size(y0)) stop "Bad sizes in rescale_1D_dp"

        call spline(x0, y0, b, c, d)

        do i = 1, size(x1)
            if (x1(i) >= minval(x0) .and. x1(i) <= maxval(x0)) then
                y1(i) = seval(x1(i), x0, y0, b, c, d)
            else
                y1(i) = 0
            end if
        enddo
    end subroutine rescale_1D_dp

    ! function: GetValAt_3D
    ! Interpolates a 3D array at an arbitrary point.
    !
    ! Uses polynomial interpolation, <polint3>, to interpolate the array.
    !
    ! Parameters:
    !   e - (3D array) the original 3D array
    !   x0a - (1D array) the positions of the first dimension of e
    !   x1a - (1D array) the positions of the second dimension of e
    !   x2a - (1D array) the positions of the third dimension of e
    !   x0 - the location of the result in the first dimension
    !   x1 - the location of the result in the second dimension
    !   x2 - the location of the result in the third dimension
    !   N - (optional) the order of the interpolation (default: 2)
    real(dp) function GetValAt_3D(e, x0a, x1a, x2a, x0, x1, x2, N) result(Z)
        use nrutils
        real(dp), intent(in   ) :: e(:,:,:)
        real(dp), intent(in   ) :: x0a(:), x1a(:), x2a(:)
        real(dp), intent(in   ) :: x0, x1, x2
        integer,  intent(in   ), optional :: N

        integer, parameter :: N0 = 2

        integer  :: i, j, k
        integer  :: i0, j0, k0
        integer  :: i1, j1, k1
        integer  :: Nt

        Nt = N0

        if (present(N)) Nt = N

        Z = 0.0_dp

        if (x0 < minval(x0a) .or. x0 > maxval(x0a)) return
        if (x1 < minval(x1a) .or. x1 > maxval(x1a)) return
        if (x2 < minval(x2a) .or. x2 > maxval(x2a)) return

        i = iminloc(abs(x0a - x0))
        j = iminloc(abs(x1a - x1))
        k = iminloc(abs(x2a - x2))

        i0 = i - Nt
        j0 = j - Nt
        k0 = k - Nt

        i1 = i + Nt
        j1 = j + Nt
        k1 = k + Nt

        if (i0 < 1) i0 = 1
        if (j0 < 1) j0 = 1
        if (k0 < 1) k0 = 1

        if (i1 > size(x0a)) i1 = size(x0a)
        if (j1 > size(x1a)) j1 = size(x1a)
        if (k1 > size(x2a)) k1 = size(x2a)

        Z = polint3(x0a(i0:i1), x1a(j0:j1), x2a(k0:k1), e(i0:i1,j0:j1,k0:k1), x0, x1, x2)
    end function GetValAt_3D

    ! function: GetValAt_2D
    ! Interpolates a 2D array at an arbitrary point.
    !
    ! Uses polynomial interpolation, <polint3>, to interpolate the array.
    !
    ! Parameters:
    !   e - (2D array) the original 3D array
    !   x0a - (1D array) the positions of the first dimension of e
    !   x1a - (1D array) the positions of the second dimension of e
    !   x0 - the location of the result in the first dimension
    !   x1 - the location of the result in the second dimension
    !   N - (optional) the order of the interpolation (default: 2)
    real(dp) function GetValAt_2D(e, x0a, x1a, x0, x1, N) result(Z)
        use nrutils
        real(dp), intent(in   ) :: e(:,:)
        real(dp), intent(in   ) :: x0a(:), x1a(:)
        real(dp), intent(in   ) :: x0, x1
        integer,  intent(in   ), optional :: N

        integer, parameter :: N0 = 2

        integer  :: i, j
        integer  :: i0, j0
        integer  :: i1, j1
        integer  :: Nt

        Nt = N0

        if (present(N)) Nt = N

        Z = 0.0_dp

        if (x0 < minval(x0a) .or. x0 > maxval(x0a)) return
        if (x1 < minval(x1a) .or. x1 > maxval(x1a)) return

        i = iminloc(abs(x0a - x0))
        j = iminloc(abs(x1a - x1))

        i0 = i - Nt
        j0 = j - Nt

        i1 = i + Nt
        j1 = j + Nt

        if (i0 < 1) i0 = 1
        if (j0 < 1) j0 = 1

        if (i1 > size(x0a)) i1 = size(x0a)
        if (j1 > size(x1a)) j1 = size(x1a)

        Z = polint2(x0a(i0:i1), x1a(j0:j1), e(i0:i1,j0:j1), x0, x1)
    end function GetValAt_2D

    ! function: GetValAt_1D
    ! Interpolates a 1D array at an arbitrary point.
    !
    ! Uses cubic spline, <spline> and <seval>, to interpolate the array.
    !
    ! Parameters:
    !   e - (1D array) the original 3D array
    !   x0 - (1D array) the positions of the first dimension of e
    !   x1 - the location of the result in the first dimension
    real(dp) function GetValAt_1D(e, x0, x1) result(Z)
        real(dp), intent(in   ) :: e(:)
        real(dp), intent(in   ) :: x0(:)
        real(dp), intent(in   ) :: x1

        real(dp) :: b(size(e)), c(size(e)), d(size(e))

        Z = 0.0_dp

        if (x1 < minval(x0) .or. x1 > maxval(x0)) return

        call spline(  x0, e, b, c, d)
        Z = seval(x1, x0, e, b, c, d)
    end function GetValAt_1D
	
	
	
	! function: GetValAt_1D_dpc
    ! Interpolates a complex 1D array at an arbitrary point.
    !
    ! Uses cubic spline, <spline> and <seval>, to interpolate the array.
    !
    ! Parameters:
    !   e - (1D array) the original 3D array
    !   x0 - (1D array) the positions of the first dimension of e
    !   x1 - the location of the result in the first dimension
    complex(dp) function GetValAt_1D_dpc(e, x0, x1) result(Z)
        complex(dp), intent(in   ) :: e(:)
        real(dp),    intent(in   ) :: x0(:)
        real(dp),    intent(in   ) :: x1

        complex(dp) :: b(size(e)), c(size(e)), d(size(e))

        Z = 0.0_dp

        if (x1 < minval(x0) .or. x1 > maxval(x0)) return

        call spline(  x0, e, b, c, d)
        Z = seval(x1, x0, e, b, c, d)
    end function GetValAt_1D_dpc
	

    ! function: polint3
    ! Interpolates a 3D array at an arbitrary point.
    !
    ! Uses polynomial interpolation, <polint2> and <polint1>, to interpolate the array.
    !
    ! The polynomial order is determined by the array size (N-1).
    !
    ! Parameters:
    !   x1a - (1D array) the positions of the first dimension of e
    !   x2a - (1D array) the positions of the second dimension of e
    !   x3a - (1D array) the positions of the third dimension of e
    !   ya - (3D array) the original 3D array
    !   x1 - the location of the result in the first dimension
    !   x2 - the location of the result in the second dimension
    !   x3 - the location of the result in the third dimension
    real(dp) function polint3(x1a, x2a, x3a, ya, x1, x2, x3, dy) result(y)
        implicit none
        real(dp), intent(in   ) :: x1a(:), x2a(:), x3a(:)
        real(dp), intent(in   ) :: ya(:,:,:)
        real(dp), intent(in   ) :: x1, x2, x3
        real(dp), intent(  out), optional :: dy
        integer :: j, m
        real(dp) :: ymtmp(size(x1a))
        m = size(x1a)
        do j = 1, m
            ymtmp(j) = polint2(x2a, x3a, ya(j,:,:), x2, x3)
        end do

        if (present(dy)) then
            y = polint1(x1a,ymtmp,x1,dy)
        else
            y = polint1(x1a, ymtmp, x1)
        end if
    end function polint3

    ! function: polint2
    ! Interpolates a 2D array at an arbitrary point.
    !
    ! Uses polynomial interpolation, <polint1>, to interpolate the array.
    !
    ! The polynomial order is determined by the array size (N-1).
    !
    ! Parameters:
    !   x1a - (1D array) the positions of the first dimension of e
    !   x2a - (1D array) the positions of the second dimension of e
    !   ya - (2D array) the original 3D array
    !   x1 - the location of the result in the first dimension
    !   x2 - the location of the result in the second dimension
    real(dp) function polint2(x1a, x2a, ya, x1, x2, dy) result(y)
        implicit none
        real(dp), intent(in   ) :: x1a(:), x2a(:)
        real(dp), intent(in   ) :: ya(:,:)
        real(dp), intent(in   ) :: x1,x2
        real(dp), intent(  out), optional :: dy
        integer :: j,m
        real(dp) :: ymtmp(size(x1a))
        m = size(x1a)
        do j = 1, m
            ymtmp(j) = polint1(x2a, ya(j,:), x2, dy)
        end do
        y = polint1(x1a, ymtmp, x1, dy)
    end function polint2

    ! function: polint1
    ! Interpolates a 1D array at an arbitrary point.
    !
    ! The polynomial order is determined by the array size (N-1).
    !
    ! Parameters:
    !   x1a - (1D array) the positions of the first dimension of e
    !   ya - (1D array) the original 3D array
    !   x1 - the location of the result in the first dimension
    real(dp) function polint1(xa, ya, x, dy) result(y)
        use nrutils, only : iminloc,nrerror
        implicit none
        real(dp), intent(in   ) :: xa(:), ya(:)
        real(dp), intent(in   ) :: x
        real(dp), intent(  out), optional :: dy
        integer :: m,n,ns
        real(dp) :: c(size(xa)), d(size(xa)), den(size(xa)), ho(size(xa)), dyt
        n = size(xa)
        c = ya
        d = ya
        ho = xa - x
        ns = iminloc(abs(x - xa))
        y = ya(ns)
        ns = ns - 1
        do m = 1, n-1
            den(1:n-m) = ho(1:n-m) - ho(1+m:n)
            if (any(den(1:n-m) == 0.0)) &
                call nrerror('polint: calculation failure')
            den(1:n-m) = (c(2:n-m+1) - d(1:n-m)) / den(1:n-m)
            d(1:n-m) = ho(1+m:n) * den(1:n-m)
            c(1:n-m) = ho(1:n-m) * den(1:n-m)
            if (2*ns < n-m) then
                dyt = c(ns + 1)
            else
                dyt = d(ns)
                ns = ns - 1
            end if
            y = y + dyt
        end do
        if (present(dy)) dy = dyt
    end function polint1

    ! subroutine: bcuint
    ! Bicubic interpolation.
    !
    ! TODO: Document this.
    real(dp) function bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy1,ansy2) result(ansy)
        real(dp), intent(in   ) :: y(4),y1(4),y2(4),y12(4)
        real(dp), intent(in   ) :: x1l,x1u,x2l,x2u,x1,x2
        real(dp), intent(  out) :: ansy1,ansy2
        integer :: i
        real(dp) :: t,u
        real(dp) :: c(4,4)
        call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
        t=(x1-x1l)/(x1u-x1l)
        u=(x2-x2l)/(x2u-x2l)
        ansy=0.0
        ansy2=0.0
        ansy1=0.0
        do i=4,1,-1
                ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
                ansy2=t*ansy2+(3.0_sp*c(i,4)*u+2.0_sp*c(i,3))*u+c(i,2)
                ansy1=u*ansy1+(3.0_sp*c(4,i)*t+2.0_sp*c(3,i))*t+c(2,i)
        end do
        ansy1=ansy1/(x1u-x1l)
        ansy2=ansy2/(x2u-x2l)
    end function bcuint

    ! subroutine: bcucof
    ! Bicubic interpolation coefficients.
    !
    ! TODO: Document this.
    subroutine bcucof(y,y1,y2,y12,d1,d2,c)
        real(dp), intent(in   ) :: d1,d2
        real(dp), intent(in   ) :: y(4), y1(4), y2(4) ,y12(4)
        real(dp), intent(  out) :: c(4,4)

        real(dp) :: x(16)
        real(dp) :: wt(16,16)
        data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
                8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
                2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
                2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
                -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
                -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
                -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
        x(1:4)=y
        x(5:8)=y1*d1
        x(9:12)=y2*d2
        x(13:16)=y12*d1*d2
        x=matmul(wt,x)
        c=reshape(x,(/4,4/),order=(/2,1/))
    end subroutine bcucof

    ! subroutine: spline2_dpc
    ! Computes the cubic spline interpolant for complex arrays.
    !
    ! Parameters:
    !   x - the x array
    !   z - the tabulated function
    !   z2 - the second derivatives
    subroutine spline2_dpc(x, z, z2)
        real(dp),    intent(in   ) :: x(:)
        complex(dp), intent(in   ) :: z(:)
        complex(dp), intent(  out) :: z2(:)

        complex(dp) :: u(size(x))
        real(dp) :: sig
        complex(dp) :: p
        integer :: N, i

        N = size(x)

        z2(1) = 0.0_dp
        u(1) = 0.0_dp

        do i = 2, N-1
            sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            p = sig * z2(i-1) + 2.0_dp
            z2(i) = (sig - 1.0_dp) / p
            u(i) = (z(i+1) - z(i)) / (x(i+1) - x(i)) - (z(i) - z(i-1)) / (x(i) - x(i-1))
            u(i) = (6.0_dp * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
        end do

        z2(N) = 0.0_dp

        do i = N-1, 1, -1
            z2(i) = z2(i) * z2(i+1) + u(i)
        end do
    end subroutine spline2_dpc

    ! function: seval2_dpc
    ! Evaluates the cubic spline for complex arrays.
    !
    ! Parameters:
    !   x0 - the position to interpolate
    !   x - the x array
    !   z - the tabulated function
    !   z2 - the second derivatives
    complex(dp) function seval2_dpc(x0, x, z, z2) result(z0)
        real(dp),    intent(in) ::  x0, x(:)
        complex(dp), intent(in) ::  z(:), z2(:)

        integer :: i
        real(dp) :: h, a, b

        i = locate(x, x0)

        h = x(i+1) - x(i)
        a = (x(i+1) - x0) / h
        b = (x0 - x(i)) / h

        z0 = a * z(i) + b * z(i+1) + ((a**3 - a) * z2(i) + (b**3 - b) * z2(i+1)) * h**2 / 6.0_dp
    end function seval2_dpc

    ! subroutine: spline2_dp
    ! Computes the cubic spline interpolant for real arrays.
    !
    ! Parameters:
    !   x - the x array
    !   y - the tabulated function
    !   y2 - the second derivatives
    subroutine spline2_dp(x, y, y2)
        real(dp), intent(in   ) :: x(:)
        real(dp), intent(in   ) :: y(:)
        real(dp), intent(  out) :: y2(:)

        real(dp) :: u(size(x))
        real(dp) :: sig, p
        integer :: N, i

        N = size(x)

        y2(1) = 0.0_dp
        u(1) = 0.0_dp

        do i = 2, N-1
            sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            p = sig * y2(i-1) + 2.0_dp
            y2(i) = (sig - 1.0_dp) / p
            u(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))
            u(i) = (6.0_dp * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
        end do

        y2(N) = 0.0_dp

        do i = N-1, 1, -1
            y2(i) = y2(i) * y2(i+1) + u(i)
        end do
    end subroutine spline2_dp

    ! function: seval2_dp
    ! Evaluates the cubic spline for real arrays.
    !
    ! Parameters:
    !   x0 - the position to interpolate
    !   x - the x array
    !   y - the tabulated function
    !   y2 - the second derivatives
    real(dp) function seval2_dp(x0, x, y, y2) result(y0)
        real(dp), intent(in) ::  x0, x(:)
        real(dp), intent(in) ::  y(:), y2(:)

        integer :: i
        real(dp) :: h, a, b

        i = locate(x, x0)

        h = x(i+1) - x(i)
        a = (x(i+1) - x0) / h
        b = (x0 - x(i)) / h

        y0 = a * y(i) + b * y(i+1) + ((a**3 - a) * y2(i) + (b**3 - b) * y2(i+1)) * h**2 / 6.0_dp
    end function seval2_dp

    ! subroutine: spline_dpc
    ! Computes the cubic spline interpolant for complex arrays.
    !
    ! Parameters:
    !   x - the x array
    !   y - the tabulated function
    !   b - the interpolant
    !   c - the interpolant
    !   d - the interpolant
    subroutine spline_dpc (x, y, b, c, d)
        implicit none
        real(dp),     intent(in   ) :: x(:)
        complex(dp), intent(in   ) :: y(:)
        complex(dp), intent(  out) :: b(:), c(:), d(:)

        integer :: n, i
        complex(dp) :: t

        n = size(x)

        if ( n < 2 ) stop 'Not enough points to spline.'
        if ( n < 3 ) then
            b(1) = (y(2)-y(1))/(x(2)-x(1))
            c(1) = cmplx(0.0_dp, 0.0_dp)
            d(1) = cmplx(0.0_dp, 0.0_dp)
            b(2) = b(1)
            c(2) = cmplx(0.0_dp, 0.0_dp)
            d(2) = cmplx(0.0_dp, 0.0_dp)
            return
        endif

        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, n-1
            d(i) = x(i+1) - x(i)
            b(i) = 2.0_dp*(d(i-1) + d(i))
            c(i+1) = (y(i+1) - y(i))/d(i)
            c(i) = c(i+1) - c(i)
        end do

        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = cmplx(0.0_dp, 0.0_dp)
        c(n) = cmplx(0.0_dp, 0.0_dp)

        if ( n /= 3 ) then
            c(1) = c(3) / (x(4) - x(2)) - c(2) / (x(3) - x(1))
            c(n) = c(n-1) / (x(n) - x(n-2)) - c(n-2) / (x(n-1) - x(n-3))
            c(1) =  c(1) * d(1  )**2 / (x(4) - x(1  ))
            c(n) = -c(n) * d(n-1)**2 / (x(n) - x(n-3))
        endif

        do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
        end do

        c(n) = c(n)/b(n)
        do i = n-1, 1, -1
            c(i) = (c(i) - d(i) * c(i+1)) / b(i)
        end do

        b(n) = (y(n) - y(n-1))/d(n-1) + d(n-1)*(c(n-1) + 2.0_dp * c(n))
        do i = 1, n-1
            b(i) = (y(i+1) - y(i)) / d(i) - d(i) * (c(i+1) + 2.0_dp * c(i))
            d(i) = (c(i+1) - c(i)) / d(i)
            c(i) = 3.0_dp * c(i)
        end do
        c(n) = 3.0_sp * c(n)
        d(n) = d(n-1)
    end subroutine spline_dpc

    ! function: seval_dpc
    ! Evaluates the cubic spline for complex arrays.
    !
    ! Parameters:
    !   x0 - the position to interpolate
    !   x - the x array
    !   y - the tabulated function
    !   b - the interpolant
    !   c - the interpolant
    !   d - the interpolant
    complex(dp) function seval_dpc(u, x, y, b, c, d)
        implicit none

        real(dp), intent(in) ::  u, x(:)
        complex(dp), intent(in) ::  y(:), b(:), c(:), d(:)

        integer :: i
        real(dp) :: dx

        i = locate(x, u)

        ! evaluate spline
        dx = u - x(i)
        seval_dpc = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    end function seval_dpc

    ! subroutine: spline_dp
    ! Computes the cubic spline interpolant for real arrays.
    !
    ! Parameters:
    !   x - the x array
    !   y - the tabulated function
    !   b - the interpolant
    !   c - the interpolant
    !   d - the interpolant
    subroutine spline_dp(x, y, b, c, d)
        implicit none
        real(dp), intent(in   ) :: x(:)
        real(dp), intent(in   ) :: y(:)
        real(dp), intent(  out) :: b(:), c(:), d(:)

        integer :: n, i
        real(dp) :: t

        n = size(x)

        if ( n < 2 ) stop 'Not enough points to spline.'
        if ( n < 3 ) then
            b(1) = (y(2)-y(1))/(x(2)-x(1))
            c(1) = 0.0
            d(1) = 0.0
            b(2) = b(1)
            c(2) = 0.0
            d(2) = 0.0
            return
        endif

        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, n-1
            d(i) = x(i+1) - x(i)
            b(i) = 2.*(d(i-1) + d(i))
            c(i+1) = (y(i+1) - y(i))/d(i)
            c(i) = c(i+1) - c(i)
        end do

        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0

        if ( n /= 3 ) then
            c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(1) = c(1)*d(1)**2/(x(4)-x(1))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        endif

        do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
        end do

        c(n) = c(n)/b(n)
        do i = n-1, 1, -1
            c(i) = (c(i) - d(i) * c(i+1)) / b(i)
        end do

        b(n) = (y(n) - y(n-1))/d(n-1) + d(n-1)*(c(n-1) + 2.0 * c(n))
        do i = 1, n-1
            b(i) = (y(i+1) - y(i)) / d(i) - d(i) * (c(i+1) + 2.0 * c(i))
            d(i) = (c(i+1) - c(i)) / d(i)
            c(i) = 3.*c(i)
        end do
        c(n) = 3.0 * c(n)
        d(n) = d(n-1)
    end subroutine spline_dp

    ! function: seval_dp
    ! Evaluates the cubic spline for real arrays.
    !
    ! Parameters:
    !   x0 - the position to interpolate
    !   x - the x array
    !   y - the tabulated function
    !   b - the interpolant
    !   c - the interpolant
    !   d - the interpolant
    real(dp) function seval_dp(u, x, y, b, c, d)
        implicit none

        real(dp), intent(in) ::  u, x(:)
        real(dp), intent(in) ::  y(:), b(:), c(:), d(:)

        integer :: i
        real(dp) :: dx

        i = locate(x, u)

        ! evaluate spline
        dx = u - x(i)
        seval_dp = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    end function seval_dp


    ! function: locate
    ! Locates the position in the array braketing a position.
    !
    ! Returns the position, i, in x where x(i) and x(i+1)
    ! braket u.  If a subsequent call is in the same interval
    ! nosearching is done.  Uses a binary search.
    !
    ! Parameters:
    !   x - the x array
    !   u - the position to braket
    integer function locate(x, u)
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

        locate = i
    end function locate
end module spliner
