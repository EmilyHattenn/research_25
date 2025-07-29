! module: strings
module strings
    use types
    use intlength
    implicit none

    ! interface: int2str
    ! Converts integer numbers to strings.
    !
    ! Implemented in:
    !   <int2str1>, <int2str2>
    interface int2str
        module procedure int2str1, int2str2
    end interface

    ! interface: bool2str
    ! Converts logicals to strings.
    !
    ! Implemented in:
    !   <bool2str0>, <bool2str1>
    interface bool2str
        module procedure bool2str0, bool2str1
    end interface

    ! interface: real2str
    ! Converts real numbers, single or double precision, to strings.
    !
    ! Implemented in:
    !   <dbl2str1>, <dbl2str2>, <dbl2str3>, <dbl2str4>, <dbl2str5>, <dbl2str6>,
    !   <sgl2str1>, <sgl2str2>, <sgl2str3>, <sgl2str4>
    interface real2str
        module procedure dbl2str1, dbl2str2, dbl2str3, dbl2str4, dbl2str5, dbl2str6
        module procedure sgl2str1, sgl2str2, sgl2str3, sgl2str4, sgl2str5
    end interface

    ! interface: dbl2str
    ! Converts double precision numbers to strings.
    !
    ! Implemented in:
    !   <dbl2str1>, <dbl2str2>, <dbl2str3>, <dbl2str4>, <dbl2str5>, <dbl2str6>
    interface dbl2str
        module procedure dbl2str1, dbl2str2, dbl2str3, dbl2str4, dbl2str5, dbl2str6
    end interface

    ! interface: sgl2str
    ! Converts single presicion numbers to strings.
    !
    ! Implemented in:
    !   <sgl2str1>, <sgl2str2>, <sgl2str3>, <sgl2str4>
    interface sgl2str
        module procedure sgl2str1, sgl2str2, sgl2str3, sgl2str4, sgl2str5
    end interface

    ! interface: wordwrap
    ! Wraps text for output.
    !
    ! Implemented in:
    !   <wordwrap1>, <wordwrap2>
    interface wordwrap
        module procedure wordwrap1, wordwrap2
    end interface

    private :: CalcIntLength, CalcDblLength
contains
    ! function: CalcLines
    ! Calculates the number of lines for word wrapping.
    !
    ! Parameters:
    !   str - the string to wrap
    !   n - the line length
    !
    ! Returns:
    !   The number of lines that <wordwrap1> and <wordwrap2> will produce.
    integer pure function CalcLines(str, n) result(M)
        character(len=*), intent(in) :: str
        integer,          intent(in) :: n

        integer :: i, j

        if (len(str) < n) then
            M = 1
            return
        end if

        M = 1
        i = 0

        do
            if ((i+n) >= len_trim(str)) exit

            j = index(str(max(i,1):min(i+n, len_trim(str))), " ", .true.)

            if (j == 0) j = n

            M = M + 1
            i = i + j
        end do
    end function CalcLines

    ! function: wordwrap1
    ! Wraps text lines for output.
    !
    ! Parameters:
    !   str - the string to wrap
    !   n - the maximum length of a line
    !
    ! Returns:
    !   An array of wrapped lines.
    pure function wordwrap1(str, n) result(W)
        character(len=*), intent(in) :: str
        integer,          intent(in) :: n
        character(len=n+1) :: W(CalcLines(str,n))

        integer :: i, j, k
        character(len=n+1) :: tmp

        W=repeat(" ", n+1)

        i = 1

        do k = 1, size(W)
            tmp = str(i:min(i+n+1, len_trim(str)))
            j = index(tmp, " ", .true.)

            if (j == 0) then
                W(k) = tmp
                j = n+1
            else
                W(k) = tmp(:j)
            end if

            i = i + j
        end do

    end function wordwrap1

    ! function: wordwrap2
    ! Calls <wordwrap1> with a line length of 80.
    pure function wordwrap2(str) result(W)
        character(len=*), intent(in) :: str
        character(len=80) :: W(CalcLines(str,80))
        W = wordwrap1(str,80)
    end function wordwrap2

    ! function: trimb
    ! Trims spaces off both ends of a string.
    pure function trimb(X) result(Y)
        character(len=*), intent(in) :: X
        character(len=len_trim(adjustl(X))) :: Y
        Y = trim(adjustl(X))
    end function

    ! subroutine: writestr
    ! Cannot write a function that uses an internal write.
    !
    ! Fortran95, in its infinite wisdom, prevents recursive
    ! writes.  Even when the inner writes are all internal.
    ! This works around the deficiency by creating a temporary
    ! string for the external write.
    !
    ! Note:
    !   I understand that this has been fixed in F2K.
    !
    ! Parameters:
    !   str - the string to write
    subroutine writestr(str)
        character(len=*), intent(in) :: str
        write(*,"(A)") str
    end subroutine

    ! function: int2str1
    ! Outputs a minimally sized string representation of an integer.
    !
    ! Parameters:
    !   i - the integer
    !
    ! Returns:
    !   The string representation on i.
    pure function int2str1(i) result(str)
        integer, intent(in) :: i
        character(len=CalcIntLength(i)) :: str

        call int2str_sub(i, str)
    end function

    ! subroutine: int2str_sub
    ! Creates a string from an integer without using write.
    !
    ! This allows the various int2str methods to be used in
    ! directly in write statements.
    !
    ! Parameters:
    !   i - the integer
    !   str - the string to store i in
    pure subroutine int2str_sub(i, str, zero_fill, overflow)
        integer,          intent(in ) :: i
        character(len=*), intent(out) :: str
        logical, optional, intent(in ) :: zero_fill
        logical, optional, intent(out) :: overflow

        integer :: n1, n2, k

        if (CalcIntLength(i) > len(str)) then
            str = repeat("*", len(str))
            if (present(overflow)) overflow = .true.
            return
        end if

        n1 = abs(i)

        str = ""

        if (present(zero_fill)) then
            if (zero_fill) str = repeat("0", len(str))
        end if

        do k = len(str), 1, -1
            n2 = (n1 / 10)
            str(k:k) = char(n1 - n2 * 10 + ichar('0'))
            n1 = n2
            if (n1 == 0) exit
        end do

        if (i < 0 .and. k > 1) str((k-1):(k-1)) = '-'
    end subroutine

    ! function: int2str2
    ! Outputs a fixed width string representation of an integer.
    !
    ! Parameters:
    !   i - the integer
    !   n - the string length
    !
    ! Returns:
    !   The string representation on i.
    pure function int2str2(i, n) result(str)
        integer, intent(in) :: i, n
        character(len=(n)) :: str

        call int2str_sub(i, str)
    end function

    ! function: dbl2str1
    ! Converts a double to a string using scientific format.
    !
    ! Parameters:
    !   x - the real value
    !   p - the precision
    !   e - the exponent length
    pure function dbl2str1(x,p,e) result(str)
        real(dp), intent(in) :: x
        integer,  intent(in) :: p
        integer,  intent(in) :: e
        character(len=CalcDblLength(x, p, e, "ES")) :: str

        str = dbl2str6(x, "ES", p, e)
    end function

    ! function: dbl2str1
    ! Converts a double to a string using scientific format.
    !
    ! Parameters:
    !   x - the real value
    !   p - the precision
    !   e - the exponent length
    pure function dbl2str6(x,frmt,p,e) result(str)
        character(len=*), intent(in) :: frmt
        real(dp),         intent(in) :: x
        integer,          intent(in) :: p
        integer,          intent(in) :: e
        character(len=CalcDblLength(x, p, e, frmt)) :: str

        integer :: i, m, g, n
        real(dp) :: x0
        logical :: overflow

        overflow = .false.

        g = floor(log10(abs(x)))
        x0 = x / 10.0_dp**g
        i = int(x0)
        m = abs(nint((x0 - i) * 10.0_dp**p))

        n = 1
        select case(frmt)
            case("ES")
                str = int2str(i); n = len_trim(str) + 1
                str(n:n) = "."; n = n + 1
                call int2str_sub(m, str(n:(p+n-1)), zero_fill=.true., overflow=overflow); n = n + p
                str(n:n+1) = "E+"; if (g < 0) str(n:n+1) = "E-"; n = n + 2
                call int2str_sub(abs(g), str(n:), zero_fill=.true., overflow=overflow)
            case("EN")
                g = floor(log10(abs(x))/3)*3
                x0 = x / 10.0_dp**g
                i = int(x0)
                m = abs(nint((x0 - i) * 10.0_dp**p))
                str = int2str(i); n = len_trim(str) + 1
                str(n:n) = "."; n = n + 1
                call int2str_sub(m, str(n:(p+n-1)), zero_fill=.true., overflow=overflow); n = n + p
                str(n:n+1) = "E+"; if (g < 0) str(n:n+1) = "E-"; n = n + 2
                call int2str_sub(abs(g), str(n:), zero_fill=.true., overflow=overflow)
            case("E","D")
                if (x < 0.0_dp) then
                    str = "-0."; n = 4
                else
                    str = "0."; n = 3
                end if
                m = abs(nint(x0 * 10.0_dp**(p-1)))
                call int2str_sub(m, str(n:(p+n-1)), zero_fill=.true., overflow=overflow); n = n + p
                if (frmt == "D") then
                    str(n:n+1) = "D+"; if (g < 0) str(n:n+1) = "D-"; n = n + 2
                else
                    str(n:n+1) = "E+"; if (g < 0) str(n:n+1) = "E-"; n = n + 2
                end if
                call int2str_sub(abs(g+1), str(n:), zero_fill=.true., overflow=overflow)
            case default
                overflow = .true.
        end select

        if (overflow) str = repeat("*", len(str))
    end function

    ! function: dbl2str2
    ! Converts a double to a string using the specified format.
    !
    ! Parameters:
    !   x - the real value
    !   frmt - the format specifier
    !   n - the length of string to return
    function dbl2str2(x,frmt,n) result(str)
        real(dp), intent(in) :: x
        integer,  intent(in) :: n
        character(len=*), intent(in) :: frmt
        character(len=n) :: str
        write(str,'('//frmt//')') x
    end function

    ! function: dbl2str3
    ! Converts a double array to a string using the specified format.
    !
    ! Parameters:
    !   x - the real array
    !   frmt - the format specifier
    !   n - the length of each element
    function dbl2str3(x,frmt,n) result(str)
        real(dp), intent(in) :: x(:)
        integer,  intent(in) :: n
        character(len=*), intent(in) :: frmt
        character(len=n*size(x)) :: str
        write(str,"(" // int2str(size(x)) // "(" // frmt // "))") x
    end function

    ! function: dbl2str4
    ! Converts a double to a string using <dbl2str1>(x, 5, 3).
    !
    ! Parameters:
    !   x - the real value
    pure function dbl2str4(x) result(str)
        real(dp), intent(in) :: x
        character(len=CalcDblLength(x, 5, 3, "ES")) :: str
        str = dbl2str1(x, 5, 3)
    end function

    ! function: dbl2str5
    ! Converts a double matrix to a string using the specified format.
    !
    ! Parameters:
    !   x - the real matrix
    !   frmt - the format specifier
    !   n - the length of each element
    !
    ! Returns:
    !   A string array representing the matrix
    function dbl2str5(x,frmt,n) result(str)
        real(dp),         intent(in) :: x(:,:)
        character(len=*), intent(in) :: frmt
        integer,          intent(in) :: n
        character(len=n*size(x,1)) :: str(size(x,2))

        integer :: i

        do i = 1, size(x,2)
            str(i) = dbl2str3(x(:,i), frmt, n)
        end do
    end function dbl2str5

    ! function: cmplx2str
    ! Converts a complex to a string using <dbl2str4>.
    !
    ! Parameters:
    !   z - the complex value
    pure function cmplx2str(z) result(str)
        complex(dp), intent(in) :: z
        character(len=(5+3+5)*2+4) :: str
        str = "(" // trimb(dbl2str4(real(z))) // ", " // trimb(dbl2str4(aimag(z))) // ")"
    end function

    ! function: sgl2str1
    ! Converts a single to a string using scientific format.
    !
    ! Parameters:
    !   x - the real value
    !   p - the precision
    !   e - the exponent length
    pure function sgl2str1(x,p,e) result(str)
        real(sp), intent(in) :: x
        integer,  intent(in) :: p,e
        character(len=p+e+5) :: str

        str = dbl2str1(real(x,dp), p, e)
    end function

    ! function: sgl2str2
    ! Converts a single to a string using the specified format.
    !
    ! Parameters:
    !   x - the real value
    !   frmt - the format specifier
    !   n - the length of string to return
    function sgl2str2(x,frmt,n) result(str)
        real(sp), intent(in) :: x
        integer,  intent(in) :: n
        character(len=*), intent(in) :: frmt
        character(len=n) :: str
        write(str,frmt) x
    end function

    ! function: sgl2str3
    ! Converts a single array to a string using the specified format.
    !
    ! Parameters:
    !   x - the real array
    !   frmt - the format specifier
    !   n - the length of each element
    function sgl2str3(x,frmt,n) result(str)
        real(sp), intent(in) :: x(:)
        integer,  intent(in) :: n
        character(len=*), intent(in) :: frmt
        character(len=n*size(x)) :: str
        write(str,"(" // int2str(size(x)) // "(" // frmt // "))") x
    end function

    ! function: sgl2str4
    ! Converts a single to a string using <sgl2str1>(x, 5, 3).
    !
    ! Parameters:
    !   x - the real value
    pure function sgl2str4(x) result(str)
        real(sp), intent(in) :: x
        character(len=5+3+5) :: str
        str = sgl2str1(x, 5, 3)
    end function

    ! function: sgl2str5
    ! Converts a single matrix to a string using the specified format.
    !
    ! Parameters:
    !   x - the real matrix
    !   frmt - the format specifier
    !   n - the length of each element
    !
    ! Returns:
    !   A string array representing the matrix
    function sgl2str5(x,frmt,n) result(str)
        real(sp),         intent(in) :: x(:,:)
        character(len=*), intent(in) :: frmt
        integer,          intent(in) :: n
        character(len=n*size(x,1)) :: str(size(x,2))

        integer :: i

        do i = 1, size(x,2)
            str(i) = sgl2str3(x(:,i), frmt, n)
        end do
    end function sgl2str5

    ! function: bool2str0
    ! Converts a logical to a fixed length string.
    !
    ! Parameters:
    !   x - a logical value
    !   n - the length of the string
    pure function bool2str0(x,n) result(str)
        logical, intent(in) :: x
        integer, intent(in) :: n
        character(len=n) :: str
        str = "F"
        if (x) str = "T"
        str=adjustr(str)
    end function bool2str0

    ! function: bool2str1
    ! Converts a logical to a fixed string.
    !
    ! Parameters:
    !   x - a logical value
    pure function bool2str1(x) result(str)
        logical, intent(in) :: x
        character(len=1) :: str
        str = "F"
        if (x) str = "T"
    end function bool2str1

    ! function: str2str
    function str2str(x,n) result(str)
        character(len=*), intent(in) :: x
        integer,  intent(in) :: n
        character(len=n) :: str
        write(str,"(A)") x
        str=adjustr(str)
    end function

    ! function: toupper
    ! Converts a ASCII string into all upper case.
    function toupper(string) result(upper)
        character(len=*), intent(in) :: string
        character(len=len(string)) :: upper

        integer :: j

        do j = 1,len(string)
            if(string(j:j) >= "a" .and. string(j:j) <= "z") then
                upper(j:j) = achar(iachar(string(j:j)) - 32)
            else
                upper(j:j) = string(j:j)
            end if
        end do
    end function toupper

    ! function: tolower
    ! Converts a ASCII string into all lower case.
    function tolower(string) result(lower)
        character(len=*), intent(in) :: string
        character(len=len(string)) :: lower

        integer :: j

        do j = 1,len(string)
            if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
                lower(j:j) = achar(iachar(string(j:j)) + 32)
            else
                lower(j:j) = string(j:j)
            end if
        end do
    end function tolower

end module strings
