! module: fileio
! This module wraps many file input and ourput functions.
module fileio
    use strings
    use constants
    implicit none

    ! const: BUFSIZE
    ! Default buffer size for many reading routines.
    integer :: BUFSIZE = 1024*1024

    ! const: DEFAULT_ERROR_UNIT
    ! The default unit for <constants::stderr>
    integer, parameter :: DEFAULT_ERROR_UNIT = 0
    ! const: DEFAULT_INPUT_UNIT
    ! The default unit for <constants::stdin>
    integer, parameter :: DEFAULT_INPUT_UNIT = 5
    ! const: DEFAULT_OUTPUT_UNIT
    ! The default unit for <constants::stdout>
    integer, parameter :: DEFAULT_OUTPUT_UNIT = 6

    ! const: NUMBER_OF_PRECONNECTED_UNITS
    ! PRIVATE: The default number of preconnected unit numbers.
    integer, parameter, private :: NUMBER_OF_PRECONNECTED_UNITS = 3
    ! const: PRECONNECTED_UNITS
    ! PRIVATE: The default unit numbers of the preconnected units.
    integer, parameter, private :: PRECONNECTED_UNITS(NUMBER_OF_PRECONNECTED_UNITS) = (/ 0, 5, 6 /)
    ! const: MAX_UNIT_NUMBER
    ! PRIVATE: Maximum unit number to try to allocate.
    integer, parameter, private :: MAX_UNIT_NUMBER = 1000


    interface IsOpen
        module procedure IsOpen_u, IsOpen_fn
    end interface

    interface readmatrix
        module procedure readmatrix_dp
    end interface

contains

    ! subroutine: ReadIniTagStr
    ! Reads a tag from a section of a INI fromatted file.
    !
    ! | [sec]
    ! | tag=
    !
    ! Parameters:
    !   file - the filename of the INI file
    !   sec - the section ([sec]) to look for the tag in
    !   tag - the tag (tag=str) to read
    !   str - the return value of tag
    !   err - (optional) error return
    !   def - (optional) default value
    subroutine ReadIniTagStr(file, sec, tag, str, err, def)
        use logger
        character(len=*),               intent(in   ) :: file
        character(len=*),               intent(in   ) :: sec
        character(len=*),               intent(in   ) :: tag
        character(len=*),               intent(  out) :: str
        integer,            optional,   intent(  out) :: err
        character(len=*),   optional,   intent(in   ) :: def

        integer                 :: U, err0
        character(len=BUFSIZE)      :: line

        err0 = 0

        call debug2("Opening INI file, " // trim(file))
        
        U = new_unit()    

        open(U, file=file, form="FORMATTED", action="READ", status="OLD")

        call debug2("Looking for section, " // "[" // trim(sec) // "]")
        do
            read(U, *, iostat = err0) line
            if (line == "[" // trim(sec) // "]" .or. err0 /= 0) exit                      
            if (line(1:1) == "[") call debug3(trim(line), __FILE__, __LINE__)   
        end do

        if (err0 /= 0) then
            close(U)
            if (present(err)) then
                call debug2("Can't find section, " // trim(sec) // ", of ini file, " // trim(file) // ".", &
                    __FILE__, __LINE__)
                if (present(def)) then
                    str = def
                end if
                err = err0
                return
            end if

            call error("Can't find section, " // trim(sec) // ", of ini file, " // trim(file) // ".", &
                1, __FILE__, __LINE__)
        end if

        call debug2("Looking for tag, " // trim(tag))
        do
            read(U, "(A)", iostat = err0) line
            call debug3(trim(line), __FILE__, __LINE__)
            if (line(1:1) == "[") err0 = -1
            if (err0 /= 0) exit
            if (line(1:len_trim(tag)+1) == tag // "=") exit
        end do

        close(U)

        if (err0 /= 0) then
            if (present(err)) then
                call debug2("Can't find tag, " // trim(tag) // ", of section, " // trim(sec) // ", of ini file, " &
                    // trim(file) // ".", &
                    __FILE__, __LINE__)
                if (present(def)) then
                    str = def
                end if
                err = err0
                return
            end if

            call error("Can't find tag, " // trim(tag) // ", of section, " // trim(sec) // ", of ini file, " &
                // trim(file) // ".", &
                1, __FILE__, __LINE__)
        end if

        str = line(len_trim(tag)+2:len_trim(line))
    end subroutine ReadIniTagStr

    ! function: IsOpen_u
    ! Quick inquiry function to determine if a unit number is open.
    !
    ! Parameters:
    !   u - the unit number
    logical function IsOpen_u(u)
        integer,          intent(in) :: u
        inquire(UNIT=u, OPENED=IsOpen_u)
    end function

    ! function: IsOpen_fn
    ! Quick inquiry function to determine if a filename is open.
    !
    ! Parameters:
    !   fn - the filename
    logical function IsOpen_fn(fn)
        character(len=*), intent(in) :: fn
        inquire(FILE=fn, OPENED=IsOpen_fn)
    end function

    ! function: GetFileParam
    ! Gets the next numeric value from a file, skipping lines beginning with "#".
    !
    ! Parameters:
    !   unit - the unit number to read
    !   iostat - (optional) value of the iostat of the read
    !
    ! Returns:
    !   The first numeric value at the beginning of the line.
    real(dp) function GetFileParam(unit, iostat) result(val)
        integer,           intent(in   ) :: unit
        integer, optional, intent(  out) :: iostat

        character(len=256) :: line

        if (present(iostat)) then
            call GetNextLine(unit, line, iostat)
            if (iostat == 0) read(line, *, iostat=iostat) val
        else
            call GetNextLine(unit, line)
            read(line, *) val
        end if
    end function

    ! subroutine: GetNextLine
    ! Skips lines in file that begin with a "#".
    !
    ! Parameters:
    !   unit - the unit to read from
    !   line - the next uncommented line
    !   iostat - (optional) value of the iostat of the read
    subroutine GetNextLine(unit, line, iostat)
        integer,           intent(in   ) :: unit
        character(len=*),  intent(  out) :: line
        integer, optional, intent(  out) :: iostat

        do
            if (present(iostat)) then
                read(unit, "(A)", iostat=iostat) line
                if (iostat /= 0) return
            else
                read(unit, "(A)") line
            end if
            if (line(1:1) /= "#") exit
        end do
    end subroutine

    ! subroutine: Write1DArray
    ! Writes a 1D array to a file or <constants::stdout>
    !
    ! Parameters:
    !   X - the array to write
    !   fn - the filename to write to.  If fn is <constants::stdout>
    !       write to "*".
    subroutine Write1DArray(X, fn)
        real(dp),         intent(in) :: X(:)
        character(len=*), intent(in) :: fn

        integer :: i

        integer :: U


        if (fn == stdout) then
            do i=1,size(X)
                write(*, efrmt) X(i)
            end do
        else
            U = new_unit()
            open(U, file=fn, form='FORMATTED', action="WRITE")
            do i=1,size(X)
                write(U, efrmt) X(i)
            end do
            close(U)
        end if
    end subroutine

    ! subroutine: readmatrix_dp
    ! Reads a text matrix from a unit.
    !
    ! This function allocates an arbitrary sized
    ! matrix to hold the contents of a file.  The
    ! size does not need to be specified beforehand.
    !
    ! Parameters:
    !   U - the unit to read from
    !   M - a pointer to hold the matrix
    subroutine readmatrix_dp(U, M)
        integer, intent(in   ) :: U
        real(dp), pointer      :: M(:,:)

        integer :: err, i, N1, N2
        character(len=BUFSIZE) :: A
        character(len=50) :: frmt
        real(dp), pointer :: temp(:,:)

        nullify(M)

        ! Read in first line.
        frmt = "(A"//adjustl(int2str(len(A)))//")"
        read(U,frmt) A

        N1 = 2
        do i = 1, len(A)/2
            allocate(M(N1,1))
            read(A, *, iostat=err) M
            if (err /= 0) exit
            N1 = N1 + 1
            deallocate(M)
        end do

        N1 = N1-1

        deallocate(M)
        allocate(M(N1,N1))

        read(A, *) M(1:N1,1)

        N2 = 2

        do
            read(U, *, iostat=err) M(1:N1,N2)
            if (err /= 0) exit
            if (N2 == size(M,2)) then
                temp => M
                nullify(M)
                allocate(M(N1,2*N2))
                M = temp
            end if
            N2 = N2 + 1
        end do

        N2 = N2-1

        if (N2 /= size(M,2)) then
            temp => M
            nullify(M)
            allocate(M(N1,N2))
            M = temp
        end if
    end subroutine readmatrix_dp

    ! subroutine: GetTextMatrixSize
    ! Reads the size of a text matrix file.
    !
    ! This has been superceded by <readmatrix_dp>
    !
    ! Parameters:
    !   fn - the filename to read
    !   N1 - returns the size of the first dimension
    !   N2 - returns the size of the second dimension
    subroutine GetTextMatrixSize(fn, N1, N2)
        character(len=*), intent(in   ) :: fn
        integer,          intent(inout) :: N1
        integer,          intent(inout) :: N2

        integer :: err, i
        real(dp), allocatable :: x(:)
        character(len=BUFSIZE) :: A
        character(len=50) :: frmt

        integer :: U

        U = new_unit()

        open(U, file=fn, action="READ", status="OLD", recl=len(A))

        frmt = "(A"//adjustl(int2str(len(A)))//")"
        read(U,frmt) A

        N1 = 2
        do i = 1, len(A)/2
            allocate(x(N1))
            read(A, *, iostat=err) x
            if (err /= 0) exit
            N1 = N1 + 1
            deallocate(x)
        end do

        N1 = N1-1
        N2 = 1

        do
            read(U,frmt, advance='yes', iostat=err) A
            if (err /= 0) exit
            N2 = N2 + 1
        end do

        close(U)
        if (allocated(x)) deallocate(x)
    end subroutine

    ! subroutine: Read2DMatrix
    ! Reads a text matrix file.
    !
    ! This has been superceded by <readmatrix_dp>
    !
    ! Parameters:
    !   X - returns the matrix
    !   fn - the filename to read
    subroutine Read2DMatrix(X, fn)
        real(dp),         intent(  out) :: X(:,:)
        character(len=*), intent(in   ) :: fn

        integer :: i,j

        integer :: U

        U = new_unit()

        open(U, file=fn, action="READ", status="OLD", recl=BUFSIZE)
        read(U,*) ((X(i,j), i = 1,size(X,1)), j = 1,size(X,2))
        close(U)
    end subroutine Read2DMatrix

    ! subroutine: Write2DArray
    ! Writes a matrix to a text file.
    !
    ! Parameters:
    !   A - the matrix
    !   fn - the filename to write
    subroutine Write2DArray(A, fn)
        implicit none
        real(dp),         intent(in) :: A(:,:)
        character(len=*), intent(in) :: fn

        character(len=50) :: frmt
        integer :: i, j

        integer :: U

        U = new_unit()

        frmt = "(SP," // int2str(size(A,1)) // "(ES23.15E3,1X)" // ')'

        if (fn == stdout) then
            do  j = 1, size(A,2)
                write(*, frmt) (A(i,j), i = 1, size(A,1))
            end do
        else
            open(U, file=fn, form="FORMATTED", recl=26*size(A,1))
            do  j = 1, size(A,2)
                write(U, frmt) (A(i,j), i = 1, size(A,1))
            end do
            close(U)
        end if
    end subroutine

    ! subroutine: Write2DIgorArray
    ! Attempts to write a matix with dimensions as understood bu Igor.
    !
    ! Parameters:
    !   A - the matix
    !   X - the x coordinates
    !   Y - the y coordinates
    !   fn - the filename to write
    subroutine Write2DIgorArray(A, X, Y, fn)
        implicit none
        real(dp),         intent(in) :: A(:,:)
        real(dp),         intent(in) :: X(:)
        real(dp),         intent(in) :: Y(:)
        character(len=*), intent(in) :: fn

        character(len=50) :: frmt1, frmt2
        character(len=*), parameter :: frmt_sub = "(ES23.15E3,1X)"
        integer :: i, j

        integer :: U

        U = new_unit()

        frmt1 = "(SP,"    // int2str(size(A,1)+1) // frmt_sub // ")"
        frmt2 = "(SP,24X" // int2str(size(A,1)  ) // frmt_sub // ")"

        if (fn == stdout) then
            write(*, frmt2) X
            do  j = 1, size(A,2)
                write(*, frmt1) Y(j), (A(i,j), i = 1, size(A,1))
            end do
        else
            open(U,file=fn, form="FORMATTED", recl=26*size(A,1))
            write(1,frmt2) X
            do  j = 1, size(A,2)
                write(U,frmt1) Y(j), (A(i,j), i = 1, size(A,1))
            end do
            close(U)
        end if
    end subroutine

    ! function: GetStdoutUnit
    ! A wrapper function to provide special handling of <constants::stdout>.
    !
    ! Unformatted writing to standard output is problematic.  Several
    ! compiler/OS combinations have difficulty with it.
    !
    ! This function provides a specially opened file for the Intel ifort
    ! compiler for BINARY formatted files.  With this, ifort can pipe large
    ! files from one process to another.
    !
    ! Parameters:
    !   binmode - are we performing an unformatted write?
    !
    ! Returns:
    !   A unit for standard out.
    !
    ! TODO:
    !   Determine if this works on non-Linux systems.
    integer function GetStdoutUnit(binmode) result(U)
        logical, intent(in) :: binmode

#ifdef __INTEL_COMPILER
        logical :: opened
        integer :: oldu

        U = DEFAULT_OUTPUT_UNIT

        if (binmode) then
            inquire(file=stdout, opened=opened, number=oldu)
            if (opened) close(oldu)
            U = new_unit()
            open(U, file=stdout, action="WRITE", form="BINARY")
        end if
#else
        U = DEFAULT_OUTPUT_UNIT
#endif
    end function GetStdoutUnit

    ! function: GetStdinUnit
    ! A wrapper function to provide special handling of <constants::stdin>.
    !
    ! Unformatted reading from standard input is problematic.  Several
    ! compiler/OS combinations have difficulty with it.
    !
    ! This function provides a specially opened file for the Intel ifort
    ! compiler for BINARY formatted files.  With this, ifort can pipe large
    ! files from one process to another.
    !
    ! Parameters:
    !   binmode - are we performing an unformatted write?
    !
    ! Returns:
    !   A unit for standard in.
    !
    ! TODO:
    !   Determine if this works on non-Linux systems.
    integer function GetStdinUnit(binmode) result(U)
        logical, intent(in) :: binmode

#ifdef __INTEL_COMPILER
        logical :: opened
        integer :: oldu

        U = DEFAULT_INPUT_UNIT

        if (binmode) then
            inquire(file=stdin, opened=opened, number=oldu)
            if (opened) close(oldu)
            U = new_unit()
            open(U, file=stdin, action="READ", form="BINARY")
        end if
#else
        U = DEFAULT_INPUT_UNIT
#endif
    end function GetStdinUnit

    ! function: new_unit
    ! Returns a unit number of a unit that exists and is not connected.
    !
    ! Searched for a unit between 11 and <MAX_UNIT_NUMBER> that is
    ! availible for use.  Eliminates the need to hardcode unit numbers.
    !
    ! Returns:
    !   A newly minted unit number.
    integer function new_unit() result(U)
        logical :: exists, opened
        integer :: ios

        do U = 11, MAX_UNIT_NUMBER
            inquire (unit = U, exist = exists, opened = opened, iostat = ios)
            if (exists .and. .not. opened .and. ios == 0) return
        end do

        U = -1
    end function new_unit


end module fileio
