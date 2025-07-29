! module: typespace
module typespace
    use constants
    use fileio
    use logger
    use helpers
    implicit none


#include "config.f"

    ! Struct: space
    ! The spatial grid structure
    type ss
        private                         ! Allows us to change the internals with impunity
        ! integer:  Dimensions
        integer  :: Dims                ! Dimensionality of grid (1=1D, 2=2D, 3=3D)
        ! integer:  Nx
        integer  :: Nx                  ! X-window # of pixels
        ! integer:  Ny
        integer  :: Ny                  ! Y-window # of pixels [Ny=1 for 1D]      
        ! integer:  Nz
        integer  :: Nz                  ! Z-window # of pixels [Nz=1 for 1D and 2D] 
        ! double:   dx
        real(dp) :: dx                  ! Size of x-pixel (m)
        ! double:   dy
        real(dp) :: dy                  ! Size of y-pixel (m), [dy=1.0 for 1D]
        ! double:   dz
        real(dp) :: dz                  ! Size of z-pixel (m), [dz=1.0 for 1D and 2D]   
        ! double:   epsr
        real(dp) :: epsr                ! Background dielectric constant  
    end type ss

    private initialize_field

contains
    ! subroutine: readspaceparams_sub
    ! Reads the space parameters(ss) from a unit.
    !
    ! parameters:
    ! u - the unit to read from.
    ! space - The space type, <ss>, to read into.
    subroutine readspaceparams_sub(u, space)
        implicit none
        integer,  intent(in ) :: u
        type(ss), intent(out) :: space

        integer :: err
        character(len=25) :: host

        space%dims  = GetFileParam(u)
        space%Nx    = GetFileParam(u)
        space%Ny    = GetFileParam(u)
        space%Nz    = GetFileParam(u)
        space%dx    = GetFileParam(u)
        space%dy    = GetFileParam(u)
        space%dz    = GetFileParam(u)
        space%epsr  = GetFileParam(u)
    end subroutine readspaceparams_sub


    ! subroutine: readspaceparams
    ! Reads the space parameters(ss) from a file.
    !
    ! paramters:
    ! cmd - The filename to read.
    ! space - The space type, <ss>, to read into.
    subroutine ReadSpaceParams(cmd, space)
        implicit none
        character(len=*), intent(in ) :: cmd
        type(ss),         intent(out) :: space

        integer :: U

        U = new_unit()
        open(unit=U, file=cmd, status='old', action='read')
        call readspaceparams_sub(U, space)
        close(U)
        call dumpspace(space)
    end subroutine ReadSpaceParams


    ! subroutine: writespaceparams_sub
    ! Writes space parameters(fs) to a file.  Now with decorations.
    subroutine WriteSpaceParams_sub(u, space)
        implicit none
        integer,  intent(in) :: u
        type(ss), intent(in) :: space

        write(u, '(I25,A)') space%dims, " : Number of Dimensions"
        write(u, '(I25,A)') space%Nx,   " : Number of space X points."
        write(u, '(I25,A)') space%Ny,   " : Number of space Y points."
        write(u, '(I25,A)') space%Ny,   " : Number of space Z points."
        write(u,  pfrmtA  ) space%dx,   " : The With of the X pixel. (m)"
        write(u,  pfrmtA  ) space%dy,   " : The With of the X pixel. (m)"
        write(u,  pfrmtA  ) space%dz,   " : The With of the X pixel. (m)"
        write(u,  pfrmtA  ) space%epsr, " : The relative background dielectric constant"
    end subroutine writespaceparams_sub


    ! subroutine: writespaceparams
    ! Writes space parameters(ss) to a file.  Now with decorations.
    subroutine writespaceparams(cmd, space)
        implicit none
        character(len=*), intent(in) :: cmd
        type(ss),         intent(in) :: space

        integer :: U

        U = new_unit()
        open (unit=U, file=cmd, action="WRITE", form="FORMATTED")

        call writespaceparams_sub(U, space)

        close(U)
    end subroutine writespaceparams


    ! subroutine: dumpspace
    subroutine dumpspace(params, level)
        type(ss), intent(in) :: params
        integer, optional, intent(in) :: level

        if (present(level)) then
            if (GetLogLevel() >= level) call writespaceparams_sub(0, params)
        else
            if (GetLogLevel() >= LOGVERBOSE) call writespaceparams_sub(0, params)
        end if
    end subroutine


    !---------------------------------------------------------------
    ! function: GetNx
    ! Returns the number of points in the x direction
    pure integer function GetNx(space)
        type(ss), intent(in) :: space
 
       GetNx = space%Nx
    end function GetNx


    ! function: GetNy
    ! Returns the number of points in the y direction
    pure integer function GetNy(space)
        type(ss), intent(in) :: space

        GetNy = space%Ny
    end function GetNy

    ! function: GetNz
    ! Returns the number of points in the z direction
    pure integer function GetNz(space)
        type(ss), intent(in) :: space

        GetNz = space%Nz
    end function GetNz

    !---------------------------------------------------------------
    ! function: GetDx
    ! Returns the stepsize for the various dimentions.
    pure real(dp) function GetDx(space)
        type(ss), intent(in) :: space
        if(space%Nx==1) then
            GetDx = 1.0_dp
        else
            GetDx = space%dx
        endif
    end function GetDx


    ! function: GetDy
    pure real(dp)  function GetDy(space)
        type(ss), intent(in) :: space
        if(space%Ny==1) then
            GetDy = 1.0_dp
        else
            GetDy = space%dy
        endif    
    end function GetDy


    ! function: GetDz
    pure real(dp)  function GetDz(space)
        type(ss), intent(in) :: space

        if (space%Nz == 1) then
            GetDz= 1.0_dp
        else
            GetDz = space%dz
        end if
    end function GetDz

    ! function: GetEpsr
    pure real(dp)  function GetEpsr(space)
        type(ss), intent(in) :: space

        GetEpsr = space%Epsr
    end function GetEpsr

    !---------------------------------------------------------------

    ! subroutine: SetNx
    ! Sets the number of points in the x dimention.
    subroutine SetNx(space, N)
        type(ss), intent(inout) :: space
        integer,  intent(in   ) :: N

        space%Nx = N
    end subroutine SetNx

    ! subroutine: SetNy
    ! Sets the number of points in the y dimention.
    subroutine SetNy(space, N)
        type(ss), intent(inout) :: space
        integer,  intent(in   ) :: N

        space%Ny = N
    end subroutine SetNy
    
    ! subroutine: SetNz
    ! Sets the number of points in the z dimention.
    subroutine SetNz(space, N)
        type(ss), intent(inout) :: space
        integer,  intent(in   ) :: N

        space%Nz = N
    end subroutine SetNz

    !---------------------------------------------------------------

    ! subroutine: SetDx
    subroutine SetDx(space, dl)
        type(ss), intent(inout) :: space
        real(dp), intent(in   ) :: dl

        space%dx = dl
    end subroutine SetDx

    ! subroutine: SetDy
    subroutine SetDy(space, dl)
        type(ss), intent(inout) :: space
        real(dp), intent(in   ) :: dl

        space%dy = dl
    end subroutine SetDy
    
    ! subroutine: SetDz
    subroutine SetDz(space, dl)
        type(ss), intent(inout) :: space
        real(dp), intent(in   ) :: dl

        space%dz = dl
    end subroutine SetDz

    !---------------------------------------------------------------

    ! function: GetXWidth
    ! Returns width of the window.
    pure real(dp) function GetXWidth(space)
        type(ss), intent(in) :: space

        GetXWidth = space%dx * (space%Nx-1)
    end function GetXWidth


    ! function: GetYWidth
    ! Returns width of the window.
    pure real(dp) function GetYWidth(space)
        type(ss), intent(in) :: space

        GetYWidth = space%dy * (space%Ny-1)
    end function GetYWidth

    ! function: GetZWidth
    ! Returns width of the window.
    pure real(dp) function GetZWidth(space)
        type(ss), intent(in) :: space

        GetZWidth = space%dz * (space%Nz-1)
    end function GetZWidth

    !---------------------------------------------------------------

    ! function: GetXArray
    ! Returns an array of the x positions.
    function GetXArray(space)
        type(ss), intent(in) :: space

        real(dp) :: GetXArray(space%Nx)

        if(space%Nx == 1) GetXArray = 0.0_dp
        if(space%Nx /= 1) GetXArray = GetSpaceArray(space%Nx, GetXWidth(space))

    end function GetXArray

    ! function: GetYArray
    ! Returns an array of the y positions.
    function GetYArray(space)
        type(ss), intent(in) :: space

        real(dp) :: GetYArray(space%Ny)

        if(space%Ny == 1) GetYArray = 0.0_dp
        if(space%Ny /= 1) GetYArray = GetSpaceArray(space%Ny, GetYWidth(space))
    end function GetYArray

    ! function: GetZArray
    ! Returns an array of the z positions.
    function GetZArray(space)
        type(ss), intent(in) :: space

        real(dp) :: GetZArray(space%Nz)

        if(space%Nz == 1) GetZArray = 0.0_dp
        if(space%Nz /= 1) GetZArray = GetSpaceArray(space%Nz, GetZWidth(space))
    end function GetZArray

    !---------------------------------------------------------------
    ! function: GetKxArray
    ! Returns the arrays for the conjugate coordinate system
    function GetKxArray(space)
        type(ss), intent(in) :: space
        real(dp) :: GetKxArray(space%Nx)

        if(space%Nx == 1) then
            GetKxArray = 0.0_dp
        else
            GetKxArray = GetKArray(GetNx(space), GetXWidth(space))
        end if
    end function GetKxArray

    ! function: GetKyArray
    ! Returns the arrays for the conjugate coordinate system
    function GetKyArray(space)
        type(ss), intent(in) :: space
        real(dp) :: GetKyArray(space%Ny)

        if(space%Ny == 1) then
            GetKyArray = 0.0_dp
        else
            GetKyArray = GetKArray(GetNy(space), GetYWidth(space))
        end if
    end function GetKyArray

    ! function: GetKzArray
    ! Returns the arrays for the conjugate coordinate system
    function GetKzArray(space)
        type(ss), intent(in) :: space
        real(dp) :: GetKzArray(space%Nz)

        if(space%Nz == 1) then
            GetKzArray = 0.0_dp
        else
            GetKzArray = GetKArray(GetNz(space), GetZWidth(space))
        end if
    end function GetKzArray

    !---------------------------------------------------------------

    ! function: GetdQx
    ! Returns the differential for the conjugate coordinate system
    real(dp) function GetDQx(space) result(dq)
        type(ss),    intent(in) :: space

        dq = 0.0_dp
        dq = twopi / GetXWidth(space)
    end function

    ! function: GetdQy
    ! Returns the differential for the conjugate coordinate system
    real(dp) function GetDQy(space) result(dq)
        type(ss),    intent(in) :: space

        dq = 0.0_dp
        dq = twopi / GetYWidth(space)
    end function

    ! function: GetdQz
    ! Returns the differential for the conjugate coordinate system
    real(dp) function GetDQz(space) result(dq)
        type(ss),    intent(in) :: space

        dq = 0.0_dp
        dq = twopi / GetZWidth(space)
    end function
    
    !---------------------------------------------------------------

    ! function: GetDVol
    ! Returns the volume element
    real(dp) function GetDVol(space) result(dVol)
        type(ss),    intent(in) :: space

        dVol = 0.0_dp
        dVol = GetDx(space) * GetDy(space) * GetDz(space)
    end function

    ! function: GetDQVol
    ! Returns the volume element
    real(dp) function GetDQVol(space) result(dQVol)
        type(ss),    intent(in) :: space

        dQVol = 0.0_dp
        dQVol = GetDQx(space) * GetDQy(space) * GetDQz(space)
    end function

    !---------------------------------------------------------------


    ! subroutine: writespace
    !
    ! Writes the space structure <ss> and the space to a single file.
    !
    ! Paramters:
    !   e - the complex field.
    !   space - the space structure <ss>
    !   binmode - binary mode?
    !   fnout - the filename to write.
    subroutine writefield(fnout, e, space, binmode, single, fnspace)
        implicit none
        character(len=*), intent(in) :: fnout
        complex(dp),      intent(in) :: e(:,:,:)
        type(ss),         intent(in) :: space
        logical,          intent(in) :: binmode, single
        character(len=*), intent(in), optional :: fnspace

        integer :: U

        U = 1

        if (fnout == stdout .or. fnout == "-") then
            U = GetStdoutUnit(binmode)
        else
            if (binmode) then
                open(unit=U, file=fnout, action="WRITE", form="UNFORMATTED")
            else
                open(unit=U, file=fnout, action="WRITE", form="FORMATTED")
            end if
        end if

        if (single) then
            if (binmode) then
                call unformatted_write_space(U, space)
            else
                call writespaceparams_sub(U, space)
            end if
        else
            if (present(fnspace)) call writespaceparams(fnspace, space)
        end if

        call writefield_to_unit(U, e, binmode)

        if (fnout == stdout .or. fnout == "-") then
            if (U /= 6) close(U)
        else
            close(U)
        end if
    end subroutine


    ! subroutine: readspace_only
    !
    ! Reads the space structure, <ss>, from a field file.
    !
    ! Paramters:
    !   fnin - the filename to read.
    !   space - the field space structure <ss>
    !   binmode - binary mode?
    !   single - Single file?
    !   fnspace - Optional space parameters filename.
    subroutine readspace_only(fnin, space, binmode, single, fnspace)
        implicit none
        character(len=*),         intent(in   ) :: fnin
        type(ss),                 intent(inout) :: space
        logical,                  intent(in   ) :: binmode, single
        character(len=*),         intent(in   ), optional :: fnspace

        integer :: U

        U = new_unit();

        if (fnin == stdin .or. fnin == "-") then
            U = GetStdinUnit(binmode)
        else
            if (binmode) then
                open(unit=U, file=fnin, action="READ", form="UNFORMATTED")
            else
                open(unit=U, file=fnin, action="READ", form="FORMATTED")
            end if
        end if

        if (single) then
            if (binmode) then
                call unformatted_read_space(U, space)
            else
                call readspaceparams_sub(U, space)
            end if
        else
            if (present(fnspace)) call readspaceparams(fnspace, space)
        end if

        if (fnin == stdin .or. fnin == "-") then
            if (U /= 5) close(U)
        else
            close(U)
        end if
    end subroutine readspace_only



    ! subroutine: readfield
    !
    ! Reads the space structure, <ss>, and the field.
    !
    ! Paramters:
    !   fnin - the filename to write.
    !   e - the complex field.
    !   space - the field space structure <ss>
    !   binmode - binary mode?
    !   single - Single file?
    !   fnspace - Optional space parameters filename.
#ifdef HAVE_ALLOC_DUMMY
    subroutine readfield(fnin, e, space, binmode, single, fnspace)
        implicit none
        character(len=*),         intent(in   ) :: fnin
        complex(dp), allocatable, intent(inout) :: e(:,:,:)
        type(ss),                 intent(inout) :: space
        logical,                  intent(in   ) :: binmode, single
        character(len=*),         intent(in   ), optional :: fnspace

        integer :: U

        U = new_unit();

        if (fnin == stdin .or. fnin == "-") then
            U = GetStdinUnit(binmode)
        else
            if (binmode) then
                open(unit=U, file=fnin, action="READ", form="UNFORMATTED")
            else
                open(unit=U, file=fnin, action="READ", form="FORMATTED")
            end if
        end if

        if (single) then
            if (binmode) then
                call unformatted_read_space(U, space)
            else
                call readspaceparams_sub(U, space)
            end if
        else
            if (present(fnspace)) call readspaceparams(fnspace, space)
        end if

        if (allocated(e)) then
            if (size(e,1) /= GetNx(space) .or. size(e,2) /= GetNy(space) .or. size(e,3) /= GetNz(space)) then
                deallocate(e)
                allocate(e(GetNx(space), GetNy(space), GetNz(space)))
                call initialize_field(e)
            end if
        else
	        allocate(e(GetNx(space), GetNy(space), GetNz(space)))
            call initialize_field(e)
        end if

        call readfield_from_unit(U, e, binmode)

        if (fnin == stdin .or. fnin == "-") then
            if (U /= 5) close(U)
        else
            close(U)
        end if
    end subroutine readfield
#else
    subroutine readfield(fnin, e, space, binmode, single, fnspace)
        implicit none
        character(len=*),         intent(in   ) :: fnin
        complex(dp), pointer :: e(:,:,:)
        type(ss),                 intent(inout) :: space
        logical,                  intent(in   ) :: binmode, single
        character(len=*),         intent(in   ), optional :: fnspace

        integer :: U

        U = new_unit();

        if (fnin == stdin .or. fnin == "-") then
            U = GetStdinUnit(binmode)
        else
            if (binmode) then
                open(unit=U, file=fnin, action="READ", form="UNFORMATTED")
            else
                open(unit=U, file=fnin, action="READ", form="FORMATTED")
            end if
        end if

        if (single) then
            if (binmode) then
                call unformatted_read_space(U, space)
            else
                call readspaceparams_sub(U, space)
            end if
        else
            if (present(fnspace)) call readspaceparams(fnspace, space)
        end if

        if (associated(e)) then
            if (size(e,1) /= GetNx(space) .or. size(e,2) /= GetNy(space) .or. size(e,3) /= GetNz(space)) then
                deallocate(e)
                allocate(e(GetNx(space), GetNy(space), GetNz(space)))
            end if
        else
            allocate(e(GetNx(space), GetNy(space), GetNz(space)))
        endif
        
        call initialize_field(e)

        if (fnin == stdin .or. fnin == "-") then
            if (U /= 5) close(U)
        else
            close(U)
        end if
    end subroutine readfield
#endif


    ! subroutine: readfield_from_unit
    subroutine readfield_from_unit(u, e, binmode)
        integer,     intent(in   ) :: u
        complex(dp), intent(  out) :: e(:,:,:)
        logical,     intent(in   ) :: binmode

        integer :: i, j, k
        real(dp) :: re, im

        if (binmode) then
            call unformatted_read_e(u, e)
        else
            do k = 1, size(e,3)
                do j = 1, size(e,2)
                    do i = 1, size(e,1)
                        read(u, *) re, im
                        e(i,j,k) = cmplx(re, im)
                    end do
                end do
            end do
        end if
    end subroutine readfield_from_unit


    ! subroutine: writefield_to_unit
    subroutine writefield_to_unit(u, e, binmode)
        integer,     intent(in   ) :: u
        complex(dp), intent(in   ) :: e(:,:,:)
        logical,     intent(in   ) :: binmode

        integer :: i, j, k

        if (binmode) then
            call unformatted_write_e(u, e)
        else
            do k = 1, size(e,3)
                do j = 1, size(e,2)
                    do i = 1, size(e,1)
                        write(u, efrmt2x) real(e(i,j,k)), aimag(e(i,j,k))
                    end do
                end do
            end do
        end if
    end subroutine writefield_to_unit


    ! subroutine: initialize_field
    subroutine initialize_field(e)
        complex(dp), intent(inout) :: e(:,:,:)
#ifdef _OPENMP
        integer :: i,j,k

        !$omp parallel do
        do k = 1, size(e,3)
            do j = 1, size(e,2)
                do i = 1, size(e,1)
                    e(i,j,k) = 0.0_dp
                end do
            end do
        end do
        !$omp end parallel do
#else
        e = cmplx(0.0_dp, 0.0_dp)
#endif
    end subroutine initialize_field


    ! subroutine: unformatted_write_space
    subroutine unformatted_write_space(U, space)
        integer,     intent(in) :: U
        type(ss),    intent(in) :: space

#ifdef __INTEL_COMPILER
        character(len=32) :: frm

        inquire(U, form=frm)
        select case(frm)
            case ("BINARY")
                write(U) int(sizeof(space)), space, int(sizeof(space))
            case default
                write(U) space
        end select
#else
        write(U) space
#endif
    end subroutine unformatted_write_space


    ! subroutine: unformatted_write_e
    subroutine unformatted_write_e(U, e)
        integer,     intent(in) :: U
        complex(dp), intent(in) :: e(:,:,:)

#ifdef __INTEL_COMPILER
        character(len=32) :: frm

        inquire(U, form=frm)
        select case(frm)
            case ("BINARY")
                write(U) int(sizeof(e(1,1,1))*size(e)), e, int(sizeof(e(1,1,1))*size(e))
            case default
                write(U) e
        end select
#else
        write(U) e
#endif
    end subroutine unformatted_write_e


    ! subroutine: unformatted_read_space
    subroutine unformatted_read_space(U, space)
        integer,     intent(in ) :: U
        type(ss),    intent(out) :: space

#ifdef __INTEL_COMPILER
        character(len=32)   :: frm
        integer             :: N1, N2

        inquire(U, form=frm)
        select case(frm)
            case ("BINARY")
                read(U) N1, space, N2
                if (N1 /= N2) &
                    call error("Error reading binary space structure, guard bytes do not match (" &
                        // int2str(int(N1)) // ", " // int2str(int(N2)) // ").", 1, __FILE__, __LINE__)
                if (int(N1) /= int(sizeof(space))) &
                    call error("Error reading binary space structure, size (" // int2str(int(sizeof(space))) // ") " &
                        // "doesn't match guard bytes (" // int2str(int(N1)) // ").", 1, __FILE__, __LINE__)
            case default
                read(U) space
        end select
#else
        read(U) space
#endif
    end subroutine unformatted_read_space


    ! subroutine: unformatted_read_e
    subroutine unformatted_read_e(U, e)
        integer,     intent(in ) :: U
        complex(dp), intent(out) :: e(:,:,:)

#ifdef __INTEL_COMPILER
        character(len=32)   :: frm
        integer             :: N1, N2

        inquire(U, form=frm)
        select case(frm)
            case ("BINARY")
                read(U) N1, e, N2
                if (N1 /= N2) &
                    call error("Error reading binary field, guard bytes do not match (" &
                        // int2str(int(N1)) // ", " // int2str(int(N2)) // ").", 1, __FILE__, __LINE__)
                if (int(N1) /= int(sizeof(e(1,1,1))*size(e))) &
                    call error("Error reading binary field, size (" // int2str(int(sizeof(e(1,1,1))*size(e))) // ") " &
                        // "doesn't match guard bytes (" // int2str(int(N1)) // ").", 1, __FILE__, __LINE__)
            case default
                read(U) e
        end select
#else
        read(U) e
#endif
    end subroutine unformatted_read_e

end module typespace
