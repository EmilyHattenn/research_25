! module: typetime
module typetime
    use constants
    use fileio
    use logger
    use helpers
    implicit none

#include "config.f"

    ! Struct: time
    ! The time grid structure
    type ts
        private                         ! Allows us to change the internals with impunity
        ! double:   t
        real(dp) :: t                   ! Current time (s)
        ! double:   tf
        real(dp) :: tf                  ! Final time for the simulation (s)
        ! double:   dt
        real(dp) :: dt                  ! Size of t-pixel (s)        
        ! integer:  n
        integer  :: n                   ! Current time index
    end type ts

    private initialize_field

contains
    ! subroutine: readtimeparams_sub
    ! Reads the time parameters(ts) from a unit.
    !
    ! parameters:
    ! u - the unit to read from.
    ! time - The time type, <ts>, to read into.
    subroutine readtimeparams_sub(u, time)
        implicit none
        integer,  intent(in ) :: u
        type(ts), intent(out) :: time

        integer :: err
        character(len=25) :: host

        time%t    = GetFileParam(u)
        time%tf   = GetFileParam(u)
        time%dt   = GetFileParam(u)        
        time%n    = GetFileParam(u)
    end subroutine readtimeparams_sub

    !---------------------------------------------------------------
    ! subroutine: readtimeparams
    ! Reads the time parameters(ss) from a file.
    !
    ! paramters:
    ! cmd - The filename to read.
    ! time - The time type, <ts>, to read into.
    subroutine ReadTimeParams(cmd, time)
        implicit none
        character(len=*), intent(in ) :: cmd
        type(ts),         intent(out) :: time

        integer :: U

        U = new_unit()
        open(unit=U, file=cmd, status='old', action='read')
        call readtimeparams_sub(U, time)
        close(U)
        call dumptime(time)
    end subroutine ReadTimeParams

    !---------------------------------------------------------------
    ! subroutine: writetimeparams_sub
    ! Writes time parameters(fs) to a file.  Now with decorations.
    subroutine WriteTimeParams_sub(u, time)
        implicit none
        integer,  intent(in) :: u
        type(ts), intent(in) :: time

        write(u,  pfrmtA  ) time%t ,       " : Current time of simulation. (s)"
        write(u,  pfrmtA  ) time%tf,       " : Final time of simulation. (s)"
        write(u,  pfrmtA  ) time%dt,       " : Time pixel size [dt]. (s)"
        write(u, '(I25,A)') time%n ,       " : Current time index."
    end subroutine writetimeparams_sub


    ! subroutine: writetimeparams
    ! Writes time parameters(ss) to a file.  Now with decorations.
    subroutine writetimeparams(cmd, time)
        implicit none
        character(len=*), intent(in) :: cmd
        type(ts),         intent(in) :: time

        integer :: U

        U = new_unit()
        open (unit=U, file=cmd, action="WRITE", form="FORMATTED")

        call writetimeparams_sub(U, time)

        close(U)
    end subroutine writetimeparams

    !---------------------------------------------------------------
    ! subroutine: dumptime
    subroutine dumptime(params, level)
        type(ts), intent(in) :: params
        integer, optional, intent(in) :: level

        if (present(level)) then
            if (GetLogLevel() >= level) call writetimeparams_sub(0, params)
        else
            if (GetLogLevel() >= LOGVERBOSE) call writetimeparams_sub(0, params)
        end if
    end subroutine

    !---------------------------------------------------------------
    ! function: GetT
    ! Returns the current time of the simulation.
    pure real(dp) function GetT(time)
        type(ts), intent(in) :: time
        GetT = time%t
    end function GetT   
 
    ! function: GetTf
    ! Returns the final time of the simulation.
    pure real(dp) function GetTf(time)
        type(ts), intent(in) :: time
        GetTf = time%tf
    end function GetTf      
 
    ! function: GetDf
    ! Returns the time pixel of the simulation.
    pure real(dp) function GetDt(time)
        type(ts), intent(in) :: time
        GetDt = time%dt
    end function GetDt
    
    ! function: GetN
    ! Returns the current time index
    pure integer function GetN(time)
        type(ts), intent(in) :: time
       GetN = time%n
    end function GetN

    !---------------------------------------------------------------
    ! subroutine: SetT
    subroutine SetT(time, t)
        type(ts), intent(inout) :: time
        real(dp), intent(in   ) :: t

        time%t = t
    end subroutine SetT 
    
    ! subroutine: SetTf
    subroutine SetTf(time, tf)
        type(ts), intent(inout) :: time
        real(dp), intent(in   ) :: tf

        time%tf = tf
    end subroutine SetTf    

    ! subroutine: SetDt
    subroutine SetDt(time, dt)
        type(ts), intent(inout) :: time
        real(dp), intent(in   ) :: dt

        time%dt = dt
    end subroutine SetDt  
    
    ! subroutine: SetN
    ! Sets the current time index.
    subroutine SetN(time, n)
        type(ts), intent(inout) :: time
        integer,  intent(in   ) :: n

        time%N = n
    end subroutine SetN

    !---------------------------------------------------------------

    ! function: GetNt
    ! Returns width of the time window left based on current time
    pure integer function CalcNt(time) 
        type(ts), intent(in) :: time
        CalcNt = int((time%tf - time%t) / time%dt) 
    end function CalcNt
    

    ! subroutine: UpdateT
    subroutine UpdateT(time, dt)
        type(ts), intent(inout) :: time
        real(dp), intent(in   ) :: dt

        time%t = time%t + dt
    end subroutine UpdateT 

    ! subroutine: UpdateN
    ! Updates the current time index.
    subroutine UpdateN(time, dn)
        type(ts), intent(inout) :: time
        integer,  intent(in   ) :: dn

        time%N = time%N + dn
    end subroutine UpdateN


    !---------------------------------------------------------------

    ! function: GetTArray
    ! Returns an array of the time.
    function GetTArray(time) result(t)
        type(ts), intent(in)            :: time
        real(dp), allocatable           :: t(:)
        integer :: Nt, i
        
        Nt = CalcNt(time)
        allocate(t(Nt)) 

        do i=1, Nt
            if(Nt == 1) then 
                t(i) = 0.0_dp
            else 
                t(i) = time%t + (i - 1) * time%dt
            end if
        end do

    end function GetTArray

    !---------------------------------------------------------------
    ! function: GetOmegaArray
    ! Returns the arrays for the conjugate coordinate system
    function GetOmegaArray(time) result(w)
        type(ts), intent(in)        :: time
        real(dp), allocatable       :: w(:)
        integer                     :: Nt, n
        real(dp)                    :: Tw

        Nt = CalcNt(time)
        Tw = Nt * time%dt
        allocate(w(Nt))

        do n=1, Nt
            if(n <= Nt/2) then
                w(n) = twopi * (n-1) / Tw
            else
                w(n) = twopi * (n - Nt - 1) / Tw
            end if
        end do
    end function GetOmegaArray

    !---------------------------------------------------------------

    ! function: GetdOmega
    ! Returns the differential for the conjugate coordinate system
    real(dp) function GetdOmega(time) 
        type(ts), intent(in) :: time

        GetdOmega = twopi / (CalcNt(time) * time%dt)
    end function

    !---------------------------------------------------------------

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


end module typetime
