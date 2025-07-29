! module: typepulse
module typepulse
    use constants
    use fileio
    use logger
    use helpers
    implicit none


#include "config.f"

    ! Struct: pulse
    ! The pulse parameters
    type ps
        private                         ! Allows us to change the internals with impunity
        ! integer:  lambda
        real(dp) :: lambda              ! Vacuum wavelength of the field (m)
        ! integer:  Amp
        real(dp) :: Amp                 ! Peak pulse amplitude (V/m) 
        ! double:   Tw
        real(dp) :: Tw                  ! Pulsewidth (s)
        ! double:   Tp
        real(dp) :: Tp                  ! Time the pulse peak passes through the origin (s)
        ! double:   chirp
        real(dp) :: chirp               ! Pulse chirp (rad/s^2)  
        ! integer:  pol
        integer  :: pol                 ! Polarization Index (optional)
    end type ps

contains

    !--------------------------------------------------------------- File I/O ---------------------------------------------------------------!

    ! subroutine: readpulseparams_sub
    ! Reads the pulse parameters(ps) from a unit.
    !
    ! parameters:
    ! u - the unit to read from.
    ! pulse - The pulse type, <ps>, to read into.
    subroutine readpulseparams_sub(u, pulse)
        implicit none
        integer,  intent(in ) :: u
        type(ps), intent(out) :: pulse

        integer :: err
        character(len=25) :: host

        pulse%lambda = GetFileParam(u)
        pulse%amp    = GetFileParam(u)
        pulse%tw     = GetFileParam(u)
        pulse%tp     = GetFileParam(u)
        pulse%chirp  = GetFileParam(u)
    end subroutine readpulseparams_sub

    
    ! subroutine: ReadPulseParams
    ! Reads the pulse parameters(ss) from a file.
    !
    ! paramters:
    ! cmd - The filename to read.
    ! pulse - The pulse type, <ps>, to read into.
    subroutine ReadPulseParams(cmd, pulse)
        implicit none
        character(len=*), intent(in ) :: cmd
        type(ps),         intent(out) :: pulse

        integer :: U

        U = new_unit()
        open(unit=U, file=cmd, status='old', action='read')
        call readpulseparams_sub(U, pulse)
        close(U)
        call dumppulse(pulse)
    end subroutine ReadPulseParams


    !---------------------------------------------------------------
    ! subroutine: writepulseparams_sub
    ! Writes pulse parameters(ps) to a file.  Now with decorations.
    subroutine writepulseparams_sub(u, pulse)
        implicit none
        integer,  intent(in) :: u
        type(ps), intent(in) :: pulse

        write(u,  pfrmtA  ) pulse%lambda,   " : The pulse wavelength. (m)"
        write(u,  pfrmtA  ) pulse%amp,      " : The pulse amplitude. (V/m)"
        write(u,  pfrmtA  ) pulse%tw,       " : The pulsewidth. (s)"
        write(u,  pfrmtA  ) pulse%tp,       " : The time the pulse crosses the origin. (s)"
        write(u,  pfrmtA  ) pulse%chirp,    " : The pulse chirp constant. (rad/s^2)"
        
    end subroutine writepulseparams_sub


    ! subroutine: WritePulseParams
    ! Writes pulse parameters(ps) to a file.  Now with decorations.
    subroutine WritePulseParams(cmd, pulse)
        implicit none
        character(len=*), intent(in) :: cmd
        type(ps),         intent(in) :: pulse

        integer :: U

        U = new_unit()
        open (unit=U, file=cmd, action="WRITE", form="FORMATTED")

        call writepulseparams_sub(U, pulse)

        close(U)
    end subroutine WritePulseParams


    ! subroutine: dumppulse
    subroutine dumppulse(params, level)
        type(ps), intent(in) :: params
        integer, optional, intent(in) :: level

        if (present(level)) then
            if (GetLogLevel() >= level) call writepulseparams_sub(0, params)
        else
            if (GetLogLevel() >= LOGVERBOSE) call writepulseparams_sub(0, params)
        end if
    end subroutine


    !--------------------------------------------------------------- Get Pulse Properties ---------------------------------------------------------------!
    ! function: GetLambda
    ! Returns the Lambda Value from Pulse object
    pure real(dp) function GetLambda(pulse)
        type(ps), intent(in) :: pulse
 
       GetLambda = pulse%lambda
    end function GetLambda


    ! function: GetAmp
    ! Returns the Amp Value from Pulse object
    pure real(dp) function GetAmp(pulse)
        type(ps), intent(in) :: pulse

        GetAmp = pulse%Amp
    end function GetAmp

    ! function: GetTw
    ! Returns the Tw Value from Pulse object
    pure real(dp) function GetTw(pulse)
        type(ps), intent(in) :: pulse

        GetTw = pulse%Tw
    end function GetTw

    ! function: GetTp
    ! Returns the Tp Value from Pulse object
    pure real(dp) function GetTp(pulse)
        type(ps), intent(in) :: pulse

        GetTp = pulse%Tp
    end function GetTp

    ! function: GetChirp
    ! Returns the chrip Value from Pulse object
    pure real(dp) function GetChirp(pulse)
        type(ps), intent(in) :: pulse

        GetChirp = pulse%chirp
    end function GetChirp

    ! function: GetPol
    ! Returns the pol Value from Pulse object
    pure integer function GetPol(pulse)
        type(ps), intent(in) :: pulse

        GetPol = pulse%pol
    end function GetPol

    !--------------------------------------------------------------- Set Pulse Attributes ---------------------------------------------------------------!

    ! subroutine: SetLambda
    ! Sets vale of lambda
    subroutine SetLambda(pulse, lambda)
        type(ps), intent(inout) :: pulse
        integer,  intent(in   ) :: lambda

        pulse%lambda = lambda
    end subroutine SetLambda

    ! subroutine: SetAmp
    ! Sets Tp value
    subroutine SetAmp(pulse, Amp)
        type(ps), intent(inout) :: pulse
        real(dp), intent(in   ) :: Amp

        pulse%Amp = Amp
    end subroutine SetAmp
    
    ! subroutine: SetTw
    ! Sets PulseWidth Value (Tw)
    subroutine SetTw(pulse, Tw)
        type(ps), intent(inout) :: pulse
        real(dp),  intent(in  ) :: Tw

        pulse%Tw = Tw
    end subroutine SetTw

    ! subroutine: SetTp
    ! Sets 
    subroutine SetTp(pulse, Tp)
        type(ps), intent(inout) :: pulse
        real(dp),  intent(in  ) :: Tp

        pulse%Tp = Tp
    end subroutine SetTp

    ! subroutine: SetChirp
    ! Sets chirp value
    subroutine SetChirp(pulse, chirp)
        type(ps),  intent(inout) :: pulse
        real(dp),  intent(in   ) :: chirp

        pulse%chirp = chirp
    end subroutine SetChirp

    ! subroutine: SetPol
    ! Sets pol value
    subroutine SetPol(pulse, pol)
        type(ps), intent(inout) :: pulse
        integer,  intent(in   ) :: pol

        pulse%pol = pol
    end subroutine SetPol

    !--------------------------------------------------------------- Temporal Properties ---------------------------------------------------------------!
    ! function: CalcK0
    ! Returns Wave Number
    pure real(dp) function CalcK0(pulse)
        type(ps), intent(in) :: pulse
        CalcK0 = twopi / pulse%lambda
    end function


    ! function: 
    ! Returns 
    pure real(dp) function CalcFreq0(pulse)
        type(ps), intent(in) :: pulse
        CalcFreq0 = c0 / pulse%lambda
    end function

    ! function: CalcOmega0
    ! Returns Carrier Frequency
    pure real(dp) function CalcOmega0(pulse)
        type(ps), intent(in) :: pulse
        CalcOmega0 = twopi * c0 / pulse%lambda
    end function

!--------------------------------------------------------------- Spectural Properties ---------------------------------------------------------------!
    ! INCLUDE????

    ! function: CalcTau
    ! Returns 
    pure real(dp) function CalcTau(pulse)
        type(ps), intent(in) :: pulse

        CalcTau = pulse%Tw / (2.0_dp * sqrt(log(2.0_dp)))
    end function

    ! function: CalcDeltaOmega
    ! Returns Fourier Width
    pure real(dp) function CalcDeltaOmega(pulse)
        type(ps), intent(in) :: pulse
    
        CalcDeltaOmega = 0.44_dp / CalcTau(pulse)
    end function

    ! function: 
    ! Returns 
    pure real(dp) function CalcTime_BandWidth(pulse)
        type(ps), intent(in) :: pulse
        CalcTime_BandWidth = CalcDeltaOmega(pulse) * CalcTau(pulse)
    end function

!--------------------------------------------------------------- Spatial Properties ---------------------------------------------------------------!

    ! function: CalcRayleigh
    ! Returns 
    pure real(dp) function CalcRayleigh(pulse)
        type(ps), intent(in) :: pulse
        real(dp) :: w0

        w0 = CalcOmega0(pulse)
        CalcRayleigh = pi * w0**2 / pulse%lambda
    end function

    ! function: CalcCurvature
    ! Returns 
    pure real(dp) function CalcCurvature(pulse, x)
        type(ps), intent(in) :: pulse
        real(dp), intent(in) :: x                      ! Axis in which pulse should propogate
        real(dp) :: xR

        xR = CalcRayleigh(pulse)
        if (x == 0.0_dp) then
            CalcCurvature = huge(1.0_dp)
        else
            CalcCurvature = x * (1.0_dp + (xR/x)**2)
        end if
    end function

    ! function: CalcGouyPhase
    ! Returns 
    pure real(dp) function CalcGouyPhase(pulse, x)
        type(ps), intent(in) :: pulse
        real(dp), intent(in) :: x                      ! Axis in which pulse should propogate

        CalcGouyPhase = atan(x / CalcRayleigh(pulse))
    end function

    !--------------------------------------------------------------- Other Functions ---------------------------------------------------------------!

    ! function: PulseFieldXT
    ! Returns 
    pure complex(dp) function PulseFieldXT(x, t, pulse) result(E)
        implicit none
        type(ps), intent(in) :: pulse
        real(dp), intent(in) :: x, t

        ! Local variables
        real(dp) :: Tw, Tp, omega0, chirp, delay, phase
        complex(dp) :: envelope

        Tw     = pulse%Tw
        Tp     = pulse%Tp
        omega0 = CalcOmega0(pulse)
        chirp  = pulse%chirp

        ! Effective delay: time in pulse frame (retarded time)
        delay = t - Tp - x / c0

        ! Gaussian envelope
        envelope = pulse%Amp * exp( - (delay**2) / (2.0_dp * Tw**2) )

        ! Total phase (carrier + chirp)
        phase = omega0 * delay + chirp * delay**2

        E = envelope * cmplx(cos(phase), sin(phase), kind=dp)
    end function PulseFieldXT


end module typepulse
