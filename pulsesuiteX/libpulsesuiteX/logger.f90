! module: logger
! Provides standard output routines.
! Also controls the amount of output.
module logger
    implicit none

    ! const: LOGERROR
    ! Error loglevel.
    ! Problems that have no fix.
    ! i.e. If the file doesn't exist, there is _no way_ to open it.
    integer, parameter :: LOGERROR   = 0
    ! const: LOGWARN
    ! Loglevel for warnings to the user.
    ! Recoverable problems that may change the
    ! expected output of a program go here.
    integer, parameter :: LOGWARN    = 1
    ! const: LOGSTD
    ! The default loglevel.
    ! Common user information is output at this level.
    integer, parameter :: LOGSTD     = 2
    ! const: LOGVERBOSE
    ! Verbose loglevel.
    ! Use this for information that may be of interest for the user.
    integer, parameter :: LOGVERBOSE = 3
    ! const: LOGDEBUG
    ! Debug loglevel.
    ! Information the user should never need, but that a developer might.
    ! Shouldn't be used in inner loops.
    integer, parameter :: LOGDEBUG   = 4
    ! const: LOGDEBUG2
    ! Debug2 loglevel.
    ! More debugging information.
    ! Can be used in inner loops.
    integer, parameter :: LOGDEBUG2  = 5
    ! const: LOGDEBUG3
    ! Debug3 loglevel.
    ! A whole lot of details.
    ! Inner loops should debugging should use this level.
    integer, parameter :: LOGDEBUG3  = 6

    ! const: LogNames
    ! The loglevels in string format , if needed.
    character(len=7), parameter :: LogNames(-1:7) = (/ "NONE   ", &
                                                       "ERROR  ", &
                                                       "WARNING", &
                                                       "STD    ", &
                                                       "VERBOSE", &
                                                       "DEBUG  ", &
                                                       "DEBUG2 ", &
                                                       "DEBUG3 ", &
                                                       "ALL    " /)

    ! variable: loglevel
    ! PRIVATE: Holds the current log level.
    integer, private :: loglevel = LOGSTD

    ! No one should call these functions but us.
    private :: KillItDefault, LogSub
contains

    ! subroutine: error
    ! Writes an error string to stderr and exits with value n.
    !
    ! This subroutine displays an error message to the user and exits
    ! the process with a nonzero return value.  If the file and line parameter
    ! are _both_ specified, they are used in generating the error message.
    !
    ! This subroutine works at the <LOGERROR> log level.  At lower levels no message
    ! is printed, but the program will always exit.
    !
    ! Parameters:
    !   str - the message to display
    !   n - (optional) the exit value (default: 1)
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine error(str, n, file, line)
        character(len=*), intent(in) :: str
        integer,          intent(in), optional :: n
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGERROR)

        if (present(n)) then
            call KillIt(n)
        else
            call KillIt(1)
        end if
    end subroutine

    ! subroutine: assert
    ! If the test id false, writes an error string to stderr and exits with value of 16.
    !
    ! This subroutine works at the <LOGERROR> log level.  At lower levels no message
    ! is printed, but the program will always exit.
    !
    ! Parameters:
    !   test - logical function that should be true
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine assertt(test, str, file, line)
        logical,          intent(in) :: test
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (test) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGERROR)
        call KillIt(16)
    end subroutine

    ! subroutine: warning
    ! Issue a warning to stderr and continue the program.
    !
    ! Use this subroutine to warn the user of recoverable problems.
    ! i.e. If the program had to use default values for a parameter, the
    ! user should be warned.
    !
    ! This subroutine works at the <LOGWARN> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine warning(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGWARN) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGWARN)
    end subroutine

    ! subroutine: stdlev
    ! This is the default log level.
    ! Use this for standard monitoring statements.
    !
    ! This subroutine works at the <LOGSTD> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine stdlev(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGSTD) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        if (loglevel >= LOGSTD) write(0,"(A)") trim(msg)
    end subroutine

    ! subroutine: verbose
    ! Use this for messages a user might need, but not normally.
    !
    ! This subroutine works at the <LOGVERBOSE> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine verbose(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGVERBOSE) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGVERBOSE)
    end subroutine

    ! subroutine: debug
    ! This is for debugging information.
    !
    ! Users not should need to see this information unless there's a problem.
    !
    ! This subroutine works at the <LOGDEBUG> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine debug(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGDEBUG) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGDEBUG)
    end subroutine

    ! subroutine: debug2
    ! Use this for information that really cluters the screen, inner loops and such.
    !
    ! This subroutine works at the <LOGDEBUG2> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine debug2(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGDEBUG2) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGDEBUG2)
    end subroutine

    ! subroutine: debug3
    ! Use this for information that really, really cluters the screen, inner, inner loops and such.
    !
    ! This subroutine works at the <LOGDEBUG3> log level.  At lower levels no message
    ! is printed.
    !
    ! Parameters:
    !   str - the message to display
    !   file - (optional) file name for error location (i.e. the __FILE__ macro)
    !   line - (optional) line number for error location (i.e. the __LINE__ macro)
    subroutine debug3(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer,          intent(in), optional :: line

        character(len=255) :: msg

        if (loglevel < LOGDEBUG3) return

        msg = str
        if (present(line) .and. present(file)) msg = MakeMsg(str, file, line)

        call LogSub(msg, LOGDEBUG3)
    end subroutine

    ! function: GetLogLevel
    ! Returns the current log level.
    integer function GetLogLevel() result(level)
        level = loglevel
    end function

    ! subroutine: IncrLogLevel
    ! Increases the log level.
    subroutine IncrLogLevel()
        loglevel = loglevel + 1
    end subroutine

    ! subroutine: DecrLogLevel
    ! Decreases to log level.
    subroutine DecrLogLevel()
        loglevel = loglevel - 1
    end subroutine

    ! subroutine: SetLogLevel
    ! Fixes the log level to a particular value.
    !
    ! Parameters:
    !   N - the new loglevel
    subroutine SetLogLevel(N)
        integer, intent(in) :: N
        loglevel = N
    end subroutine

    ! function: GetLogLevelStr
    ! Returns a human readable string of the current log level or the specified level.
    !
    ! Parameters:
    !   level - (optional) the loglevel to get the string representation of
    function GetLogLevelStr(level) result(str)
        integer, intent(in), optional :: level
        character(len=len(LogNames(1))) str

        if (present(level)) then
            str = LogNames(min(max(   level, lbound(LogNames,1)), ubound(LogNames,1)))
        else
            str = LogNames(min(max(loglevel, lbound(LogNames,1)), ubound(LogNames,1)))
        end if
    end function

    ! subroutine: LogSub
    ! Performs the actual output of the log message.
    subroutine LogSub(str, level)
        character(len=*), intent(in) :: str
        integer,          intent(in) :: level
        write(0,"(3A)") GetLogLevelStr(level), ": ", trim(str)
    end subroutine

    ! function: MakeMsg
    ! Builds the message when file and line are specified.
    !
    ! Parameters:
    !   str - the message to display
    !   file - file name for error location (i.e. the __FILE__ macro)
    !   line - line number for error location (i.e. the __LINE__ macro)
    function MakeMsg(str, file, line) result(msg)
        character(len=*), intent(in) :: str, file
        integer,          intent(in) :: line

        character(len=len_trim(file)+4+len_trim(str)+6) :: msg
        character(len=8) :: ln

        write(ln,"(I8)") line

        msg = trim(file) // ":" // trim(adjustl(ln)) // " : " // trim(str)
    end function MakeMsg

    ! subroutine: KillIt
    ! Our exit function.
    ! This wraps the various exit functions.
    subroutine KillIt(n)
#if defined(__INTEL_COMPILER)
        use IFPOSIX                         ! Needed for exit subroutine.
        use IFCORE                          ! tracebackqq
#endif
#if defined(NAG)
        use F90_unix                        ! Needed for exit subroutine.
#endif

        integer, intent(in) :: n            ! Always included.

#if defined(__INTEL_COMPILER)
        if (GetLogLevel() >= LOGDEBUG2) call tracebackqq(user_exit_code=abs(n))
#endif
#if defined(XLF)
        call exit_(n)
#endif
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(NAG)
        call exit(n)
#endif
        call KillItDefault(n)               ! Default exit routine.
    end subroutine KillIt

    ! subroutine: KillItDefault
    ! This is supposed to cause a run time error.
    !
    subroutine KillItDefault(n)
        integer, intent(in) :: n
        real x
        x = 0.0
        write(0,'(A,I3,A)') 'Trying to create runtime error for braindead fortran: ',n, x, 1.0/x
        stop 1
    end subroutine KillItDefault

end module logger
