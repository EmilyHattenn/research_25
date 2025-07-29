program pstd3d
    use types
    use fftw
    use constants
    use helpers
    use typespace
    use typetime
    use typepulse
    use rhopj
    use usefulsubs

    implicit none

    type(ss) :: space
    type(ts) :: time
    type(ps) :: pulse
    
    
    !  Read in the space structure for the simulation
    call ReadSpaceParams("params/space.params", space)
    
    !  Read in the time structure for the simulation
    call ReadTimeParams("params/time.params", time)   
   
    !  Read in the time structure for the simulation
    call ReadPulseParams("params/pulse.params", pulse)    
    
    !  Initialize the source calculations scheme
    call InitializeSources(space, time, pulse)
    
    !  Do the actual work of propagation
    call PSTD_3D_Propagator(space, time, pulse)
    
    !  Program Complete
contains

    subroutine PSTD_3D_Propagator(space, time, pulse)
        type(ss), intent(in   ) :: space
        type(ts), intent(inout) :: time
        type(ps), intent(in   ) :: pulse
                        
        complex(dp) :: Ex(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Ey(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Ez(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Bx(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: By(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Bz(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Jx(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Jy(GetNx(space), GetNy(space), GetNz(space)) 
        complex(dp) :: Jz(GetNx(space), GetNy(space), GetNz(space)) 
                    
        integer     :: n, n1, Nt
        real(dp)    :: v, dt

        n1 = GetN(time)
        Nt = CalcNt(time)
        v  = c0 / sqrt(GetEpsr(space))
        dt = GetDt(time)
    
        call InitializeFields(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz)
        
        ! Do the time loop
        do n=n1, Nt
            ! Update the user
            print*, n1, n, Nt
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! UPDATE THE E-FIELDS
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !  Update E-field by TFSF pulse injection
            call UpdateTFSC(space, time, pulse, GetAmp(pulse), Ey)
            
            !  Update E-field by Maxwell-Ampere's Law
            call UpdateE3D(space, time, Bx, By, Bz, Jx, Jy, Jz, Ex, Ey, Ez)
            
            !  Update time by one-half time-step (n+1/2)
            call UpdateT(time, dt/2d0)
            
            !  Calculate new source terms for time n+1
            call CalcJ(space, time, Ex, Ey, Ez, Jx, Jy, Jz)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! UPDATE THE B-FIELDS
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
            !  Update B-field by TFSF pulse injection
            call UpdateTFSC(space, time, pulse, GetAmp(pulse) / v, Bz)
            
            !  Update B-field by Faraday's Law
            call UpdateB3D(space, time, Ex, Ey, Ez, Bx, By, Bz)
            
            !  Update time by one-half time-step
            call UpdateT(time, dt/2d0)
            
            !  Update current time index by one (n-> n+1)
            call UpdateN(time, n+1)
            
        end do
        
        ! WriteAllFields
    end subroutine
    

    subroutine UpdateE3D(space, time, Bx, By, Bz, Jx, Jy, Jz, Ex, Ey, Ez)
        type(ss),    intent(in   ) :: space
        type(ts),    intent(in   ) :: time
        complex(dp), intent(in   ) :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
        complex(dp), intent(in   ) :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
        complex(dp), intent(inout) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        
        integer  :: i, j, k
        real(dp) :: qx(GetNx(space)), qy(GetNy(space)), qz(GetNz(space))
        real(dp) :: v2, dt
        
        v2 = c0**2 * GetEpsr(space)
        dt = GetDt(time)
        
        !  Get Q-Space Arrays    
        qx = GetKxArray(space)
        qy = GetKxArray(space)
        qz = GetKyArray(space)

        ! Update Bx
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                Ex(:,j,k) = Ex(:,j,k) + ii * (qy(j) * Bz(:,j,k) - qz(k) * By(:,j,k)) * v2 * dt - mu0 * Jx(:,i,j) * v2 * dt
            end do
        end do
        !$omp end parallel do

        ! Update By
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                Ey(:,j,k) = Ey(:,j,k) + ii * (qz(k) * Bx(:,j,k) - qx(:) * Bz(:,j,k)) * v2 * dt - mu0 * Jy(:,i,j) * v2 * dt
            end do
        end do
        !$omp end parallel do       

        ! Update Bz
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                Ez(:,j,k) = Ez(:,j,k) + ii * (qx(:) * By(:,j,k) - qy(j) * Bx(:,j,k)) * v2 * dt - mu0 * Jz(:,i,j) * v2  * dt
            end do
        end do
        !$omp end parallel do 
        
    end subroutine

    
    subroutine UpdateB3D(space, time, Ex, Ey, Ez, Bx, By, Bz)
        type(ss),    intent(in   ) :: space
        type(ts),    intent(in   ) :: time
        complex(dp), intent(in   ) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        complex(dp), intent(inout) :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
        
        integer  :: i, j, k
        real(dp) :: qx(GetNx(space)), qy(GetNy(space)), qz(GetNz(space))
        real(dp) :: dt
        
        dt = GetDt(time)
        
        !  Get Q-Space Arrays    
        qx = GetKxArray(space)
        qy = GetKxArray(space)
        qz = GetKyArray(space)

        ! Update Bx
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                Bx(:,j,k) = Bx(:,j,k) - ii * (qy(j) * Ez(:,j,k) - qz(k) * Ey(:,j,k)) * dt
            end do
        end do
        !$omp end parallel do

        ! Update By
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                By(:,j,k) = By(:,j,k) - ii * (qz(k) * Ex(:,j,k) - qx(:) * Ez(:,j,k)) * dt
            end do
        end do
        !$omp end parallel do       

        ! Update Bz
        !$omp parallel do private(k,j)
        do k=1, GetNz(space)
            do j=1, GetNy(space)
                Bz(:,j,k) = Bz(:,j,k) - ii * (qx(:) * Ey(:,j,k) - qy(j) * Ex(:,j,k)) * dt
            end do
        end do
        !$omp end parallel do 
        
    end subroutine    
    
    subroutine UpdateTFSC(space, time, pulse, Emax, E)
        type(ss),    intent(in   ) :: space
        type(ts),    intent(in   ) :: time
        type(ps),    intent(in   ) :: pulse
        real(dp),    intent(in   ) :: Emax
        complex(dp), intent(inout) :: E(:,:,:)
        complex(dp)                :: dE(:,:,:)
        
        
        !  Do all the real space field creation here!!!
        dE = 0d0
        
        call fft(dE)
        
        E = E + dE
    end subroutine    
    
    
    subroutine InitializeFields(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz)
        complex(dp), intent(inout) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        complex(dp), intent(inout) :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
        complex(dp), intent(inout) :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)

        Ex = 0d0
        Ey = 0d0
        Ez = 0d0
        Bx = 0d0
        By = 0d0
        Bz = 0d0        
        Jx = 0d0
        Jy = 0d0
        Jz = 0d0                
    end subroutine

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
    write(*,"(A)")"Report bugs to Jeremy Gulley <jgulley@furman.edu>"
end subroutine help

! subroutine: version
! Ouputs the version information.
! Uses CVS RCSfile and Revision tags to generate file name and version.
subroutine version
    use pscommandline
end subroutine
