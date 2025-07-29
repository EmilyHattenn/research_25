! module: fft
!
! A collection of routines to perform FFTs on various dimentionallities of
! arrays.
module fftw
    use types
    implicit none

#include "config.f"
#include <fftw3.f>

    interface fft
        module procedure fft_1D, fft_2D, fft_3D
    end interface

    interface ifft
        module procedure ifft_1D, ifft_2D, ifft_3D
    end interface

    interface fftc
        module procedure fftc_1D, fftc_2D, fftc_3D
    end interface

    interface ifftc
        module procedure ifftc_1D, ifftc_2D, ifftc_3D
    end interface

    interface nyquist
        module procedure nyquist_1D, nyquist_2D, nyquist_3D
    end interface

    interface fftw_initialize
        module procedure fftw_initialize_2D, fftw_initialize_3D
    end interface

    interface
        INTEGER FUNCTION OMP_GET_MAX_THREADS()
        end function
    end interface


    real(dp), private, allocatable :: HT(:,:)
    real(dp), private, allocatable :: a(:)
    character(len=256), private :: J0zeros_file = '/home/jgulley/Desktop/pulsesuite/J0zeros.dat'

contains

    ! function: fft_1D
    !
    ! Computes a one dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fft_1D(Z)
        complex(dp), intent(inout) :: Z(:)

        integer(i8b) :: plan
        complex(dp) :: T(size(Z))

        T = Z

        call dfftw_plan_dft_1d(plan,size(T),T,T,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute(plan)
        call dfftw_destroy_plan(plan)

        Z = T
    end subroutine fft_1D

    ! function: ifft_1D
    !
    ! Computes a one dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifft_1D(Z)
        complex(dp), intent(inout) :: Z(:)

        integer(i8b) :: plan
        complex(dp) :: T(size(Z))

        T = Z

        call dfftw_plan_dft_1d(plan,size(T),T,T,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute(plan)
        call dfftw_destroy_plan(plan)

        Z = T
        Z = Z / size(Z)
    end subroutine ifft_1D

    ! function: fftc_1D
    !
    ! Computes a centered, one dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fftc_1D(Z)
        complex(dp), intent(inout) :: Z(:)

        call nyquist_1D(Z)

        call fft_1D(Z)

        call nyquist_1D(Z)
    end subroutine fftc_1D

    ! function: ifftc_1D
    !
    ! Computes a centered, one dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifftc_1D(Z)
        complex(dp), intent(inout) :: Z(:)

        call nyquist_1D(Z)

        call ifft_1D(Z)

        call nyquist_1D(Z)
    end subroutine ifftc_1D

    ! subroutine: nyquist_1D
    !
    ! Multiplies a 1D array by its Nyquist frequencies.
    !
    ! Paramters:
    !   Z - The array.
    subroutine nyquist_1D(Z)
        complex(dp), intent(inout) :: Z(:)

        integer :: i

        do i = 1, size(Z)
            if (mod(i, 2) == 0) Z(i) = -Z(i)
        end do

    end subroutine nyquist_1D

    ! function: fft_2D
    !
    ! Computes a two dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fft_2D(Z)
        complex(dp), intent(inout) :: Z(:,:)

        integer(i8b) :: plan

!         call fftw_initialize_2D(Z)
        call dfftw_plan_dft_2d(plan, size(Z,1), size(Z,2), Z, Z, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan)
!         call dfftw_destroy_plan(plan)
    end subroutine fft_2D

    ! function: ifft_2D
    !
    ! Computes a two dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifft_2D(Z)
        complex(dp), intent(inout) :: Z(:,:)

        integer(i8b) :: plan

!         call fftw_initialize_2D(Z)
        call dfftw_plan_dft_2d(plan, size(Z,1), size(Z,2), Z, Z, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan)
!         call dfftw_destroy_plan(plan)

        Z = Z  / size(Z)
    end subroutine ifft_2D

    ! function: fftc_2D
    !
    ! Computes a centered, two dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fftc_2D(Z)
        complex(dp), intent(inout) :: Z(:,:)

        call nyquist_2D(Z)
        call fft_2D(Z)
        call nyquist_2D(Z)
    end subroutine fftc_2D

    ! function: ifftc_2D
    !
    ! Computes a centered, two dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifftc_2D(Z)
        complex(dp), intent(inout) :: Z(:,:)

        call nyquist_2D(Z)
        call ifft_2D(Z)
        call nyquist_2D(Z)
    end subroutine ifftc_2D

    ! subroutine: nyquist_2D
    !
    ! Multiplies a 2D array by its Nyquist frequencies.
    !
    ! Paramters:
    !   Z - The array.
    subroutine nyquist_2D(Z)
        complex(dp), intent(inout) :: Z(:,:)

        integer :: i, j

        do j = 1, size(Z,2)
            do i = 1, size(Z,1)
                if (mod(i+j, 2) == 0) Z(i,j) = -Z(i,j)
            end do
        end do

    end subroutine nyquist_2D

    ! subroutine: fftw_initialize_2D
    subroutine fftw_initialize_2D(Z)
        complex(dp), intent(in) :: Z(:,:)

        logical, save :: initialized = .false.

        if (.not.initialized) then
            call fftw_initialize_sub_2D(size(Z,1),size(Z,2))
            initialized = .true.
        end if
    end subroutine fftw_initialize_2D

    ! subroutine: fftw_initialize_sub_2D
    subroutine fftw_initialize_sub_2D(i, j)
        integer, intent(in) :: i, j

        integer(i8b) :: plan
        integer :: isuccess
        complex(dp) :: T(i, j)

#ifdef HAVE_FC_LIBFFTW3_THREADS
        call dfftw_init_threads
#ifdef _OPENMP
        call dfftw_plan_with_nthreads(omp_get_max_threads())
#else
        call dfftw_plan_with_nthreads(HAVE_FC_LIBFFTW3_THREADS)
#endif
#endif
        call dfftw_import_system_wisdom(isuccess)
!         if (isuccess == 0) write(0,*) "FFTW: Problem loading the system wisdom."
        call dfftw_plan_dft_2d(plan, i, j, T, T, FFTW_FORWARD, FFTW_MEASURE)
        call dfftw_destroy_plan(plan)
        call dfftw_plan_dft_2d(plan, i, j, T, T, FFTW_BACKWARD, FFTW_MEASURE)
        call dfftw_destroy_plan(plan)

    end subroutine fftw_initialize_sub_2D

    ! function: fft_3D
    !
    ! Computes a three dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fft_3D(Z)
        complex(dp), intent(inout) :: Z(:,:,:)

        integer(i8b) :: plan

#ifndef XLF
        call fftw_initialize_3D(Z)
#endif
        call dfftw_plan_dft_3d(plan, size(Z,1), size(Z,2), size(Z,3), Z,Z, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan)
#ifndef XLF
        call dfftw_destroy_plan(plan)
#endif
    end subroutine fft_3D

    ! function: ifft_3D
    !
    ! Computes a three dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifft_3D(Z)
        complex(dp), intent(inout) :: Z(:,:,:)

        integer(i8b) :: plan

#ifndef XLF
        call fftw_initialize_3D(Z)
#endif
        call dfftw_plan_dft_3d(plan, size(Z,1), size(Z,2), size(Z,3), Z,Z, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan)
#ifndef XLF
        call dfftw_destroy_plan(plan)
#endif

        Z = Z  / size(Z)
    end subroutine ifft_3D

    ! function: fftc_3D
    !
    ! Computes a centered, three dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine fftc_3D(Z)
        complex(dp), intent(inout) :: Z(:,:,:)

        call nyquist_3D(Z)
        call fft_3D(Z)
        call nyquist_3D(Z)
    end subroutine fftc_3D

    ! function: ifftc_3D
    !
    ! Computes a centered, three dimensional, in-place FFT.
    !
    ! Parameters:
    !   Z - The array to be transformed
    subroutine ifftc_3D(Z)
        complex(dp), intent(inout) :: Z(:,:,:)

        call nyquist_3D(Z)
        call ifft_3D(Z)
        call nyquist_3D(Z)
    end subroutine ifftc_3D

    ! subroutine: nyquist_3D
    !
    ! Multiplies a 3D array by its Nyquist frequencies.
    !
    ! Paramters:
    !   Z - The array.
    subroutine nyquist_3D(Z)
        complex(dp), intent(inout) :: Z(:,:,:)

        integer :: i, j, k

        do k = 1, size(Z,3)
            do j = 1, size(Z,2)
                do i = 1, size(Z,1)
                    if (mod(i+j+k, 2) == 0) Z(i,j,k) = -Z(i,j,k)
                end do
            end do
        end do

    end subroutine nyquist_3D

    ! subroutine: fftw_initialize_3D
    subroutine fftw_initialize_3D(Z)
        complex(dp), intent(in) :: Z(:,:,:)

        logical, save :: initialized = .false.

        if (.not.initialized) then
            call fftw_initialize_sub_3D(size(Z,1),size(Z,2),size(Z,3))
            initialized = .true.
        end if
    end subroutine fftw_initialize_3D

    ! subroutine: fftw_initialize_sub_3D
    subroutine fftw_initialize_sub_3D(i, j, k)
        integer, intent(in) :: i, j, k

        integer(i8b) :: plan
        integer :: isuccess
        complex(dp) :: T(i, j, k)

#ifdef HAVE_FC_LIBFFTW3_THREADS
        call dfftw_init_threads
#ifdef _OPENMP
        call dfftw_plan_with_nthreads(omp_get_max_threads())
#else
        call dfftw_plan_with_nthreads(HAVE_FC_LIBFFTW3_THREADS)
#endif
#endif
!         write(0,*) "Initializing FFTW"
        call dfftw_import_system_wisdom(isuccess)
!         if (isuccess == 0) write(0,*) "FFTW: Problem loading the system wisdom."
        call dfftw_plan_dft_3d(plan, i, j, k, T, T, FFTW_FORWARD, FFTW_MEASURE)
        call dfftw_destroy_plan(plan)
        call dfftw_plan_dft_3d(plan, i, j, k, T, T, FFTW_BACKWARD, FFTW_MEASURE)
        call dfftw_destroy_plan(plan)

    end subroutine fftw_initialize_sub_3D

    ! subroutine: fft_field
    subroutine fft_field(e)
        complex(dp), intent(inout) :: e(:,:,:)

        integer :: i,j,k
        integer :: Nx, Ny, Nt

        Nx = size(e,1)
        Ny = size(e,2)
        Nt = size(e,3)

        do k = 1, Nt
            call fft(e(:,:,k))
        end do

        do j = 1, Ny
            do i = 1, Nx
                call ifft(e(i,j,:))
            end do
        end do

        e = e * Nt
    end subroutine fft_field

    ! subroutine: ifft_field
    subroutine ifft_field(e)
        complex(dp), intent(inout) :: e(:,:,:)

        integer :: i,j,k
        integer :: Nx, Ny, Nt

        Nx = size(e,1)
        Ny = size(e,2)
        Nt = size(e,3)

        do k = 1, Nt
            call ifft(e(:,:,k))
        end do

        do j = 1, Ny
            do i = 1, Nx
                call fft(e(i,j,:))
            end do
        end do

        e = e / Nt
    end subroutine ifft_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The rest of this module adds mixing Hankel transforms with FFTW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Transform(Z)
        complex(dp), intent(inout) :: Z(:,:,:)
        integer :: j, k

        if(size(Z,1) == 1) then
            if(.not.allocated(HT)) call CreateHT(size(Z,2))
            !$omp parallel do private(k)
            do k=1, size(Z,3)
                call HankelTransform(Z(1,:,k))
            end do
            !$omp end parallel do
            do j=1, size(Z,2)
                call FFT(Z(1,j,:))
            end do
        else
            call FFT(Z)
        end if       
    end subroutine


    subroutine iTransform(Z)
        complex(dp), intent(inout) :: Z(:,:,:)
        integer :: j, k

        if(size(Z,1) == 1) then
            if(.not.allocated(HT)) call CreateHT(size(Z,2))
            !$omp parallel do private(k)
            do k=1, size(Z,3)
                call HankelTransform(Z(1,:,k))
            end do
            !$omp end parallel do
            do j=1, size(Z,2)
                call iFFT(Z(1,j,:))
            end do
        else
            call iFFT(Z)
        end if       
    end subroutine


     subroutine HankelTransform(f)
         complex(dp), intent(inout) :: f(:)
 
         if(size(f) .ne. size(HT,2)) then
             print*, "ERROR:    Size[E(r)] =/ Size[HT] in fftw.F90"
             stop
         end if
 
         f=matmul(HT,f)

     end subroutine



     subroutine CreateHT(Nr)
         !!  USE IFLPORT FOR ifort COMPILER
         !use IFLPORT
         use typefield
         integer, intent(in) :: Nr
         integer             :: n,m
         real(dp)            :: J1
 
         if(allocated(HT)) return
         allocate(HT(Nr,Nr))
         allocate (a(Nr+1))
         a=GetJ0zeros(size(a))

         !! For ifort compiler
         !do n=1, Nr
         !    J1 = DBESJ1(a(n))
         !    do m=1, Nr
         !        HT(m,n) = (2.0_dp/a(Nr+1)) * DBESJ0(a(m)*a(n)/a(Nr+1)) / J1**2
         !    enddo
         !enddo

         !! For gfortran compiler
         do n=1, Nr
             J1 = Bessel_JN(1, a(n))
             do m=1, Nr
                 HT(m,n) = (2.0_dp/a(Nr+1)) * Bessel_JN(0, (a(m)*a(n)/a(Nr+1)) / J1**2)
             enddo
         enddo

     end subroutine

 

end module fftw
