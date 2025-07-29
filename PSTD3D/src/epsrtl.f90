module epsrtl

    use types       ! Defines data types (double precision = (dp), etc.)
    use constants   ! Defines common math/physics constants (pi, e0, me0, hbar, etc.)
    use fftw        ! Contains functions for fast Fourier transform
    use helpers
    use usefulsubs
    use qwoptics
    use coulomb
    implicit none

    integer,  private :: Nw = 500
    real(dp), private :: dw = 2.35e15 / 500d0 * 2
    real(dp), private :: R0 = 1d-9
    real(dp), private :: g  = 1d0 / 1d-12
    real(dp), private :: kB = 1.38064852d-23
    real(dp), private :: epsb = 3.011**2
    real(dp), private :: dcv
    real(dp), private :: n00
    real(dp), private :: kf

contains

    subroutine GetEpsrLEpsrT(n1D, dcv0, Te, Th, me, mh, Eg, ky)
        real(dp), intent(in   ) :: n1D, dcv0, Te, Th, me, mh, Eg, ky(:)
        real(dp), allocatable   :: qy(:)
        integer                 :: Nq, n
        real(dp)                :: dq
        character(len=1)        :: B

        Nq = (size(ky) - 1)*100 + 1
        allocate(qy(Nq))
        dq = (ky(2)-Ky(1))/50d0

        do n=1, Nq
            qy(n) = - (Nq-1)/2d0*dq + (n-1)*dq
        end do

        n00 = n1D
        kf  = n00 / 2d0
        dcv = dcv0

        B="E"
        call ZeroT_L(B,me, qy, kf)
        B="H"
        call ZeroT_L(B,mh, qy, kf)

        call ZeroT_T(me, mh, Eg, dcv0, qy, kf)

!stop
!        call RecordEpsrL_T0(me, mh, qy(:))
!stop
!        call RecordEpsrL(Te, Th, me, mh, qy(:))



!        call RecordEpsrT(Te, Th, me, mh, Eg, qy(:))

        open(unit=922,file='dataQW/Wire/qw.dat')
        write(922,*) "Nq", Nq
        write(922,*) "Nw", Nw*2+1
        write(922,*) "ky(1)", qy(1)
        write(922,*) "dky"  , qy(2)-qy(1)
        write(922,*) "w(1)", -Nw*dw
        write(922,*) "dw"  , dw
        close(922)

        stop
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    subroutine RecordEpsrT(Te, Th, me, mh, Eg, ky)
        real(dp), intent(in   ) :: Te, Th, me, mh, Eg, ky(:)
        real(dp)                :: epsR(size(ky),-Nw:Nw), epsI(size(ky),-Nw:Nw)
        real(dp)                :: Ek(size(ky)), Ekq(size(ky))
        real(dp)                :: dk, ww, T, m, a, b
        complex(dp)             :: tmp
        integer                 :: q, w

        a = 1d0 / eps0 / epsb * dcv**2 / pi / R0**2
        b = 0d0

        epsR = 0d0
        epsI = 0d0
        dk   = ky(2)-ky(1)
        g    = 1d12

        Ek = Eng(mh,ky)

        do w= -Nw, Nw
print*, "T w", w
            ww   =  w*dw
            
            !$omp parallel do private(q, Ekq, tmp)
            do q=1, size(ky)
                Ekq  =  Eng(me,ky+ky(q)) + Eg
                tmp  =  PiT(ky(q), ww, me, mh, Te, Th, dk, Ek, Ekq)

                epsR(q,w) = 1d0 - a *  real(tmp)
                epsI(q,w) =     - a * aimag(tmp)
            end do 
            !$omp end parallel do
        end do

        print*, "min val EpsT real", "max val EpsT real"
        print*, minval(epsR), maxval(epsR)
        print*, "min val EpsT imag", "max val EpsT imag"
        print*, minval(epsI), maxval(epsI)

        open(unit=922,file='dataQW/Wire/EpsT.dat')
        do w= -Nw, Nw
            do q=1, size(ky)
                write(922,*) epsR(q,w), epsI(q,w) 
            end do
        end do
        close(922)

open(unit=937,file="chi.0.w.dat")
do w=-Nw,Nw
    write(937,*) w, epsR(Floor(1141d0/2d0+2),w)
end do
 close(937)

    end subroutine


    function PiT(q, w, me, mh, Te, Th, dk, Ek, Ekq)
        real(dp), intent(in) :: q, w, me, mh, Te, Th, dk, Ek(:), Ekq(:)
        complex(dp) :: PiT
        real(dp) :: a

        a =  2d0 / pi * dk

        g = 2.35d15 / 1000d0

        PiT =  a * sum((1d0 - ff0(Ek(:),Th,mh)- ff0(Ekq(:),Te,me)) * &
                       ( + 1d0 / ( hbar*w - Ekq(:)-Ek(:) + ii*hbar*g ) &
                         - 0d0 / ( hbar*w + Ekq(:)+Ek(:) - ii*hbar*g ) ) )
                                             
    end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine RecordEpsrL(Te, Th, me, mh, ky)
        real(dp), intent(in   ) :: Te, Th, me, mh, ky(:)
        complex(dp)             :: eps(size(ky),-Nw:Nw)
        complex(dp)             :: PiE(size(ky),-Nw:Nw)
        complex(dp)             :: PiH(size(ky),-Nw:Nw)

        real(dp)                :: tmpR(size(ky),-Nw:Nw), tmpI(size(ky),-Nw:Nw)
        real(dp)                :: Ek(size(ky)), Ekq(size(ky))
        real(dp)                :: Vc(size(ky))
        real(dp)                :: dk, ww, T, m
        real(dp)                :: PILqw, Qq
        complex(dp)             :: tmp
        integer                 :: q, w


        eps = 1d0
        PiE = 0d0
        PiH = 0d0

        dk   = ky(2)-ky(1)

        do q=1, size(ky)
            Vc(q) = (e0**2 / twopi / eps0 / epsb) * K03( max(abs(ky(q)*R0), 1d-10) )
        end do

        ! For electrons
        T  = Te
        m  = me
        Ek = Eng(m,ky)
        do w= -Nw, Nw
            print*, "L e w", w
            ww   =  w*dw
            !$omp parallel do private(q, Ekq)
            do q=1, size(ky)
                Ekq       =  Eng(m,ky+ky(q))
                PiE(q,w)  =  PiL(ky(q), ww, m, T, dk, Ek, Ekq)
            end do 
            !$omp end parallel do
        end do


        !call QqGq(ky,size(ky),dk,dw,1d0-EpsR,-EpsI,'e')


        ! For Holes
        T  = Th
        m  = mh
        Ek = Eng(m,ky)
        do w= -Nw, Nw
            ww   =  w*dw
            print*, "L h w", w
            !$omp parallel do private(q, Ekq)
            do q=1, size(ky)
                Ekq       =  Eng(m,ky+ky(q))
                PiH(q,w)  =  PiL(ky(q), ww, m, T, dk, Ek, Ekq)
            end do 
            !$omp end parallel do
        end do

        !call QqGq(ky,size(ky),dk,dw,1d0-tmpR,-tmpI,'h')

        !$omp parallel do private(ww, q)
        do w= -Nw, Nw
            do q=1, size(ky)
                eps(q,w) = 1d0 - Vc(q) * PiE(q,w) - Vc(q) * PiH(q,w)
            end do 
        end do 
        !$omp end parallel do

        !call QqGq(ky,size(ky),dk,dw,epsR,epsI,'eh')

        print*, "min val EpsL real", "max val EpsL real"
        print*, minval(real(eps)), maxval(real(eps))
        print*, "min val EpsL imag", "max val EpsL imag"
        print*, minval(aimag(eps)), maxval(aimag(eps))


        open(unit=922,file='dataQW/Wire/EpsL.dat')
        do w= -Nw, Nw
            do q=1, size(ky)
                write(922,*) real(eps(q,w)), aimag(eps(q,w))
            end do
        end do
        close(922)



    end subroutine





    function PiL(q, w, m, T, dk, Ek, Ekq)
        real(dp), intent(in) :: q, w, m, T, dk, Ek(:), Ekq(:)
        complex(dp)          :: PiL

        PiL = 0d0

        g = 0.01*e0*1d-3 / hbar



        PiL =  2d0/pi * dk * (sum((ff0(Ek(:),T,m)- ff0(Ekq(:),T,m)) * &
                                       ( + 1d0 / ( hbar*w - (Ekq(:)-Ek(:)) + ii*hbar*g &
                                         - 0d0 / ( hbar*w + (Ekq(:)-Ek(:)) - ii*hbar*g) ) ) ) )
                                             
    end function



    subroutine RecordEpsrL_T0(me, ky)
        real(dp), intent(in) :: me, ky(:)
        complex(dp)             :: eps(size(ky),-Nw:Nw)
        real(dp)               :: PiE(size(ky),-Nw:Nw)
        real(dp)               :: PiH(size(ky),-Nw:Nw)

        real(dp)                :: tmpR(size(ky),-Nw:Nw), tmpI(size(ky),-Nw:Nw)
        real(dp)                :: Ek(size(ky)), Ekq(size(ky))
        real(dp)                :: Vc(size(ky))
        real(dp)                :: dk, ww, T, m
        real(dp)                :: PILqw, Qq
        complex(dp)             :: tmp
        integer                 :: q, w, Nk, k

        Nk = size(ky)

        eps = 1d0
        PiE = 0d0
        PiH = 0d0
        dk  = ky(2)-ky(1)

        do q=1, size(ky)
            Vc(q) = (e0**2 / twopi / eps0 / epsb) * K03( max(abs(ky(q)*R0), 1d-10) )
        end do

        
        do w=-Nw, Nw
            do q=-Nk, Nk
                do k=-Nk, Nk
                    
                end do
            end do
        end do
            
            



!        ! For electrons
!        T  = 0d0
!        m  = me
!        Ek = Eng(m,ky)
!        do w= -Nw, Nw
!            print*, "L e w", w
!            ww   =  w*dw
!            !$omp parallel do private(q, Ekq)
!            do q=1, size(ky)
!                Ekq       =  Eng(m,ky+ky(q))
!                PiE(q,w)  =  PiL_T0(ky(q), ww, m, T, dk, Ek, Ekq)
!            end do 
!            !$omp end parallel do
!        end do


        !call QqGq(ky,size(ky),dk,dw,1d0-EpsR,-EpsI,'e')


        ! For Holes
!        T  = 0d0
!        m  = mh
!        Ek = Eng(m,ky)
!        do w= -Nw, Nw
!            ww   =  w*dw
!            print*, "L h w", w
!            !$omp parallel do private(q, Ekq)
!            do q=1, size(ky)
!                Ekq       =  Eng(m,ky+ky(q))!
!                PiH(q,w)  =  PiL_T0(ky(q), ww, m, T, dk, Ek, Ekq)
!            end do 
!            !$omp end parallel do
!        end do

        !call QqGq(ky,size(ky),dk,dw,1d0-tmpR,-tmpI,'h')

        !$omp parallel do private(ww, q)
        do w= -Nw, Nw
            do q=1, size(ky)
                eps(q,w) = 1d0 - Vc(q) * PiE(q,w) - Vc(q) * PiH(q,w)
            end do 
        end do 
        !$omp end parallel do

        !call QqGq(ky,size(ky),dk,dw,epsR,epsI,'eh')

        print*, "min val EpsL real", "max val EpsL real"
        print*, minval(real(eps)), maxval(real(eps))
        print*, "min val EpsL imag", "max val EpsL imag"
        print*, minval(aimag(eps)), maxval(aimag(eps))



        open(unit=922,file='dataQW/Wire/EpsL.dat')
        do w= -Nw, Nw
print*, "writing", w
            do q=1, size(ky)
                write(922,*) real(eps(q,w)), aimag(eps(q,w))
            end do
        end do
        close(922)

    end subroutine




    function PiL_T0(q, w, m, T, dk, Ek, Ekq)
        !use IFLPORT

        real(dp), intent(in) :: q, w, m, T, dk, Ek(:), Ekq(:)
        complex(dp)          :: PiL_T0
        real(dp)             :: q_inv

        PiL_T0 = 0d0
        q_inv     = q / (q**2 + dk**2)

        g = 0.01*e0*1d-3 / hbar


        PiL_T0 =  - m / pi / hbar**2 * q_inv * &
                  ( &
                      atanhc(2*hbar*kf*q / (ii*g*m + hbar*q**2 - m*w) ) + &
                      atanhc(2*hbar*kf*q / (ii*g*m + hbar*q**2 + m*w) )   &
                  )
                                             
    end function


    elemental function atanhc(x)
        complex(dp), intent(in) :: x
        complex(dp)             :: atanhc

        atanhc = 0.5d0 * log( (1d0+x)/(1d0 - x) )
    end function


    subroutine QqGq(ky,Nk,dk,dw,EpsR,EpsI,eh)
        real(dp),         intent(in) :: ky(:)
        integer,          intent(in) :: Nk
        real(dp),         intent(in) :: dk, dw, EpsR(:,:), EpsI(:,:)
        character(len=*), intent(in) :: eh
        integer                      :: q, k, w
        real(dp)                     :: Omega(Nk), Gam(Nk), depsRdw(-Nw:Nw)
        real(dp)                     :: tmp

        Omega = 0d0
        Gam   = 0d0

        dEpsRdw = 0d0
        tmp     = 0d0

!        do q=Nk/2+10, Nk/2+10
         do q=1, Nk

            dEpsRdw(:) = ( cshift(EpsR(q,:),1) - cshift(EpsR(q,:),-1) ) / 2d0 / dw

            tmp     = 0d0

            do w=-Nw, Nw

!print*, q, w, epsR(q,w)   !   This doesn't produce what I saw on 2D graphs earlier!!!!!


                if( 1d0 / (abs(EpsR(q,w)) + 0d-3) > tmp) then

                    tmp      = 1d0 / (abs(EpsR(q,w)) + 0d-3)

                    Omega(q) = w*dw
                    Gam(q)   = epsI(q,w) / dEpsRdw(w)

                end if 
            end do
        end do
!stop
        open(unit=739,file='dataQW/Wire/Omega_qp.'//eh//'.dat')
        do q=1, Nk
	    !write(739,*) ky(q), maxval(abs(epsR(q,:))), 1d0/ (minval(abs(epsR(q,:))) + 1d-10)
            write(739,*) ky(q), Omega(q), Gam(q)
        end do
        close(739)
    end subroutine



    subroutine ZeroT_L(B,m, qy, kf)
        Character(len=1), intent(In) :: B
        real(dp), intent(in) :: m, qy(:), kf
        real(dp)             :: Pi1(size(qy),-Nw:Nw)
        real(dp)             :: Pi2(size(qy),-Nw:Nw)
        complex(dp)          :: eps(size(qy),-Nw:Nw)
        real(dp)             :: q, dq, k, Aq
        integer              :: qq, ww
        complex(dp)          :: hw, xqw
        real(dp)             :: Vc(size(qy))

        dq = qy(2)-qy(1)

        Vc=0d0
        do qq=1, size(qy)
            Vc(qq) = (e0**2 / twopi / eps0 / epsb) * K03( max(abs(qy(qq)*R0), dq*R0) )
        end do

        !  Calculate Pi_1_L and Pi_2_L and eps
        !$omp parallel do private(ww, qq, hw, q, Aq, k)
        do ww=-Nw, Nw
            do qq=1, size(qy)

                hw  = hbar * ww*dw + ii * 1d-4 * e0
                q   = qy(qq)
                Aq  = m / hbar**2 * q / (q**2 + dq**2)
                xqw = (Eng(m,kf-q)-Eng(m,kf)-hw) * (Eng(m,kf-q)-Eng(m,kf)+hw) /   &
                      (Eng(m,kf+q)-Eng(m,kf)-hw) / (Eng(m,kf+q)-Eng(m,kf)+hw)
                k   = real(Aq * (hw - Eng(m,q)))

                Pi1(qq,ww) = Aq / pi * real(log( xqw ))

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Aq  = m / hbar**2 * q / (q**2 + dq**2)
                xqw = (-Eng(m,kf-q)+Eng(m,kf)-hw) * (-Eng(m,kf-q)+Eng(m,kf)+hw) /   &
                      (-Eng(m,kf+q)+Eng(m,kf)-hw) / (-Eng(m,kf+q)+Eng(m,kf)+hw)
                k   = real(Aq * (hw + Eng(m,q)))

                Pi1(qq,ww) = Pi1(qq,ww) + Aq / pi * real(log( xqw ))                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                Pi2(qq,ww) = - Aq * (fT0(k,kf) - fT0(k+q,kf))

                !Pi2(qq,ww) = Aq / pi * aimag(log( xqw ))

                eps(qq,ww) = - Vc(qq) * (Pi1(qq,ww) + ii*Pi2(qq,ww)) 

            end do
         end do
        !$omp end parallel do

        open(unit=922,file='dataQW/Wire/ChiL.'//B//'.dat')
        do ww= -Nw, Nw
            print*, "writing", ww
            do qq=1, size(qy)
                write(922,*) real(eps(qq,ww)), aimag(eps(qq,ww))
            end do
        end do
        close(922)
                

    end subroutine



    subroutine ZeroT_T(me, mh, Egap, dcv, qy, kf)
        real(dp), intent(in) :: me, mh, Egap, dcv, qy(:), kf
        real(dp)             :: Pi1(size(qy),-Nw:Nw)
        real(dp)             :: Pi2(size(qy),-Nw:Nw)
        complex(dp)          :: Pi3(size(qy),-Nw:Nw)
        real(dp)             :: q, dq, k, Aq
        integer              :: qq, ww, kk
        complex(dp)          :: hw, a, d
        real(dp)             :: Vc, b, c

        dq = qy(2)-qy(1)

        Vc = dcv**2 / eps0 / epsb / pi / R0**2
        b  = hbar**2 / 2d0 / me
        c  = hbar**2 / 2d0 / mh

        Pi1 = 0d0
        Pi2 = 0d0
        Pi3 = 0d0

        !  Calculate Pi_1_L and Pi_2_L and eps
        !$omp parallel do private(ww, qq, hw, q, a,d, k, kk)
        do ww=-Nw, Nw
print*, ww
            do qq=1, size(qy)

                q  = qy(qq)                
                a  = hbar * ww*dw - Egap + ii * 1d-3 * e0
                d  =  sqrt( a*(b+c) + b*c*q**2 )

                Pi1(qq,ww) = real(- Vc * atanJG((+kf*(b+c) + b*q)/d) / d &
                                  + Vc * atanJG((-kf*(b+c) + b*q)/d) / d &
                                  - Vc * atanJG((+kf*(b+c) - c*q)/d) / d &
                                  + Vc * atanJG((-kf*(b+c) - c*q)/d) / d )

                Pi2(qq,ww) = q**2 * c0**2 / epsb / ((ww*dw)**2 + 9*dw**2)

                do kk=1, size(qy)
                    k = qy(kk)
                    Pi3(qq,ww) = Pi3(qq,ww) + Vc*dq * (1d0 - fT0(k,kf)-fT0(k+q,kf)) * &
                                             (+1d0/(hbar*ww*dw - Egap - Eng(me,k+q) - Eng(mh,k) + ii*e0*5d-3) &
                                              -1d0/(hbar*ww*dw + Egap + Eng(me,k+q) + Eng(mh,k) + ii*e0*5d-3) )
                end do
            end do
         end do
        !$omp end parallel do

        open(unit=922,file='dataQW/Wire/ChiT.dat')
        do ww= -Nw, Nw
            print*, "writing", ww
            do qq=1, size(qy)
                write(922,*) Pi2(qq,ww), -real(Pi3(qq,ww))
            end do
        end do
        close(922)
    end subroutine



    elemental function atanJG(z)
        complex(dp), intent(in) :: z
        complex(dp)             :: atanJG

        atanJG = log((ii-z)/(ii+z)) / (2*ii)
    end function



    elemental function Eng(m,k)
        real(dp), intent(in) :: m, k
        real(dp) :: Eng
        Eng =  hbar**2 * k**2 / 2d0 / m
    end function


    elemental function fT0(k,kf)
        real(dp), intent(in) :: k, kf
        real(dp)             :: fT0
        fT0 = 1d0 - theta(abs(k)-kf)
    end function


    elemental Function ff0(E,T,m)
        real(dp), intent(in) :: E, T, m
        real(dp) :: ff0

        ff0 = n00 * sqrt(hbar**2 / twopi / m / kB / T) * Exp(-E/kB/T)
    end function

end module
