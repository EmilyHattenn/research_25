! module: types
! Kinds lifted from numerical recipes nrtype.f90
module types
    ! const: i8b
    ! 64-bit integer
    integer, parameter :: i8b = selected_int_kind(18)

    ! const: i4b
    ! 32-bit integer
    integer, parameter :: i4b = selected_int_kind(9)

    ! const: i2b
    ! 16-bit integer
    integer, parameter :: i2b = selected_int_kind(4)

    ! const: i1b
    ! 8-bit integer
    integer, parameter :: i1b = selected_int_kind(2)

    ! const: sp
    ! Single precision real
    integer, parameter :: sp = kind(1.0)

    ! const: dp
    ! Double precision real
    integer, parameter :: dp = kind(1.0d0)

    ! const: qp
    ! Quad precision real (?)
    integer, parameter :: qp = 16
end module types
