module kind_mod
implicit none
private
public :: dp,bad_real,long_int,i32,i64,bad_int
integer, parameter :: dp = selected_real_kind(15, 307)
! integer, parameter :: dp = selected_real_kind(32) ! to get quadruple precision
real(kind=dp), parameter :: bad_real = -999.0_dp
integer, parameter :: i32=selected_int_kind(9), long_int = selected_int_kind(15), i64 = long_int, &
                      bad_int = -999
end module kind_mod
