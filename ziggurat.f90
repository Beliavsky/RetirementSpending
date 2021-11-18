! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001
module ziggurat_mod
implicit none
private
integer,  parameter       ::  dp=selected_real_kind(12, 60)
real(kind=dp), parameter  ::  m1=2147483648.0_dp,   m2=2147483648.0_dp, &
                              half=0.5_dp, ve=0.003949659822581572_dp, &
                              vn=0.00991256303526217_dp
real(kind=dp)             ::  dn=3.442619855899_dp, tn=3.442619855899_dp, &
                              de=7.697117470131487_dp, &
                              te=7.697117470131487_dp
real(kind=dp)             ::  q
integer,  save            ::  iz, jz, jsr=123456789, kn(0:127), ke(0:255), hz
real(kind=dp), save       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
logical,  save            ::  initialized=.false.
public                    ::  zigset, shr3, uni, rnor, rexp, rnor_vec_inline, jsrset
interface uni
   module procedure uni_scalar, uni_vec, uni_vec_given_seed
end interface uni
interface rnor
   module procedure rnor_scalar,rnor_vec
end interface rnor
contains
subroutine jsrset(jsrseed)
integer, intent(in) :: jsrseed
!  set the seed
jsr = jsrseed
end subroutine jsrset
!
subroutine zigset(jsrseed)
integer, intent(in) :: jsrseed
integer             :: i
!  set the seed
jsr = jsrseed
!  tables for rnor
q = vn*exp(half*dn*dn)
kn(0) = (dn/q)*m1
kn(1) = 0
wn(0) = q/m1
wn(127) = dn/m1
fn(0) = 1.0_dp
fn(127) = exp(-half*dn*dn)
do  i = 126, 1, -1
   dn = sqrt(-2.0_dp * log(vn/dn + exp(-half*dn*dn)))
   kn(i+1) = (dn/tn)*m1
   tn = dn
   fn(i) = exp(-half*dn*dn)
   wn(i) = dn/m1
end do
!  tables for rexp
q = ve*exp(de)
ke(0) = (de/q)*m2
ke(1) = 0
we(0) = q/m2
we(255) = de/m2
fe(0) = 1.0_dp
fe(255) = exp(-de)
do  i = 254, 1, -1
   de = -log(ve/de + exp(-de))
   ke(i+1) = m2 * (de/te)
   te = de
   fe(i) = exp(-de)
   we(i) = de/m2
end do
initialized = .true.
end subroutine zigset

!  generate random 32-bit integers
function shr3() result(ival)
integer  ::  ival
jz = jsr
jsr = ieor(jsr, ishft(jsr,  13))
jsr = ieor(jsr, ishft(jsr, -17))
jsr = ieor(jsr, ishft(jsr,   5))
ival = jz + jsr
end function shr3

!  generate uniformly distributed random numbers
function uni_scalar() result(ran)
real(kind=dp)  ::  ran
ran = half + 0.2328306e-9_dp * shr3()
end function uni_scalar
!
function uni_vec(n) result(ran)
integer, intent(in) :: n
real(kind=dp)       :: ran(n)
integer             :: i
do i=1,n
   ran(i) = half + 0.2328306e-9_dp * shr3()   
end do
end function uni_vec
!
function uni_vec_given_seed(n,seed) result(ran)
! return n uniform variates
integer, intent(in) :: n
integer, intent(in) :: seed
real(kind=dp)       :: ran(n)
integer             :: i
jsr = seed
do i=1,n
   ran(i) = half + 0.2328306e-9_dp * shr3()   
end do
end function uni_vec_given_seed
!
function rnor_scalar() result(fn_val)
!  generate random normal variate 
real(kind=dp)             ::  fn_val
real(kind=dp), parameter  ::  r = 3.442620_dp
real(kind=dp)             ::  x, y
if (.not. initialized) call zigset(jsr)
hz = shr3()
iz = iand(hz, 127)
if (abs(hz) < kn(iz)) then
   fn_val = hz * wn(iz)
else
   do
      if (iz == 0) then
         do
            x = -0.2904764_dp* log(uni())
            y = -log(uni())
            if (y+y >= x*x) exit
         end do
         fn_val = r+x
         if (hz <= 0) fn_val = -fn_val
         return
      end if
      x = hz * wn(iz)
      if (fn(iz) + uni()*(fn(iz-1)-fn(iz)) < exp(-half*x*x)) then
         fn_val = x
         return
      end if
      hz = shr3()
      iz = iand(hz, 127)
      if (abs(hz) < kn(iz)) then
         fn_val = hz * wn(iz)
         return
      end if
   end do
end if
end function rnor_scalar
!
function rnor_vec_inline(n) result(fn_val)
!  generate n random normal variates
integer      , intent(in) ::  n
real(kind=dp)             ::  fn_val(n)
real(kind=dp), parameter  ::  r = 3.442620_dp
real(kind=dp)             ::  x, y
integer                   ::  i
if (.not. initialized) call zigset(jsr)
loop_variate: do i=1,n
hz = shr3()
iz = iand(hz, 127)
if (abs(hz) < kn(iz)) then
   fn_val(i) = hz * wn(iz)
else
   do
      if (iz == 0) then
         do
            x = -0.2904764_dp* log(uni())
            y = -log(uni())
            if (y+y >= x*x) exit
         end do
         fn_val(i) = r+x
         if (hz <= 0) fn_val(i) = -fn_val(i)
         cycle loop_variate
      end if
      x = hz * wn(iz)
      if (fn(iz) + uni()*(fn(iz-1)-fn(iz)) < exp(-half*x*x)) then
         fn_val(i) = x
         cycle loop_variate
      end if
      hz = shr3()
      iz = iand(hz, 127)
      if (abs(hz) < kn(iz)) then
         fn_val(i) = hz * wn(iz)
         cycle loop_variate
      end if
   end do
end if
end do loop_variate
end function rnor_vec_inline
!
function rnor_vec(n) result(ran)
! generate n random normal variates
integer, intent(in) :: n
real(kind=dp)       :: ran(n)
integer             :: i
do i=1,n
   ran(i) = rnor_scalar()
end do
end function rnor_vec
!
function rexp() result(fn_val)
!  generate random exponential variate
real(kind=dp)  ::  fn_val
real(kind=dp)  ::  x
if (.not. initialized) call zigset(jsr)
jz = shr3()
iz = iand(jz, 255)
if (abs(jz) < ke(iz)) then
   fn_val = abs(jz) * we(iz)
   return
end if
do
   if (iz == 0) then
      fn_val = 7.69711 - log(uni())
      return
   end if
   x = abs(jz) * we(iz)
   if (fe(iz) + uni()*(fe(iz-1) - fe(iz)) < exp(-x)) then
      fn_val = x
      return
   end if
   jz = shr3()
   iz = iand(jz, 255)
   if (abs(jz) < ke(iz)) then
      fn_val = abs(jz) * we(iz)
      return
   end if
end do
end function rexp
end module ziggurat_mod
