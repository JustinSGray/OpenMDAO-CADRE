subroutine computeP(n, w, T, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, w, T
  !f2py intent(out) P
  !f2py depend(n) w, T, P

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  w(3,n), T(3,n)

  !Output
  double precision, intent(out) ::  P(3,n)

  !Working
  integer i, k
  double precision V, a, b, I0

  V = 4.0
  a = 4.9e-4
  b = 4.5e2
  I0 = 0.017

  do i=1,n
     do k=1,3
        P(k,i) = V * (a * w(k,i) + b * T(k,i))**2 + V * I0
     end do
  end do

end subroutine computeP



subroutine computeJacobianP(n, w, T, dP_dw, dP_dT)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, w, T
  !f2py intent(out) dP_dw, dP_dT
  !f2py depend(n) w, T, dP_dw, dP_dT

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  w(3,n), T(3,n)

  !Output
  double precision, intent(out) ::  dP_dw(n,3), dP_dT(n,3)

  !Working
  integer i, k
  double precision V, a, b, I0

  V = 4.0
  a = 4.9e-4
  b = 4.5e2
  I0 = 0.017

  do i=1,n
     do k=1,3
        dP_dw(i,k) = 2 * V * a * (a * w(k,i) + b * T(k,i))
        dP_dT(i,k) = 2 * V * b * (a * w(k,i) + b * T(k,i))
     end do
  end do

end subroutine computeJacobianP
