subroutine getConstants(alpha)

  implicit none

  double precision, intent(out) ::  alpha
  double precision pi, c, Gr, Ll, f, k, SNR, T
  
  pi = 2*acos(0.0)
  c = 299792458
  Gr = 10**(12.9/10.0)
  Ll = 10**(-2.0/10.0)
  f = 437e6
  k = 1.3806503e-23
  SNR = 10**(5.0/10.0)
  T = 500.0

  alpha = c**2 * Gr * Ll / 16.0 / pi**2 / f**2 / k / SNR / T / 1e6

end subroutine getConstants



subroutine computeDr(n, P, Gt, S, LOS, Dr)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, P, Gt, S, LOS
  !f2py intent(out) Dr
  !f2py depend(n) P, Gt, S, LOS, Dr

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P(n), Gt(n), S(n), LOS(n)

  !Output
  double precision, intent(out) ::  Dr(n)

  !Working
  integer i
  double precision alpha, S2

  call getConstants(alpha)

  do i=1,n
     if (abs(S(i)) .gt. 1e-10) then
        S2 = S(i) * 1e3
     else
        S2 = 1e-10
     end if
     Dr(i) = alpha * P(i) * Gt(i) * LOS(i) / S2**2
  end do

end subroutine computeDr



subroutine computeJacobianDr(n, P, Gt, S, LOS, dD_dP, dD_dGt, dD_dS, dD_dLOS)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, P, Gt, S, LOS
  !f2py intent(out) dD_dP, dD_dGt, dD_dS, dD_dLOS
  !f2py depend(n) P, Gt, S, LOS, dD_dP, dD_dGt, dD_dS, dD_dLOS

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P(n), Gt(n), S(n), LOS(n)

  !Output
  double precision, intent(out) ::  dD_dP(n), dD_dGt(n), dD_dS(n), dD_dLOS(n)

  !Working
  integer i
  double precision alpha, S2

  call getConstants(alpha)

  do i=1,n
     if (abs(S(i)) .gt. 1e-10) then
        S2 = S(i) * 1e3
     else
        S2 = 1e-10
     end if
     !Dr(i) = alpha * P(i) * Gt(i) * LOS(i) / S2**2
     dD_dP(i) = alpha * Gt(i) * LOS(i) / S2**2
     dD_dGt(i) = alpha * P(i) * LOS(i) / S2**2
     dD_dS(i) = - 2 * alpha * P(i) * Gt(i) * LOS(i) / S2**3 * 1e3
     dD_dLOS(i) = alpha * P(i) * Gt(i) / S2**2
  end do

end subroutine computeJacobianDr
