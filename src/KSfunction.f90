subroutine computeKS(n, rho, g, KS)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, rho, g
  !f2py intent(out) KS
  !f2py depend(n) g

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  rho, g(n)

  !Output
  double precision, intent(out) ::  KS

  !Working
  integer i, imax
  double precision sum

  imax = maxloc(g, dim=1)

  sum = 0
  do i=1,n
     sum = sum + exp(rho * (g(i) - g(imax)))
  end do

  KS = g(imax) + 1.0/rho * log(sum)

end subroutine computeKS




subroutine computeJacobianKS(n, rho, g, dKS_dg)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, rho, g
  !f2py intent(out) dKS_dg
  !f2py depend(n) g, dKS_dg

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  rho, g(n)

  !Output
  double precision, intent(out) ::  dKS_dg(n)

  !Working
  integer i, imax
  double precision sum
  double precision dsum_dg(n), dKS_dsum

  imax = maxloc(g, dim=1)

  sum = 0
  do i=1,n
     sum = sum + exp(rho * (g(i) - g(imax)))
     dsum_dg(i) = rho * exp(rho * (g(i) - g(imax)))
  end do

  !KS = g(imax) + 1.0/rho * log(sum)
  dKS_dsum = 1.0 / rho / sum
  dKS_dg(:) = dKS_dsum * dsum_dg(:)

end subroutine computeJacobianKS
