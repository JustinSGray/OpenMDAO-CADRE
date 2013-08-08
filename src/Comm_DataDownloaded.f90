subroutine getf(dat, u, s)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), u(1)

  !Output
  double precision, intent(out) ::  s(1)

  s(1) = dat(1)

end subroutine getf




subroutine getdfdy(dat, u, dfdu)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), u(1)

  !Output
  double precision, intent(out) ::  dfdu(1,1)

  dfdu(1,1) = 0.0

end subroutine getdfdy



subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), y(1)

  !Output
  double precision, intent(out) ::  dfdx(1,1)

  dfdx(1,1) = 1.0

end subroutine getdfdx


