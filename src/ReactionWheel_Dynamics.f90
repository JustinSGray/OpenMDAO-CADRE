subroutine getf(dat, u, f)

  implicit none

  !Input
  double precision, intent(in) ::  dat(6), u(3)

  !Output
  double precision, intent(out) ::  f(3)

  !Working
  double precision wx(3,3)

  wx(1,:) = (/ dble(0), -dat(3),  dat(2) /)
  wx(2,:) = (/  dat(3), dble(0), -dat(1) /)
  wx(3,:) = (/ -dat(2),  dat(1), dble(0) /)

  f(:) = -dat(4:6) - matmul(wx, u)

end subroutine getf




subroutine getdfdy(dat, u, dfdu)

  implicit none

  !Input
  double precision, intent(in) ::  dat(6), u(3)

  !Output
  double precision, intent(out) ::  dfdu(3,3)

  !Working
  double precision wx(3,3)

  wx(1,:) = (/ dble(0), -dat(3),  dat(2) /)
  wx(2,:) = (/  dat(3), dble(0), -dat(1) /)
  wx(3,:) = (/ -dat(2),  dat(1), dble(0) /)

  !f(:) = -dat(4:6) - matmul(wx, u)
  dfdu = -wx

end subroutine getdfdy



subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(6), y(3)

  !Output
  double precision, intent(out) ::  dfdx(3,6)

  !Working
  double precision dwx_dw(3,3,3)
  integer k

  dwx_dw(1,:,1) = (/  0,  0,  0 /)
  dwx_dw(2,:,1) = (/  0,  0, -1 /)
  dwx_dw(3,:,1) = (/  0,  1,  0 /)

  dwx_dw(1,:,2) = (/  0,  0,  1 /)
  dwx_dw(2,:,2) = (/  0,  0,  0 /)
  dwx_dw(3,:,2) = (/ -1,  0,  0 /)

  dwx_dw(1,:,3) = (/  0, -1,  0 /)
  dwx_dw(2,:,3) = (/  1,  0,  0 /)
  dwx_dw(3,:,3) = (/  0,  0,  0 /)

  !f(:) = -dat(4:6) - matmul(wx, u)
  dfdx(:,:) = 0.0
  do k=1,3
     dfdx(:,k) = -matmul(dwx_dw(:,:,k), y)
     dfdx(k,k+3) = -1.0
  end do

end subroutine getdfdx
