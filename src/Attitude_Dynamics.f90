subroutine normalize(n, y)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n
  !f2py intent(inout) y
  !f2py depend(n) y

  !Input
  integer, intent(in) ::  n

  !Output
  double precision, intent(inout) ::  y(4,n)

  !Working
  integer i
  double precision norm

  do i=1,n
     norm = dot_product(y(:,i), y(:,i))**0.5
     if (norm .gt. 1e-10) then
        y(:,i) = y(:,i)/norm
     end if
  end do

end subroutine normalize



subroutine getf(dat, y, f)

  implicit none

  !Input
  double precision, intent(in) ::  dat(3), y(4)

  !Output
  double precision, intent(out) ::  f(4)

  !Working
  double precision w1, w2, w3
  double precision A(4,4)

  w1 = dat(1)
  w2 = dat(2)
  w3 = dat(3)

  A(1,:) = (/ dble(0),     -w1,     -w2,     -w3 /)
  A(2,:) = (/      w1, dble(0),      w3,     -w2 /)
  A(3,:) = (/      w2,     -w3, dble(0),      w1 /)
  A(4,:) = (/      w3,      w2,     -w1, dble(0) /)

  f = 0.5 * matmul(A,y)

end subroutine getf



subroutine getdfdy(dat, y, dfdy)

  implicit none

  !Input
  double precision, intent(in) ::  dat(3), y(4)

  !Output
  double precision, intent(out) ::  dfdy(4,4)

  !Working
  double precision w1, w2, w3
  double precision A(4,4)

  w1 = dat(1)
  w2 = dat(2)
  w3 = dat(3)

  A(1,:) = (/ dble(0),     -w1,     -w2,     -w3 /)
  A(2,:) = (/      w1, dble(0),      w3,     -w2 /)
  A(3,:) = (/      w2,     -w3, dble(0),      w1 /)
  A(4,:) = (/      w3,      w2,     -w1, dble(0) /)

  dfdy = 0.5 * A

end subroutine getdfdy



subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(3), y(4)

  !Output
  double precision, intent(out) ::  dfdx(4,3)

  !Working
  double precision dA_dw(4,4,3)
  integer k

  dA_dw(1,:,1) = (/  0, -1,  0,  0 /)
  dA_dw(2,:,1) = (/  1,  0,  0,  0 /)
  dA_dw(3,:,1) = (/  0,  0,  0,  1 /)
  dA_dw(4,:,1) = (/  0,  0, -1,  0 /)

  dA_dw(1,:,2) = (/  0,  0, -1,  0 /)
  dA_dw(2,:,2) = (/  0,  0,  0, -1 /)
  dA_dw(3,:,2) = (/  1,  0,  0,  0 /)
  dA_dw(4,:,2) = (/  0,  1,  0,  0 /)

  dA_dw(1,:,3) = (/  0,  0,  0, -1 /)
  dA_dw(2,:,3) = (/  0,  0,  1,  0 /)
  dA_dw(3,:,3) = (/  0, -1,  0,  0 /)
  dA_dw(4,:,3) = (/  1,  0,  0,  0 /)

  do k=1,3
     dfdx(:,k) = 0.5 * matmul(dA_dw(:,:,k),y)
  end do

end subroutine getdfdx

