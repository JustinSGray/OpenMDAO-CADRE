subroutine getInertiaMtx(J)

  implicit none

  double precision, intent(out) ::  J(3,3)

  J(1,:) = (/ 0.018,   0.0,   0.0 /)
  J(2,:) = (/   0.0, 0.018,   0.0 /)
  J(3,:) = (/   0.0,   0.0, 0.006 /)

end subroutine getInertiaMtx




subroutine computeT(n, w, wdot, T)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, w, wdot
  !f2py intent(out) T
  !f2py depend(n) w, wdot, T

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  w(3,n), wdot(3,n)

  !Output
  double precision, intent(out) ::  T(3,n)

  !Working
  double precision J(3,3), wx(3,3)
  integer i

  call getInertiaMtx(J)

  do i=1,n
     wx(1,:) = (/ dble(0), -w(3,i),  w(2,i) /)
     wx(2,:) = (/  w(3,i), dble(0), -w(1,i) /)
     wx(3,:) = (/ -w(2,i),  w(1,i), dble(0) /)
     T(:,i) = matmul(J,wdot(:,i)) + matmul(wx, matmul(J,w(:,i)))
  end do

end subroutine computeT




subroutine computeJacobianT(n, w, wdot, dT_dw, dT_dwdot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, w, wdot
  !f2py intent(out) dT_dw, dT_dwdot
  !f2py depend(n) w, wdot, dT_dw, dT_dwdot

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  w(3,n), wdot(3,n)

  !Output
  double precision, intent(out) ::  dT_dw(n,3,3), dT_dwdot(n,3,3)

  !Working
  double precision J(3,3), wx(3,3), dwx_dw(3,3,3)
  integer i, k

  call getInertiaMtx(J)

  dwx_dw(1,:,1) = (/  0,  0,  0 /)
  dwx_dw(2,:,1) = (/  0,  0, -1 /)
  dwx_dw(3,:,1) = (/  0,  1,  0 /)

  dwx_dw(1,:,2) = (/  0,  0,  1 /)
  dwx_dw(2,:,2) = (/  0,  0,  0 /)
  dwx_dw(3,:,2) = (/ -1,  0,  0 /)

  dwx_dw(1,:,3) = (/  0, -1,  0 /)
  dwx_dw(2,:,3) = (/  1,  0,  0 /)
  dwx_dw(3,:,3) = (/  0,  0,  0 /)

  do i=1,n
     wx(1,:) = (/ dble(0), -w(3,i),  w(2,i) /)
     wx(2,:) = (/  w(3,i), dble(0), -w(1,i) /)
     wx(3,:) = (/ -w(2,i),  w(1,i), dble(0) /)
     
     !T(:,i) = matmul(J,wdot(:,i)) + matmul(wx, matmul(J,w(:,i)))
     dT_dwdot(i,:,:) = J
     dT_dw(i,:,:) = matmul(wx, J)
     do k=1,3
        dT_dw(i,:,k) = dT_dw(i,:,k) + matmul(dwx_dw(:,:,k), matmul(J,w(:,i)))
     end do
  end do

end subroutine computeJacobianT
