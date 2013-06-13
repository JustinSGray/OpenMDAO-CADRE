subroutine computeqE(n, t, q_E)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, t
  !f2py intent(out) q_E
  !f2py depend(n) t, q_E

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  t(n)

  !Output
  double precision, intent(out) ::  q_E(4,n)

  !Working
  double precision pi, theta
  integer i
  
  pi = 2*acos(0.0)
  
  q_E(:,:) = 0.0
  do i=1,n
     theta = pi * t(i) / 3600.0 / 24.0
     q_E(1,i) = cos(theta)
     q_E(4,i) = -sin(theta)
  end do

end subroutine computeqE



subroutine computeJacobianqE(n, t, dq_dt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, t
  !f2py intent(out) dq_dt
  !f2py depend(n) t, dq_dt

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  t(n)

  !Output
  double precision, intent(out) ::  dq_dt(n,4)

  !Working
  double precision pi, theta, dtheta_dt
  integer i
  
  pi = 2*acos(0.0)

  dq_dt(:,:) = 0.0
  do i=1,n
     theta = pi * t(i) / 3600.0 / 24.0
     dtheta_dt = pi / 3600.0 / 24.0
     dq_dt(i,1) = -sin(theta) * dtheta_dt
     dq_dt(i,4) = -cos(theta) * dtheta_dt
  end do

end subroutine computeJacobianqE
