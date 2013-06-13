subroutine computeqA(n, t, q_A)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, t
  !f2py intent(out) q_A
  !f2py depend(n) q_A

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  t

  !Output
  double precision, intent(out) ::  q_A(4,n)

  !Working
  double precision rt2
  integer i
  
  rt2 = 2**0.5
  
  q_A(:,:) = 0.0
  do i=1,n
     q_A(1,i) = cos(t/2)
     q_A(2,i) = sin(t/2) / rt2
     q_A(3,i) = - sin(t/2) / rt2
     q_A(4,i) = 0.0
  end do

end subroutine computeqA




subroutine computeJacobianqA(t, dq_dt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) t
  !f2py intent(out) dq_dt

  !Input
  double precision, intent(in) ::  t

  !Output
  double precision, intent(out) ::  dq_dt(4)

  !Working
  double precision rt2
  
  rt2 = 2**0.5
  
  dq_dt(1) = - sin(t/2) / 2
  dq_dt(2) = cos(t/2) / rt2 / 2
  dq_dt(3) = - cos(t/2) / rt2 / 2
  dq_dt(4) = 0.0

end subroutine computeJacobianqA
