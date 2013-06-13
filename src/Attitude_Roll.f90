subroutine computeOBR(n, g, O)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, g
  !f2py intent(out) O
  !f2py depend(n) g, O

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  g(n)

  !Output
  double precision, intent(out) ::  O(3,3,n)

  !Working
  integer i

  do i=1,n
     O(1,:,i) = (/  cos(g(i)) ,  sin(g(i)) , dble(0) /)
     O(2,:,i) = (/ -sin(g(i)) ,  cos(g(i)) , dble(0) /)
     O(3,:,i) = (/   dble(0)  ,   dble(0)  , dble(1) /)
  end do

end subroutine computeOBR




subroutine computeJacobianOBR(n, g, dO_dg)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, g
  !f2py intent(out) dO_dg
  !f2py depend(n) g, dO_dg

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  g(n)

  !Output
  double precision, intent(out) ::  dO_dg(n,3,3)

  !Working
  integer i

  do i=1,n
     dO_dg(i,1,:) = (/ -sin(g(i)) ,  cos(g(i)) , dble(0) /)
     dO_dg(i,2,:) = (/ -cos(g(i)) , -sin(g(i)) , dble(0) /)
     dO_dg(i,3,:) = (/   dble(0)  ,   dble(0)  , dble(0) /)
  end do

end subroutine computeJacobianOBR
