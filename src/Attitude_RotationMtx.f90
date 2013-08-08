subroutine computeOBI(n, O1, O2, O)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, O1, O2
  !f2py intent(out) O
  !f2py depend(n) O1, O2, O

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  O1(3,3,n), O2(3,3,n)

  !Output
  double precision, intent(out) ::  O(3,3,n)

  !Working
  integer i

  do i=1,n
     O(:,:,i) = matmul(O1(:,:,i), O2(:,:,i))
  end do

end subroutine computeOBI
