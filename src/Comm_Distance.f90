subroutine computeD(n, r, d)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r
  !f2py intent(out) d
  !f2py depend(n) r, d

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r(3,n)

  !Output
  double precision, intent(out) ::  d(n)

  !Working
  integer i

  do i=1,n
     d(i) = dot_product(r(:,i), r(:,i))**0.5
  end do

end subroutine computeD



subroutine computeJacobianD(n, r, J)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r
  !f2py intent(out) J
  !f2py depend(n) r, J

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r(3,n)

  !Output
  double precision, intent(out) ::  J(n,3)

  !Working
  integer i
  double precision norm

  do i=1,n
     norm = dot_product(r(:,i), r(:,i))**0.5
     if (norm .gt. 1e-10) then
        J(i,:) = r(:,i) / norm
     else
        J(i,:) = 0
     end if
  end do

end subroutine computeJacobianD
