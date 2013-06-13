subroutine computeFwdV(n, v, f)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, v
  !f2py intent(out) f
  !f2py depend(n) v, f

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  v(3,n)

  !Output
  double precision, intent(out) ::  f(n)

  !Working
  integer i
  double precision norm

  do i=1,n
     norm = dot_product(v(:,i), v(:,i))**0.5
     if (norm .lt. 1e-15) then
        norm = 1e-5
     end if
     f(i) = - v(3,i) / norm
  end do

end subroutine computeFwdV



subroutine computeJacobianFwdV(n, v, df_dv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, v
  !f2py intent(out) df_dv
  !f2py depend(n) v, df_dv

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  v(3,n)

  !Output
  double precision, intent(out) ::  df_dv(n,3)

  !Working
  integer i
  double precision norm, dnorm_dv(3)
  
  do i=1,n
     norm = dot_product(v(:,i), v(:,i))**0.5
     if (norm .lt. 1e-15) then
        norm = 1e-5
     end if
     dnorm_dv(:) = v(:,i) / norm
     !f(i) = - v(3,i) / norm
     df_dv(i,:) = v(3,i) / norm**2 * dnorm_dv(:)
     df_dv(i,3) = df_dv(i,3) - 1.0 / norm
  end do

end subroutine computeJacobianFwdV
