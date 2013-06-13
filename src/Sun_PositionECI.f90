subroutine computePosition(n, t, r_e2s_I)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, t
  !f2py intent(out) r_e2s_I
  !f2py depend(n) t
  !f2py depend(n) r_e2s_I

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  t(n)

  !Output
  double precision, intent(out) ::  r_e2s_I(3,n)

  !Working
  integer i
  double precision pi, d2r, L, g, lambda, eps

  pi = 2*acos(0.0)
  d2r = pi/180.0
  do i=1,n
     L = d2r*280.460 + d2r*0.9856474*t(i)
     g = d2r*357.528 + d2r*0.9856003*t(i)
     lambda = L + d2r*1.914666*sin(g) + d2r*0.01999464*sin(2*g)
     eps = d2r*23.439 - d2r*3.56e-7*t(i)
     r_e2s_I(1,i) = cos(lambda)
     r_e2s_I(2,i) = sin(lambda)*cos(eps)
     r_e2s_I(3,i) = sin(lambda)*sin(eps)
  end do

end subroutine computePosition



subroutine computePositionJacobian(n, nJ, t, Ja, Ji, Jj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, nJ, t
  !f2py intent(out) Ja, Ji, Jj
  !f2py depend(n) t
  !f2py depend(nJ) Ja, Ji, Jj

  !Input
  integer, intent(in) ::  n, nJ
  double precision, intent(in) ::  t(n)

  !Output
  double precision, intent(out) ::  Ja(nJ)
  integer, intent(out) ::  Ji(nJ), Jj(nJ)

  !Working
  integer k, i, iJ
  double precision pi, d2r
  double precision L, g, lambda, eps
  double precision dL_dt, dg_dt, dlambda_dt, deps_dt
  double precision dr_dt(3)

  pi = 2*acos(0.0)
  d2r = pi/180.0
  do i=1,n
     L = d2r*280.460 + d2r*0.9856474*t(i)
     g = d2r*357.528 + d2r*0.9856003*t(i)
     lambda = L + d2r*1.914666*sin(g) + d2r*0.01999464*sin(2*g)
     eps = d2r*23.439 - d2r*3.56e-7*t(i)

     dL_dt = d2r*0.9856474
     dg_dt = d2r*0.9856003
     dlambda_dt = dL_dt + d2r*1.914666*cos(g)*dg_dt + d2r*0.01999464*cos(2*g)*2*dg_dt
     deps_dt = -d2r*3.56e-7
     
     dr_dt(1) = -sin(lambda)*dlambda_dt
     dr_dt(2) = cos(lambda)*cos(eps)*dlambda_dt - sin(lambda)*sin(eps)*deps_dt
     dr_dt(3) = cos(lambda)*sin(eps)*dlambda_dt + sin(lambda)*cos(eps)*deps_dt

     do k=1,3
        iJ = (i-1)*3 + k
        Ja(iJ) = dr_dt(k)
        Ji(iJ) = iJ - 1
        Jj(iJ) = i - 1
     end do
  end do

end subroutine computePositionJacobian
