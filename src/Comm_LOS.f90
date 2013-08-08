subroutine computeLOS(n, r_b2g_I, r_e2g_I, LOS)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r_b2g_I, r_e2g_I
  !f2py intent(out) LOS
  !f2py depend(n) r_b2g_I, r_e2g_I, LOS

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r_b2g_I(3,n), r_e2g_I(3,n)

  !Output
  double precision, intent(out) ::  LOS(n)

  !Working
  integer i
  double precision Re, Rb, proj, x

  Re = 6378.137
  Rb = 100.0

  do i=1,n
     proj = dot_product(r_b2g_I(:,i), r_e2g_I(:,i)) / Re
     if (proj .gt. 0) then
        LOS(i) = 0.0
     else if (proj .lt. -Rb) then
        LOS(i) = 1.0
     else
        x = (proj - 0) / (-Rb - 0)
        LOS(i) = 3*x**2 - 2*x**3
     end if
  end do

end subroutine computeLOS




subroutine computeJacobianLOS(n, r_b2g_I, r_e2g_I, dLOS_drb, dLOS_dre)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r_b2g_I, r_e2g_I
  !f2py intent(out) dLOS_drb, dLOS_dre
  !f2py depend(n) r_b2g_I, r_e2g_I, dLOS_drb, dLOS_dre

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r_b2g_I(3,n), r_e2g_I(3,n)

  !Output
  double precision, intent(out) ::  dLOS_drb(n,3), dLOS_dre(n,3)

  !Working
  integer i
  double precision Re, Rb, proj, x
  double precision dproj_drb(3), dproj_dre(3), dLOS_dx, dx_dproj

  Re = 6378.137
  Rb = 10.0

  do i=1,n
     proj = dot_product(r_b2g_I(:,i), r_e2g_I(:,i)) / Re
     dproj_drb(:) = r_e2g_I(:,i)
     dproj_dre(:) = r_b2g_I(:,i)
     if (proj .gt. 0) then
        dLOS_drb(i,:) = 0.0
        dLOS_dre(i,:) = 0.0
     else if (proj .lt. -Rb) then
        dLOS_drb(i,:) = 0.0
        dLOS_dre(i,:) = 0.0
     else
        x = (proj - 0) / (-Rb - 0)
        dx_dproj = -1.0 / Rb
        dLOS_dx = 6*x - 6*x**2
        dLOS_drb(i,:) = dLOS_dx * dx_dproj * dproj_drb(:)
        dLOS_dre(i,:) = dLOS_dx * dx_dproj * dproj_dre(:)
     end if
  end do

end subroutine computeJacobianLOS
