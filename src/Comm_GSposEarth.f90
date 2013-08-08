subroutine computeGS(n, lon, lat, alt, r_e2g_E)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, lon, lat, alt
  !f2py intent(out) r_e2g_E
  !f2py depend(n) r_e2g_E

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  lon, lat, alt

  !Output
  double precision, intent(out) ::  r_e2g_E(3,n)

  !Working
  double precision Re, pi, d2r

  Re = 6378.137
  pi = 2*acos(0.0)
  d2r = pi/180.0

  r_e2g_E(1,:) = (Re + alt) * cos(d2r*lat) * cos(d2r*lon)
  r_e2g_E(2,:) = (Re + alt) * cos(d2r*lat) * sin(d2r*lon)
  r_e2g_E(3,:) = (Re + alt) * sin(d2r*lat)

end subroutine computeGS



subroutine computeJacobianGS(lon, lat, alt, dr_dlon, dr_dlat, dr_dalt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) lon, lat, alt
  !f2py intent(out) dr_dlon, dr_dlat, dr_dalt

  !Input
  double precision, intent(in) ::  lon, lat, alt

  !Output
  double precision, intent(out) ::  dr_dlon(3), dr_dlat(3), dr_dalt(3)

  !Working
  double precision Re, pi, d2r

  Re = 6378.137
  pi = 2*acos(0.0)
  d2r = pi/180.0

  !r_e2g_E(1,:) = (Re + alt) * cos(d2r*lat) * cos(d2r*lon)
  dr_dlon(1) = - d2r * (Re + alt) * cos(d2r*lat) * sin(d2r*lon)
  dr_dlat(1) = - d2r * (Re + alt) * sin(d2r*lat) * cos(d2r*lon)
  dr_dalt(1) = cos(d2r*lat) * cos(d2r*lon)

  !r_e2g_E(2,:) = (Re + alt) * cos(d2r*lat) * sin(d2r*lon)
  dr_dlon(2) =   d2r * (Re + alt) * cos(d2r*lat) * cos(d2r*lon)
  dr_dlat(2) = - d2r * (Re + alt) * sin(d2r*lat) * sin(d2r*lon)
  dr_dalt(2) = cos(d2r*lat) * sin(d2r*lon)

  !r_e2g_E(3,:) = (Re + alt) * sin(d2r*lat)
  dr_dlon(3) = 0.0
  dr_dlat(3) = d2r * (Re + alt) * cos(d2r*lat)
  dr_dalt(3) = sin(d2r*lat)

end subroutine computeJacobianGS
