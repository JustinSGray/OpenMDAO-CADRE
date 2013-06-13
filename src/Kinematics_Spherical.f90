subroutine computePositionSpherical(n, v, azimuth, elevation)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, v
  !f2py intent(out) azimuth, elevation
  !f2py depend(n), v, azimuth, elevation

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  v(3,n)

  !Output
  double precision, intent(out) ::  azimuth(n), elevation(n)

  !Working
  integer i
  double precision x, y, z, r
  double precision arctan

  do i=1,n
     x = v(1,i)
     y = v(2,i)
     z = v(3,i)
     r = (x**2 + y**2 + z**2)**0.5
     if (r .lt. 1e-15) then
        r = 1e-5
     end if

     azimuth(i) = arctan(x, y)
     elevation(i) = acos(z/r)
  end do

end subroutine computePositionSpherical



subroutine computePositionSphericalJacobian(n, nJ, v, &
     Ja1, Ji1, Jj1, Ja2, Ji2, Jj2)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, nJ, v
  !f2py intent(out) Ja1, Ji1, Jj1, Ja2, Ji2, Jj2
  !f2py depend(n) v
  !f2py depend(nJ) Ja1, Ji1, Jj1, Ja2, Ji2, Jj2

  !Input
  integer, intent(in) ::  n, nJ
  double precision, intent(in) ::  v(3,n)

  !Output
  double precision, intent(out) ::  Ja1(nJ)
  integer, intent(out) ::  Ji1(nJ), Jj1(nJ)
  double precision, intent(out) ::  Ja2(nJ)
  integer, intent(out) ::  Ji2(nJ), Jj2(nJ)

  !Working
  integer i, k, iJ
  double precision a, e, da_dr(3), de_dr(3), x, y, z, r
  double precision arctan
  
  do i=1,n
     x = v(1,i)
     y = v(2,i)
     z = v(3,i)
     r = (x**2 + y**2 + z**2)**0.5
     if (r .lt. 1e-15) then
        r = 1e-5
     end if

     a = arctan(x, y)
     e = acos(z/r)

     if (e .lt. 1e-15) then
        e = 1e-5
     end if
     if (e .gt. (2*acos(0.0) - 1e-15)) then
        e = 2*acos(0.0) - 1e-5
     end if

     da_dr = 1.0/r * (/-sin(a)/sin(e), cos(a)/sin(e), dble(0) /)
     de_dr = 1.0/r * (/ cos(a)*cos(e), sin(a)*cos(e), -sin(e) /)

     do k=1,3
        iJ = (i-1)*3 + k
        Ja1(iJ) = da_dr(k)
        Ji1(iJ) = i - 1
        Jj1(iJ) = iJ - 1
        Ja2(iJ) = de_dr(k)
        Ji2(iJ) = i - 1
        Jj2(iJ) = iJ - 1
     end do
  end do

end subroutine computePositionSphericalJacobian



function arctan(x, y)

  implicit none

  double precision, intent(in) ::  x, y
  double precision ::  arctan
  double precision pi

  pi = 2*acos(0.0)
  if (x .eq. 0) then
     if (y .gt. 0) then
        arctan = pi/2.0
     else if (y .lt. 0) then
        arctan = 3*pi/2.0
     else
        arctan = 0.0
     end if
  else if (y .eq. 0) then
     if (x .gt. 0) then
        arctan = 0.0
     else if (x .lt. 0) then
        arctan = pi
     else
        arctan = 0.0
     end if
  else if (x .lt. 0) then
     arctan = atan(y/x) + pi
  else if (y .lt. 0) then
     arctan = atan(y/x) + 2*pi
  else if (y .gt. 0) then
     arctan = atan(y/x)
  else
     arctan = 0.0
  end if

end function arctan




subroutine fixAngles(n, azimuth0, elevation0, azimuth, elevation)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, azimuth0, elevation0
  !f2py intent(out) azimuth, elevation
  !f2py depend(n) azimuth0, elevation0, azimuth, elevation

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  azimuth0(n), elevation0(n)

  !Output
  double precision, intent(out) ::  azimuth(n), elevation(n)

  !Working
  integer i
  double precision pi

  pi = 2*acos(0.0)
  do i=1,n
     azimuth(i) = modulo(azimuth0(i), 2*pi)
     elevation(i) = modulo(elevation0(i), 2*pi)
     if (elevation(i) .gt. pi) then
        elevation(i) = 2*pi - elevation(i)
        azimuth(i) = pi + azimuth(i)
        azimuth(i) = modulo(azimuth(i), 2*pi)
     end if
  end do

end subroutine fixAngles
