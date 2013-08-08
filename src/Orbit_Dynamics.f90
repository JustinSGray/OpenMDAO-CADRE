subroutine getConstants(C1, C2, C3, C4)

  implicit none

  !Output
  double precision, intent(out) ::  C1, C2, C3, C4

  !Working
  double precision mu, Re, J2, J3, J4

  mu = 398600.44
  Re = 6378.137
  J2 = 1.08264e-3
  J3 = -2.51e-6
  J4 = -1.60e-6

  C1 = -mu
  C2 = -1.5*mu*J2*Re**2
  C3 = -2.5*mu*J3*Re**3
  C4 = 1.875*mu*J4*Re**4

end subroutine getConstants




subroutine getf(dat, u, s)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), u(6)

  !Output
  double precision, intent(out) ::  s(6)

  !Working
  double precision r, x, y, z
  double precision C1, C2, C3, C4
  double precision T2, T3, T4, T3z, T4z

  call getConstants(C1, C2, C3, C4)
    
  x = u(1)
  y = u(2)
  z = u(3)  
  if (abs(z) .lt. 1e-15) then
     z = 1e-5
  end if
  r = dat(1)
  r = (x**2 + y**2 + z**2)**0.5

  T2 = 1 - 5*z**2/r**2
  T3 = 3*z - 7*z**3/r**2
  T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
  T3z = 3*z - 0.6*r**2/z
  T4z = 4 - 28.0/3.0*z**2/r**2

  s(1:3) = u(4:6)
  s(4:6) = u(1:3)*(C1/r**3 + C2/r**5*T2 + C3/r**7*T3 + C4/r**7*T4)
  s(6) = s(6) + z*(C2/r**5*2 + C3/r**7*T3z + C4/r**7*T4z)

end subroutine getf




subroutine getdfdy(dat, u, dfdu)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), u(6)

  !Output
  double precision, intent(out) ::  dfdu(6,6)

  !Working
  double precision r, x, y, z, drdx, drdy, drdz
  double precision C1, C2, C3, C4
  double precision T2, T3, T4, T3z, T4z
  double precision dT2_dx, dT3_dx, dT4_dx, dT3z_dx, dT4z_dx
  double precision dT2_dy, dT3_dy, dT4_dy, dT3z_dy, dT4z_dy
  double precision dT2_dz, dT3_dz, dT4_dz, dT3z_dz, dT4z_dz
  double precision eye(3,3)

  call getConstants(C1, C2, C3, C4)
    
  x = u(1)
  y = u(2)
  z = u(3)
  if (abs(z) .lt. 1e-15) then
     z = 1e-5
  end if
  r = dat(1)
  r = (x**2 + y**2 + z**2)**0.5

  drdx = x/r
  drdy = y/r
  drdz = z/r

  T2 = 1 - 5*z**2/r**2
  T3 = 3*z - 7*z**3/r**2
  T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
  T3z = 3*z - 0.6*r**2/z
  T4z = 4 - 28.0/3.0*z**2/r**2

  dT2_dx = 10*z**2/r**3*drdx
  dT2_dy = 10*z**2/r**3*drdy
  dT2_dz = 10*z**2/r**3*drdz - 10*z/r**2

  dT3_dx = 14*z**3/r**3*drdx
  dT3_dy = 14*z**3/r**3*drdy
  dT3_dz = 14*z**3/r**3*drdz - 21*z**2/r**2 + 3

  dT4_dx = 28*z**2/r**3*drdx - 84*z**4/r**5*drdx
  dT4_dy = 28*z**2/r**3*drdy - 84*z**4/r**5*drdy
  dT4_dz = 28*z**2/r**3*drdz - 84*z**4/r**5*drdz - 28*z/r**2 + 84*z**3/r**4

  dT3z_dx = -1.2*r/z*drdx
  dT3z_dy = -1.2*r/z*drdy
  dT3z_dz = -1.2*r/z*drdz + 0.6*r**2/z**2 + 3

  dT4z_dx = 56.0/3.0*z**2/r**3*drdx
  dT4z_dy = 56.0/3.0*z**2/r**3*drdy
  dT4z_dz = 56.0/3.0*z**2/r**3*drdz - 56.0/3.0*z/r**2 

  eye(:,:) = 0.0
  eye(1,1) = 1.0
  eye(2,2) = 1.0
  eye(3,3) = 1.0

  dfdu(:,:) = 0.0
  dfdu(1:3,4:6) = dfdu(1:3,4:6) + eye

  dfdu(4:6,1:3) = dfdu(4:6,1:3) + eye*(C1/r**3 + C2/r**5*T2 + C3/r**7*T3 + C4/r**7*T4)
  dfdu(4:6,1) = dfdu(4:6,1) + drdx*u(1:3)*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
  dfdu(4:6,2) = dfdu(4:6,2) + drdy*u(1:3)*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
  dfdu(4:6,3) = dfdu(4:6,3) + drdz*u(1:3)*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
  dfdu(4:6,1) = dfdu(4:6,1) + u(1:3)*(C2/r**5*dT2_dx + C3/r**7*dT3_dx + C4/r**7*dT4_dx)
  dfdu(4:6,2) = dfdu(4:6,2) + u(1:3)*(C2/r**5*dT2_dy + C3/r**7*dT3_dy + C4/r**7*dT4_dy)
  dfdu(4:6,3) = dfdu(4:6,3) + u(1:3)*(C2/r**5*dT2_dz + C3/r**7*dT3_dz + C4/r**7*dT4_dz)
  dfdu(6,1) = dfdu(6,1) + drdx*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
  dfdu(6,2) = dfdu(6,2) + drdy*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
  dfdu(6,3) = dfdu(6,3) + drdz*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
  dfdu(6,1) = dfdu(6,1) + z*(C3/r**7*dT3z_dx + C4/r**7*dT4z_dx)
  dfdu(6,2) = dfdu(6,2) + z*(C3/r**7*dT3z_dy + C4/r**7*dT4z_dy)
  dfdu(6,3) = dfdu(6,3) + z*(C3/r**7*dT3z_dz + C4/r**7*dT4z_dz)
  dfdu(6,3) = dfdu(6,3) + (C2/r**5*2 + C3/r**7*T3z + C4/r**7*T4z)

end subroutine getdfdy



subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(1), y(5)

  !Output
  double precision, intent(out) ::  dfdx(5,1)

  dfdx(:,:) = 0.0

end subroutine getdfdx
