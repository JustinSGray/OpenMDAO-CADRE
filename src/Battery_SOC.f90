subroutine getConstants(sigma, eta, Cp, IR, T0, alpha)

  double precision, intent(out) ::  sigma, eta, Cp, IR, T0, alpha

  sigma = 1e-10
  eta = 0.99
  Cp = 2900*0.001*3600
  IR = 0.9
  T0 = 293
  alpha = log(1/1.1**5)

end subroutine getConstants




subroutine getVoc(SOC, Voc, dVoc_dSOC)

  implicit none

  double precision, intent(in) ::  SOC
  double precision, intent(out) ::  Voc, dVoc_dSOC
  double precision a(4)
  
  a = (/ 3.8360, -2.0488, 1.9938, -1.7244 /)
  
  !Voc = a(1) + a(2)*SOC + a(3)*SOC**2 + a(4)*SOC**3
  !dVoc_dSOC = a(2) + 2*a(3)*SOC + 3*a(4)*SOC**2

  !Voc = 3 + SOC
  !dVoc_dSOC = 1.0

  Voc = 3 + (exp(SOC)-1) / (exp(1.0)-1)
  dVoc_dSOC = exp(SOC) / (exp(1.0)-1)

end subroutine getVoc




subroutine computeI(n, P, T, SOC, I)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, P, T, SOC
  !f2py intent(out) I
  !f2py depend(n) P, T, SOC, I

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P(n), T(n), SOC(n)

  !Output
  double precision, intent(out) ::  I(n)

  !Working
  double precision sigma, eta, Cp, IR, T0, alpha, a(4)
  double precision Voc, dVoc_dSOC, V
  integer k

  call getConstants(sigma, eta, Cp, IR, T0, alpha)

  do k=1,n
     call getVoc(SOC(k), Voc, dVoc_dSOC)
     V = IR * Voc * (2 - exp(alpha*(T(k)-T0)/T0))
     I(k) = P(k) / V
  end do

end subroutine computeI




subroutine computeJacobianI(n, P, T, SOC, dI_dP, dI_dT, dI_dSOC)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, P, T, SOC
  !f2py intent(out) dI_dP, dI_dT, dI_dSOC
  !f2py depend(n) P, T, SOC, dI_dP, dI_dT, dI_dSOC

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P(n), T(n), SOC(n)

  !Output
  double precision, intent(out) ::  dI_dP(n), dI_dT(n), dI_dSOC(n)

  !Working
  double precision sigma, eta, Cp, IR, T0, alpha, a(4)
  double precision Voc, V, dVoc_dSOC, dV_dVoc, dV_dT
  integer k

  call getConstants(sigma, eta, Cp, IR, T0, alpha)

  do k=1,n
     call getVoc(SOC(k), Voc, dVoc_dSOC)

     V = IR * Voc * (2 - exp(alpha*(T(k)-T0)/T0))
     dV_dVoc = IR * (2 - exp(alpha*(T(k)-T0)/T0))
     dV_dT = - IR * Voc * exp(alpha*(T(k)-T0)/T0) * alpha/T0

     dI_dP(k) = 1.0 / V
     dI_dT(k) = -P(k)/V**2 * dV_dT
     dI_dSOC(k) = -P(k)/V**2 * dV_dVoc * dVoc_dSOC
  end do

end subroutine computeJacobianI




subroutine getf(dat, y, f)

  implicit none

  !Input
  double precision, intent(in) ::  dat(2), y(1)

  !Output
  double precision, intent(out) ::  f(1)

  !Working
  double precision sigma, eta, Cp, IR, T0, alpha, a(4)
  double precision SOC, P, T, Voc, dVoc_dSOC, V, I

  call getConstants(sigma, eta, Cp, IR, T0, alpha)
  
  SOC = y(1)
  P = dat(1)
  T = dat(2)
  call getVoc(SOC, Voc, dVoc_dSOC)
  V = IR * Voc * (2 - exp(alpha*(T-T0)/T0))
  I = P/V

  f(1) = -sigma/24*SOC + eta/Cp*I

end subroutine getf




subroutine getdfdy(dat, y, dfdy)

  implicit none

  !Input
  double precision, intent(in) ::  dat(2), y(1)

  !Output
  double precision, intent(out) ::  dfdy(1,1)

  !Working
  double precision sigma, eta, Cp, IR, T0, alpha, a(4)
  double precision SOC, P, T, Voc, V, I
  double precision dVoc_dSOC, dV_dSOC, dI_dSOC

  call getConstants(sigma, eta, Cp, IR, T0, alpha)

  SOC = y(1)
  P = dat(1)
  T = dat(2)
  call getVoc(SOC, Voc, dVoc_dSOC)
  V = IR * Voc * (2 - exp(alpha*(T-T0)/T0))
  I = P/V

  dV_dSOC = IR * dVoc_dSOC * (2 - exp(alpha*(T-T0)/T0))
  dI_dSOC = -P/V**2 * dV_dSOC

  dfdy(1,1) = -sigma/24 + eta/Cp*dI_dSOC

end subroutine getdfdy



subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(2), y(1)

  !Output
  double precision, intent(out) ::  dfdx(1,2)

  !Working
  double precision sigma, eta, Cp, IR, T0, alpha, a(4)
  double precision SOC, P, T, Voc, dVoc_dSOC, V, I
  double precision dV_dT, dI_dT, dI_dP

  call getConstants(sigma, eta, Cp, IR, T0, alpha)
  
  SOC = y(1)
  P = dat(1)
  T = dat(2)
  call getVoc(SOC, Voc, dVoc_dSOC)
  V = IR * Voc * (2 - exp(alpha*(T-T0)/T0))
  I = P/V

  dV_dT = - IR * Voc * exp(alpha*(T-T0)/T0) * alpha/T0
  dI_dT = - P/V**2 * dV_dT
  dI_dP = 1.0/V

  dfdx(1,1) = eta/Cp*dI_dP
  dfdx(1,2) = eta/Cp*dI_dT

end subroutine getdfdx

