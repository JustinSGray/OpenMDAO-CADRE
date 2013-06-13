subroutine bracketedSecant(n, C1, C2, I1, I2, I3)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, C1, C2, I1, I2
  !f2py intent(out) I3
  !f2py depend(n) C1, C2, I1, I2, I3

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  C1(12,n), C2(12,n), I1(12,n), I2(12,n)

  !Output
  double precision, intent(out) ::  I3(12,n)

  !Working
  integer i, k

  do k=1,12
     do i=1,n
        if (abs(C1(k,i) - C2(k,i)) .gt. 1e-10) then
           I3(k,i) = I1(k,i) - C1(k,i) * (I2(k,i) - I1(k,i)) / (C2(k,i) - C1(k,i))
           if ((I3(k,i) .lt. min(I1(k,i), I2(k,i))) .or. (I3(k,i) .gt. max(I1(k,i), I2(k,i)))) then
              I3(k,i) = 0.5*I1(k,i) + 0.5*I2(k,i)
           end if
        else
           I3(k,i) = 0.5*I1(k,i) + 0.5*I2(k,i)
        end if
     end do
  end do

end subroutine bracketedSecant




subroutine chooseBracket(n, C1, C2, C3, I3, I1, I2)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, C1, C2, C3, I3
  !f2py intent(inout) I1, I2
  !f2py depend(n) C1, C2, C3, I3, I1, I2

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  C1(12,n), C2(12,n), C3(12,n), I3(12,n)

  !Output
  double precision, intent(inout) ::  I1(12,n), I2(12,n)

  !Working
  integer i, k

  do k=1,12
     do i=1,n
        if (C1(k,i)*C3(k,i) .lt. 0) then
           I2(k,i) = I3(k,i)
        else
           I1(k,i) = I3(k,i)
        end if
     end do
  end do

end subroutine chooseBracket



subroutine getConstants(n, k, q, Voc, V0, AT, Isat)

  implicit none

  !Output
  double precision, intent(out) ::  n, k, q, Voc, V0, AT, Isat

  n = 1.35
  k = 1.38065e-23
  q = 1.60218e-19
  Voc = 2.68/3.0
  V0 = -0.6
  AT = 1.834e-3
  Isat = 2.809e-12

end subroutine getConstants
  



subroutine compute(n, I, LOS, T, A, P_sol)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, I, LOS, T, A
  !f2py intent(out) P_sol
  !f2py depend(n) LOS, T, A, P_sol

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  I(12), LOS(n), T(5,n), A(7,12,n)

  !Output
  double precision, intent(out) ::  P_sol(n)

  !Working
  integer j, c, p, id
  double precision Isc, Vt
  double precision nn, k, q, Voc, V0, AT, Isat
  double precision arg, den

  call getConstants(nn, k, q, Voc, V0, AT, Isat)

  P_sol(:) = 0.0
  do j=1,n
     do p=1,12
        if (p .lt. 5) then
           id = 5
        else
           id = mod(p,4) + 1
        end if
        do c=1,7
           Isc = 0.453*LOS(j)*A(c,p,j)/AT
           Vt = nn*k*T(id,j)/q
           arg = (Isc+Isat-I(p))/Isat
           den = I(p) - Isc - V0*Isat/Vt
           if (I(p) .le. Isc) then
              P_sol(j) = P_sol(j) + I(p)*Vt*log(arg)
           else
              P_sol(j) = P_sol(j) + I(p)*(V0**2*Isat/Vt/den + V0)
           end if
        end do
     end do
  end do

end subroutine compute




subroutine computeDerivatives(n, I, LOS, T, A, dP_dI, dP_dL, dP_dT, dP_dA)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, I, LOS, T, A
  !f2py intent(out) dP_dI, dP_dL, dP_dT, dP_dA
  !f2py depend(n) LOS, T, A, dP_dI, dP_dL, dP_dT, dP_dA

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  I(12), LOS(n), T(5,n), A(7,12,n)

  !Output
  double precision, intent(out) ::  dP_dI(n,12), dP_dL(n), dP_dT(n,5), dP_dA(n,7,12)

  !Working
  integer j, c, p, id
  double precision Isc, Vt, dIsc_dL, dIsc_dA, dVt_dT
  double precision nn, k, q, Voc, V0, AT, Isat
  double precision arg, darg_dI, darg_dL, darg_dA, darg_dT
  double precision den, dden_dI, dden_dL, dden_dA, dden_dT

  call getConstants(nn, k, q, Voc, V0, AT, Isat)

  dP_dI(:,:) = 0.0
  dP_dT(:,:) = 0.0
  dP_dA(:,:,:) = 0.0
  do j=1,n
     do p=1,12
        if (p .lt. 5) then
           id = 5
        else
           id = mod(p,4) + 1
        end if
        do c=1,7
           Isc = 0.453*LOS(j)*A(c,p,j)/AT
           dIsc_dL = 0.453*A(c,p,j)/AT
           dIsc_dA = 0.453*LOS(j)/AT

           Vt = nn*k*T(id,j)/q
           dVt_dT = nn*k/q

           arg = (Isc+Isat-I(p))/Isat
           darg_dI = -1.0/Isat
           darg_dL = dIsc_dL/Isat
           darg_dA = dIsc_dA/Isat
           darg_dT = 0.0

           den = I(p) - Isc - V0*Isat/Vt
           dden_dI = 1.0
           dden_dL = -dIsc_dL
           dden_dA = -dIsc_dA
           dden_dT = V0*Isat/Vt**2*dVt_dT
           if (I(p) .le. Isc) then
              !P_sol(j) = P_sol(j) + I(p)*Vt*log(arg)
              dP_dI(j,p) = dP_dI(j,p) + Vt*log(arg) + I(p)*Vt/arg*darg_dI
              dP_dL(j) = dP_dL(j) + I(p)*Vt/arg*darg_dL
              dP_dA(j,c,p) = dP_dA(j,c,p) + I(p)*Vt/arg*darg_dA
              dP_dT(j,id) = dP_dT(j,id) + I(p)*dVt_dT*log(arg)
           else
              !P_sol(j) = P_sol(j) + I(p)*(V0**2*Isat/Vt/den + V0)
              dP_dI(j,p) = dP_dI(j,p) + V0**2*Isat/Vt/den + V0 - I(p)*V0**2*Isat/Vt/den**2*dden_dI
              dP_dL(j) = dP_dL(j) - I(p)*V0**2*Isat/Vt/den**2*dden_dL
              dP_dT(j,id) = dP_dT(j,id) - I(p)*V0**2*Isat/Vt**2/den*dVt_dT - I(p)*V0**2*Isat/Vt/den**2*dden_dT
              dP_dA(j,c,p) = dP_dA(j,c,p) - I(p)*V0**2*Isat/Vt/den**2*dden_dA
           end if
        end do
     end do
  end do

end subroutine computeDerivatives
