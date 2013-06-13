subroutine fixTemps(n, T0, T)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, T0
  !f2py intent(out) T
  !f2py depend(n) T0, T

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  T0(5,n)

  !Output
  double precision, intent(out) ::  T(5,n)

  !Working
  integer i, k

  do i=1,n
     do k=1,5
        T(k,i) = T0(k,i)
        if (T(k,i) .lt. 0) then
           T(k,i) = 0
        end if           
     end do
  end do

end subroutine fixTemps



subroutine getConstants(m_f, cp_f, m_b, cp_b, A_T, &
     alpha_c, alpha_r, eps_c, eps_r, q_sol, K)

  implicit none

  !Output
  double precision, intent(out) ::   m_f, cp_f, m_b, cp_b, A_T 
  double precision, intent(out) ::   alpha_c, alpha_r, eps_c, eps_r
  double precision, intent(out) ::   q_sol, K

  m_f = 0.4
  m_b = 2.0
  cp_f = 0.6e3
  cp_b = 2.0e3
  A_T = 2.66e-3

  alpha_c = 0.9
  alpha_r = 0.2
  eps_c = 0.87
  eps_r = 0.88

  q_sol = 1360.0
  K = 5.6704e-8

end subroutine getConstants



subroutine getf(dat, y, f)

  implicit none

  !Input
  double precision, intent(in) ::  dat(170), y(5)

  !Output
  double precision, intent(out) ::  f(5)

  !Working
  double precision m_f, cp_f, m_b, cp_b, A_T 
  double precision alpha_c, alpha_r, eps_c, eps_r
  double precision q_sol, K
  double precision m, cp, alpha, eps
  double precision A_exp, w, LOS, Pcomm
  integer c, p, id, idat

  call getConstants(m_f, cp_f, m_b, cp_b, A_T, &
     alpha_c, alpha_r, eps_c, eps_r, q_sol, K)

  f(:) = 0.0
  LOS = dat(169)
  Pcomm = dat(170)
  do p=1,12
     if (p .lt. 5) then
        id = 5
        m = m_b
        cp = cp_b
     else
        id = mod(p,4) + 1
        m = m_f
        cp = cp_f
     end if
     do c=1,7
        idat = (p-1)*7 + c
        A_exp = dat(idat)
        w = dat(84 + idat)
        alpha = alpha_c*w + alpha_r*(1-w)
        eps = eps_c*w + eps_r*(1-w)
        f(id) = f(id) + alpha * q_sol * A_exp * LOS / m / cp
        f(id) = f(id) - eps * k * A_T * y(id)**4 / m / cp
     end do
  end do
  f(5) = f(5) + 4.0 * Pcomm / m_b / cp_b

end subroutine getf




subroutine getdfdy(dat, y, dfdy)

  implicit none

  !Input
  double precision, intent(in) ::  dat(170), y(5)

  !Output
  double precision, intent(out) ::  dfdy(5,5)

  !Working
  double precision m_f, cp_f, m_b, cp_b, A_T 
  double precision alpha_c, alpha_r, eps_c, eps_r
  double precision q_sol, K
  double precision m, cp, alpha, eps
  double precision A_exp, w, LOS
  integer c, p, id, idat

  call getConstants(m_f, cp_f, m_b, cp_b, A_T, &
     alpha_c, alpha_r, eps_c, eps_r, q_sol, K)

  dfdy(:,:) = 0.0
  do p=1,12
     if (p .lt. 5) then
        id = 5
        m = m_b
        cp = cp_b
     else
        id = mod(p,4) + 1
        m = m_f
        cp = cp_f
     end if
     do c=1,7
        idat = (p-1)*7 + c
        A_exp = dat(idat)
        w = dat(84 + idat)
        LOS = dat(169)
        alpha = alpha_c*w + alpha_r*(1-w)
        eps = eps_c*w + eps_r*(1-w)
        dfdy(id,id) = dfdy(id,id) - 4 * eps * k * A_T * y(id)**3 / m / cp
     end do
  end do

end subroutine getdfdy




subroutine getdfdx(dat, y, dfdx)

  implicit none

  !Input
  double precision, intent(in) ::  dat(170), y(5)

  !Output
  double precision, intent(out) ::  dfdx(5,170)

  !Working
  double precision m_f, cp_f, m_b, cp_b, A_T 
  double precision alpha_c, alpha_r, eps_c, eps_r
  double precision q_sol, K
  double precision m, cp, alpha, eps
  double precision A_exp, w, LOS, Pcomm
  double precision dalpha_dw, deps_dw
  integer c, p, id, idat, iA, iw

  call getConstants(m_f, cp_f, m_b, cp_b, A_T, &
     alpha_c, alpha_r, eps_c, eps_r, q_sol, K)

  dfdx(:,:) = 0.0
  LOS = dat(169)
  Pcomm = dat(170)
  do p=1,12
     if (p .lt. 5) then
        id = 5
        m = m_b
        cp = cp_b
     else
        id = mod(p,4) + 1
        m = m_f
        cp = cp_f
     end if
     do c=1,7
        idat = (p-1)*7 + c
        iA = idat
        iw = 84 + idat
        A_exp = dat(iA)
        w = dat(iw)
        alpha = alpha_c*w + alpha_r*(1-w)
        eps = eps_c*w + eps_r*(1-w)

        dalpha_dw = alpha_c - alpha_r
        deps_dw = eps_c - eps_r

        !f(id) = f(id) + alpha * q_sol * A_exp * LOS / m / cp
        !f(id) = f(id) - eps * k * A_T * y(id)**4 / m / cp
        dfdx(id,iA) = dfdx(id,iA) + alpha * q_sol * LOS / m / cp
        dfdx(id,iw) = dfdx(id,iw) + dalpha_dw * q_sol * A_exp * LOS / m / cp
        dfdx(id,iw) = dfdx(id,iw) - deps_dw * k * A_T * y(id)**4 / m / cp
        dfdx(id,169) = dfdx(id,169) + alpha * q_sol * A_exp / m / cp
     end do
  end do
  dfdx(5,170) = dfdx(5,170) + 4.0 / m_b / cp_b

end subroutine getdfdx
