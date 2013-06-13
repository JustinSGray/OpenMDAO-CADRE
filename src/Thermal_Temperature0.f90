subroutine getDat(nc, np, n, ny, m_f, cp_f, m_b, cp_b, A_T, &
     alpha_c, alpha_r, eps_c, eps_r, LOS, cellInstd, A_exp, Dat)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nc, np, n, ny, m_f, cp_f, m_b, cp_b, A_T, alpha_c, alpha_r, eps_c, eps_r, LOS, cellInstd, A_exp
  !f2py intent(out) Dat
  !f2py depend(n) LOS
  !f2py depend(nc,np) cellInstd
  !f2py depend(nc,np,n) A_exp
  !f2py depend(ny) Dat

  !Input
  integer, intent(in) ::  nc, np, n, ny
  double precision, intent(in) ::  m_f, cp_f, m_b, cp_b, A_T 
  double precision, intent(in) ::  alpha_c, alpha_r, eps_c, eps_r
  double precision, intent(in) ::  LOS(n), cellInstd(nc,np), A_exp(nc,np,n)

  !Output
  double precision, intent(out) ::  Dat(ny,6)

  !Working
  double precision m, cp, w, alpha, eps, C1, C2
  double precision K, q_sol
  integer c, p, i, id, iy

  K = 5.67051e-8
  q_sol = 1360.0

  Dat(:,:) = 0.0
  do p=1,np
     if (p .lt. 5) then
        id = 5
        m = m_b
        cp = cp_b
     else
        id = mod(p,4) + 1
        m = m_f
        cp = cp_f
     end if
     do c=1,nc
        w = cellInstd(c,p)
        alpha = alpha_c*w + alpha_r*(1-w)
        eps = eps_c*w + eps_r*(1-w)
        C1 = q_sol*alpha/m/cp
        C2 = k*eps/m/cp
        do i=1,n
           iy = 5*(i-1) + id
           Dat(iy,6) = Dat(iy,6) + C1*A_exp(c,p,i)*LOS(i)
           Dat(iy,id) = Dat(iy,id) - C2*A_T
        end do
     end do
  end do

end subroutine getDat




subroutine getf(dat, u, s)

  implicit none

  !Input
  double precision, intent(in) ::  dat(5,6), u(5)

  !Output
  double precision, intent(out) ::  s(5)

  !Working
  integer i

  do i=1,5
     s(i) = dat(i,6) + dat(i,i)*u(i)**4
  end do

end subroutine getf




subroutine getdfdy(dat, u, dfdu)

  implicit none

  !Input
  double precision, intent(in) ::  dat(5,6), u(5)

  !Output
  double precision, intent(out) ::  dfdu(5,5)

  !Working
  integer i

  dfdu(:,:) = 0.0
  do i=1,5
     dfdu(i,i) = 4*dat(i,i)*u(i)**3
  end do

end subroutine getdfdy
