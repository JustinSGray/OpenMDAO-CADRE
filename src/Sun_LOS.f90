subroutine computeLOS(n, r1, r2, r_e2b_I, r_e2s_I, LOS)

  implicit none

  !Fortran-pyton interface directives
  !f2py intent(in) n, r1, r2, r_e2b_I, r_e2s_I
  !f2py intent(out) LOS
  !f2py depend(n) r_e2b_I, r_e2s_I, LOS

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r1, r2
  double precision, intent(in) ::  r_e2b_I(6,n), r_e2s_I(3,n)

  !Output
  double precision, intent(out) ::  LOS(n)

  !Working
  integer i
  double precision r_b(3), r_s(3), dot, dist, Bx(3,3), cross(3), x

  do i=1,n
     r_b = r_e2b_I(:3,i)
     r_s = r_e2s_I(:3,i)
     call crossMatrix(r_b, Bx)
     dot = dot_product(r_b,r_s)
     cross = matmul(Bx,r_s)
     dist = dot_product(cross,cross)**0.5
     if (dot .ge. 0) then
        LOS(i) = 1.0
     else if (dist .le. r1) then
        LOS(i) = 0.0
     else if (dist .ge. r2) then
        LOS(i) = 1.0
     else
        x = (dist-r1)/(r2-r1)
        LOS(i) = 3*x**2 - 2*x**3
     end if
  end do

end subroutine computeLOS



subroutine computeLOSJacobian(n, nJ, r1, r2, r_e2b_I, r_e2s_I, &
     Jab, Jib, Jjb, Jas, Jis, Jjs)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, nJ, r1, r2, r_e2b_I, r_e2s_I
  !f2py intent(out) Jab, Jib, Jjb, Jas, Jis, Jjs
  !f2py depend(n) r_e2b_I, r_e2s_I
  !f2py depend(nJ) Jab, Jib, Jjb, Jas, Jis, Jjs

  !Input
  integer, intent(in) ::  n, nJ
  double precision, intent(in) ::  r1, r2
  double precision, intent(in) ::  r_e2b_I(6,n), r_e2s_I(3,n)

  !Output
  double precision, intent(out) ::  Jab(nJ)
  integer, intent(out) ::  Jib(nJ), Jjb(nJ)
  double precision, intent(out) ::  Jas(nJ)
  integer, intent(out) ::  Jis(nJ), Jjs(nJ)

  !Working
  integer i, k, iJ
  double precision r_b(3), r_s(3), dot, dist, Bx(3,3), Sx(3,3), cross(3), x
  double precision dx_ddist, ddist_dcross(3), dcross_drb(3,3), dcross_drs(3,3)
  double precision LOS, dLOS_dx, dLOS_drb(3), dLOS_drs(3)

  do i=1,n
     r_b = r_e2b_I(:3,i)
     r_s = r_e2s_I(:3,i)
     call crossMatrix( r_b, Bx)
     call crossMatrix(-r_s, Sx)
     dot = dot_product(r_b,r_s)
     cross = matmul(Bx,r_s)
     dist = dot_product(cross,cross)**0.5
     if (dot .ge. 0) then
        dLOS_drb(:) = 0.0
        dLOS_drs(:) = 0.0        
     else if (dist .le. r1) then
        dLOS_drb(:) = 0.0
        dLOS_drs(:) = 0.0      
     else if (dist .ge. r2) then
        dLOS_drb(:) = 0.0
        dLOS_drs(:) = 0.0  
     else
        x = (dist-r1)/(r2-r1)
        LOS = 3*x**2 - 2*x**3
        ddist_dcross = cross/dist
        dcross_drb = Sx
        dcross_drs = Bx
        dx_ddist = 1.0/(r2-r1)
        dLOS_dx = 6*x - 6*x**2
        dLOS_drb = dLOS_dx*dx_ddist*matmul(ddist_dcross,dcross_drb)
        dLOS_drs = dLOS_dx*dx_ddist*matmul(ddist_dcross,dcross_drs)
     end if
     do k=1,3
        iJ = (i-1)*3 + k
        Jab(iJ) = dLOS_drb(k)
        Jib(iJ) = i - 1
        Jjb(iJ) = (i-1)*6 + k - 1
        Jas(iJ) = dLOS_drs(k)
        Jis(iJ) = i - 1
        Jjs(iJ) = (i-1)*3 + k - 1
     end do
  end do

end subroutine computeLOSJacobian



subroutine crossMatrix(v, M)

  implicit none

  !Input
  double precision, intent(in) ::  v(3)

  !Output
  double precision, intent(out) ::  M(3,3)

  M(1,:) = (/ dble(0), -v(3), v(2) /)
  M(2,:) = (/ v(3), dble(0), -v(1) /)
  M(3,:) = (/ -v(2), v(1), dble(0) /)

end subroutine crossMatrix
