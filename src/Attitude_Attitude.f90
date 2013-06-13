subroutine computeORI(n, r_e2b_I, O_RI)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r_e2b_I
  !f2py intent(out) O_RI
  !f2py depend(n) r_e2b_I, O_RI

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r_e2b_I(6,n)

  !Output
  double precision, intent(out) ::  O_RI(3,3,n)

  !Working
  integer i
  double precision normr, normv, r(3), v(3), vx(3,3)
  double precision iB(3), jB(3)

  O_RI(:,:,:) = 0.0
  do i=1,n
     r = r_e2b_I(1:3,i)
     v = r_e2b_I(4:6,i)

     normr = sqrt(dot_product(r,r))
     normv = sqrt(dot_product(v,v))
     if (normr .lt. 1e-10) then
        normr = 1e-10
     end if
     if (normv .lt. 1e-10) then
        normv = 1e-10
     end if

     r = r / normr
     v = v / normv

     vx(1,:) = (/ dble(0),  -v(3) ,   v(2)  /)
     vx(2,:) = (/   v(3) , dble(0),  -v(1)  /)
     vx(3,:) = (/  -v(2) ,   v(1) , dble(0) /)

     iB =  matmul(vx, r)
     jB = -matmul(vx, iB)

     O_RI(1,:,i) = iB
     O_RI(2,:,i) = jB
     O_RI(3,:,i) = -v
  end do  

end subroutine computeORI



subroutine computeJacobianORI(n, r_e2b_I, dO_dr)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, r_e2b_I
  !f2py intent(out) dO_dr
  !f2py depend(n) r_e2b_I, dO_dr

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  r_e2b_I(6,n)

  !Output
  double precision, intent(out) ::  dO_dr(n,3,3,6)

  !Working
  integer i, k
  double precision normr, normv, r(3), v(3), vx(3,3), dvx_dv(3,3,3)
  double precision iB(3), jB(3)
  double precision diB_dr(3,3), diB_dv(3,3)
  double precision djB_diB(3,3), djB_dv(3,3)
  double precision dr_dr(3,3), dv_dv(3,3)

  dO_dr(:,:,:,:) = 0.0
  do i=1,n
     r = r_e2b_I(1:3,i)
     v = r_e2b_I(4:6,i)

     normr = sqrt(dot_product(r,r))
     normv = sqrt(dot_product(v,v))
     if (normr .lt. 1e-10) then
        normr = 1e-10
     end if
     if (normv .lt. 1e-10) then
        normv = 1e-10
     end if

     r = r / normr     
     v = v / normv
     dr_dr(:,:) = 0.0
     dv_dv(:,:) = 0.0
     do k=1,3
        dr_dr(k,k) = dr_dr(k,k) + 1.0 / normr
        dv_dv(k,k) = dv_dv(k,k) + 1.0 / normv
        dr_dr(:,k) = dr_dr(:,k) - r_e2b_I(1:3,i) / normr**2 * r_e2b_I(k,i) / normr
        dv_dv(:,k) = dv_dv(:,k) - r_e2b_I(4:6,i) / normv**2 * r_e2b_I(3+k,i)/ normv
     end do     

     vx(1,:) = (/ dble(0),  -v(3) ,   v(2)  /)
     vx(2,:) = (/   v(3) , dble(0),  -v(1)  /)
     vx(3,:) = (/  -v(2) ,   v(1) , dble(0) /)

     dvx_dv(1,:,1) = (/  0,  0,  0 /)
     dvx_dv(2,:,1) = (/  0,  0, -1 /)
     dvx_dv(3,:,1) = (/  0,  1,  0 /)
     
     dvx_dv(1,:,2) = (/  0,  0,  1 /)
     dvx_dv(2,:,2) = (/  0,  0,  0 /)
     dvx_dv(3,:,2) = (/ -1,  0,  0 /)
     
     dvx_dv(1,:,3) = (/  0, -1,  0 /)
     dvx_dv(2,:,3) = (/  1,  0,  0 /)
     dvx_dv(3,:,3) = (/  0,  0,  0 /)

     iB =  matmul(vx, r)
     jB = -matmul(vx, iB)

     diB_dr = vx
     diB_dv(:,1) = matmul(dvx_dv(:,:,1),r)
     diB_dv(:,2) = matmul(dvx_dv(:,:,2),r)
     diB_dv(:,3) = matmul(dvx_dv(:,:,3),r)

     djB_diB = -vx
     djB_dv(:,1) = -matmul(dvx_dv(:,:,1),iB)
     djB_dv(:,2) = -matmul(dvx_dv(:,:,2),iB)
     djB_dv(:,3) = -matmul(dvx_dv(:,:,3),iB)

     !O_RI(1,:,i) = iB
     dO_dr(i,1,:,1:3) = matmul(diB_dr , dr_dr)
     dO_dr(i,1,:,4:6) = matmul(diB_dv , dv_dv)

     !O_RI(2,:,i) = jB
     dO_dr(i,2,:,1:3) = matmul(matmul(djB_diB, diB_dr) , dr_dr)
     dO_dr(i,2,:,4:6) = matmul(matmul(djB_diB, diB_dv) + djB_dv , dv_dv)

     !O_RI(3,:,i) = -v
     dO_dr(i,3,:,4:6) = -dv_dv
  end do

end subroutine computeJacobianORI
