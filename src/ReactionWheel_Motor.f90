subroutine computeT(n, T_RW, w_B, h_RW, T_m)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, T_RW, w_B, h_RW
  !f2py intent(out) T_m
  !f2py depend(n) T_RW, w_B, h_RW, T_m

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  T_RW(3,n), w_B(3,n), h_RW(3,n)

  !Output
  double precision, intent(out) ::  T_m(3,n)

  !Working
  integer i
  double precision w_Bx(3,3)

  do i=1,n
     w_Bx(1,:) = (/   dble(0) , -w_B(3,i) ,  w_B(2,i) /)
     w_Bx(2,:) = (/  w_B(3,i) ,   dble(0) , -w_B(1,i) /)
     w_Bx(3,:) = (/ -w_B(2,i) ,  w_B(1,i) ,   dble(0) /)

     T_m(:,i) = -T_RW(:,i) - matmul(w_Bx , h_RW(:,i))
  end do

end subroutine computeT



subroutine computeJacobianT(n, T_RW, w_B, h_RW, dT_dTm, dT_dwb, dT_dh)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, T_RW, w_B, h_RW
  !f2py intent(out) dT_dTm, dT_dwb, dT_dh
  !f2py depend(n) T_RW, w_B, h_RW, dT_dTm, dT_dwb, dT_dh

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  T_RW(3,n), w_B(3,n), h_RW(3,n)

  !Output
  double precision, intent(out) ::  dT_dTm(n,3,3), dT_dwb(n,3,3), dT_dh(n,3,3)

  !Working
  integer i, k
  double precision w_Bx(3,3), dwx_dwb(3,3,3)

  dT_dTm(:,:,:) = 0.0
  dT_dwb(:,:,:) = 0.0
  dT_dh(:,:,:) = 0.0
  do i=1,n
     w_Bx(1,:) = (/   dble(0) , -w_B(3,i) ,  w_B(2,i) /)
     w_Bx(2,:) = (/  w_B(3,i) ,   dble(0) , -w_B(1,i) /)
     w_Bx(3,:) = (/ -w_B(2,i) ,  w_B(1,i) ,   dble(0) /)

     dwx_dwb(1,:,1) = (/  0 ,  0 ,  0 /)
     dwx_dwb(2,:,1) = (/  0 ,  0 , -1 /)
     dwx_dwb(3,:,1) = (/  0 ,  1 ,  0 /)

     dwx_dwb(1,:,2) = (/  0 ,  0 ,  1 /)
     dwx_dwb(2,:,2) = (/  0 ,  0 ,  0 /)
     dwx_dwb(3,:,2) = (/ -1 ,  0 ,  0 /)

     dwx_dwb(1,:,3) = (/  0 , -1 ,  0 /)
     dwx_dwb(2,:,3) = (/  1 ,  0 ,  0 /)
     dwx_dwb(3,:,3) = (/  0 ,  0 ,  0 /)

     !T_m(:,i) = -T_RW(:,i) - matmul(w_Bx , h_RW(:,i))
     do k=1,3
        dT_dTm(i,k,k) = - 1.0
        dT_dwb(i,:,k) = - matmul(dwx_dwb(:,:,k) , h_RW(:,i))
     end do
     dT_dh(i,:,:) = - w_Bx
  end do

end subroutine computeJacobianT
