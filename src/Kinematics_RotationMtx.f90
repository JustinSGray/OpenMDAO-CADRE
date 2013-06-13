subroutine computeRotationMtx(n, q, O_21)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, q
  !f2py intent(out) O_21
  !f2py depend(n) q
  !f2py depend(n) O_21

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  q(4,n)

  !Output
  double precision, intent(out) ::  O_21(3,3,n)

  !Working
  integer i
  double precision A(4,3), B(4,3)

  do i=1,n
     A(1,:) = (/ q(1,i),-q(4,i), q(3,i) /)
     A(2,:) = (/ q(4,i), q(1,i),-q(2,i) /)
     A(3,:) = (/-q(3,i), q(2,i), q(1,i) /)
     A(4,:) = (/ q(2,i), q(3,i), q(4,i) /)

     B(1,:) = (/ q(1,i), q(4,i),-q(3,i) /)
     B(2,:) = (/-q(4,i), q(1,i), q(2,i) /)
     B(3,:) = (/ q(3,i),-q(2,i), q(1,i) /)
     B(4,:) = (/ q(2,i), q(3,i), q(4,i) /)

     O_21(:,:,i) = matmul(transpose(A),B)
  end do

end subroutine computeRotationMtx



subroutine computeRotMtxJacobian(n, q, J)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, q
  !f2py intent(out) J
  !f2py depend(n) q, J

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  q(4,n)

  !Output
  double precision, intent(out) ::  J(n,3,3,4)

  !Working
  integer i, k
  double precision A(4,3), B(4,3), dA_dq(4,3,4), dB_dq(4,3,4)

  dA_dq(1,:,1) = (/ 1, 0, 0 /)
  dA_dq(2,:,1) = (/ 0, 1, 0 /)
  dA_dq(3,:,1) = (/ 0, 0, 1 /)
  dA_dq(4,:,1) = (/ 0, 0, 0 /)
  
  dA_dq(1,:,2) = (/ 0, 0, 0 /)
  dA_dq(2,:,2) = (/ 0, 0,-1 /)
  dA_dq(3,:,2) = (/ 0, 1, 0 /)
  dA_dq(4,:,2) = (/ 1, 0, 0 /)
  
  dA_dq(1,:,3) = (/ 0, 0, 1 /)
  dA_dq(2,:,3) = (/ 0, 0, 0 /)
  dA_dq(3,:,3) = (/-1, 0, 0 /)
  dA_dq(4,:,3) = (/ 0, 1, 0 /)
  
  dA_dq(1,:,4) = (/ 0,-1, 0 /)
  dA_dq(2,:,4) = (/ 1, 0, 0 /)
  dA_dq(3,:,4) = (/ 0, 0, 0 /)
  dA_dq(4,:,4) = (/ 0, 0, 1 /)


  dB_dq(1,:,1) = (/ 1, 0, 0 /)
  dB_dq(2,:,1) = (/ 0, 1, 0 /)
  dB_dq(3,:,1) = (/ 0, 0, 1 /)
  dB_dq(4,:,1) = (/ 0, 0, 0 /)
  
  dB_dq(1,:,2) = (/ 0, 0, 0 /)
  dB_dq(2,:,2) = (/ 0, 0, 1 /)
  dB_dq(3,:,2) = (/ 0,-1, 0 /)
  dB_dq(4,:,2) = (/ 1, 0, 0 /)
  
  dB_dq(1,:,3) = (/ 0, 0,-1 /)
  dB_dq(2,:,3) = (/ 0, 0, 0 /)
  dB_dq(3,:,3) = (/ 1, 0, 0 /)
  dB_dq(4,:,3) = (/ 0, 1, 0 /)
  
  dB_dq(1,:,4) = (/ 0, 1, 0 /)
  dB_dq(2,:,4) = (/-1, 0, 0 /)
  dB_dq(3,:,4) = (/ 0, 0, 0 /)
  dB_dq(4,:,4) = (/ 0, 0, 1 /)

  do i=1,n
     A(1,:) = (/ q(1,i),-q(4,i), q(3,i) /)
     A(2,:) = (/ q(4,i), q(1,i),-q(2,i) /)
     A(3,:) = (/-q(3,i), q(2,i), q(1,i) /)
     A(4,:) = (/ q(2,i), q(3,i), q(4,i) /)

     B(1,:) = (/ q(1,i), q(4,i),-q(3,i) /)
     B(2,:) = (/-q(4,i), q(1,i), q(2,i) /)
     B(3,:) = (/ q(3,i),-q(2,i), q(1,i) /)
     B(4,:) = (/ q(2,i), q(3,i), q(4,i) /)

     do k=1,4
        J(i,:,:,k) = matmul(transpose(dA_dq(:,:,k)),B) + &
             matmul(transpose(A),dB_dq(:,:,k))
     end do
  end do

end subroutine computeRotMtxJacobian
  
