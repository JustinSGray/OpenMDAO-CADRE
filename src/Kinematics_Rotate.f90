subroutine computePositionRotd(n, v1, O_21, v2)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, v1, O_21
  !f2py intent(out) v2
  !f2py depend(n) v1, O_21, v2

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  v1(3,n), O_21(3,3,n)

  !Output
  double precision, intent(out) ::  v2(3,n)

  !Working
  integer i

  do i=1,n
     v2(:,i) = matmul(O_21(:,:,i),v1(:,i))
  end do

end subroutine computePositionRotd



subroutine computePositionRotdJacobian(n, v1, O_21, J1, J2)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, v1, O_21
  !f2py intent(out) J1, J2
  !f2py depend(n) v1, O_21, J1, J2

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  v1(3,n), O_21(3,3,n)

  !Output
  double precision, intent(out) ::  J1(n,3,3,3), J2(n,3,3)

  !Working
  integer i, j, k, u, v, iJ1, iJ2
  double precision eye(3,3)

  eye(:,:) = 0.0
  eye(1,1) = 1.0
  eye(2,2) = 1.0
  eye(3,3) = 1.0
  do i=1,n
     !v2(:,i) = matmul(O_21(:,:,i),v1(:,i))
     do k=1,3
        do u=1,3
           do v=1,3
              J1(i,k,u,v) = eye(k,u)*v1(v,i)
           end do
        end do
        do j=1,3
           J2(i,k,j) = O_21(k,j,i)
        end do
     end do
  end do  

end subroutine computePositionRotdJacobian
