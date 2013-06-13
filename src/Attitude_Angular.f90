subroutine computew(n, O, Odot, w)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, O, Odot
  !f2py intent(out) w
  !f2py depend(n) O, Odot, w

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  O(3,3,n), Odot(3,3,n)

  !Output
  double precision, intent(out) ::  w(3,n)

  !Working
  integer i

  do i=1,n
     w(1,i) = dot_product(Odot(3,:,i) , O(2,:,i))
     w(2,i) = dot_product(Odot(1,:,i) , O(3,:,i))
     w(3,i) = dot_product(Odot(2,:,i) , O(1,:,i))
  end do

end subroutine computew




subroutine computeJacobianw(n, O, Odot, dw_dO, dw_dOdot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, O, Odot
  !f2py intent(out) dw_dO, dw_dOdot
  !f2py depend(n) O, Odot, dw_dO, dw_dOdot

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  O(3,3,n), Odot(3,3,n)

  !Output
  double precision, intent(out) ::  dw_dO(n,3,3,3), dw_dOdot(n,3,3,3)

  !Working
  integer i

  dw_dO(:,:,:,:) = 0.0
  dw_dOdot(:,:,:,:) = 0.0
  do i=1,n
     !w(1,i) = Odot(3,:,i) * O(2,:,i)
     dw_dOdot(i,1,3,:) = O(2,:,i)
     dw_dO(i,1,2,:) = Odot(3,:,i)

     !w(2,i) = Odot(1,:,i) * O(3,:,i)
     dw_dOdot(i,2,1,:) = O(3,:,i)
     dw_dO(i,2,3,:) = Odot(1,:,i)

     !w(3,i) = Odot(2,:,i) * O(1,:,i)
     dw_dOdot(i,3,2,:) = O(1,:,i)
     dw_dO(i,3,1,:) = Odot(2,:,i)
  end do

end subroutine computeJacobianw
