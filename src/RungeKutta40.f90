subroutine computeResiduals(ny, nv, nd, h, y0, y, dat, R)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ny, nv, nd, h, y0, y, dat
  !f2py intent(out) R
  !f2py depend(nv) y0
  !f2py depend(ny) y
  !f2py depend(ny,nd) dat
  !f2py depend(ny) R

  !Input
  integer, intent(in) ::  ny, nv, nd
  double precision, intent(in) ::  h, y0(nv), y(ny), dat(ny,nd)

  !Output
  double precision, intent(out) ::  R(ny)

  !Working
  integer k, k1, k2, n
  double precision a(nv), b(nv), c(nv), d(nv)

  n = ny/nv

  R(1:nv) = y(1:nv) - y0(:)

  do k=1,n-1
     k1 = (k-1)*nv + 1
     k2 = k*nv
     call getf(dat(k1:k2,:), y(k1:k2), a)
     call getf(dat(k1:k2,:), y(k1:k2) + h/2*a, b)
     call getf(dat(k1:k2,:), y(k1:k2) + h/2*b, c)
     call getf(dat(k1:k2,:), y(k1:k2) + h*c, d)
     R(nv+k1:nv+k2) = y(nv+k1:nv+k2) - y(k1:k2) - h/6*(a + 2*b + 2*c + d)
  end do

end subroutine computeResiduals



subroutine computeJacobian(ny, nv, nJ, nd, h, y, dat, Ja, Ji, Jj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ny, nv, nJ, nd, h, y, dat
  !f2py intent(out) Ja, Ji, Jj
  !f2py depend(ny) y
  !f2py depend(ny,nd) dat
  !f2py depend(nJ) Ja, Ji, Jj

  !Input
  integer, intent(in) ::  ny, nv, nJ, nd
  double precision, intent(in) ::  h, y(ny), dat(ny,nd)
  
  !Output
  double precision, intent(out) ::  Ja(nJ)
  integer, intent(out) ::  Ji(nJ), Jj(nJ)

  !Working
  integer i, j, k, k1, k2, iJ, n
  double precision a(nv), b(nv), c(nv), d(nv)
  double precision dR_dy(nv,nv), eye(nv,nv)
  double precision da_dy(nv,nv), db_dy(nv,nv), dc_dy(nv,nv), dd_dy(nv,nv)
  double precision df_dy(nv,nv), dg_dy(nv,nv), dh_dy(nv,nv), di_dy(nv,nv)

  n = ny/nv

  eye(:,:) = 0.0
  do k=1,nv
     eye(k,k) = 1.0
  end do

  Ja(1:ny) = 1.0
  do k=1,ny
     Ji(k) = k - 1
     Jj(k) = k - 1
  end do

  do k=1,n-1
     k1 = (k-1)*nv + 1
     k2 = k*nv
     call    getf(dat(k1:k2,:), y(k1:k2)        , a)
     call    getf(dat(k1:k2,:), y(k1:k2) + h/2*a, b)
     call    getf(dat(k1:k2,:), y(k1:k2) + h/2*b, c)
     call    getf(dat(k1:k2,:), y(k1:k2) + h*c  , d)
     call getdfdy(dat(k1:k2,:), y(k1:k2)        , df_dy)
     call getdfdy(dat(k1:k2,:), y(k1:k2) + h/2*a, dg_dy)
     call getdfdy(dat(k1:k2,:), y(k1:k2) + h/2*b, dh_dy)
     call getdfdy(dat(k1:k2,:), y(k1:k2) + h*c  , di_dy)
     da_dy = df_dy
     db_dy = matmul(dg_dy, eye + h/2*da_dy)
     dc_dy = matmul(dh_dy, eye + h/2*db_dy)
     dd_dy = matmul(di_dy, eye + h*dc_dy)
     dR_dy = -eye - h/6*(da_dy + 2*db_dy + 2*dc_dy + dd_dy)
     do i=1,nv
        do j=1,nv
           iJ = ny + (k-1)*nv*nv + (j-1)*nv + i
           Ja(iJ) = dR_dy(i,j)
           Ji(iJ) = k*nv + i - 1
           Jj(iJ) = (k-1)*nv + j - 1
        end do
     end do
  end do  

end subroutine computeJacobian
  



subroutine computeSolution(ny, nv, nd, h, y0, dat, y)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ny, nv, nd, h, y0, dat
  !f2py intent(out) y
  !f2py depend(nv) y0
  !f2py depend(ny,nd) dat
  !f2py depend(ny) y

  !Input
  integer, intent(in) ::  ny, nv, nd
  double precision, intent(in) ::  h, y0(nv), dat(ny,nd)

  !Output
  double precision, intent(out) ::  y(ny)

  !Working
  integer k, k1, k2, n
  double precision a(nv), b(nv), c(nv), d(nv)

  n = ny/nv

  y(1:nv) = y0(:)

  do k=1,n-1
     k1 = (k-1)*nv + 1
     k2 = k*nv
     call getf(dat(k1:k2,:), y(k1:k2), a)
     call getf(dat(k1:k2,:), y(k1:k2) + h/2*a, b)
     call getf(dat(k1:k2,:), y(k1:k2) + h/2*b, c)
     call getf(dat(k1:k2,:), y(k1:k2) + h*c, d)
     y(nv+k1:nv+k2) = y(k1:k2) + h/6*(a + 2*b + 2*c + d)
  end do

end subroutine computeSolution
