subroutine update_new ( &
                              lo, hi, &
                               u, ulo, uhi, &
                               v, vlo, vhi, &
                               u_new, nulo, nuhi, &
                               v_new, nvlo, nvhi, &
                               dt, uFactor, vFactor,F,Kf) bind(C, name="update_new")

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: ulo(3), uhi(3)
  integer, intent(in) :: vlo(3), vhi(3)
  integer, intent(in) :: nulo(3), nuhi(3), nvlo(3), nvhi(3)
  real*8, intent(in)    :: u  (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
  real*8, intent(in)    :: v  (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
  real*8, intent(inout) :: u_new( nulo(1): nuhi(1), nulo(2): nuhi(2), nulo(3): nuhi(3))
  real*8, intent(inout) :: v_new( nvlo(1): nvhi(1), nvlo(2): nvhi(2), nvlo(3): nvhi(3))
  real*8, intent(in) :: dt, F, Kf, uFactor, vFactor

  ! local variables
  integer i,j,k

  ! x-fluxes
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           u_new(i,j,k) = u(i,j,k) + uFactor * ( u(i+1,j,k) + u(i-1,j,k) + &
                                    u(i,j+1,k) + u(i,j-1,k) + &
                                    u(i,j,k-1) + u(i,j,k+1) - &
                                    6.0*u(i,j,k) ) - &
                                    dt * u(i,j,k)*v(i,j,k)*v(i,j,k) - &
                                    dt * F * (u(i,j,k) - 1.0)


           v_new(i,j,k) = v(i,j,k) + vFactor * ( v(i+1,j,k) + v(i-1,j,k) + &
                                    v(i,j+1,k) + v(i,j-1,k) + &
                                    v(i,j,k-1) + v(i,j,k+1) - &
                                    6.0*v(i,j,k) ) + &
                                    dt * u(i,j,k)*v(i,j,k)*v(i,j,k) - &
                                    dt * (F+Kf) * v(i,j,k)
        end do
     end do
  end do


end subroutine update_new

