module cylindrical_curl
   use stel_kinds, only: rprec
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: curl_cyl
contains
   subroutine curl_cyl(r, z, dphi, a, valid_a, b, valid_curl)
      ! Physical components: final array index 1=R, 2=phi, 3=Z.
      ! Uniform grid; first/last toroidal planes are duplicates.
      ! Only centered stencils with all seven A samples valid are used.
      real(rprec), intent(in) :: r(:), z(:), dphi, a(:,:,:,:)
      logical, intent(in) :: valid_a(:,:,:)
      real(rprec), intent(out) :: b(:,:,:,:)
      logical, intent(out) :: valid_curl(:,:,:)
      integer :: nr, np, nz, nu, i, j, k, jm, jp
      real(rprec) :: dr, dz

      nr = size(r); np = size(a,2); nz = size(z); nu = np-1
      if (nr < 3 .or. nz < 3 .or. nu < 3) error stop 'Grid too small for curl'
      if (any(shape(a) /= [nr,np,nz,3])) error stop 'A shape mismatch'
      if (any(shape(b) /= shape(a))) error stop 'B shape mismatch'
      if (any(shape(valid_a) /= [nr,np,nz])) error stop 'A mask shape mismatch'
      if (any(shape(valid_curl) /= [nr,np,nz])) error stop 'Curl mask shape mismatch'
      if (any(r <= 0.0_rprec) .or. dphi <= 0.0_rprec) error stop 'Invalid R or dphi'
      if (any(r(2:) <= r(:nr-1))) error stop 'R must increase'
      if (any(z(2:) <= z(:nz-1))) error stop 'Z must increase'

      b = ieee_value(0.0_rprec, ieee_quiet_nan)
      valid_curl = .false.
      do k = 2, nz-1
         dz = z(k+1)-z(k-1)
         do j = 1, nu
            jm = modulo(j-2,nu)+1
            jp = modulo(j,nu)+1
            do i = 2, nr-1
               if (.not.(valid_a(i,j,k) .and. valid_a(i-1,j,k) .and. &
                   valid_a(i+1,j,k) .and. valid_a(i,jm,k) .and. &
                   valid_a(i,jp,k) .and. valid_a(i,j,k-1) .and. &
                   valid_a(i,j,k+1))) cycle
               dr = r(i+1)-r(i-1)
               b(i,j,k,1) = (a(i,jp,k,3)-a(i,jm,k,3))/(2*dphi*r(i)) &
                          -(a(i,j,k+1,2)-a(i,j,k-1,2))/dz
               b(i,j,k,2) = (a(i,j,k+1,1)-a(i,j,k-1,1))/dz &
                          -(a(i+1,j,k,3)-a(i-1,j,k,3))/dr
               b(i,j,k,3) = (r(i+1)*a(i+1,j,k,2)-r(i-1)*a(i-1,j,k,2))/(dr*r(i)) &
                          -(a(i,jp,k,1)-a(i,jm,k,1))/(2*dphi*r(i))
               valid_curl(i,j,k) = .true.
            end do
         end do
      end do
      b(:,np,:,:) = b(:,1,:,:)
      valid_curl(:,np,:) = valid_curl(:,1,:)
   end subroutine curl_cyl
end module cylindrical_curl
