program test_curl
   use stel_kinds, only: rprec
   use cylindrical_curl, only: curl_cyl
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(rprec) :: coarse, fine, ratio
   call run_case(17, coarse)
   call run_case(33, fine)
   ratio = coarse/fine
   write(*,'(a,es12.4)') 'Coarse maximum error: ', coarse
   write(*,'(a,es12.4)') 'Fine maximum error:   ', fine
   write(*,'(a,f8.4)') 'Refinement ratio:     ', ratio
   if (abs(ratio-4.0_rprec) > 0.1_rprec) error stop 'Expected second-order convergence'
   write(*,'(a)') 'PASS: analytic field, periodic endpoint, invalid masks, bounds.'
contains
   subroutine run_case(np, max_error)
      integer, intent(in) :: np
      real(rprec), intent(out) :: max_error
      integer, parameter :: nr=9, nz=9
      real(rprec), parameter :: twopi=8.0_rprec*atan(1.0_rprec)
      real(rprec), parameter :: ca=0.07_rprec, cb=0.03_rprec, cc=0.11_rprec, b0=2.0_rprec
      real(rprec) :: r(nr),z(nz),p(np),a(nr,np,nz,3),b(nr,np,nz,3)
      real(rprec) :: want(3),dphi,angle
      logical :: good(nr,np,nz),valid(nr,np,nz)
      integer :: i,j,k,nu
      nu = np-1
      dphi = twopi/(5*real(nu,rprec))
      r = [(4.0_rprec+0.1_rprec*(i-1),i=1,nr)]
      z = [(-0.4_rprec+0.1_rprec*(k-1),k=1,nz)]
      p = [(0.13_rprec+dphi*(j-1),j=1,np)]
      do k=1,nz
         do j=1,np
            angle=5*p(j)
            do i=1,nr
               a(i,j,k,1)=cc*z(k)*sin(angle)
               a(i,j,k,2)=0.5_rprec*b0*r(i)+ca*r(i)*z(k)
               a(i,j,k,3)=cb*r(i)**2*cos(angle)
            end do
         end do
      end do
      a(:,np,:,:)=a(:,1,:,:)
      good=.true.
      call curl_cyl(r,z,dphi,a,good,b,valid)
      max_error=0
      do k=2,nz-1
         do j=1,nu
            angle=5*p(j)
            do i=2,nr-1
               if (.not.valid(i,j,k)) error stop 'Interior unexpectedly invalid'
               want=[-5*cb*r(i)*sin(angle)-ca*r(i), &
                     cc*sin(angle)-2*cb*r(i)*cos(angle), &
                     b0+2*ca*z(k)-5*cc*z(k)*cos(angle)/r(i)]
               max_error=max(max_error,maxval(abs(b(i,j,k,:)-want)))
            end do
         end do
      end do
      if (any(valid(1,:,:)) .or. any(valid(nr,:,:))) error stop 'R edge valid'
      if (any(valid(:,:,1)) .or. any(valid(:,:,nz))) error stop 'Z edge valid'
      if (any(valid(:,np,:) .neqv. valid(:,1,:))) error stop 'Endpoint mask mismatch'
      if (maxval(abs(b(2:nr-1,np,2:nz-1,:)-b(2:nr-1,1,2:nz-1,:))) > &
          epsilon(1.0_rprec)) error stop 'Endpoint field mismatch'
      good(4,1,4)=.false.
      call curl_cyl(r,z,dphi,a,good,b,valid)
      if (valid(4,1,4) .or. valid(3,1,4) .or. valid(5,1,4) .or. &
          valid(4,2,4) .or. valid(4,nu,4) .or. valid(4,1,3) .or. &
          valid(4,1,5) .or. valid(4,np,4)) error stop 'Invalid stencil accepted'
      if (.not.all(ieee_is_nan(b(4,1,4,:)))) error stop 'Invalid field must be NaN'
   end subroutine run_case
end program test_curl
