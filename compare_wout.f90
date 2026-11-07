program compare_wout
   use stel_kinds, only: rprec
   use read_wout_mod, only: read_wout_file, nfp, ns, mpol, ntor, &
                          xm, xn, rmnc, zmns, rmns, zmnc, lasym
   use vmec_utils, only: GetBcyl, GetAcyl
   use cylindrical_curl, only: curl_cyl
   use, intrinsic :: iso_fortran_env, only: error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none

   integer :: nr=33, nphi=33, nz=33
   real(rprec) :: rmin, rmax, zmin, zmax, phi0=0.0_rprec
   real(rprec) :: s_min=1.e-4_rprec, s_max=0.98_rprec
   logical :: full_torus=.false.
   character(len=1024) :: output='comparison.csv'
   namelist /grid/ nr, nphi, nz, rmin, rmax, zmin, zmax, phi0, &
                   full_torus, s_min, s_max, output

   character(len=4096) :: wout_file, config_file, message
   integer :: i, j, k, ierr, ia, ib, iu, ios, nperiod
   integer :: na, nb, ncompare
   real(rprec), parameter :: twopi=8.0_rprec*atan(1.0_rprec)
   real(rprec) :: period, dphi, sa, ua, sb, ub, nan
   real(rprec) :: av(3), bv(3), diff(3), err2, ref2, cross, e2, b2
   real(rprec) :: max_abs, max_rel, rms, relrms, scale
   real(rprec), allocatable :: r(:), p(:), z(:), s(:,:,:)
   real(rprec), allocatable :: a(:,:,:,:), b(:,:,:,:), bc(:,:,:,:)
   logical, allocatable :: good_a(:,:,:), good_b(:,:,:), good_curl(:,:,:)

   if (command_argument_count() < 1 .or. command_argument_count() > 2) then
      write(error_unit,'(a)') 'Usage: compare_wout WOUT_FILE [grid.nml]'
      stop 1
   end if
   call get_command_argument(1, wout_file)
   if (trim(wout_file) == '--help') then
      write(*,'(a)') 'Usage: compare_wout WOUT_FILE [grid.nml]'
      write(*,'(a)') 'Defaults: 33 x 33 x 33, one field period, comparison.csv.'
      stop
   end if
   call read_wout_file(trim(wout_file), ierr)
   if (ierr /= 0) then
      write(error_unit,'(a,i0)') 'read_wout_file failed; status = ', ierr
      stop 1
   end if
   if (ns < 2 .or. nfp < 1) call die('Invalid WOUT dimensions')

   call boundary_bounds(rmin, rmax, zmin, zmax)
   if (command_argument_count() == 2) then
      call get_command_argument(2, config_file)
      open(newunit=iu, file=trim(config_file), status='old', action='read', iostat=ios, iomsg=message)
      if (ios /= 0) call die('Opening namelist: '//trim(message))
      read(iu, nml=grid, iostat=ios, iomsg=message)
      if (ios /= 0) call die('Reading namelist: '//trim(message))
      close(iu)
   end if
   if (nr < 3 .or. nz < 3 .or. nphi < 4) call die('Need nr,nz >= 3 and nphi >= 4')
   if (.not.all(ieee_is_finite([rmin,rmax,zmin,zmax,phi0,s_min,s_max]))) &
      call die('Nonfinite grid parameter')
   if (rmin <= 0 .or. rmax <= rmin .or. zmax <= zmin) call die('Invalid grid limits')
   if (s_min < 0 .or. s_max > 1 .or. s_max <= s_min) call die('Require 0 <= s_min < s_max <= 1')

   nperiod = nphi-1
   period = twopi/nfp
   if (full_torus) period = twopi
   dphi = period/nperiod
   allocate(r(nr), p(nphi), z(nz), s(nr,nphi,nz))
   allocate(a(nr,nphi,nz,3), b(nr,nphi,nz,3), bc(nr,nphi,nz,3))
   allocate(good_a(nr,nphi,nz), good_b(nr,nphi,nz), good_curl(nr,nphi,nz))
   r = [(rmin+(rmax-rmin)*real(i-1,rprec)/(nr-1), i=1,nr)]
   p = [(phi0+dphi*(j-1), j=1,nphi)]
   z = [(zmin+(zmax-zmin)*real(k-1,rprec)/(nz-1), k=1,nz)]
   nan = ieee_value(0.0_rprec, ieee_quiet_nan)
   a = nan; b = nan; s = nan
   good_a = .false.; good_b = .false.

   write(*,'(a,3(i0,1x))') 'Grid (R,phi,Z): ', nr, nphi, nz
   write(*,'(a,2es16.7)') 'R limits [m]: ', rmin, rmax
   write(*,'(a,2es16.7)') 'Z limits [m]: ', zmin, zmax
   write(*,'(a,2es16.7)') 'Phi limits [rad]: ', p(1), p(nphi)
   write(*,'(a)') 'Sampling GetBcyl and GetAcyl ...'
   do k = 1, nz
      do j = 1, nperiod
         do i = 1, nr
            ! Always initialize optional initial guesses. Do not use
            ! INTENT(OUT) dummies as initial guesses in GetAcyl_WOUT.
            sb = 0.0_rprec; ub = 0.0_rprec; bv = nan
            call GetBcyl(r(i), p(j), z(k), bv(1), bv(2), bv(3), &
                         sflx=sb, uflx=ub, info=ib)
            if (ib /= 0) cycle
            if (.not.all(ieee_is_finite([sb,ub]))) cycle
            if (sb <= s_min .or. sb > s_max) cycle
            s(i,j,k) = sb
            if (all(ieee_is_finite(bv))) then
               b(i,j,k,:) = bv
               good_b(i,j,k) = .true.
            end if
            sa = sb; ua = ub; av = nan
            call GetAcyl(r(i), p(j), z(k), av(1), av(2), av(3), &
                         sflx=sa, uflx=ua, info=ia)
            if (ia /= 0) cycle
            if (.not.all(ieee_is_finite([sa,ua]))) cycle
            if (sa <= s_min .or. sa > s_max) cycle
            if (.not.all(ieee_is_finite(av))) cycle
            a(i,j,k,:) = av
            good_a(i,j,k) = .true.
         end do
      end do
      write(*,'(a,i0,a,i0)') '  Z plane ', k, ' / ', nz
   end do
   a(:,nphi,:,:) = a(:,1,:,:); b(:,nphi,:,:) = b(:,1,:,:)
   s(:,nphi,:) = s(:,1,:)
   good_a(:,nphi,:) = good_a(:,1,:); good_b(:,nphi,:) = good_b(:,1,:)

   call curl_cyl(r, z, dphi, a, good_a, bc, good_curl)
   open(newunit=iu, file=trim(output), status='replace', action='write', iostat=ios, iomsg=message)
   if (ios /= 0) call die('Opening output: '//trim(message))
   write(iu,'(a)') 'R,phi,Z,s,Br_getBcyl,Bphi_getBcyl,Bz_getBcyl,'// &
                  'Br_curlA,Bphi_curlA,Bz_curlA,dBr,dBphi,dBz,abs_error,relative_error'
   na = count(good_a(:,1:nperiod,:)); nb = count(good_b(:,1:nperiod,:))
   ncompare = 0; err2 = 0; ref2 = 0; cross = 0; max_abs = 0; max_rel = 0
   do k = 2, nz-1
      do j = 1, nperiod
         do i = 2, nr-1
            if (.not.(good_b(i,j,k) .and. good_curl(i,j,k))) cycle
            if (.not.all(ieee_is_finite(bc(i,j,k,:)))) cycle
            bv = b(i,j,k,:); diff = bc(i,j,k,:)-bv
            e2 = sum(diff**2); b2 = sum(bv**2)
            ncompare = ncompare+1
            err2 = err2+e2; ref2 = ref2+b2
            cross = cross+sum(bc(i,j,k,:)*bv)
            max_abs = max(max_abs,sqrt(e2))
            max_rel = max(max_rel,sqrt(e2)/max(sqrt(b2),tiny(1.0_rprec)))
            write(iu,'(*(es24.15e3,:,","))') r(i),p(j),z(k),s(i,j,k), &
               bv,bc(i,j,k,:),diff,sqrt(e2),sqrt(e2)/max(sqrt(b2),tiny(1.0_rprec))
         end do
      end do
   end do
   close(iu)

   write(*,'(a,i0)') 'Valid A samples (unique phi planes): ', na
   write(*,'(a,i0)') 'Valid B samples (unique phi planes): ', nb
   write(*,'(a,i0)') 'Compared points with complete valid stencils: ', ncompare
   if (ncompare == 0) call die('No valid comparison points. Refine grid or adjust bounds/s limits.')
   if (ref2 <= tiny(1.0_rprec)) call die('Reference field norm is zero')
   rms = sqrt(err2/ncompare); relrms = sqrt(err2/ref2); scale = cross/ref2
   write(*,'(a,es15.7)') 'RMS vector error [T]:             ', rms
   write(*,'(a,es15.7)') 'Relative RMS vector error:        ', relrms
   write(*,'(a,es15.7)') 'Maximum vector error [T]:         ', max_abs
   write(*,'(a,es15.7)') 'Maximum pointwise relative error: ', max_rel
   write(*,'(a,es15.7)') 'Best-fit curl(A)/B scale:         ', scale
   write(*,'(a)') 'CSV: '//trim(output)

contains
   subroutine die(text)
      character(len=*), intent(in) :: text
      write(error_unit,'(a)') trim(text)
      stop 1
   end subroutine die

   subroutine boundary_bounds(rlo, rhi, zlo, zhi)
      ! Estimate a convenient rectangular box from the boundary Fourier
      ! series. These sampled extrema are defaults, not a rigorous bound.
      real(rprec), intent(out) :: rlo, rhi, zlo, zhi
      real(rprec) :: u, v, rr, zz
      real(rprec), allocatable :: angle(:)
      integer :: it, iv, nu, nv
      allocate(angle(size(xm)))
      nu = max(64,4*mpol); nv = max(16,4*(2*ntor+1))
      rlo = huge(1.0_rprec); rhi = -rlo; zlo = rlo; zhi = -rlo
      do iv = 0, nv-1
         v = twopi*iv/(nfp*real(nv,rprec))
         do it = 0, nu-1
            u = twopi*it/nu
            angle = xm*u-xn*v
            rr = sum(rmnc(:,ns)*cos(angle))
            zz = sum(zmns(:,ns)*sin(angle))
            if (lasym) then
               rr = rr+sum(rmns(:,ns)*sin(angle))
               zz = zz+sum(zmnc(:,ns)*cos(angle))
            end if
            rlo = min(rlo,rr); rhi = max(rhi,rr)
            zlo = min(zlo,zz); zhi = max(zhi,zz)
         end do
      end do
   end subroutine boundary_bounds
end program compare_wout
