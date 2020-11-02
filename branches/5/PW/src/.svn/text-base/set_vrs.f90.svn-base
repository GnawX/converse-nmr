!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  !--------------------------------------------------------------------
  ! set the total local potential vrs on the smooth mesh to be used in 
  ! h_psi, adding the (spin dependent) scf (H+xc) part and the sum of 
  ! all the local pseudopotential contributions.
  !
  USE kinds
  USE funct, only : dft_is_meta
  USE fft_base, only : dffts
  USE paw_gipaw, only : lambda_so, reconstruction_only
  implicit none

  integer :: nspin, nrxx
  ! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  ! input: the fft grid dimension
  real(DP) :: vrs (nrxx, nspin), vltot (nrxx), vr (nrxx, nspin), &
              kedtau(dffts%nnr,nspin), kedtaur(nrxx,nspin)
  ! output: total local potential on the smooth grid
  !         vrs=vltot+vr
  ! input: the total local pseudopotential
  ! input: the scf(H+xc) part of the local potential
  logical :: doublegrid
  ! input: true if a doublegrid is used

  integer:: is

  do is = 1, nspin
     !
     ! define the total local potential (external + scf) for each spin ...
     !
     if (is > 1 .and. nspin == 4) then
        !
        ! noncolinear case: only the first component contains vltot
        !
        vrs (:, is) = vr (:, is)
     else
        vrs (:, is) = vltot (:) + vr (:, is)
     end if
     !
     ! ... and interpolate it on the smooth mesh if necessary
     !
     if (doublegrid) call interpolate (vrs (1, is), vrs (1, is), - 1)
     if (dft_is_meta()) call interpolate(kedtaur(1,is),kedtau(1,is),-1)
  enddo
  !
  !
  ! UWG: ... and for spin-orbit coupling in GIPAW style:
  if (any(lambda_so /= 0.d0)) call set_dvrs(vrs, nrxx, nspin)
  !
  return

end subroutine set_vrs



!--------------------------------------------------------------------
subroutine set_dvrs (vrs, nrxx, nspin)
  !--------------------------------------------------------------------
  ! gradient of the total local potential vrs on the smooth mesh
  !
  USE kinds
  USE scf,           ONLY : dvrs
  USE gvect,         ONLY : nlm, ngm, g, nl
  USE gvecs,         ONLY : nls, ngms 
  USE fft_base,      ONLY : dffts
  implicit none

  integer :: nspin, nrxx
  real(dp) :: vrs (nrxx, nspin)
  real(dp), allocatable :: aux(:,:)
  integer:: is, ig, ipol

  allocate( aux(3,nrxx) )
  do is = 1, nspin
    ! UWG: check nrxx, nrx1s,.. for doublegrid !!!
    call gradient( nrxx, vrs(1,is), ngm, g, nl, aux )
    do ipol = 1, 3
      dvrs(1:nrxx,is,ipol) = aux(ipol,1:nrxx)
    enddo
  enddo
  deallocate (aux)
  ! 
  ! add correction due to spin-other orbit (SOO) contribution 
  call add_dvx (nrxx, nspin)
 
end subroutine set_dvrs



!--------------------------------------------------------------------
subroutine add_dvx (nrxx, nspin)
  !--------------------------------------------------------------------
  ! gradient of the total local potential vrs on the smooth mesh
  !
  USE kinds
  USE scf,           ONLY : dvrs
  USE gvect,         ONLY : nlm, ngm, g, nl
  USE gvecs,         ONLY : nls, ngms 
  USE fft_base,      ONLY : dffts, dfftp
  !
  USE scf,           ONLY : rho, rho_core, rhog_core
  USE funct,         ONLY : xc_spin
  USE constants,     ONLY : e2
  implicit none

  real(dp), allocatable :: vxc(:,:), aux(:,:)
  real(dp) :: etxc, vtxc
  ! 
  integer :: nspin, nrxx
  integer:: is, ig, ipol
  !
  real(dp) :: vrx (nrxx, nspin)
  real(dp) :: vx(2), vc(2), ex, ec 
  real(dp) :: rhox, arhox, zeta, sigma_soo
  integer:: ir 
  !
  sigma_soo = -3.d0
  !
!  goto 999

  !allocate ( vxc(dfftp%nnr,nspin) )  
  allocate ( vxc(nrxx,nspin), aux(3,nrxx) )  

!  call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  call v_x (rho, rho_core, rhog_core, etxc, vtxc, vxc)

  ! add the v_xc based correction to the 'effective gradient' of the potential
  do is = 1, nspin
     ! UWG: check nrxx, nrx1s,.. for doublegrid !!!
     call gradient( nrxx, vxc(1,is), ngm, g, nl, aux )
     do ipol = 1, 3
        dvrs(1:nrxx,is,ipol) = dvrs(1:nrxx,is,ipol) + sigma_soo * aux(ipol,1:nrxx)
     enddo
  enddo

  deallocate ( vxc, aux )  

  RETURN

  999 continue
  ! setup the density and calculate the exchange potential vx(:)
  DO ir = 1, dfftp%nnr
     rhox = rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir)
     arhox = ABS( rhox )
     IF ( arhox > 1.d-10) THEN
        zeta = ( rho%of_r(ir,1) - rho%of_r(ir,2) ) / arhox
        IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
        CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
        vrx(ir,:) = e2 * vx(:) 
!       vrx(ir,:) = e2 * ( vx(:) + vc(:) )
     END IF
  END DO

  ! add the correction to the effective gradient of the potential
  allocate( aux(3,nrxx) )
  do is = 1, nspin
     ! UWG: check nrxx, nrx1s,.. for doublegrid !!!
     call gradient( nrxx, vrx(1,is), ngm, g, nl, aux )
     do ipol = 1, 3
        dvrs(1:nrxx,is,ipol) = dvrs(1:nrxx,is,ipol) + sigma_soo * aux(ipol,1:nrxx)
     enddo
  enddo
  deallocate (aux)

end subroutine add_dvx
