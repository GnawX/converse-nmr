!
! Copyright (C) 2004-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  SUBROUTINE add_so_valence(ik, n, psi, p_psic)
  SUBROUTINE add_so_valence(ik, n, psi)
!  USE gipaw_module
  USE paw_gipaw
  USE parallel_include
  USE wavefunctions_module, ONLY : psic
  USE gvect,                ONLY : nlm, ngm, g
  USE cell_base,            ONLY : tpiba
  USE lsda_mod,             ONLY : current_spin, nspin
  USE klist,                ONLY : xk
  USE scf,                  ONLY : dvrs, vrs
  USE wvfct,                ONLY : igk
  USE gvecs,                ONLY : nls
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft
  !---------------------------------------------------------------
  ! add the SO valence term
  !---------------------------------------------------------------

  IMPLICIT NONE
  integer :: ik, n, ipol, ig, i
  complex(dp) :: p_psic(dffts%nnr,3) 
  complex(dp) :: psi(n)
  real(dp) :: gk, sigma
  integer :: anti_spin
  ! index for the cross product
  integer :: ind(2,3), ii, jj
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)

  ! derivative of the potential already available
  ! call set_dvrs(dvrs, vrs, dffts%nnr, nspin)

  ! compute (p+k)|psi> in real space
  p_psic(1:dffts%nnr,1:3) = (0.d0,0.d0)
  do ipol = 1, 3
    do ig = 1, n
      gk = xk(ipol,ik) + g(ipol,igk(ig))
      p_psic(nls(igk(ig)),ipol) = -gk * tpiba * psi(ig)
    enddo
    CALL invfft ('Wave', p_psic(:,ipol), dffts)
  enddo

  ! compute lambda_so ( sigma \cdot dvrs \times p_psi)
  sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

  do ipol = 1, 3
    if (lambda_so(ipol) == 0.d0) cycle
    ii = ind(1,ipol)
    jj = ind(2,ipol)
    do i = 1, dffts%nnr 
      psic(i) = psic(i) + lambda_so(ipol) * a2gp8 * sigma * &
        ( dvrs(i,current_spin,ii)*p_psic(i,jj) - &
          dvrs(i,current_spin,jj)*p_psic(i,ii) )
    enddo
  enddo
  END SUBROUTINE add_so_valence



  !---------------------------------------------------------------
  ! add the F_R^{NL} term of the SO reconstruction
  !---------------------------------------------------------------
  SUBROUTINE add_so_Fnl(lda, n, m, psi, hpsi)
!  USE gipaw_module
  USE paw_gipaw
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,        ONLY : current_spin
  USE wvfct,           ONLY : current_k, npwx, npw
  USE klist,           ONLY : xk
  USE wvfct,           ONLY : igk, g2kin, nbndx, nbnd
  USE becmod,          ONLY : calbec !, becp 
  USE paw_gipaw,       ONLY : paw_vkb, paw_nkb, paw_recon, paw_becp
  IMPLICIT NONE
  INTEGER :: lda, n, m
  INTEGER :: ibnd, ijkb0, nt, na, ikb, jkb, ih, jh
  INTEGER :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
  COMPLEX(DP), ALLOCATABLE :: ps (:,:)
  COMPLEX(DP) :: psi(lda,n), hpsi(lda,n)
  real(dp) :: sigma

  !!RETURN
  if (m > nbndx) call errore('add_so_Fnl', 'm > nbndx ???', m)
  if (m > nbnd) call errore('add_so_Fnl', 'm > nbnd ???', m)
  ALLOCATE (ps(paw_nkb,m))
  ps(:,:) = (0.D0, 0.D0)

  !      call init_paw_2(n, igk, xk(1,current_k), paw_vkb)
  ! it was: call calbec(paw_nkb, lda, n, m, paw_becp, paw_vkb, psi)

  call calbec (npw, paw_vkb, psi, paw_becp, m)

  sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

  ijkb0 = 0
  do nt = 1, ntyp
    do na = 1, nat
      if (ityp(na) .eq. nt) then
        do ih = 1, paw_recon(nt)%paw_nh
          ikb = ijkb0 + ih
          nbs1 = paw_recon(nt)%paw_indv(ih)
          l1 = paw_recon(nt)%paw_nhtol(ih)
          m1 = paw_recon(nt)%paw_nhtom(ih)
          lm1 = m1 + l1**2

          do jh = 1, paw_recon(nt)%paw_nh
            jkb = ijkb0 + jh
            nbs2 = paw_recon(nt)%paw_indv(jh)
            l2 = paw_recon(nt)%paw_nhtol(jh)
            m2 = paw_recon(nt)%paw_nhtom(jh)
            lm2 = m2 + l2**2

            if (l1 /= l2) cycle
            if (l1 == 0) cycle

            do ibnd = 1, m
              ps(ikb,ibnd) = ps(ikb,ibnd) + &
                lambda_so(1) * a2gp8 * sigma * lx(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)
            
              ps(ikb,ibnd) = ps(ikb,ibnd) + &
                lambda_so(2) * a2gp8 * sigma * ly(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)
             
              ps(ikb,ibnd) = ps(ikb,ibnd) + &
                lambda_so(3) * a2gp8 * sigma * lz(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)
            
            enddo   ! ibnd
          enddo   ! jh
        enddo   ! ih
        ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
      endif    ! ityp(na) .eq. nt
    enddo   ! na
  enddo   ! nt

  CALL ZGEMM( 'N', 'N', n, m, paw_nkb, ( 0.D0, 1.0D0 ) , paw_vkb, &
              lda, ps, paw_nkb, ( 1.D0, 0.D0 ) , hpsi, lda )

  deallocate (ps)
  END SUBROUTINE add_so_Fnl 



