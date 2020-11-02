!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! EPR-related routines
!
! References:
!   [1] Phys. Rev. Lett. 63, 245101 (2001)   (norm-conserving GIPAW)
!   [2] Phys. Status Solidi C 7, 157 (2010) 

!====================================================================
! Relativistic mass correction bare + GIPAW
! Eq.(5) of [1]
!====================================================================
SUBROUTINE rmc(s_weight, delta_g_rmc, delta_g_rmc_gipaw)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k
  USE wavefunctions_module,   ONLY : evc
  USE becmod,                 ONLY : calbec  
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : nbnd_occ, radial_integral_rmc

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: s_weight
  real(dp), intent(out) :: delta_g_rmc, delta_g_rmc_gipaw

  !-- local variables ----------------------------------------------------
  integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, ijkb0
  complex(dp) :: rmc_corr, bec_product

  ! bare term
  do ibnd = 1, nbnd_occ(current_k)
     delta_g_rmc = delta_g_rmc - s_weight * wg(ibnd,current_k) * &
                   sum(g2kin(1:npw)*conjg(evc(1:npw,ibnd))*evc(1:npw,ibnd))
  end do

  ! gipaw term
  rmc_corr = (0.d0,0.d0)
  !do ibnd = 1, nbnd
  do ibnd = 1, nbnd_occ(current_k)
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) == nt) then
              do ih = 1, paw_recon(nt)%paw_nh
                 ikb = ijkb0 + ih
                 nbs1 = paw_recon(nt)%paw_indv(ih)
                 l1 = paw_recon(nt)%paw_nhtol(ih)
                 m1 = paw_recon(nt)%paw_nhtom(ih)
                 lm1 = m1+l1**2

                 do jh = 1, paw_recon(nt)%paw_nh
                    jkb = ijkb0 + jh
                    nbs2 = paw_recon(nt)%paw_indv(jh)
                    l2 = paw_recon(nt)%paw_nhtol(jh)
                    m2 = paw_recon(nt)%paw_nhtom(jh)
                    lm2 = m2+l2**2
                      
                    if ( l1 /= l2 .or. m1 /= m2) cycle
                      
                    bec_product = paw_becp(jkb,ibnd)*conjg(paw_becp(ikb,ibnd))
                    rmc_corr = rmc_corr + bec_product * wg(ibnd,current_k) * &
                               radial_integral_rmc(nbs1,nbs2,nt)
                 enddo
              enddo
              ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
           endif
        enddo
     enddo
  enddo
  delta_g_rmc_gipaw = delta_g_rmc_gipaw - s_weight * real(rmc_corr,dp)

END SUBROUTINE rmc



!====================================================================
! Paramagnetic contribution to the induced magnetic field
! norm-conserving contrib.: Eq.(11) of [1]
!====================================================================
SUBROUTINE paramagnetic_correction_so (paramagnetic_tensor, paramagnetic_tensor_us, &
                                       g_vel_evc, u_svel_evc, ipol)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k
  USE becmod,                 ONLY : calbec  
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, paw_becp2, paw_becp3, alpha, &
                                     radial_integral_paramagnetic_so, nbnd_occ
  USE uspp,                   ONLY : okvan

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: paramagnetic_tensor(3,3)
  real(dp), intent(inout) :: paramagnetic_tensor_us(3,3)
  complex(dp), intent(in) :: g_vel_evc(npwx,nbnd,3)
  complex(dp), intent(in) :: u_svel_evc(npwx,nbnd,3)
  integer, intent(in) :: ipol ! cartesian index of u_i

  !-- local variables ----------------------------------------------------
  complex(dp) :: para_corr(3), para_corr_us(3)
  complex(dp) :: bec_product, cc
  integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, ijkb0, jpol

  do jpol = 1, 3 
     if ( jpol == ipol ) cycle
     call calbec (npw, paw_vkb, g_vel_evc(:,:,jpol), paw_becp2)
     if (okvan) call calbec (npw, paw_vkb, u_svel_evc(:,:,jpol), paw_becp3)

     para_corr = (0.d0, 0.d0)
     para_corr_us = (0.d0, 0.d0)

     !do ibnd = 1, nbnd
     do ibnd = 1, nbnd_occ(current_k)
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              if (ityp (na) .eq.nt) then
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
                         
                       if ( l1 /= l2 ) cycle
                       bec_product = conjg(paw_becp(ikb,ibnd)) * paw_becp2(jkb,ibnd)

                       cc = bec_product * radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                            * wg(ibnd,current_k) * alpha**2
                       para_corr(1) = para_corr(1) + cc * lx(lm1,lm2)
                       para_corr(2) = para_corr(2) + cc * ly(lm1,lm2)
                       para_corr(3) = para_corr(3) + cc * lz(lm1,lm2)

                       if (okvan) then
                          bec_product = conjg(paw_becp(ikb,ibnd)) * paw_becp3(jkb,ibnd)
                          cc = bec_product * radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                               * wg (ibnd,current_k) * alpha ** 2
                          para_corr_us(1) = para_corr_us(1) + cc * lx(lm1,lm2)
                          para_corr_us(2) = para_corr_us(2) + cc * ly(lm1,lm2)
                          para_corr_us(3) = para_corr_us(3) + cc * lz(lm1,lm2)
                       endif

                    enddo
                 enddo
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif
           enddo
        enddo
     enddo
     paramagnetic_tensor(1:3,jpol) = real(para_corr(1:3),dp)
     if (okvan) paramagnetic_tensor_us(1:3,jpol) = real(para_corr_us(1:3),dp)
  end do
END SUBROUTINE paramagnetic_correction_so


  
!====================================================================
! Diamagnetic contribution to the induced magnetic field
! Eq.(9) and (10) of [1]
!====================================================================
SUBROUTINE diamagnetic_correction_so (diamagnetic_tensor)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k
  USE becmod,                 ONLY : calbec  
  USE constants,              ONLY : pi
  USE parameters,             ONLY : lmaxx
  USE uspp,                   ONLY : ap
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, paw_becp2, alpha, &
                                     radial_integral_diamagnetic_so, nbnd_occ
  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: diamagnetic_tensor(3,3)

  !-- local variables ----------------------------------------------------
  integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, lm, ijkb0
  complex(dp) :: dia_corr(lmaxx**2)
  complex(dp) :: bec_product, cc_dia
  
  dia_corr = 0.0_dp
  do ibnd = 1, nbnd_occ(current_k)
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if ( ityp(na) == nt ) then
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
                    lm2  = m2 + l2**2 
                      
                    bec_product = paw_becp(jkb,ibnd) * conjg( paw_becp(ikb,ibnd) )
                     
                    !<apsi> s/non-trace-zero component:
                    ! 2/3 to separate the non-trace vanishing component,
                    ! 1/(2 c^2) from the equation (59) in PM-PRB.
                    cc_dia = 2.0_dp/3.0_dp * bec_product & 
                             * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                             * wg(ibnd,current_k) * alpha**2  ! UWG: or ???
!                             * wg(ibnd,current_k) * alpha**2 / 2.0_dp 

                    if ( l1 == l2 .and. m1 == m2 ) then
                       diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) + cc_dia
                       diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) + cc_dia
                       diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) + cc_dia
                    endif
                    
                    ! 2/3 to separate the non-trace vanishing component
                    do lm = 5, 9
                       dia_corr(lm) = dia_corr(lm) + cc_dia * ap(lm,lm1,lm2)/2.d0
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
           endif
        enddo
     enddo
  enddo
    
  !  transform in cartesian coordinates
  dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)
  diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
  diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
  diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) + dia_corr(5) * 2.0_dp
  diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) + dia_corr(9) * sqrt(3.0_dp)
  diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
  diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) - dia_corr(6) * sqrt(3.0_dp)
  diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
  diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) - dia_corr(7) * sqrt(3.0_dp)
  diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)

END SUBROUTINE diamagnetic_correction_so



!====================================================================
! Ultrasoft augmentation (L_R Q_R) contribution to the bare and
! paramagnetic current
!====================================================================
SUBROUTINE paramagnetic_correction_aug_so (paug_corr_tensor, j_bare_s)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                     current_k, ecutwfc
  USE lsda_mod,               ONLY : current_spin, nspin
  USE wavefunctions_module,   ONLY : evc
  USE becmod,                 ONLY : calbec, allocate_bec_type, deallocate_bec_type
  USE constants,              ONLY : pi
  USE parameters,             ONLY : lmaxx
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, radial_integral_paramagnetic_so, &
                                     q_gipaw, alpha, nbnd_occ
  USE fft_base,               ONLY : dffts
  USE uspp,                   ONLY : ap, qq, vkb, nkb
  USE uspp_param,             ONLY : nh, upf
  USE cell_base,              ONLY : tpiba, omega, tpiba2
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g, ngm
  USE io_global,              ONLY : stdout, ionode
#ifdef __BANDS
  USE gipaw_module,           ONLY : ibnd_start, ibnd_end
  USE mp_image_global_module, ONLY : inter_image_comm
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : inter_bgrp_comm, mpime
#endif

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: paug_corr_tensor(3,3)
  real(dp), intent(inout) :: j_bare_s(dffts%nnr,3,3,nspin)

  !-- local variables ----------------------------------------------------
  integer ::ibnd,ijkb0,nt,na,ih,ikb,jh,jkb,ig,ipol,jpol,kpol,kfac 
  integer :: nbs1,nbs2,l1,l2,m1,m2,lm1,lm2, ik
  real(DP) , allocatable  :: gk (:,:), gg(:,:)
  complex(dp) , allocatable :: pcorr_jpaug(:) 
  complex(dp) , allocatable :: ps(:,:) ! 1st part
  complex(dp) , allocatable :: dvkbj(:,:), dvkby(:,:), Lp(:,:,:) ! 2nd part
!  complex(dp) , allocatable :: LQ(:,:,:) ! 3rd part
  complex(dp) , allocatable :: paw_becp_gLQ(:,:), aux_LQ(:,:) ! 3rd part
!  complex(dp) , allocatable :: g_LQ_evc(:,:,:), becp2(:,:)
  complex(dp) , allocatable :: g_LQ_evc(:,:), becp2(:,:)
  complex(dp) :: cc, bec_product
  real(dp) :: epsi(3,3), xyz(3,3),emine_q(3),dvkb_dir(3),ffact
  DATA epsi/0.d0,-3.d0,2.d0,3.d0,0.d0,-1.d0,-2.d0,1.d0,0.d0/, &
       xyz/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/

  !calculating ps = q_{ji}<p_i|u>
  emine_q(1)=0.d0; emine_q(2)=0.d0;emine_q(3)=0.d0

  !
  ! UWG: to reduce memory increase from NC to US case 
  !      some fields have been removed,
  !      others have been reduced to a single ipol direction.
  ! 
  
  allocate( ps( nkb, nbnd ), becp2(nkb,nbnd) )
  ik = current_k

  call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
  vkb = (0.d0,0.d0)
  call init_us_2 (npw, igk, xk(:,ik), vkb)
  ! TODO: calbec_bands
  call calbec (npwx, vkb, evc, becp2, nbnd)
  ps(:,:) = (0.d0, 0.d0)

  ijkb0 = 0 
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ibnd = 1, nbnd_occ(current_k)
              do ih = 1, nh(nt) 
                 ikb = ijkb0 + ih
                 do jh = 1, nh(nt)
                    jkb = ijkb0 + jh  
                    ps(jkb, ibnd) = ps(jkb, ibnd) + qq(jh,ih,nt) * becp2(ikb,ibnd)
                 enddo ! jh
              enddo ! ih
           enddo ! nbnd
           ijkb0 = ijkb0 + nh(nt)
        endif ! ityp(na)==nt
     enddo ! nat
  enddo ! ntyp
  !
  deallocate ( becp2 )


  ! now we have ps (nkb x nbnd)
  ! calculating L|p> = eps(abc)r_b p_c |p> dir_a = eps(abc)p_c r_b|p>
  !                  = eps(abc)p_c (d|p>/dk_b) (derivative wrt k_b)
  !                  = eps(abc)(k+G)_c(d|p>/dk_b) (check the p part again later)
  allocate(dvkbj(npwx,nkb), dvkby(npwx,nkb), Lp(npwx,nkb,3))
  allocate(gk(3,npwx),gg (3,npwx))
  Lp(:,:,:) = (0.d0,0.d0)
  dvkbj(:,:) = (0.d0,0.d0); dvkby(:,:) = (0.d0,0.d0)
  !
  gk(:,:) = 0.d0; gg(:,:) = 0.d0
  do ig = 1,npw
     gk(1:3,ig)=(xk(1:3,ik)+g(1:3,igk(ig)))*tpiba
     g2kin(ig) = SUM (gk(1:3,ig)**2)
     if(g2kin (ig) < 1.0d-10) then
       gg (:, ig) = 0.d0
     else
       gg (1:3,ig) = gk (1:3,ig) / sqrt(g2kin(ig))
     endif
  enddo

  call gen_us_dj(ik,dvkbj)
  do ipol = 1,3
     dvkb_dir(:)= xyz(:,ipol)
     call gen_us_dy(ik,dvkb_dir, dvkby)
     do jpol = 1,3
        kpol = int(abs(epsi(ipol,jpol)))
        if(kpol.eq.0)cycle
        kfac = int(epsi(ipol,jpol)/kpol)
        ijkb0 =0
        do nt = 1,ntyp
           do na = 1,nat
              if (nt==ityp(na)) then
                 do ikb = 1, nh(nt)
                    ijkb0 = ijkb0+1
                    do ig = 1,npw
                       Lp(ig,ijkb0,kpol)=Lp(ig,ijkb0,kpol)+ &
                                (dvkby(ig,ijkb0) + dvkbj(ig,ijkb0) &
                                 *gg(ipol,ig))*gk(jpol,ig)*kfac
                    enddo!npw
                 enddo !ikb
              endif !ityp(na)=nt
           enddo !na
        enddo !nt
     enddo !jpol
  enddo !ipol 
  deallocate ( dvkbj, dvkby)
  deallocate ( gk, gg )

  ! now we have both ps and Lp (npwx x nkb,3)
  ! hence, we can construct LQ (npw,nbnd) = LQ|u> = L|p>q<p|u>
  !               and apply Green's function:
  ! (to save memory all the rest is done in a single ipol-loop) 
  allocate ( aux_LQ(npwx,nbnd), pcorr_jpaug(3) )
  allocate ( g_LQ_evc(npwx,nbnd), paw_becp_gLQ(paw_nkb,nbnd) )
  do ipol = 1,3
     pcorr_jpaug(:) = (0.0d0,0.d0)
     aux_LQ(:,:) = (0.d0,0.d0)  
     ! compute  aux_LQ(:,:) = LQ(:,:,ipol)
     call zgemm ('N', 'N', npwx, nbnd, nkb, &
          (1.d0,0.d0), Lp(:,:,ipol), npwx, ps, nkb, (1.d0,0.d0), aux_LQ, npwx )
     ! preconditioning and calculation of g_LQ_evc
     g_LQ_evc(:,:) = (0.d0,0.d0)
     call greenfunction(ik, aux_LQ, g_LQ_evc(:,:), (/0.d0, 0.d0, 0.d0/))
     call calbec (npwx, paw_vkb, g_LQ_evc(:,:), paw_becp_gLQ)
     !
     do ibnd = 1, nbnd_occ(current_k)
        ijkb0 = 0
        do nt = 1, ntyp
           !IF ( upf(nt)%tvanp ) THEN
           !
           do na = 1, nat
              if (ityp(na).eq.nt) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih); l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih); lm1 = m1 + l1**2
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh); l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh); lm2 = m2 + l2**2
                       if(l1 /= l2) cycle    ! li (lm1,lm2) for different l are zero 
                       bec_product = conjg(paw_becp(ikb,ibnd))*paw_becp_gLQ(jkb,ibnd)
                       cc = bec_product*radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                            * wg(ibnd,ik)* alpha ** 2
                       pcorr_jpaug(1) = pcorr_jpaug(1) + cc * lx ( lm1, lm2 )
                       pcorr_jpaug(2) = pcorr_jpaug(2) + cc * ly ( lm1, lm2 )
                       pcorr_jpaug(3) = pcorr_jpaug(3) + cc * lz ( lm1, lm2 )
                    enddo !jh
                 enddo !ih
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif !ityp(na)==nt
           enddo !nat
           !
           !ENDIF
        enddo !ntyp 
     enddo !bands        
     paug_corr_tensor(:,ipol) = REAL (pcorr_jpaug(:), dp)
     !
     ! adding this to sigma is easy coz there is no cross product..
     ffact =  -( 2.0_dp * q_gipaw * tpiba )
     call j_para(ffact,evc,g_LQ_evc(:,:),ik,emine_q, j_bare_s(:,:,ipol,current_spin))
     !
  enddo !ipol   
  !write(77,*) 'tpiba =', tpiba
  !
  deallocate ( ps, Lp )
  deallocate ( g_LQ_evc, paw_becp_gLQ )
  deallocate ( pcorr_jpaug, aux_LQ )
  !
END SUBROUTINE paramagnetic_correction_aug_so



!====================================================================
! SO contribution to g-tensor
! Not exactly Eq.(6) of [1]: here we use the total potential
!====================================================================
SUBROUTINE compute_delta_g_so (j_bare, s_maj, s_min, delta_g_so)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                     current_k, ecutwfc
  USE lsda_mod,               ONLY : current_spin
  USE constants,              ONLY : pi
  USE cell_base,              ONLY : tpiba, omega, tpiba2
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g, ngm, nl
  USE fft_base,               ONLY : dfftp
  USE io_global,              ONLY : stdout, ionode
  USE scf,                    ONLY : vltot, v, rho
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : ry2ha, alpha, gprime, filgradv
  USE mp_pools,               ONLY : intra_pool_comm
  USE mp,                     ONLY : mp_sum

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: j_bare(dfftp%nnr,3,3,nspin)
  integer, intent(in) :: s_maj, s_min
  real(dp), intent(out) :: delta_g_so(3,3)

  !-- local variables ----------------------------------------------------
  ! UWG:  real(dp), allocatable :: grad_vr(:,:), v_local(:,:)
  real(dp), allocatable :: grad_vr(:,:,:), v_local(:)
  real(dp) :: d_omega
  integer :: ispin, ipol
  
  ! calculate the gradient of the potential:
  !
  ! UWG: the rest of the subroutine has been changed in order to take 
  !      into account the correct spin indices for grad_vr !!!  
  !      (for grad_vr we need an additional spin index. On the other hand:
  !      to save memory we can remove the spin index in v_local)
  ! NOTE: this requires a change of 'gradient', see 'gradient_vr' below.
  !
  allocate ( grad_vr(dfftp%nnr,3,nspin), v_local(dfftp%nnr) )
  do ispin = 1, nspin
!    this reproduces the old, slightly wrong values:
!                v_local(:) = vltot(:) + v%of_r(:,1)
!    correct:    v_local(:) = vltot(:) + v%of_r(:,ispin)
     v_local(:) = vltot(:) + v%of_r(:,1)
     call gradient_vr(dfftp%nnr, v_local, ngm, g, nl, grad_vr(:,:,ispin) ) 
     if (trim(filgradv) /= '') then
        call write_tensor_field_xsf(filgradv, 0, ispin, grad_vr(:,:,ispin))
     endif
  enddo
  deallocate (v_local)
  grad_vr = grad_vr * ry2ha
  
  do ipol = 1, 3
     delta_g_so ( ipol, 1 ) = SUM ( &
          ( j_bare(:,2,ipol,s_maj)*grad_vr(:,3,s_maj) -j_bare(:,2,ipol,s_min)*grad_vr(:,3,s_min) ) &
          -(j_bare(:,3,ipol,s_maj)*grad_vr(:,2,s_maj) -j_bare(:,3,ipol,s_min)*grad_vr(:,2,s_min) ) )
     delta_g_so ( ipol, 2 ) = SUM ( &
          ( j_bare(:,3,ipol,s_maj)*grad_vr(:,1,s_maj) -j_bare(:,3,ipol,s_min)*grad_vr(:,1,s_min) ) &
          -(j_bare(:,1,ipol,s_maj)*grad_vr(:,3,s_maj) -j_bare(:,1,ipol,s_min)*grad_vr(:,3,s_min) ) )
     delta_g_so ( ipol, 3 ) = SUM ( &
          ( j_bare(:,1,ipol,s_maj)*grad_vr(:,2,s_maj) -j_bare(:,1,ipol,s_min)*grad_vr(:,2,s_min) ) &
          -(j_bare(:,2,ipol,s_maj)*grad_vr(:,1,s_maj) -j_bare(:,2,ipol,s_min)*grad_vr(:,1,s_min) ) )
  end do
  deallocate (grad_vr)
  
#ifdef __MPI
  call mp_sum(delta_g_so, intra_pool_comm)
#endif

  d_omega = omega / real(dfftp%nr1*dfftp%nr2*dfftp%nr3, dp)
  delta_g_so = delta_g_so * d_omega
  delta_g_so = -delta_g_so * gprime / 2.d0 * alpha * 1.0d6

END SUBROUTINE compute_delta_g_so



!====================================================================
! SOO contribution to g-tensor
! There are two ways of computing it: 
! the first (in compute_delta_g_soo) is 'a la PARATEC',
! the second (in compute_delta_g_soo2) is according Eq.(3) of [2],
! which improves/clarifies Eq.(7) of [1].
!====================================================================
SUBROUTINE compute_delta_g_soo (j_bare, s_maj, s_min, delta_g_soo)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                     current_k, ecutwfc
  USE lsda_mod,               ONLY : current_spin
  USE constants,              ONLY : pi
  USE cell_base,              ONLY : tpiba, omega, tpiba2
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g, ngm, nl
  USE scf,                    ONLY : vltot, v, rho
  USE io_global,              ONLY : stdout, ionode
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : ry2ha, alpha, gprime
  USE fft_base,               ONLY : dfftp
  USE fft_interfaces,         ONLY : fwfft
  USE mp_pools,               ONLY : intra_pool_comm
  USE mp,                     ONLY : mp_sum

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: j_bare(dfftp%nnr,3,3,nspin)
!  real(dp), intent(in) :: B_ind_r(dfftp%nnr,3,3,nspin)
  integer, intent(in) :: s_maj, s_min
  real(dp), intent(out) :: delta_g_soo(3,3) !, delta_g_soo2(3,3)

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: grad_vh(:,:), vh(:,:)
  complex(dp), allocatable :: aux1(:)
  real(dp) :: d_omega, e_hartree, charge
  integer :: ipol, jpol

  ! Paratec-way:  int dr j_up(r) \times grad v_h[n_unpaired](r) 
  allocate (grad_vh(dfftp%nnr,3), vh(dfftp%nnr,nspin), aux1(dfftp%nnr))

  aux1(:) = rho%of_r(:,s_maj) - rho%of_r(:,s_min)
  call fwfft('Dense', aux1, dfftp)

  rho%of_g(:,1) = aux1(nl(:))
  rho%of_g(:,2) = 0.d0
  vh = 0.d0

  call v_h(rho%of_g, e_hartree, charge, vh)
  call gradient_vr (dfftp%nnr, vh, ngm, g, nl, grad_vh)
  grad_vh = grad_vh * ry2ha

  deallocate (vh, aux1)
  
  ! -2*j_bare_dw(r) because of the self-interaction correction:
  ! j_bare(r) - [j_bare_up(r)-j_bare_dw(r)] = -2*j_bare_dw(r)
  do ipol = 1, 3
    delta_g_soo(ipol,1) = 2.d0 * sum( j_bare(:,2,ipol,s_min) * grad_vh(:,3) &
                                    - j_bare(:,3,ipol,s_min) * grad_vh(:,2) )
    delta_g_soo(ipol,2) = 2.d0 * sum( j_bare(:,3,ipol,s_min) * grad_vh(:,1) &
                                    - j_bare(:,1,ipol,s_min) * grad_vh(:,3) )
    delta_g_soo(ipol,3) = 2.d0 * sum( j_bare(:,1,ipol,s_min) * grad_vh(:,2) &
                                    - j_bare(:,2,ipol,s_min) * grad_vh(:,1) )
  enddo
  deallocate ( grad_vh )
  
#ifdef __MPI
  call mp_sum(delta_g_soo, intra_pool_comm)
#endif

  d_omega = omega / real(dfftp%nr1*dfftp%nr2*dfftp%nr3, dp)
  delta_g_soo = delta_g_soo * d_omega
  delta_g_soo = delta_g_soo * 2.d0 * alpha * 1.0d6

END SUBROUTINE compute_delta_g_soo


SUBROUTINE compute_delta_g_soo2 (j_bare, B_ind_r, s_maj, s_min, delta_g_soo, delta_g_soo2)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                     current_k, ecutwfc
  USE lsda_mod,               ONLY : current_spin
  USE constants,              ONLY : pi
  USE cell_base,              ONLY : tpiba, omega, tpiba2
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g, ngm, nl
  USE scf,                    ONLY : vltot, v, rho
  USE io_global,              ONLY : stdout, ionode
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : ry2ha, alpha, gprime
  USE fft_base,               ONLY : dfftp
  USE fft_interfaces,         ONLY : fwfft
  USE mp_global,              ONLY : intra_pool_comm
  USE mp,                     ONLY : mp_sum

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: j_bare(dfftp%nnr,3,3,nspin)
  real(dp), intent(in) :: B_ind_r(dfftp%nnr,3,3,nspin)
  integer, intent(in) :: s_maj, s_min
  real(dp), intent(out) :: delta_g_soo(3,3), delta_g_soo2(3,3)

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: grad_vh(:,:), vh(:,:)
  complex(dp), allocatable :: aux1(:)
  real(dp) :: d_omega, e_hartree, charge
  integer :: ipol, jpol

  ! Paratec-way:  int dr j_up(r) \times grad v_h[n_unpaired](r) 
  allocate (grad_vh(dfftp%nnr,3), vh(dfftp%nnr,nspin), aux1(dfftp%nnr))

  aux1(:) = rho%of_r(:,s_maj) - rho%of_r(:,s_min)
  call fwfft('Dense', aux1, dfftp)

  rho%of_g(:,1) = aux1(nl(:))
  rho%of_g(:,2) = 0.d0
  vh = 0.d0

  call v_h(rho%of_g, e_hartree, charge, vh)
  call gradient_vr (dfftp%nnr, vh, ngm, g, nl, grad_vh)
  grad_vh = grad_vh * ry2ha

  deallocate (vh, aux1)
  
  ! -2*j_bare_dw(r) because of the self-interaction correction:
  ! j_bare(r) - [j_bare_up(r)-j_bare_dw(r)] = -2*j_bare_dw(r)
  do ipol = 1, 3
    delta_g_soo(ipol,1) = 2.d0 * sum( j_bare(:,2,ipol,s_min) * grad_vh(3,:) &
                                    - j_bare(:,3,ipol,s_min) * grad_vh(2,:) )
    delta_g_soo(ipol,2) = 2.d0 * sum( j_bare(:,3,ipol,s_min) * grad_vh(1,:) &
                                    - j_bare(:,1,ipol,s_min) * grad_vh(3,:) )
    delta_g_soo(ipol,3) = 2.d0 * sum( j_bare(:,1,ipol,s_min) * grad_vh(2,:) &
                                    - j_bare(:,2,ipol,s_min) * grad_vh(1,:) )
  enddo
  deallocate ( grad_vh )
  
#ifdef __MPI
  call mp_sum(delta_g_soo, intra_pool_comm)
#endif

  d_omega = omega / real(dfftp%nr1*dfftp%nr2*dfftp%nr3, dp)
  delta_g_soo = delta_g_soo * d_omega
  delta_g_soo = delta_g_soo * 2.d0 * alpha * 1.0d6


  ! This is according to Eq.(7) of [1] and Eq.(3) of [2]:
  do jpol = 1, 3
    do ipol = 1, 3
      delta_g_soo2(ipol,jpol) = sum( &
          (B_ind_r(:,ipol,jpol,1) + B_ind_r(:,ipol,jpol,2)) * &
          (rho%of_r(:,s_maj) - rho%of_r(:,s_min)) )
    enddo
  enddo

#ifdef __MPI
   call mp_sum(delta_g_soo2, intra_pool_comm)
#endif
  
  d_omega = omega / real(dfftp%nr1*dfftp%nr2*dfftp%nr3, dp)
  delta_g_soo2 = delta_g_soo2 * d_omega
  delta_g_soo2 = delta_g_soo2 * 2.d0 * 1.0d6

END SUBROUTINE compute_delta_g_soo2



!====================================================================
! Print the contributions to the g-tensor
!====================================================================
SUBROUTINE print_g_tensor(delta_g_rmc, delta_g_rmc_gipaw, delta_g_so, &
                          delta_g_soo, delta_g_soo2, delta_g_so_para, &
                          delta_g_so_para_us, delta_g_so_para_aug, &
                          delta_g_so_dia)
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, tau, atm, ityp, ntyp => nsp
  USE io_global,            ONLY : stdout
  USE symme,                ONLY : symtensor, symmatrix
  USE uspp,                 ONLY : okvan
  USE gipaw_module,         ONLY : g_e, gprime, alpha, ry2ha, filfield

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: delta_g_rmc, delta_g_rmc_gipaw
  real(dp), intent(inout) :: delta_g_so(3,3)
  real(dp), intent(inout) :: delta_g_soo(3,3), delta_g_soo2(3,3)
  real(dp), intent(inout) :: delta_g_so_para(3,3)
  real(dp), intent(inout) :: delta_g_so_para_us(3,3)
  real(dp), intent(inout) :: delta_g_so_para_aug(3,3)
  real(dp), intent(inout) :: delta_g_so_dia(3,3)
  !-- local variables ---------------------------------------------------
  real(dp) :: delta_g_tot(3,3), delta_g_tot2(3,3)

  write(stdout,'(5X,''Contributions to the EPR g-tensor (in ppm): -----------------------------'')')
  write(stdout,*)

  delta_g_rmc = delta_g_rmc * alpha**2 * g_e * ry2ha * 1.0d6
  delta_g_rmc_gipaw = delta_g_rmc_gipaw * alpha**2 * g_e * ry2ha * 1.0d6

  delta_g_so_dia = delta_g_so_dia * gprime / 4.d0 * ry2ha * 1.0d6
  delta_g_so_para = delta_g_so_para * gprime / 2.d0 * ry2ha * 1.0d6 
  delta_g_so_para_us = delta_g_so_para_us * gprime / 2.d0 * ry2ha * 1.0d6 
  delta_g_so_para_aug = delta_g_so_para_aug * gprime / 2.d0 * ry2ha * 1.0d6 

  ! calculate total  
  delta_g_tot = delta_g_so + delta_g_so_para + delta_g_so_para_aug + delta_g_so_para_us &
       + delta_g_so_dia + delta_g_soo
  delta_g_tot(1,1) = delta_g_tot(1,1) +  delta_g_rmc + delta_g_rmc_gipaw 
  delta_g_tot(2,2) = delta_g_tot(2,2) +  delta_g_rmc + delta_g_rmc_gipaw 
  delta_g_tot(3,3) = delta_g_tot(3,3) +  delta_g_rmc + delta_g_rmc_gipaw

  delta_g_tot2 = delta_g_so + delta_g_so_para + delta_g_so_para_aug + delta_g_so_para_us &
       + delta_g_so_dia + delta_g_soo2
  delta_g_tot2(1,1) = delta_g_tot2(1,1) +  delta_g_rmc + delta_g_rmc_gipaw 
  delta_g_tot2(2,2) = delta_g_tot2(2,2) +  delta_g_rmc + delta_g_rmc_gipaw 
  delta_g_tot2(3,3) = delta_g_tot2(3,3) +  delta_g_rmc + delta_g_rmc_gipaw


  write(stdout,'(5X,''Relativistic mass correction bare :'',F12.2)') delta_g_rmc
  write(stdout,'(5X,''Relativistic mass correction gipaw:'',F12.2)') delta_g_rmc_gipaw
  write(stdout,*)

  write(stdout,'(5X,''Delta_g SO - bare term:'')')
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_so(:,:)
  write(stdout,*)

  write(stdout,'(5X,''Delta_g SO - diamagnetic term:'')')
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_so_dia(:,:)
  write(stdout,*)

  write(stdout,'(5X,''Delta_g SO - paragnetic term:'')')
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_so_para(:,:)
  write(stdout,*)

  if (okvan) then
    write(stdout,'(5X,''Delta_g SO - paragnetic US occ-occ term:'')')
    write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_so_para_us(:,:)
    write(stdout,*)
    write(stdout,'(5X,''Delta_g SO - paramagnetic US L_R Q_R term:'')')
    write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_so_para_aug(:,:)
    write(stdout,*)
  endif

  write(stdout,'(5X,''Delta_g SOO - a la Paratec:'')')
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_soo(:,:)
  write(stdout,*)


  if (trim(filfield) /= '') then
    ! write(stdout,'(5X,''Delta_g SOO - Eq.(7):'')')
    write(stdout,'(5X,''Delta_g SOO - via B_ind:'')')
    write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_soo2(:,:)
    write(stdout,*)
  endif

  write(stdout,'(5X,''Delta_g total -----------------------------------------------------------'')')
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_tot(:,:)
  ! call symtensor(1, delta_g_tot)
  ! use symmatrix instead of symtensor
  call symmatrix(delta_g_tot)
  write(stdout,'(5X,''Delta_g total ----------------------- symmetrized -----------------------'')')
  write (stdout,*) 
  write(stdout,'(3(5X,3(F12.2,2X)/))') delta_g_tot(:,:)
  write(stdout,'(5X,''Delta_g total ==========================================================='')')

END SUBROUTINE print_g_tensor


!
!----------------------------------------------------------------------------
SUBROUTINE gradient_vr( nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga(nrxx,3)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(:,ipol) = tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE gradient_vr
!



