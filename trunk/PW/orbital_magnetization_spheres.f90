!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE orbital_magnetization_spheres
!-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_files,             ONLY : nwordwfc, iunwfc  
  USE wvfct,                ONLY : npw, g2kin, igk, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                ONLY : xk, nks
  USE cell_base,            ONLY : tpiba2, tpiba
  USE gvect,                ONLY : nr1, nr2, nr3, nrxx, ngm, ecutwfc, g
  USE gsmooth,              ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE io_global,            ONLY : stdout
  USE uspp,                 ONLY : nkb, vkb
  USE becmod,               ONLY : becp
  USE wavefunctions_module, ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
  USE g_tensor_module,      ONLY : lambda_so, add_so_Fnl, calc_g_tensor
  USE paw,                  ONLY : paw_vkb, paw_nkb, paw_becp
  USE control_flags,        ONLY : iverbosity
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : my_pool_id, me_pool, root_pool
  USE mp,                   ONLY : mp_barrier
  USE g_tensor_module,      ONLY : init_paw_2_no_phase, init_us_2_no_phase, r_mt
  USE gsmooth,              ONLY : nr1s, nr2s, nr3s, nrxxs, ngms, nls
  !--------------------------------------------------------------------
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE uspp_param, ONLY : nh
  USE uspp,       ONLY : deeq
  USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
  USE g_tensor_module, ONLY : radial_integral_L
  USE g_tensor_module, ONLY : lx, ly, lz
  USE g_tensor_module, ONLY: pointlist, distlist, pointnum, r_mt
  !--------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------
  complex(dp), external :: ZDOTC
  integer :: ik, ipol, ig, ir, ip
  integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
  integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
  complex(dp), allocatable :: ang_mom(:,:)
  complex(dp) :: becp_product
  complex(dp), allocatable, dimension(:) :: psic
  complex(dp), allocatable, dimension(:,:) :: p_psic
  real(dp) :: gk, dist(3)
  complex(dp) :: tmpL(3)
  !--------------------------------------------------------------------

  if (paw_nkb == 0) call errore("orb_magn_spheres", "no GIPAW projectors?", 1)

  ! allocate arrays
  allocate(psic(nrxx), p_psic(nrxx,3))
  allocate(ang_mom(3, nat))
  ang_mom = 0.d0

  ! setup atomic spheres
  write(stdout,*)
  allocate(pointlist(nrxx,nat))
  allocate(distlist(3,nrxx,nat))
  allocate(pointnum(nat))
  call make_pointlists_my(r_mt)
  do na = 1, nat
    write(stdout,'(5X,''atom:'',I4,''  npoints='',I6)') na, pointnum(nat)
  enddo

  ! loop over k-points
  do ik = 1, nks
    current_k = ik
    current_spin = 1
    if (lsda) current_spin = isk(ik)
   
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(1:npw) = g2kin(1:npw) * tpiba2

    ! load wfcs from disk
    call get_buffer(evc, nwordwfc, iunwfc, ik)

    ! loop over bands
    do ibnd = 1, nbnd
      ! fft evc to real space
      psic(1:nrxx) = ( 0.D0, 0.D0 )
      psic(nls(igk(1:npw))) = evc(1:npw,ibnd)
      CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )

      ! calculate p|evc> in real space
      p_psic(1:nrxx,1:3) = (0.d0,0.d0)
      do ipol = 1, 3
        do ig = 1, npw
          gk = xk(ipol,ik) + g(ipol,igk(ig))
          p_psic(nls(igk(ig)),ipol) = -gk * tpiba * evc(ig,ibnd)
        enddo
        call cft3s( p_psic(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
      enddo

      ! calculate <L> for this band and k-point
      do na = 1, nat
        tmpL(1:3) = (0.d0, 0.d0)
        do ir = 1, pointnum(na)
          ip = pointlist(ir,na)
          dist(1:3) = distlist(1:3,ir,na)
          !!write(90,'(3(E12.4,2X))') dist
          tmpL(1) = tmpL(1) + conjg(psic(ip)) * (dist(2)*p_psic(ip,3) - dist(3)*p_psic(ip,2))
          tmpL(2) = tmpL(2) + conjg(psic(ip)) * (dist(3)*p_psic(ip,1) - dist(1)*p_psic(ip,3))
          tmpL(3) = tmpL(3) + conjg(psic(ip)) * (dist(1)*p_psic(ip,2) - dist(2)*p_psic(ip,1))

          !tmpL(1) = tmpL(1) + conjg(psic(ip)) * psic(ip)
          !tmpL(2) = tmpL(2) + conjg(psic(ip)) * psic(ip)
          !tmpL(3) = tmpL(3) + conjg(psic(ip)) * psic(ip)

          !tmpL(1) = tmpL(1) + conjg(psic(ip)) * p_psic(ip,1)
          !tmpL(2) = tmpL(2) + conjg(psic(ip)) * p_psic(ip,2)
          !tmpL(3) = tmpL(3) + conjg(psic(ip)) * p_psic(ip,3)
        enddo
        !!write(*,'(2I4,4X,6(F10.4,2X))') ik, ibnd, wg(ibnd,ik)*tmpL(1:3)/DBLE(nrxx)
        ang_mom(1:3,na) = ang_mom(1:3,na) + wg(ibnd,ik) * tmpL(1:3)/DBLE(nrxx)
      enddo ! nat

   enddo ! ibnd
  
  enddo ! ik
  ! reduction in parallel?

  ! print results
  write(stdout,*)
  write(stdout,'(5X,''Orbital magnetization in the atomic spheres:'')')
  write(stdout,'(5X,''(bare valence term)'')')
  do na = 1, nat
    write(stdout, '(5X,''atom:'',I4,''  L='',3(F10.4))') na, real(ang_mom(:,na))
    !!print*, ang_mom(1,na)
    !!print*, ang_mom(2,na)
    !!print*, ang_mom(3,na)
  enddo

  deallocate(pointlist, distlist, pointnum)
  deallocate(psic, p_psic)

  ! NOW THE RECONSTUCTION PART
  ang_mom = 0.d0

  ! loop over k-points
  do ik = 1, nks
    current_k = ik
    current_spin = 1
    if (lsda) current_spin = isk(ik)
   
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(1:npw) = g2kin(1:npw) * tpiba2

    ! load wfcs from disk
    call get_buffer(evc, nwordwfc, iunwfc, ik)

    ! setup the projectors
    call init_paw_2(npw, igk, xk(1,ik), paw_vkb)
    call ccalbec(paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)

    ! loop over atoms and projectors
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

              if ( l1 /= l2 ) cycle
              if (l1 == 0) cycle

              do ibnd = 1, nbnd
                becp_product = conjg(paw_becp(ikb,ibnd)) * &
                                     paw_becp(jkb,ibnd)
                ang_mom(1,na) = ang_mom(1,na) + lx(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_L(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                ang_mom(2,na) = ang_mom(2,na) + ly(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_L(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                ang_mom(3,na) = ang_mom(3,na) + lz(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_L(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
              enddo   ! ibnd
            enddo   ! jh
          enddo   ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif   ! if (ityp...
      enddo   ! na
    enddo  ! nt

  enddo ! ik
  ! reduction in parallel??

  ! print results
  write(stdout,*)
  write(stdout,'(5X,''Orbital magnetization in the atomic spheres:'')')
  write(stdout,'(5X,''(GIPAW term)'')')
  do na = 1, nat
    write(stdout, '(5X,''atom:'',I4,''  L='',3(F10.4))') na, real(ang_mom(:,na))
    !!print*, ang_mom(1,na)
    !!print*, ang_mom(2,na)
    !!print*, ang_mom(3,na)
  enddo

!-----------------------------------------------------------------------
END SUBROUTINE orbital_magnetization_spheres
!-----------------------------------------------------------------------









