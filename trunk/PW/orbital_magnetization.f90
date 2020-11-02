!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE orbital_magnetization
  !-----------------------------------------------------------------------
  !
  ! ... the title is self explaining
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  SAVE
  REAL(DP) :: orb_magn_LC(3)
  REAL(DP) :: orb_magn_IC(3)
  REAL(DP) :: berry_curvature(3)
  REAL(DP) :: orb_magn_tot(3)
  REAL(DP) :: delta_M_bare(3)
  REAL(DP) :: delta_M_para(3)
  REAL(DP) :: delta_M_dia(3)
  COMPLEX(dp), ALLOCATABLE :: dbecp(:,:,:)
  COMPLEX(dp), ALLOCATABLE :: paw_dbecp(:,:,:)

  !real(dp), parameter :: q_gipaw = 0.02_dp
  real(dp) :: delta_k

  LOGICAL :: torbmagn 
  CHARACTER(80) :: dudk
  
CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE calc_orbital_magnetization
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_files,             ONLY : nwordwfc, iunwfc  
  USE wvfct,                ONLY : npw, g2kin, igk, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                ONLY : xk, nks, lgauss
  USE gsmooth,              ONLY : nrxxs
  USE cell_base,            ONLY : tpiba2, tpiba
  USE gvect,                ONLY : nr1, nr2, nr3, ngm, ecutwfc, g
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunsat, nwordatwfc
  USE uspp,                 ONLY : nkb, vkb
  USE becmod,               ONLY : becp
  USE wavefunctions_module, ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
  USE g_tensor_module,      ONLY : lambda_so, add_so_Fnl, calc_g_tensor
  USE nmr_mod,              ONLY : m_0_atom, m_0, calc_chemical_shift
  USE paw,                  ONLY : paw_vkb
  USE ener,                 ONLY : ef
  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE control_flags,        ONLY : iverbosity
  USE paw,                  ONLY : paw_nkb, paw_becp
  USE buffers,              ONLY : get_buffer
  USE scf,                  ONLY : dvrs, vrs
  USE mp_global,            ONLY : my_pool_id, me_pool, root_pool, mpime
  USE mp,                   ONLY : mp_barrier
  USE g_tensor_module,      ONLY : init_paw_2_no_phase, init_us_2_no_phase
  USE g_tensor_module,      ONLY : q_gipaw, nmr_shift_core
  implicit none
  complex(dp), external :: ZDOTC
  integer, parameter :: iundudk1 = 75, iundudk2 = 76, iundudk3 = 77
  real(dp), parameter :: rydtohar = 0.5d0
  complex(dp), allocatable :: dudk_bra(:,:), dudk_ket(:,:), hpsi(:)
  complex(dp), allocatable :: vkb_save(:,:), aux(:,:)
  complex(dp) :: braket
  real(dp) :: kp_berry(3), kp_M_LC(3), kp_M_IC(3), tmp1(3), tmp2(3)
  integer :: ik, ibnd, jbnd, kk, ii, jj, occ
  real(dp) :: tmp(3), emin, emax
  ! index for the cross product
  integer :: ind(2,3)
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)

  ! set delta_k
  !delta_k = q_gipaw/2.d0/tpiba
  delta_k = q_gipaw/tpiba

  ! calculation inside atom centered spheres
  if (trim(dudk) == 'null') then
    call orbital_magnetization_spheres
    return
  endif

  ! compute the covariant derivatives
  call compute_dudk(dudk)

  ! allocate memory
  allocate (dudk_bra(npwx,nbnd), dudk_ket(npwx,nbnd), hpsi(npwx))

  ! zero the accumulators
#ifdef __PARA
  CALL mp_barrier()
#endif
  write(stdout,*)
  write(stdout,'(5X,''Computing the orbital magnetization (bohr mag/cell):'')')

  berry_curvature = 0.d0
  orb_magn_LC = 0.d0
  orb_magn_IC = 0.d0
  delta_M_bare = 0.d0
  delta_M_para = 0.d0
  delta_M_dia = 0.d0

  ! allocate the derivatives of the projectors
  allocate(becp(nkb,nbnd), dbecp(nkb,nbnd,3), paw_dbecp(paw_nkb,nbnd,3))
  allocate(vkb_save(npwx,nkb), aux(nkb,nbnd))   

  ! set the external potential
  call set_dvrs(dvrs, vrs, nrxxs, nspin)

  ! loop over k-points
  do ik = 1, nks
    call find_nbnd_occ(ik, occ, emin, emax)
    if (lgauss) occ = nbnd

    ! setup the hamiltonian
    current_k = ik
    current_spin = 1
    if (lsda) current_spin = isk(ik)
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(1:npw) = g2kin(1:npw) * tpiba2
    call get_buffer(evc, nwordwfc, iunwfc, ik)
    if (nkb > 0) then
      call init_us_2(npw, igk, xk(1,ik), vkb)
      call ccalbec(nkb, npwx, npw, nbnd, becp, vkb, evc)
    endif
    if (lda_plus_u) call davcio(swfcatom, nwordatwfc, iunsat, ik, -1)

    ! compute the diamagnetic terms
    call init_paw_2(npw, igk, xk(1,ik), paw_vkb)
    call ccalbec(paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)
    if (any(m_0 /= 0.d0))       call calc_delta_M_dia_nmr
    if (any(lambda_so /= 0.d0)) call calc_delta_M_dia_so

    call compute_dbecp
    call compute_paw_dbecp

    ! loop over the magnetization directions
    do kk =  1, 3
      ii = ind(1,kk)
      jj = ind(2,kk)

      ! read the bra and the ket
      call davcio(dudk_bra, 2*nwordwfc, iundudk1 + ii - 1, ik, -1)
      call davcio(dudk_ket, 2*nwordwfc, iundudk1 + jj - 1, ik, -1)

      ! compute the orbital magnetization
      kp_berry(kk) = 0.d0
      kp_M_IC(kk) = 0.d0
      kp_M_LC(kk) = 0.d0
      do ibnd = 1, occ
        ! IC term and Berry curvature
        braket = zdotc(npw, dudk_bra(1,ibnd), 1, dudk_ket(1,ibnd), 1)
        kp_berry(kk) = kp_berry(kk) + 2.d0*wg(ibnd,ik)*imag(braket)
        berry_curvature(kk) = berry_curvature(kk) + &
                              2.d0*wg(ibnd,ik)*imag(braket)
        kp_M_IC(kk) = kp_M_IC(kk) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
        orb_magn_IC(kk) = orb_magn_IC(kk) + &
                         2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)

        ! this is the LC term
        call h_psi(npwx, npw, 1, dudk_ket(1:npwx,ibnd), hpsi)
        braket = zdotc(npw, dudk_bra(1,ibnd), 1, hpsi, 1)
        kp_M_LC(kk) = kp_M_LC(kk) + wg(ibnd,ik)*imag(braket)
        orb_magn_LC(kk) = orb_magn_LC(kk) + wg(ibnd,ik)*imag(braket)

        call h_psi(npwx, npw, 1, dudk_bra(1:npwx,ibnd), hpsi)
        braket = zdotc(npw, dudk_ket(1,ibnd), 1, hpsi, 1)
        kp_M_LC(kk) = kp_M_LC(kk) - wg(ibnd,ik)*imag(braket)
        orb_magn_LC(kk) = orb_magn_LC(kk) - wg(ibnd,ik)*imag(braket)

      enddo

      ! compute the GIPAW corrections
      call calc_delta_M_bare
      if (any(lambda_so /= 0.d0)) call calc_delta_M_para_so
      if (any(m_0 /= 0.d0))       call calc_delta_M_para_nmr
    enddo ! kk
#ifdef __PARA
    call reduce(3, kp_berry)
    call reduce(3, kp_M_LC)
    call reduce(3, kp_M_IC)
#endif
    if (me_pool == root_pool) then
      write(*,'(''BC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_berry
      write(*,'(''LC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_M_LC*rydtohar
      write(*,'(''IC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_M_IC*rydtohar
      !!write(*,'(5X,''k-point:'',I4,4X,''pool:'',I4,4X,3(F10.4),/,5X,(9F12.6))') &
      !!      ik, my_pool_id+1, xk(:,ik), kp_berry, kp_M_LC*rydtohar, kp_M_IC*rydtohar
    endif
  enddo ! ik
  write(stdout,*)

#ifdef __PARA
  call reduce(3, orb_magn_LC)
  call reduce(3, orb_magn_IC)
  call reduce(3, berry_curvature)
  ! no reduction for delta_M_bare and delta_M_para and delta_M_dia
#endif

#ifdef __PARA
  call poolreduce(3, orb_magn_LC)
  call poolreduce(3, orb_magn_IC)
  call poolreduce(3, berry_curvature)
  call poolreduce(3, delta_M_bare)
  call poolreduce(3, delta_M_para)
  call poolreduce(3, delta_M_dia)
#endif

  ! close files
  close(unit=iundudk1, status='keep')
  close(unit=iundudk2, status='keep')
  close(unit=iundudk3, status='keep')

  ! convert to hartree and sum up all terms
  ! <\uwg> conversion no more needed
  !ef = ef * rydtohar
  !orb_magn_LC = orb_magn_LC * rydtohar
  !orb_magn_IC = orb_magn_IC * rydtohar
  !delta_M_bare = delta_M_bare * rydtohar
  !delta_M_para = delta_M_para * rydtohar
  !delta_M_dia = delta_M_dia * rydtohar
  
  orb_magn_tot = orb_magn_LC + orb_magn_IC + &
                 delta_M_bare + delta_M_dia + delta_M_para

  ! print results
  write(stdout,'(5X,''SPIN-ORBIT:'')')
  write(stdout,'(5X,''lambda_so          = '',3(F14.6))') lambda_so
  write(stdout,*)

  write(stdout,'(5X,''NUCLEAR DIPOLE ON ATOM'',I4,'':'')') m_0_atom
  write(stdout,'(5X,''m_0                = '',3(F14.6))') m_0
  write(stdout,*)

  write(stdout,'(5X,''Berry curvature    = '',3(F14.6))') berry_curvature
  write(stdout,'(5X,''Fermi energy       = '',F14.6,'' rydberg'')') ef
  write(stdout,*)

  if (iverbosity > 0) then
    write(stdout,'(5X,''(without Berry curvature term)'')')
    write(stdout,'(5X,''M_LC               = '',3(F14.6))') orb_magn_LC
    write(stdout,'(5X,''M_IC               = '',3(F14.6))') orb_magn_IC
    write(stdout,'(5X,''Delta_M_bare       = '',3(F14.6))') delta_M_bare
    !!print*, delta_M_bare
    write(stdout,'(5X,''Delta_M_para       = '',3(F14.6))') delta_M_para
    !!print*, delta_M_para
    write(stdout,'(5X,''Delta_M_dia        = '',3(F14.6))') delta_M_dia
    !!print*, delta_M_dia
    write(stdout,'(5X,''M_tot              = '',3(F14.6))') orb_magn_tot
    write(stdout,*)
    write(stdout,'(5X,''(with Berry curvature term)'')')
  endif
  orb_magn_LC = orb_magn_LC - ef*berry_curvature
  orb_magn_IC = orb_magn_IC - ef*berry_curvature
  orb_magn_tot = orb_magn_tot - 2.d0*ef*berry_curvature
  write(stdout,'(5X,''M_LC               = '',3(F14.6))') orb_magn_LC
  write(stdout,'(5X,''M_IC               = '',3(F14.6))') orb_magn_IC
  write(stdout,'(5X,''Delta_M            = '',3(F14.6))') &
        delta_M_bare + delta_M_para + delta_M_dia
  write(stdout,'(5X,''M_tot              = '',3(F14.6))') orb_magn_tot
  
  ! free memory
  deallocate( dudk_bra, dudk_ket, hpsi, becp )

  ! close files
  close(unit=iundudk1)
  close(unit=iundudk2)
  close(unit=iundudk3)

  ! go on, reporting the g-tensor
  if (any(lambda_so /= 0.d0)) call calc_g_tensor(orb_magn_tot)

  ! go on, reporting the chemical shift
  if (any(m_0 /= 0.d0)) call calc_chemical_shift(orb_magn_tot, nmr_shift_core)

  CONTAINS
    !------------------------------------------------------------------
    ! derivative of the beta's
    !------------------------------------------------------------------
    SUBROUTINE compute_dbecp
    USE cell_base, ONLY : tpiba
    USE g_tensor_module, ONLY : init_us_2_no_phase
    implicit none
    integer :: ipol, sig
    real(dp) :: kq(3)
    if (nkb == 0) return
    vkb_save = vkb
    do ipol = 1, 3
      dbecp(:,:,ipol) = (0.d0,0.d0)
      do sig = -1, 1, 2
        kq(:) = xk(:,ik)
        kq(ipol) = kq(ipol) + sig * delta_k
        call init_us_2_no_phase(npw, igk, kq, vkb)
        call ccalbec(nkb, npwx, npw, nbnd, aux, vkb, evc)
        dbecp(:,:,ipol) = dbecp(:,:,ipol) + &
                          0.5d0*sig/(delta_k*tpiba) * aux(:,:)
      enddo
    enddo
    vkb = vkb_save
    END SUBROUTINE compute_dbecp


    !------------------------------------------------------------------
    ! GIPAW correction (delta_M_bare)
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_bare
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq, indv, nhtol
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: nbs1, nbs2, l1, l2
    if (nkb == 0) return
    tmp = (0.d0,0.d0)
    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if ( ityp(na) == nt ) then
          do jh = 1, nh(nt)
            jkb = ijkb0 + jh
            nbs1 = indv(jh,nt)
            l1 = nhtol(jh,nt)
            do ih = 1, nh(nt)
              ikb = ijkb0 + ih
              nbs2 = indv(ih,nt)
              l2 = nhtol(ih,nt)
              if (l1 /= l2) cycle
              do ibnd = 1, occ
                becp_product = conjg(dbecp(ikb,ibnd,ii)) * dbecp(jkb,ibnd,jj)
                tmp = tmp + wg(ibnd,ik) * deeq(ih,jh,na,current_spin) * &
                            becp_product
              enddo
            enddo
          enddo
          ijkb0 = ijkb0 + nh(nt)
        endif
      enddo
    enddo
    !!PRINT*, mpime, kk, -2.d0*tmp
    ! check the sign and real or imag!!
    delta_M_bare(kk) = delta_M_bare(kk) - 2.d0*imag(tmp)
    END SUBROUTINE calc_delta_M_bare

      
    !------------------------------------------------------------------
    ! derivative of GIPAW projectors
    !------------------------------------------------------------------
    SUBROUTINE compute_paw_dbecp
    USE cell_base, ONLY : tpiba
    USE g_tensor_module, ONLY : init_paw_2_no_phase
    implicit none
    integer :: ipol, sig
    real(dp) :: kq(3)
    if (paw_nkb == 0) return
    do ipol = 1, 3
      paw_dbecp(:,:,ipol) = (0.d0,0.d0)
      do sig = -1, 1, 2
        kq(:) = xk(:,ik)
        kq(ipol) = kq(ipol) + sig * delta_k  
        !!!!call init_paw_2(npw, igk, kq, paw_vkb)
        call init_paw_2_no_phase(npw, igk, kq, paw_vkb)
        call ccalbec(paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)
        paw_dbecp(:,:,ipol) = paw_dbecp(:,:,ipol) + &
                              0.5d0*sig/(delta_k*tpiba) * paw_becp(:,:)
      enddo
    enddo
    END SUBROUTINE compute_paw_dbecp


    !------------------------------------------------------------------
    ! GIPAW correction (Delta_M_para), SO case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_para_so
    USE constants,  ONLY : e2
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq  
    USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE g_tensor_module, ONLY : a2gp4, a2gp8, radial_integral_paramagnetic_so
    USE g_tensor_module, ONLY : lx, ly, lz
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
    real(dp) :: sigma
    if (paw_nkb == 0) return  
    sigma = 1.d0; if (current_spin == 2) sigma = -1.d0
    tmp = (0.d0,0.d0)
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
              do ibnd = 1, occ
                becp_product = conjg(paw_dbecp(ikb,ibnd,ii)) * &
                                     paw_dbecp(jkb,ibnd,jj)
                tmp = tmp + lambda_so(1) * a2gp8 * &
                      sigma * lx(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp + lambda_so(2) * a2gp8 * &
                      sigma * ly(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp + lambda_so(3) * a2gp8 * &
                      sigma * lz(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
              enddo   ! ibnd
            enddo   ! jh
          enddo   ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif   ! if (ityp...
      enddo   ! na
    enddo  ! nt
    !!PRINT*, kk, tmp
    ! check the sign and real or imag!!
    delta_M_para(kk) = delta_M_para(kk) - 2.d0*imag(tmp)
    END SUBROUTINE calc_delta_M_para_so


    !------------------------------------------------------------------
    ! GIPAW correction (delta_M_dia), SO case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_dia_so
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE parameters, ONLY : lmaxx
    USE constants,  ONLY : pi
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq, ap
    USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE g_tensor_module, ONLY : gprime, a2gp4, a2gp8, radial_integral_diamagnetic_so
    USE g_tensor_module, ONLY : lx, ly, lz, alpha
    IMPLICIT NONE
    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, nrc, ijkb0
    complex(dp) :: bec_product
    complex(dp), allocatable :: dia_corr(:)
    real(dp) :: diamagnetic_tensor(3,3), sigma

    if (paw_nkb == 0) return  
    sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

    allocate(dia_corr(lmaxx**2))
    dia_corr = 0.0_dp
    diamagnetic_tensor = 0.d0

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
              lm2 = m2 + l2**2

              do ibnd = 1, nbnd
                bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
                !!!PRINT*, ikb,jkb,bec_product
                !<apsi> s/non-trace-zero component
                ! 2/3 to separate the non-trace vanishing component
                ! 1/(2c^2) from the equation (59) in PM-PRB
                if ( l1 == l2 .AND. m1 == m2 ) then
                  diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                  diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                  diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                endif

                ! 2/3 to separate the non-trace vanishing component
                do lm = 5, 9
                  dia_corr(lm) = dia_corr(lm) + bec_product / 3.0_dp &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * ap(lm,lm1,lm2) * wg(ibnd,ik) * alpha**2.d0
                enddo

              enddo  ! ibnd
            enddo  ! jh
          enddo  ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif
      enddo  ! na
    enddo  ! nt
    
    !  transform in cartesian coordinates
    dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)
    diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
      + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
      - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
      + dia_corr(5) * 2.0_dp
    diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) &
      +  dia_corr(9) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
    diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) &
      - dia_corr(6) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
    diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) &
      - dia_corr(7) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)
    deallocate (dia_corr)

    !!!WRITE(*,'(''dia:'',/,3(3(E14.4,2X),/))') diamagnetic_tensor
    !!! </uwg>  (gprime/4.d0) correction factor added, compare GIPAW  
    delta_M_dia = delta_M_dia + (gprime/4.d0)*0.5d0*sigma &
                                 *matmul(diamagnetic_tensor, lambda_so)
   
    END SUBROUTINE calc_delta_M_dia_so


    !------------------------------------------------------------------
    ! GIPAW correction (Delta_M_para), NMR case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_para_nmr
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq
    USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE nmr_mod,    ONLY : m_0, m_0_atom, fine_struct
    USE g_tensor_module, ONLY : radial_integral_paramagnetic_nmr
    USE g_tensor_module, ONLY : lx, ly, lz
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
    !!RETURN
    if (paw_nkb == 0) return
    tmp = (0.d0,0.d0)
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
              !!if (l1 == 0) cycle
              if (na /= m_0_atom) cycle
              do ibnd = 1, occ
                becp_product = conjg(paw_dbecp(ikb,ibnd,ii)) * &
                                     paw_dbecp(jkb,ibnd,jj)

                tmp = tmp - fine_struct**2 * &
                      m_0(1) * lx(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_nmr(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp - fine_struct**2 * &
                      m_0(2) * ly(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_nmr(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp - fine_struct**2 * &
                      m_0(3)* lz(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_nmr(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
              enddo   ! ibnd
            enddo   ! jh
          enddo   ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif   ! if (ityp...
      enddo   ! na
    enddo  ! nt
    !!PRINT*, kk, tmp
    ! check the sign and real or imag!!
    delta_M_para(kk) = delta_M_para(kk) - 2.d0*imag(tmp)/sqrt(2.d0)
    END SUBROUTINE calc_delta_M_para_nmr


    !------------------------------------------------------------------
    ! GIPAW correction (Delta_M_dia), NMR case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_dia_nmr
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE constants,  ONLY : pi
    USE parameters, ONLY : lmaxx
    USE uspp,       ONLY : deeq, ap  
    USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE nmr_mod,    ONLY : m_0, m_0_atom, fine_struct
    USE g_tensor_module, ONLY : radial_integral_diamagnetic_nmr
    USE g_tensor_module, ONLY : lx, ly, lz, alpha
    implicit none
    complex(dp) :: tmp, bec_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2, lm
    complex(dp), allocatable :: dia_corr(:)
    real(dp) :: diamagnetic_tensor(3,3)
   
    allocate (dia_corr(lmaxx**2))
    diamagnetic_tensor = 0.d0
    dia_corr = 0.0_dp

    if (paw_nkb == 0) return  
    tmp = (0.d0,0.d0)
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
              lm2=m2+l2**2 
              if (na /= m_0_atom) cycle
              do ibnd = 1, nbnd
                bec_product = paw_becp(jkb,ibnd) * CONJG( paw_becp(ikb,ibnd) )
                !<apsi> s/non-trace-zero component
                ! 2/3 to separate the non-trace vanishing component
                ! 1/(2c^2) from the equation (59) in PM-PRB
                IF ( l1 == l2 .AND. m1 == m2 ) THEN
                  diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_nmr(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                  diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_nmr(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                  diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_nmr(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                END IF
                ! 2/3 to separate the non-trace vanishing component
                do lm = 5, 9
                   dia_corr(lm) =  dia_corr(lm) + bec_product / 3.0_dp &
                              * radial_integral_diamagnetic_nmr(nbs1,nbs2,nt) &
                              * ap(lm,lm1,lm2) * wg(ibnd,ik) * alpha ** 2 &
                              / 2.0_dp
                enddo
              enddo  ! ibnd
            enddo  ! jh
          enddo  ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif  
      enddo  ! na
    enddo  ! bt
    
    !
    !  transform in cartesian coordinates
    !
    
    dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)
    
    diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
         + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
         - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
         + dia_corr(5) * 2.0_dp
    diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) &
         +  dia_corr(9) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
    diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) &
         - dia_corr(6) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
    diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) &
         - dia_corr(7) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)
    deallocate(dia_corr)

    delta_M_dia = delta_M_dia - sqrt(2.d0)*matmul(diamagnetic_tensor, m_0)

    END SUBROUTINE calc_delta_M_dia_nmr


  END SUBROUTINE calc_orbital_magnetization

!-----------------------------------------------------------------------
END MODULE orbital_magnetization
!-----------------------------------------------------------------------

