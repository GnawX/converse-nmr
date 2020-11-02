
!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!-----------------------------------------------------------------------
MODULE g_tensor_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for GIPAW calculations
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : a0_to_cm => bohr_radius_cm
  USE parameters, ONLY: npk, ntypx, lmaxx, nbrx
  IMPLICIT NONE
  SAVE
 
  ! strenght of spin-orbit
  REAL(DP) :: lambda_so(3)
 
  ! alpha
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp
  ! g_e
  REAL(DP), PARAMETER :: g_e = 2.0023192778_dp
  ! g'
  REAL(DP), PARAMETER :: gprime = 2.0046385556_dp
  ! (alpha^2 g')/4
  REAL(DP), PARAMETER :: a2gp4 = alpha*alpha*gprime/4_dp
  ! (alpha^2 g')/8
  REAL(DP), PARAMETER :: a2gp8 = alpha*alpha*gprime/8_dp

  ! avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142e23_dp

  ! format for a rank-2 tensor
  CHARACTER(*), PARAMETER :: tens_fmt = '(3(5X,3(F14.4,2X)/))'
  
  CHARACTER(256) :: file_reconstruction ( ntypx )
  LOGICAL :: read_recon_in_paratec_fmt
  LOGICAL :: broken_time_reversal
  REAL(dp) :: q_conv_thr, q_gipaw
  REAL(dp) :: rc(ntypx,0:lmaxx)
  !!COMPLEX(dp), ALLOCATABLE :: paw_becp2 ( :, : )
  REAL(dp), ALLOCATABLE, DIMENSION ( :, : ) :: lx, ly, lz
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic_nmr(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic_nmr(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_rmc(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_L(:,:,:)
  
  REAL(dp) :: r_mt ! muffin-tin radius
  INTEGER, ALLOCATABLE :: pointlist(:,:), pointnum(:)
  REAL(dp), ALLOCATABLE :: distlist(:,:,:)

  REAL(dp) :: nmr_shift_core(ntypx)

  LOGICAL, ALLOCATABLE :: vloc_present ( : )
  REAL(dp), ALLOCATABLE :: gipaw_ae_vloc ( :, : ), gipaw_ps_vloc ( :, : )
  
  !!PUBLIC :: read_recon_paratec
  !<apsi>

CONTAINS
  
  !-----------------------------------------------------------------------
  ! GIPAW setup
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_setup
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE ions_base,     ONLY : tau, nat, ntyp => nsp
    USE atom,          ONLY : nlcc, mesh, msh, r, rab
    USE wvfct,         ONLY : nbnd, nbndx, et, wg, npwx
    USE lsda_mod,      ONLY : nspin, lsda
    USE scf,           ONLY : vr, vrs, vltot, rho, rho_core
    USE gvect,         ONLY : nrxx, ngm
    USE gsmooth,       ONLY : doublegrid
    USE klist,         ONLY : xk, degauss, ngauss, nks, nelec
    USE constants,     ONLY : degspin, pi
    USE paw,           ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp, &
                              read_recon, read_recon_paratec
    IMPLICIT none
    integer :: ik, nt, ibnd
    logical :: nlcc_any
    real(dp) :: emin, emax
    
    !<apsi>
    integer :: il, lm, l, m, lm1, lm2, m1, m2, abs_m1, abs_m2
    integer :: sign_m1, sign_m2, il1, il2, l1, l2, j, kkpsi, nrc
    real(dp) :: alpha_lm, beta_lm
    integer, allocatable :: lm2l ( : ), lm2m ( : )
    integer :: kkpsi_max
    real(dp), allocatable :: work(:), kinetic_aephi(:), kinetic_psphi(:)
    real(dp), allocatable :: aephi_dvloc_dr(:), psphi_dvloc_dr(:)
    integer :: core_orb
    real(dp) :: occupation, integral
    real(dp) :: mysum1 ( 3, lmaxx ) !TMPTMPTMP
    real(dp) :: mysum2 ( 3, 1:lmaxx ) !TMPTMPTMP
    logical :: vloc_set
    !</apsi>
    
    call start_clock ('gipaw_setup')
    
    ! initialize pseudopotentials
    call init_us_1

    !<apsi>
    
    ! Read in qe format
    IF ( .NOT. ALLOCATED ( paw_recon ) ) THEN
       ALLOCATE ( paw_recon(ntyp) )
       paw_recon(:)%gipaw_data_in_upf_file = .FALSE.
    END IF
    DO nt = 1, ntyp
       IF ( read_recon_in_paratec_fmt ) THEN
          
          ! Read in paratec format
          CALL read_recon_paratec ( file_reconstruction(nt), nt, &
               paw_recon(nt), vloc_set )
          IF ( .NOT. vloc_set ) THEN
             CALL errore ( "gipaw_setup", &
                  "no local potential set in read_recon_paratec", -1 )
             !stop
          END IF
          
       ELSE
          CALL read_recon ( file_reconstruction(nt), nt, paw_recon(nt) )
       END IF
    END DO
    
    ! initialize paw
    do nt = 1, ntyp
       do il = 1, paw_recon(nt)%paw_nbeta
          IF ( paw_recon(nt)%psphi(il)%label%rc < -0.99_dp ) THEN
             rc(nt,paw_recon(nt)%psphi(il)%label%l) = 1.6_dp
             rc(nt,paw_recon(nt)%aephi(il)%label%l) = 1.6_dp
             paw_recon(nt)%psphi(il)%label%rc &
                  = rc(nt,paw_recon(nt)%psphi(il)%label%l)
             paw_recon(nt)%aephi(il)%label%rc &
                  = rc(nt,paw_recon(nt)%aephi(il)%label%l)
          ELSE
             rc(nt,paw_recon(nt)%psphi(il)%label%l) &
                  = paw_recon(nt)%psphi(il)%label%rc
             rc(nt,paw_recon(nt)%aephi(il)%label%l) &
                  = paw_recon(nt)%aephi(il)%label%rc
          END IF
       enddo
    enddo
    
    call init_paw_1()
    allocate ( paw_vkb(npwx,paw_nkb) )
    allocate ( paw_becp(paw_nkb,nbndx) )
    
    !<apsi>
    allocate ( radial_integral_diamagnetic_so(2*nbrx,2*nbrx,ntypx) )
    allocate ( radial_integral_paramagnetic_so(2*nbrx,2*nbrx,ntypx) )
    allocate ( radial_integral_rmc(2*nbrx,2*nbrx,ntypx) )
    radial_integral_diamagnetic_so = 0.0_dp
    radial_integral_paramagnetic_so = 0.0_dp
    radial_integral_rmc = 0.0_dp

    allocate ( radial_integral_diamagnetic_nmr(2*nbrx,2*nbrx,ntypx) )
    allocate ( radial_integral_paramagnetic_nmr(2*nbrx,2*nbrx,ntypx) )
    allocate ( radial_integral_L(2*nbrx,2*nbrx,ntypx) )
    radial_integral_diamagnetic_nmr = 0.0_dp
    radial_integral_paramagnetic_nmr = 0.0_dp
    radial_integral_L = 0.0_dp
  
    do nt=1, ntyp
       
       do il1=1, paw_recon(nt)%paw_nbeta
          if (il1 > 2*nbrx) STOP 'AHIA'
          l1 = paw_recon(nt)%psphi(il1)%label%l
          kkpsi = paw_recon(nt)%aephi(il1)%kkpsi
          
          nrc = paw_recon(nt)%psphi(il1)%label%nrc
          
          allocate ( work(kkpsi) )

          do il2 = 1, paw_recon(nt)%paw_nbeta
             if (il2 > 2*nbrx) STOP 'AHIA'
             l2 = paw_recon(nt)%psphi(il2)%label%l
             
             IF ( l1 /= l2 ) CYCLE
             
             !
             ! NMR shielding, diamagnetic
             !
             do j = 1, nrc
                work(j) = (paw_recon(nt)%aephi(il1)%psi(j)*paw_recon(nt)%aephi(il2)%psi(j)-&
                     paw_recon(nt)%psphi(il1)%psi(j)*paw_recon(nt)%psphi(il2)%psi(j))/r(j,nt)
             enddo
             
             CALL simpson( nrc, work, rab(:nrc,nt), &
                  radial_integral_diamagnetic_nmr(il1,il2,nt) )
             if (il1 >= il2) &
               write(stdout,'(''DIA (NMR) :'',I3,4X,2I3,2X,F14.4)') nt, l1, l2, &
                     radial_integral_diamagnetic_nmr(il1,il2,nt)
             
             !
             ! NMR shielding, paramagnetic
             !
             do j = 1, nrc
                work(j) = &
                     ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / r(j,nt) ** 3
             end do
             call simpson( nrc, work, rab(:,nt), &
                  radial_integral_paramagnetic_nmr(il1,il2,nt) )
             if (il1 >= il2) &
               write(stdout,'(''PARA (NMR):'',I3,4X,2I3,2X,F14.4)') nt, l1, l2, &
                     radial_integral_paramagnetic_nmr(il1,il2,nt)

             !
             ! Angular moment (AE)
             !
             do j = 1, nrc
                work(j) = &
                     ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) )
             end do        
             call simpson( nrc, work, rab(:,nt), radial_integral_L(il1,il2,nt) )
             if (il1 >= il2) &
               write(stdout,'(''ANG. MOM. :'',I3,4X,2I3,2X,F14.4)') nt, l1, l2, &
                     radial_integral_L(il1,il2,nt)

             ! Calculate the radial integral only if the radial potential
             !    is present
             IF ( .NOT. paw_recon(nt)%vloc_present ) CYCLE
             if ( all(lambda_so == 0.d0)) CYCLE
 
             !
             ! g tensor, relativistic mass correction
             !
             ALLOCATE ( kinetic_aephi ( kkpsi ), kinetic_psphi ( kkpsi ) )
             CALL radial_kinetic_energy ( l2, r(:nrc,nt), &
                  paw_recon(nt)%aephi(il2)%psi(:nrc), kinetic_aephi(:nrc) )
             CALL radial_kinetic_energy ( l2, r(:nrc,nt), &
                  paw_recon(nt)%psphi(il2)%psi(:nrc), kinetic_psphi(:nrc) )
             
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * kinetic_aephi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * kinetic_psphi(j) )
             end do
             DEALLOCATE ( kinetic_aephi, kinetic_psphi )
             
             CALL simpson ( nrc, work, rab(:,nt), &
                  radial_integral_rmc(il1,il2,nt) )
             
             ALLOCATE ( aephi_dvloc_dr ( nrc ), psphi_dvloc_dr ( nrc ) )
             
             CALL radial_derivative ( r(:nrc,nt), &
                  paw_recon(nt)%gipaw_ae_vloc(:nrc), &
                  aephi_dvloc_dr(:nrc) )
             CALL radial_derivative ( r(:nrc,nt), &
                  paw_recon(nt)%gipaw_ps_vloc(:nrc), &
                  psphi_dvloc_dr ( :nrc ) )
             
             !
             ! g tensor, diamagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     * r(j,nt)
             end do
             
             call simpson( nrc, work, rab(:,nt), &
                  radial_integral_diamagnetic_so(il1,il2,nt) )
             if (il1 >= il2) &
               write(stdout,'(''DIA (SO)  :'',I3,4X,2I3,2X,F14.4)') nt, l1, l2, &
                     radial_integral_diamagnetic_so(il1,il2,nt)
             
             !
             ! g tensor, paramagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / r(j,nt)
             end do
             
             call simpson( nrc,work,rab(:,nt), &
                  radial_integral_paramagnetic_so(il1,il2,nt) )
             if (il1 >= il2) &
               write(stdout,'(''PARA (SO) :'',I3,4X,2I3,2X,F14.4)') nt, l1, l2, &
                     radial_integral_paramagnetic_so(il1,il2,nt)

             DEALLOCATE ( aephi_dvloc_dr, psphi_dvloc_dr )
             
          enddo
          
          DEALLOCATE ( work )
          
       enddo
    enddo
    
    !</apsi>
    
    ! terms for paramagnetic
    
    !<apsi>
    allocate ( lx ( lmaxx**2, lmaxx**2 ) )
    allocate ( ly ( lmaxx**2, lmaxx**2 ) )
    allocate ( lz ( lmaxx**2, lmaxx**2 ) )
    
    allocate ( lm2l ( lmaxx**2 ), lm2m ( lmaxx**2 ) )
    
    lm = 0
    do l = 0, lmaxx - 1
       do m = 0, l
          lm = lm + 1
          lm2l ( lm ) = l
          lm2m ( lm ) = m
          if ( m /= 0 ) then
             lm = lm + 1
             lm2l ( lm ) = l
             lm2m ( lm ) = - m
          end if
       end do
    end do
    
    lx = 0.0_dp
    ly = 0.0_dp
    lz = 0.0_dp
    do lm2 = 1, lmaxx**2
       do lm1 = 1, lmaxx**2
          ! same L?
          if (lm2l(lm1) /= lm2l(lm2)) cycle

          l = lm2l(lm1)
          m1 = lm2m(lm1)
          m2 = lm2m(lm2)

          ! L_x, L_y
          if ( m2 == 0 ) then
             if ( m1 == -1 ) then
                lx ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             else if ( m1 == +1 ) then
                ly ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             end if
          else if ( m1 == 0 ) then
             if ( m2 == -1 ) then
                lx ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             else if ( m2 == +1 ) then
                ly ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             end if
          else
             abs_m1 = abs ( m1 )
             abs_m2 = abs ( m2 )
             sign_m1 = sign ( 1, m1 )
             sign_m2 = sign ( 1, m2 )
             alpha_lm = sqrt(real(l*(l+1)-abs_m2*(abs_m2+1),dp))
             beta_lm  = sqrt(real(l*(l+1)-abs_m2*(abs_m2-1),dp))
             if ( abs_m1 == abs_m2 + 1 ) then
                lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * alpha_lm / 4.0_dp
                ly ( lm1, lm2 ) = ( sign_m2 + sign_m1 ) * alpha_lm / 4.0_dp &
                     / sign_m2
             else if ( abs_m1 == abs_m2 - 1 ) then
                lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * beta_lm / 4.0_dp
                ly ( lm1, lm2 ) =-( sign_m2 + sign_m1 ) * beta_lm / 4.0_dp &
                     / sign_m2
             end if
          end if
          
          ! L_z
          if ( m1 == - m2 ) then
             lz ( lm1, lm2 ) = - m2
          end if
          
       end do
    end do
    
    deallocate ( lm2l, lm2m )

    ! Compute the shift due to core orbitals
    DO nt = 1, ntyp
       IF ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) CYCLE
       ALLOCATE ( work(mesh(nt)) )
       nmr_shift_core(nt) = 0.0
       DO core_orb = 1, paw_recon(nt)%gipaw_ncore_orbital
          DO j = 1, SIZE(work)
             work(j) = paw_recon(nt)%gipaw_core_orbital(j,core_orb) ** 2 &
                  / r(j,nt)
          END DO
          CALL simpson( SIZE(work), work, rab(:,nt), &
               integral )
          occupation = 2 * ( &
               2 * paw_recon(nt)%gipaw_core_orbital_l(core_orb) + 1 )
          nmr_shift_core(nt) = nmr_shift_core(nt) + occupation * integral
       END DO
       DEALLOCATE ( work )
       nmr_shift_core(nt) = nmr_shift_core(nt) * 17.75045395
    END DO

    call stop_clock('gipaw_setup')

    write(stdout, '(5X,''total number of PAW projectors: '',I3)') paw_nkb

  CONTAINS
    
    SUBROUTINE radial_kinetic_energy ( l, rdata, ydata, kin_ydata )
      USE splinelib
      
      INTEGER, INTENT ( IN ) :: l
      REAL(dp), INTENT ( IN ) :: rdata ( : ), ydata ( : )
      REAL(dp), INTENT ( OUT ) :: kin_ydata ( : )
      
      REAL(dp) :: d1
      
      d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
      CALL spline ( rdata, ydata, 0.0_dp, d1, kin_ydata )
      
      kin_ydata = - kin_ydata + l*(l+1) * ydata / rdata ** 2
      
    END SUBROUTINE radial_kinetic_energy
    
    SUBROUTINE radial_derivative ( rdata, ydata, dydata_dr )
      USE splinelib
      
      ! ydata passed as y * r
      REAL(dp), INTENT ( IN ) :: rdata ( : ), ydata ( : )
      REAL(dp), INTENT ( OUT ) :: dydata_dr ( : )
      
      INTEGER :: j
      REAL(dp) :: d1, tab_d2y ( SIZE ( ydata ) )
      
      d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
      CALL spline ( rdata, ydata, 0.0_dp, d1, tab_d2y )
      
      DO j = 1, SIZE ( ydata )
         dydata_dr ( j ) = &
              ( splint_deriv ( rdata, ydata, tab_d2y, rdata ( j ) ) &
              - ydata ( j ) / rdata ( j ) ) / rdata ( j )
      END DO
      
    END SUBROUTINE radial_derivative
    
  END SUBROUTINE gipaw_setup



  !----------------------------------------------------------------------
  SUBROUTINE init_us_2_no_phase (npw_, igk_, q_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE constants,  ONLY : tpi
  USE gvect,      ONLY : eigts1, eigts2, eigts3, ig1, ig2, ig3, g
  USE wvfct,      ONLY : npw, npwx, igk
  USE us,         ONLY : nqx, dq, tab, tab_d2y, spline_ps
  USE splinelib
  USE uspp,       ONLY : nkb, vkb, nhtol, nhtolm, indv
  USE uspp_param, ONLY : lmaxkb, nbeta, nhm, nh
  !
  implicit none
  !
  integer :: npw_, igk_ (npw_)
  ! input: number of PW's
  ! input: indices of q+G
  real(DP) :: q_(3)
  ! input: q vector
  complex(DP) :: vkb_ (npwx, nkb)
  ! output: beta functions
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, l, lm, na, nt, nb, ih, jkb

  real(DP) :: px, ux, vx, wx, arg
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk(:)

  real(DP), allocatable :: xdata(:)
  integer :: iq
  !
  if (lmaxkb.lt.0) return
  call start_clock ('init_us_2')
  allocate (vkb1( npw_,nhm))    
  allocate (  sk( npw_))    
  allocate (  qg( npw_))    
  allocate (  vq( npw_))    
  allocate ( ylm( npw_, (lmaxkb + 1) **2))    
  allocate (  gk( 3, npw_))    
  !
  do ig = 1, npw_
     gk (1,ig) = q_(1) + g(1, igk_(ig) )
     gk (2,ig) = q_(2) + g(2, igk_(ig) )
     gk (3,ig) = q_(3) + g(3, igk_(ig) )
     qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  call ylmr2 ((lmaxkb+1)**2, npw_, gk, qg, ylm)
  !
  ! set now qg=|q+G| in atomic units
  !
  do ig = 1, npw_
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  jkb = 0
  do nt = 1, ntyp
     ! calculate beta in G-space using an interpolation table
     do nb = 1, nbeta (nt)
        do ig = 1, npw_
           if (spline_ps) then
             vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
           else
             px = qg (ig) / dq - int (qg (ig) / dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = INT( qg (ig) / dq ) + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             vq (ig) = tab (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                       tab (i1, nb, nt) * px * vx * wx / 2.d0 - &
                       tab (i2, nb, nt) * px * ux * wx / 2.d0 + &
                       tab (i3, nb, nt) * px * ux * vx / 6.d0
           endif
        enddo
        ! add spherical harmonic part
        do ih = 1, nh (nt)
           if (nb.eq.indv (ih, nt) ) then
              l = nhtol (ih, nt)
              lm =nhtolm (ih, nt)
              do ig = 1, npw_
                 vkb1 (ig,ih) = ylm (ig, lm) * vq (ig)
              enddo
           endif
        enddo
     enddo
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if (ityp (na) .eq.nt) then
           arg = (q_(1) * tau (1, na) + &
                  q_(2) * tau (2, na) + &
                  q_(3) * tau (3, na) ) * tpi
           phase = CMPLX (cos (arg), - sin (arg) )
           do ig = 1, npw_
              sk (ig) = eigts1 (ig1(igk_(ig)), na) * &
                        eigts2 (ig2(igk_(ig)), na) * &
                        eigts3 (ig3(igk_(ig)), na)
           enddo
           do ih = 1, nh (nt)
              jkb = jkb + 1
              pref = (0.d0, -1.d0) **nhtol (ih, nt) !* phase
              do ig = 1, npw_
                 vkb_(ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)
  deallocate (vkb1)

  call stop_clock ('init_us_2')
  return
  END SUBROUTINE init_us_2_no_phase



  !----------------------------------------------------------------------
  SUBROUTINE init_paw_2_no_phase (npw_, igk_, q_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates paw_beta functions (paw projectors), with
  !   structure factor, for all atoms, in reciprocal space
  !
  USE kinds ,     ONLY : dp
  USE constants , ONLY : tpi
  USE wvfct ,     ONLY : npwx
  USE cell_base , ONLY : tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE gvect ,     ONLY : eigts1, eigts2, eigts3, g, ig1, ig2, ig3 
  USE us,         ONLY : dq
  USE paw,        ONLY : paw_nkb, paw_recon, paw_lmaxkb
  USE us,         ONLY : nqx, dq, spline_ps
  USE splinelib
  !
  implicit none
  !
  integer :: npw_, igk_ (npw_)
  ! input: number of PW's
  ! input: indices of q+G
  real(DP) :: q_(3)
  ! input: q vector
  complex(DP) :: vkb_ (npwx, paw_nkb)
  ! output: beta functions
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, l, lm, na, nt, nb, ih, jkb

  real(DP) :: px, ux, vx, wx, arg
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk(:)

  real(DP), allocatable :: xdata(:)
  integer :: iq
  !
  !
  if (paw_lmaxkb.lt.0) return
  call start_clock ('init_paw_2')
  allocate (  sk( npw_))    
  allocate (  qg( npw_))    
  allocate (  vq( npw_))    
  allocate ( ylm( npw_, (paw_lmaxkb + 1) **2))    
  allocate (  gk( 3, npw_))    
  !

  do ig = 1, npw_
     gk (1,ig) = q_(1) + g(1, igk_(ig) )
     gk (2,ig) = q_(2) + g(2, igk_(ig) )
     gk (3,ig) = q_(3) + g(3, igk_(ig) )
     qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  call ylmr2 ((paw_lmaxkb+1)**2, npw_, gk, qg, ylm)
  !
  ! set now qg=|q+G| in atomic units
  !
  do ig = 1, npw_
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  jkb = 0
  do nt = 1, ntyp
     allocate ( vkb1(npw_,paw_recon(nt)%paw_nh) )
     
     ! calculate beta in G-space using an interpolation table
     do nb = 1, paw_recon(nt)%paw_nbeta
        do ig = 1, npw_
           if (spline_ps) then
             vq(ig) = splint ( xdata, paw_recon(nt)%paw_tab(:,nb), &
                  paw_recon(nt)%paw_tab_d2y(:,nb), qg(ig) )
           else
             px = qg (ig) / dq - int (qg (ig) / dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = qg (ig) / dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             vq (ig) = paw_recon(nt)%paw_tab(i0,nb) * ux * vx * wx / 6.d0 + &
                       paw_recon(nt)%paw_tab(i1,nb) * px * vx * wx / 2.d0 - &
                       paw_recon(nt)%paw_tab(i2,nb) * px * ux * wx / 2.d0 + &
                       paw_recon(nt)%paw_tab(i3,nb) * px * ux * vx / 6.d0
           endif
        enddo
        ! add spherical harmonic part
        do ih = 1, paw_recon(nt)%paw_nh
           if ( nb == paw_recon(nt)%paw_indv(ih) ) then
              l = paw_recon(nt)%paw_nhtol(ih)
              lm = l * l + paw_recon(nt)%paw_nhtom(ih)
              do ig = 1, npw_
                 vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
              enddo 
           endif
        enddo
     enddo
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if (ityp (na) .eq.nt) then
           arg = (q_(1) * tau (1, na) + &
                  q_(2) * tau (2, na) + &
                  q_(3) * tau (3, na) ) * tpi
           phase = CMPLX (cos (arg), - sin (arg) )
           do ig = 1, npw_
              sk (ig) = eigts1 (ig1(igk_(ig)), na) * &
                        eigts2 (ig2(igk_(ig)), na) * &
                        eigts3 (ig3(igk_(ig)), na)
           enddo
           do ih = 1, paw_recon(nt)%paw_nh
              jkb = jkb + 1
              pref = (0.d0, -1.d0) ** paw_recon(nt)%paw_nhtol(ih) ! * phase
              do ig = 1, npw_
                 vkb_(ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
              enddo
              
           enddo
        endif
 
     enddo
     
     deallocate (vkb1)
  enddo
  
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)

  call stop_clock ('init_paw_2')
  return
  END SUBROUTINE init_paw_2_no_phase

  

  !---------------------------------------------------------------
  ! add the SO valence term
  !---------------------------------------------------------------
  SUBROUTINE add_so_valence(ik, n, alpha, psi, p_psic)
  USE wavefunctions_module, ONLY : psic
  USE gvect,                ONLY : nlm, ngm, g
  USE cell_base,            ONLY : tpiba
  USE lsda_mod,             ONLY : current_spin
  USE klist,                ONLY : xk
  USE scf,                  ONLY : dvrs
  USE wvfct,                ONLY : igk
  USE gsmooth,              ONLY : nls, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs
  IMPLICIT NONE
  integer :: ik, n, ipol, ig, i
  complex(dp) :: p_psic(nrxxs,3)
  complex(dp) :: psi(n)
  real(dp) :: gk, sigma, alpha
  ! index for the cross product
  integer :: ind(2,3), ii, jj
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)

  !!RETURN
  ! compute (p+k)|psi> in real space
  p_psic(1:nrxxs,1:3) = (0.d0,0.d0)
  do ipol = 1, 3
    do ig = 1, n
      gk = xk(ipol,ik) + g(ipol,igk(ig))
      p_psic(nls(igk(ig)),ipol) = -gk * tpiba * psi(ig)
    enddo
    call cft3s( p_psic(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
  enddo

  ! compute lambda_so ( sigma \cdot dvrs \times p_psi)
  sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

  do ipol = 1, 3
    if (lambda_so(ipol) == 0.d0) cycle
    ii = ind(1,ipol)
    jj = ind(2,ipol)
    do i = 1, nrxxs
      psic(i) = psic(i) + alpha * lambda_so(ipol) * a2gp8 * sigma * ( &
        dvrs(i,current_spin,ii)*p_psic(i,jj) - &
        dvrs(i,current_spin,jj)*p_psic(i,ii) )
    enddo
  enddo
  END SUBROUTINE add_so_valence



  !---------------------------------------------------------------
  ! add the F_R^{NL} term of the SO reconstruction
  !---------------------------------------------------------------
  SUBROUTINE add_so_Fnl(lda, n, m, alpha, psi, hpsi)
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : current_spin
  USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
  USE wvfct,           ONLY : current_k
  USE klist,           ONLY : xk
  USE wvfct,    ONLY : igk, g2kin, nbndx, nbnd
  USE gsmooth,  ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE paw,      ONLY : paw_becp
  IMPLICIT NONE
  INTEGER :: lda, n, m
  INTEGER :: ibnd, ijkb0, nt, na, ikb, jkb, ih, jh
  INTEGER :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
  COMPLEX(DP), ALLOCATABLE :: ps (:,:)
  COMPLEX(DP) :: psi(lda,n), hpsi(lda,n)
  real(dp) :: sigma, alpha

  !!RETURN
  if (m > nbndx) call errore('add_so_Fnl', 'm > nbndx ???', m)
  !!if (m > nbnd) call errore('add_so_Fnl', 'm > nbnd ???', m)
  ALLOCATE (ps(paw_nkb,m))
  ps(:,:) = (0.D0, 0.D0)

  call init_paw_2(n, igk, xk(1,current_k), paw_vkb)
  call ccalbec(paw_nkb, lda, n, m, paw_becp, paw_vkb, psi)

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
              ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(1) * a2gp8 * &
                sigma * lx(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)

              ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(2) * a2gp8 * &
                sigma * ly(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)

              ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(3) * a2gp8 * &
                sigma * lz(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp(jkb,ibnd)
            enddo   ! ibnd
          enddo   ! jh
        enddo   ! ih
        ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
      endif    ! ityp(na) .eq. nt
    enddo   ! na
  enddo   ! nt
  ps = ps * alpha

  CALL ZGEMM( 'N', 'N', n, m, paw_nkb, ( 0.D0, 1.0D0 ) , paw_vkb, &
              lda, ps, paw_nkb, ( 1.D0, 0.D0 ) , hpsi, lda )

  deallocate (ps)
  END SUBROUTINE add_so_Fnl 



  !====================================================================
  ! GIPAW contribution to the RMC term
  !====================================================================
  SUBROUTINE relativistic_mass_correction (ik, rmc_gipaw)
  USE atom,              ONLY : r, rab
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE paw,               ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, &
                                paw_becp
  USE wvfct,             ONLY : wg, nbnd
  IMPLICIT NONE
  real(dp), intent(inout):: rmc_gipaw
  integer :: ik, l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, lm, nrc, ijkb0
  complex(dp) :: rmc_corr
  complex(dp) :: bec_product

  rmc_corr = 0.0_dp

  do ibnd = 1, nbnd
    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if (ityp (na) == nt) then
          do ih = 1, paw_recon(nt)%paw_nh
            ikb = ijkb0 + ih
            nbs1=paw_recon(nt)%paw_indv(ih)
            l1=paw_recon(nt)%paw_nhtol(ih)
            m1=paw_recon(nt)%paw_nhtom(ih)
            lm1=m1+l1**2
            do jh = 1, paw_recon(nt)%paw_nh
              jkb = ijkb0 + jh
              nbs2=paw_recon(nt)%paw_indv(jh)
              l2=paw_recon(nt)%paw_nhtol(jh)
              m2=paw_recon(nt)%paw_nhtom(jh)
              lm2=m2+l2**2
              if ( l1 /= l2 ) cycle
              if ( m1 /= m2 ) cycle
              bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
              rmc_corr = rmc_corr + bec_product &
                       * radial_integral_rmc(nbs1,nbs2,nt) * wg(ibnd,ik)
            enddo
          enddo
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif
      enddo
    enddo
  enddo
  rmc_gipaw = -real(rmc_corr,dp)
  END SUBROUTINE relativistic_mass_correction

  



  !-----------------------------------------------------------------------
  SUBROUTINE calc_g_tensor(orb_magn_tot)
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_files,             ONLY : nwordwfc, iunwfc  
  USE wvfct,                ONLY : npw, g2kin, igk, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                ONLY : xk, nks 
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, ecutwfc, g
  USE io_global,            ONLY : stdout
  USE wavefunctions_module, ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE control_flags,        ONLY : iverbosity
  USE buffers,              ONLY : get_buffer  
  USE paw,                  ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, & 
                                   paw_becp
  implicit none
  real(dp), parameter :: rydtohar = 0.5d0
  real(dp) :: orb_magn_tot(3)
  real(dp) :: delta_rmc, delta_rmc_gipaw, lambda_mod, delta_g_orb_magn(3)
  real(dp) :: tmp, delta_g_total(3)
  integer :: ik, ibnd, ipol, s_weight

  !--------------------------------------------------
  ! Relativistic Mass Correction and diamagnetic term
  !--------------------------------------------------
  delta_rmc = 0.d0
  delta_rmc_gipaw = 0.d0

  do ik = 1, nks
    call get_buffer(evc, nwordwfc, iunwfc, ik)

    current_k = ik
    if (lsda) current_spin = isk(ik)
    s_weight = 1; if (current_spin == 2) s_weight = -1

    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(1:npw) = g2kin(1:npw) * tpiba2
    do ibnd = 1, nbnd
      delta_rmc = delta_rmc - s_weight &
            * wg(ibnd,ik) * sum(g2kin(1:npw) &
            * conjg(evc(1:npw,ibnd)) * evc(1:npw,ibnd))
    end do

    call init_paw_2(npw, igk, xk(1,ik), paw_vkb)
    call ccalbec(paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)

    call relativistic_mass_correction (ik, tmp)
    delta_rmc_gipaw = delta_rmc_gipaw + s_weight * tmp

  enddo
#ifdef __PARA
  call reduce(1, delta_rmc)
  ! no reduction for delta_rmc_gipaw
  call poolreduce(1, delta_rmc)
  call poolreduce(1, delta_rmc_gipaw)
#endif

  delta_rmc = delta_rmc * rydtohar * alpha**2.d0 * g_e * 1d6
  delta_rmc_gipaw = delta_rmc_gipaw * rydtohar * alpha**2.d0 * g_e * 1d6

  ! report results
  lambda_mod = sqrt(sum(lambda_so(:)**2.d0))
  delta_g_orb_magn = orb_magn_tot/lambda_mod * 1d6
  write(stdout,*)
  write(stdout,'(5X,''SPIN-ORBIT:'')')
  write(stdout,'(5X,''lambda_so          = '',3(F14.6))') lambda_so
  write(stdout,'(5X,''Contributions to the g-tensor (ppm):'')')
  write(stdout,'(5X,''delta_g RMC        = '',F14.4)') delta_rmc
  write(stdout,'(5X,''delta_g RMC(GIPAW) = '',F14.4)') delta_rmc_gipaw 
  write(stdout,'(5X,''delta_g SO         = '',3(F14.4))') delta_g_orb_magn

  ! put all terms together
  !
  ! delta_g_total = delta_g_orb_magn + delta_g_dia
  ! delta_g_dia still allready in delta_g_orb_magn included
  !
  delta_g_total = delta_g_orb_magn
  do ipol = 1, 3
    if (lambda_so(ipol) == 0.d0) cycle
    delta_g_total(ipol) = delta_g_total(ipol) + delta_rmc+delta_rmc_gipaw
  enddo
  write(stdout,*)
  write(stdout,'(5X,''delta_g total      = '',3(F14.4))') delta_g_total

  END SUBROUTINE calc_g_tensor


!-----------------------------------------------------------------------
END MODULE g_tensor_module
!-----------------------------------------------------------------------



