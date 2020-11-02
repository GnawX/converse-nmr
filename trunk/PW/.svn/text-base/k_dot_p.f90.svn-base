!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"

!-----------------------------------------------------------------
subroutine cg_psi (lda, n, m, psi, h_diag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, only : DP
  implicit none

  integer :: lda, n, m
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors

  complex(DP) :: psi (lda, m)
  ! inp/out: the vector to be preconditioned

  real(DP) :: h_diag (lda, m)
  ! input: the preconditioning vector

  integer :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  do k = 1, m
     do i = 1, n
        psi (i, k) = psi (i, k) * h_diag (i, k)
     enddo
  enddo
  return
end subroutine cg_psi

!----------------------------------------------------------------------
subroutine cgsolve_all (ch_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, alpha_pv)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 ( h - e + Q ) * dpsi = d0psi                      (1)
  !
  !                 where h is a complex hermitean matrix, e is a real sca
  !                 dpsi and d0psi are complex vectors
  !
  !     on input:
  !                 h_psi    EXTERNAL  name of a subroutine:
  !                          h_psi(ndim,psi,psip)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
  USE kinds, only : DP
  implicit none
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik      ! input: the k point

  real(DP) :: &
             e(nbnd), & ! input: the actual eigenvalue
             anorm,   & ! output: the norm of the error in the solution
             h_diag(ndmx,nbnd), & ! input: an estimate of ( H - \epsilon )
             ethr       ! input: the required precision

  complex(DP) :: &
             dpsi (ndmx, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx, nbnd)   ! input: the known term

  logical :: conv_root ! output: if true the root is converged
  external ch_psi, &   ! input: the routine computing ch_psi
           cg_psi      ! input: the routine computing cg_psi
  !
  !  here the local variables
  !
  integer, parameter :: maxter = 200
  ! the maximum number of iterations
  integer :: iter, ibnd, lbnd
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged

  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(DP) ::  dcgamma, dclambda, ZDOTC
  !  the ratio between rho
  !  step length
  !  the scalar product
  real(DP), allocatable :: rho (:), rhoold (:), eu (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  real(dp) :: alpha_pv
  call start_clock ('cgsolve')
  allocate ( g(ndmx,nbnd), t(ndmx,nbnd), h(ndmx,nbnd), hold(ndmx ,nbnd) )    
  allocate (a(nbnd), c(nbnd))    
  allocate (conv ( nbnd))    
  allocate (rho(nbnd),rhoold(nbnd))    
  allocate (eu (  nbnd))    
  !      WRITE( stdout,*) g,t,h,hold

  kter_eff = 0.d0
  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo 
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  do iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     if (iter == 1) then
        call ch_psi (ndim, dpsi, g, e, ik, nbnd, alpha_pv)
        do ibnd = 1, nbnd
           call ZAXPY (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
        enddo
     endif
     !
     !    compute preconditioned residual vector and convergence check
     !
     lbnd = 0
     do ibnd = 1, nbnd
!!!        if (conv (ibnd) .eq.0) then
           lbnd = lbnd+1
           call ZCOPY (ndim, g (1, ibnd), 1, h (1, ibnd), 1)
           call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
           rho(lbnd) = ZDOTC (ndim, h(1,ibnd), 1, g(1,ibnd), 1)
!!!        endif
     enddo
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
#ifdef __PARA
     call reduce (lbnd, rho )
#endif
     !!!print*, '----------------------------------------------'
     do ibnd = nbnd, 1, -1
!!!        if (conv(ibnd).eq.0) then
           rho(ibnd)=rho(lbnd)
           lbnd = lbnd -1
           anorm = sqrt (rho (ibnd) )
           if (anorm.lt.ethr) conv (ibnd) = 1
           !!!write(*,'(3I10,E12.4)') iter, ibnd,conv(ibnd), anorm
!!!        endif
     enddo
!
     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root.and. (conv (ibnd) .eq.1)
     enddo
     if (conv_root) goto 100
     !
     !        compute the step direction h. Conjugate it to previous step
     !
     lbnd = 0
     do ibnd = 1, nbnd
!!!        if (conv (ibnd) .eq.0) then
!
!          change sign to h 
!
           call DSCAL (2 * ndim, - 1.d0, h (1, ibnd), 1)
           if (iter.ne.1) then
              dcgamma = rho (ibnd) / rhoold (ibnd)
              call ZAXPY (ndim, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
           endif

!
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h 
!
           lbnd = lbnd+1
           call ZCOPY (ndim, h (1, ibnd), 1, hold (1, lbnd), 1)
           eu (lbnd) = e (ibnd)
!!!        endif
     enddo
     !
     !        compute t = A*h
     !
     call ch_psi (ndim, hold, t, eu, ik, lbnd, alpha_pv)
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     lbnd=0
     do ibnd = 1, nbnd
!!!        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           a(lbnd) = ZDOTC (ndim, h(1,ibnd), 1, g(1,ibnd), 1)
           c(lbnd) = ZDOTC (ndim, h(1,ibnd), 1, t(1,lbnd), 1)
!!!        end if
     end do
#ifdef __PARA
     call reduce (lbnd, a)
     call reduce (lbnd, c)
#endif
     lbnd=0
     do ibnd = 1, nbnd
!!!        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           dclambda = CMPLX ( - a(lbnd) / c(lbnd), 0.d0)
           !
           !    move to new position
           !
           call ZAXPY (ndim, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
           !
           !    update to get the gradient
           !
           !g=g+lam
           call ZAXPY (ndim, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
           !
           !    save current (now old) h and rho for later use
           ! 
           call ZCOPY (ndim, h(1,ibnd), 1, hold(1,ibnd), 1)
           rhoold (ibnd) = rho (ibnd)
!!!        endif
     enddo
  enddo
100 continue
  kter = kter_eff
  deallocate (eu)
  deallocate (rho, rhoold)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)

  call stop_clock ('cgsolve')
  return
end subroutine cgsolve_all


!-----------------------------------------------------------------------
subroutine ch_psi_all (n, h, ah, e, ik, m, alpha_pv)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE pwcom, only: npwx
  USE becmod
  USE kinds, only : DP
  USE wavefunctions_module, ONLY : evc
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  real(DP) :: e (m), alpha_pv
  ! input: the eigenvalue

  complex(DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h
  
  call start_clock ('ch_psi')
  !!allocate (ps  ( nbnd , m))    
  allocate (ps  ( m , m))    
  allocate (hpsi( npwx , m))    
  allocate (spsi( npwx , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  call h_psi (npwx, n, m, h, hpsi)
  call s_psi (npwx, n, m, h, spsi)

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  !
  !   Here we compute the projector in the valence band
  !
  ps (:,:) = (0.d0, 0.d0)
  !call ZGEMM ('C', 'N', m, m, n, (1.d0, 0.d0) , evc, &
  !            npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  call ZGEMM ('C', 'N', m, m, n, (1.d0, 0.d0) , evc, &
              npwx, spsi, npwx, (0.d0, 0.d0) , ps, m)
  ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
  call reduce (2 * m * m, ps)
#endif
  hpsi (:,:) = (0.d0, 0.d0)
  !call ZGEMM ('N', 'N', n, m, nbnd, (1.d0, 0.d0) , evc, &
  !            npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  call ZGEMM ('N', 'N', n, m, m, (1.d0, 0.d0) , evc, &
              npwx, ps, m, (1.d0, 0.d0) , hpsi, npwx)
  spsi(:,:) = hpsi(:,:)
  !
  !    And apply S again
  !
  !!!call ccalbec (nkb, npwx, n, m, becp, vkb, hpsi)

  !!!call s_psi (npwx, n, m, hpsi, spsi)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
  enddo

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all



!-----------------------------------------------------------------------
SUBROUTINE greenfunction(ik, psi, nbnd_occ, alpha_pv, g_psi, ethr_q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the Green function operator
  ! ... (formula here)
  ! ... We use Hartree atomic units; since G is the inverse of an
  ! ... energy: G => G / ryd_to_hartree
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout  
  USE becmod,                      ONLY : becp
  USE wavefunctions_module,        ONLY : evc
  USE noncollin_module,            ONLY : npol
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE wvfct,                       ONLY : wg
  USE klist,                       ONLY : lgauss, degauss, ngauss
  USE pwcom
  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)  ! psi is changed on output!!!
  COMPLEX(DP), INTENT(OUT) :: g_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  real(dp), parameter :: conv_threshold = 1d-18
  complex(dp), allocatable :: ps(:,:), work (:)
  real(dp), allocatable :: h_diag (:,:), eprec (:)
  real(dp) :: anorm, thresh, gk(3), dxk(3)
  integer :: ibnd, jbnd, ig, lter, nbnd_occ
  logical :: conv_root
  real(dp) :: alpha_pv
  real(dp) :: wg1, w0g, wgp, wwg, deltae, theta
  real(dp) :: ethr_q
  complex(dp), external :: ZDOTC
  external ch_psi_all, cg_psi
  real(DP), external :: w0gauss, wgauss

  ! start clock
  call start_clock ('greenf')
  current_k = ik
  if (lsda) current_spin = isk(ik)

  ! allocate memory
  allocate (work(npwx), ps(nbnd,nbnd), h_diag(npwx,nbnd), &
            eprec(nbnd), becp(nkb,nbnd))

  !====================================================================
  ! apply -Q_k to the r.h.s.
  !====================================================================
  ! project on <evc|: ps(i,j) = <evc(i)|psi(j)>
  ps = (0.d0,0.d0)
  if (lgauss) then
    ! metallic case
    CALL ZGEMM('C', 'N', nbnd, nbnd_occ, npw, &
               (1.d0,0.d0), evc(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
               ps(1,1), nbnd)
    do ibnd = 1, nbnd_occ
      wg1 = wgauss ((ef-et(ibnd,ik)) / degauss, ngauss)
      w0g = w0gauss((ef-et(ibnd,ik)) / degauss, ngauss) / degauss
      do jbnd = 1, nbnd
        wgp = wgauss ( (ef - et (jbnd, ik) ) / degauss, ngauss)
        deltae = et (jbnd, ik) - et (ibnd, ik)
        theta = wgauss (deltae / degauss, 0)
        wwg = wg1 * (1.d0 - theta) + wgp * theta
        if (jbnd <= nbnd_occ) then
          if (abs (deltae) > 1.0d-5) then
            wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
          else
            wwg = wwg - alpha_pv * theta * w0g
          endif
        endif
        ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
      enddo
      call DSCAL (2*npw, wg1, psi(1,ibnd), 1)
    enddo
  
  else
    ! insulators
    CALL ZGEMM('C', 'N', nbnd_occ, nbnd_occ, npw, &
               (1.d0,0.d0), evc(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
               ps(1,1), nbnd)
  endif
#ifdef __PARA
  call reduce (2 * nbnd * nbnd_occ, ps(1,1))
#endif

  !! this is the case with overlap (ultrasoft)
  !! g_psi is used as work space to store S|evc>
  !!
  !!CALL ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, evc)
  !!CALL s_psi (npwx, npw, nbnd_occ(ik), evc, g_psi)
  !! |psi> = -(|psi> - S|evc><evc|psi>)
  !!
  !!CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
  !!     (1.d0,0.d0), g_psi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
  !!     psi(1,1), npwx )

  ! |psi> = -(1 - |evc><evc|) |psi>
  CALL ZGEMM('N', 'N', npw, nbnd_occ, nbnd_occ, &
             (1.d0,0.d0), evc(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
             psi(1,1), npwx)


  !====================================================================
  ! solve the linear system (apply G_{k+q})
  !====================================================================
  ! convergence treshold
  !thresh = sqrt(conv_threshold)   ! sqrt(of that of PARATEC)
  thresh = sqrt(ethr_q)   ! sqrt(of that of PARATEC)
 
  ! use the hamiltonian at k
  do ig = 1, npw
    gk(1) = (xk(1,ik) + g(1,igk(ig))) * tpiba
    gk(2) = (xk(2,ik) + g(2,igk(ig))) * tpiba
    gk(3) = (xk(3,ik) + g(3,igk(ig))) * tpiba
    g2kin (ig) = gk(1)**2 + gk(2)**2 + gk(3)**2
  enddo

  ! preconditioning of the linear system
  do ibnd = 1, nbnd_occ
     do ig = 1, npw
        work (ig) = g2kin (ig) * evc (ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0 * ZDOTC (npw, evc (1, ibnd), 1, work, 1)
  enddo
#ifdef __PARA
  call reduce (nbnd_occ, eprec)
#endif
  do ibnd = 1, nbnd_occ
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     enddo
  enddo

  call init_us_2(npw, igk, xk(1,ik), vkb)
  call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, psi)
    
  ! initial guess
  !!g_psi(:,:) = (0.d0, 0.d0)
  g_psi(:,:) = -psi(:,:)

  ! solve linear system  
  conv_root = .true.
  call cgsolve_all (ch_psi_all, cg_psi, et(1,ik), psi, g_psi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ, alpha_pv)

  !! debug  
  !!write(stdout, '(5X,''cgsolve_all converged in '',I3,'' iterations'')') &
  !!      lter

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4, &
       & " linter: root not converged ",e10.3)') ik, anorm

  ! convert to Hartree
  g_psi(:,:) = g_psi(:,:) / ryd_to_hartree

  call flush_unit( stdout )
  call stop_clock('greenf')
 
  ! free memory
  deallocate (work, h_diag, eprec, ps, becp)

END SUBROUTINE greenfunction




!-----------------------------------------------------------------------
SUBROUTINE apply_vel(psi, vel_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  ! ...   v = p + dV^{NL}_{k,k}/dk = dH_k/dk
  ! ...
  ! ... Here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, current_k
  USE lsda_mod,             ONLY : current_spin, isk
  USE becmod,               ONLY : becp
  USE cell_base,            ONLY : tpiba
  USE g_tensor_module,      ONLY : init_us_2_no_phase
  USE g_tensor_module,      ONLY : q_gipaw
  USE pwcom

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  !real(dp), parameter :: q_gipaw = 0.02d0
  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
  real(dp) :: dk, oldk(3)
  integer :: isign

  call start_clock('apply_vel')
  vel_psi = (0.d0,0.d0)

  ! set dk (= delta_k ?)
  !dk = q_gipaw/2.d0/tpiba
  dk = q_gipaw/tpiba
  
  ! allocate temporary arrays, save old NL-potential
  allocate(aux(npwx,nbnd), vkb_save(npwx,nkb))
  vkb_save = vkb
  oldk(:) = xk(:,ik)

  current_k = ik
  if (lsda) current_spin = isk(ik)
  call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)

  !====================================================================
  ! compute (1/2|dk|) ( H_{k+dk} |psi> - H_{k-dk}|psi> )
  !====================================================================
  allocate(becp(nkb,nbnd))
  do isign = -1,1,2
    xk(ipol,ik) = oldk(ipol) + isign * dk     ! k + dk

    ! compute <\beta(k \pm dk)| and project on |psi>
    call g2_kin(ik)
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    aux = (0.d0,0.d0)
    call h_psi(npwx, npw, nbnd, psi, aux)
    vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
  enddo
  deallocate(becp)

  ! restore NL-potential at k
  xk(:,ik) = oldk(:)
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel

