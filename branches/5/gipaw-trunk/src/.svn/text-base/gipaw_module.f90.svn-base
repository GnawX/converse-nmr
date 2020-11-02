!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! there are still some problems with the dimensions in calbec_bands...
#ifdef __BANDS
#error "Band parallelization has not yet been ported to this version!"
#endif

!-----------------------------------------------------------------------
MODULE gipaw_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for GIPAW calculations
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : a0_to_cm => bohr_radius_cm
  USE parameters, ONLY : npk, ntypx, lmaxx
  USE paw_gipaw ! now in PW/src
  
  IMPLICIT NONE
  SAVE

  !!INTEGER, PARAMETER:: nbrx=2*(lmaxx+1)**2  ! max number of beta functions
 
  ! Physical constants:
  ! Avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142e23_dp
 
  ! rydberg to Hartree
  REAL(DP), PARAMETER :: ry2ha = 0.5_DP

  ! number of occupied bands at each k-point
  INTEGER :: nbnd_occ(npk)
  
  ! alpha shift of the projector on the valence wfcs
  REAL(DP) :: alpha_pv
  
  ! eigenvalues and eigenfunctions at k+q
  REAL(DP), ALLOCATABLE :: etq(:,:)
  COMPLEX(DP), ALLOCATABLE :: evq(:,:)

  ! convergence threshold for diagonalizationa and greenfunction
  REAL(DP) :: conv_threshold
  
  ! q for the perturbation (in bohrradius^{-1})
  REAL(DP) :: q_gipaw
  
  ! q for the EFG
  REAL(DP) :: q_efg ( ntypx )

  ! restart mode: 'from_scratch' or 'restart'
  CHARACTER(80) :: restart_mode

  ! disk_io: 'high': restart possible from any k+q  (slow)
  !          'normal': restart from any k-point  (default)
  !          'low':  no restart information written (fast)
  !          'none': nothing writting to disk, wfc in RAM!
  CHARACTER(80) :: disk_io

  ! verbosity
  INTEGER :: iverbosity
 
  ! diagonalization method
  INTEGER :: isolve
 
  ! job: nmr, g_tensor, efg, hyperfine
  CHARACTER(80) :: job
 
  ! core-relax method to calculate change of XC
  INTEGER :: core_relax_method = 0
  REAL(dp) :: ncore_orbitals( ntypx )  ! UWG: number of PAW core orbitals to be recalculated  
 
  ! format for a rank-2 tensor
  CHARACTER(*), PARAMETER :: tens_fmt = '(3(5X,3(F14.4,2X)/))'
  
  ! for plotting the induced current and induced field
  CHARACTER(80) :: filcurr, filfield, filgradv, filnics

  ! macroscopic shape for the NMR
  LOGICAL :: use_nmr_macroscopic_shape
  REAL(DP) :: nmr_macroscopic_shape(3,3)

  ! contribution to NMR chemical shift due to core contribution
  REAL(dp) :: nmr_shift_core(ntypx)

  ! parameters for hyperfine interaction
  CHARACTER(10) :: hfi_input_unit
  CHARACTER(10) :: hfi_output_unit
  REAL(dp) :: hfi_nuclear_g_factor(ntypx)
  LOGICAL :: hfi_via_reconstruction_only = .false. !UWG: speed up for large systems 
 

  REAL(dp) :: r_rand ! EMINE: determines the randomization range used in
                     ! compute_u_kq routine. read from input. change it for debugging only

  ! parameters for paratec compatibility of the pseudos
  !CHARACTER(256) :: file_reconstruction ( ntypx )
  !LOGICAL :: read_recon_in_paratec_fmt  
  !LOGICAL :: pawproj( ntypx ) ! EMINE: if true paw projectors will be used instead of GIPAW ones
  ! UWG: parameters for scf spin-orbit coupling in GIPAW style
  !LOGICAL ::  reconstruction_only
  !LOGICAL ::  align_spin
  !REAL(DP) :: lambda_so(3) ! strength of spin-orbit

  INTEGER, PARAMETER :: iumagres = 57

#ifdef __BANDS
  INTEGER :: ibnd_start, ibnd_end
#endif
!-----------------------------------------------------------------------
END MODULE gipaw_module
!-----------------------------------------------------------------------
