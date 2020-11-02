!
! Copyright (C) 2003-2014 the Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! TODO: If Fortran95 then allocatable else pointer


MODULE paw_gipaw

  USE kinds,        ONLY: dp
  USE parameters,   ONLY: npsx
  USE radial_grids, ONLY: ndmx
  USE parameters,   ONLY: ntypx, lmaxx

  IMPLICIT NONE
  SAVE
 
  !====================================================================
  ! These parameters are needed for the paw-gipaw projectors
  !====================================================================
  type wfc_label
     integer :: na          ! Atom number
     integer :: nt          ! type
     integer :: n           ! chi index
     integer :: l           ! l
     integer :: m           ! m
     integer :: nrc         ! index of core radius on mesh
     integer :: nrs         ! index of inner core radius (where "f_step" starts)
     real(dp) :: rc         ! paw core radius
  end type wfc_label

  type at_wfc
     type(wfc_label)          :: label
     integer                  :: kkpsi
     real(dp), pointer        :: psi(:)
  end type at_wfc


  ! fine structure constant alpha (c = 1/alpha)                                                                                           
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp
  ! g_e, g' and so on...                                                                                                                  
  REAL(DP), PARAMETER :: g_e = 2.0023192778_dp
  REAL(DP), PARAMETER :: gprime = 2.d0 * (g_e - 1.d0)
  REAL(DP), PARAMETER :: a2gp4 = alpha*alpha*gprime/4_dp   ! (alpha^2 g')/4                                                               
  REAL(DP), PARAMETER :: a2gp8 = alpha*alpha*gprime/8_dp   ! (alpha^2 g')/8   

  ! max number of projectors per atom
  integer, parameter :: nbrx = 20     ! TODO: dynamic alloc, paw_nkm
  logical :: pawproj(ntypx) = .false. ! TODO: from input

  ! UWG: parameters for scf spin-orbit coupling in GIPAW style
  !logical :: broken_time_reversal
  logical ::  reconstruction_only
  logical ::  align_spin
  real(dp) :: lambda_so(3) ! strength of spin-orbit
  complex(dp), allocatable :: deeq_so(:,:,:,:)

  real(dp) :: rc(ntypx,0:lmaxx)
  complex(dp), allocatable :: paw_becp2(:,:), paw_becp3(:,:)
  real(dp), allocatable, dimension (:,:) :: lx, ly, lz

  ! maximum dimension of radial integrals, ...      
  integer :: paw_nkm
  real(dp), allocatable :: radial_integral_paramagnetic(:,:,:)
  real(dp), allocatable :: radial_integral_diamagnetic(:,:,:)
  real(dp), allocatable :: radial_integral_paramagnetic_so(:,:,:)
  real(dp), allocatable :: radial_integral_diamagnetic_so(:,:,:)
  real(dp), allocatable :: radial_integral_rmc(:,:,:)

  ! parameters for paratec compatibility of the pseudos
  character(256) :: file_reconstruction(ntypx)
  logical :: read_recon_in_paratec_fmt  
  
  integer :: paw_nkb         ! total number of beta functions, with st.fact.
  integer :: paw_nqxq        ! size of interpolation table
  integer :: paw_lmaxkb      ! max angular momentum
  integer :: paw_lmaxq       ! max angular momentum + 1 for Q functions
  integer :: paw_nqx         ! number of interpolation points
  
  ! TODO: use the becp derived type
  complex(dp), allocatable :: paw_vkb(:,:) ! all beta functions in reciprocal space
  complex(dp), allocatable :: paw_becp(:,:)!  products of wavefunctions and proj
  
  type paw_recon_type
     real(dp) :: paw_betar(ndmx,nbrx)  ! radial beta_{mu} functions
     integer :: paw_nh                 ! number of beta functions per atomic type
     integer :: paw_nbeta              ! number of beta functions
     integer :: paw_kkbeta             ! point where the beta are zero
     ! TODO: remove nbrx and allocate it
     integer :: paw_lll(nbrx)          ! angular momentum of the beta function
     integer, pointer :: paw_indv(:)   ! correspondence of betas atomic <-> soli
     integer, pointer :: paw_nhtol(:)  ! correspondence n <-> angular momentum
     integer, pointer :: paw_nhtom(:)  ! correspondence n <-> magnetic angular m
     integer, pointer :: paw_nl(:)     ! number of projectors for each 
     integer, pointer :: paw_iltonh(:,:) ! corresp l, num <--> n for each type
     
     real(dp), pointer :: paw_tab(:,:)       ! interpolation table for pps
     real(dp), pointer :: paw_tab_d2y(:,:)   ! for cubic splines
     
     type(at_wfc), pointer :: aephi(:)
     type(at_wfc), pointer :: psphi(:)
     
     logical :: vloc_present
     real(dp), pointer :: gipaw_ae_vloc(:)
     real(dp), pointer :: gipaw_ps_vloc(:)
     
     logical :: gipaw_data_in_upf_file
     
     integer :: gipaw_ncore_orbital
     real(dp), pointer :: gipaw_core_orbital(:,:)
     integer, pointer :: gipaw_core_orbital_l(:)
  end type paw_recon_type
  
  type(paw_recon_type), allocatable :: paw_recon(:)
  
CONTAINS

  !==================================================================
  SUBROUTINE paw_wfc_init(phi)
  !==================================================================
  ! 
  ! Initialize default values for labe end kkpsi
  !
  implicit none
  type(at_wfc), intent(inout) :: phi(:)
       
  phi(:)%label%na  = 0
  phi(:)%label%nt  = 0
  phi(:)%label%n   = 0
  phi(:)%label%l   = -99
  phi(:)%label%m   = -99
  phi(:)%label%nrc = 0
  phi(:)%kkpsi     = 0
    
  END SUBROUTINE paw_wfc_init


  !==================================================================
  SUBROUTINE set_paw_upf(is, upf)
  !==================================================================
  !
  ! Interface between the UPF pseudo type and the internal representation
  ! of the PAW-related variables
  !
  USE ions_base,       ONLY : nsp
  USE atom,            ONLY : rgrid
  USE pseudo_types
  implicit none

  integer, intent(in) :: is
  type(pseudo_upf) :: upf
  integer :: nb
 
  if ( upf%has_gipaw ) then
     if ( .not. allocated ( paw_recon ) ) then !cg
        allocate ( paw_recon(nsp) )
        paw_recon(:)%gipaw_data_in_upf_file = .false.
     end if

     paw_recon(is)%paw_nbeta = upf%gipaw_wfs_nchannels
     paw_recon(is)%vloc_present = .true.
     paw_recon(is)%gipaw_data_in_upf_file = .true.
     
     paw_recon(is)%gipaw_ncore_orbital = upf%gipaw_ncore_orbitals
     allocate ( paw_recon(is)%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
     allocate ( paw_recon(is)%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
     paw_recon(is)%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals) &
          = upf%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals)
     paw_recon(is)%gipaw_core_orbital_l(:upf%gipaw_ncore_orbitals) &
          = upf%gipaw_core_orbital_l(:upf%gipaw_ncore_orbitals)
     
     allocate ( paw_recon(is)%gipaw_ae_vloc(upf%mesh) )
     allocate ( paw_recon(is)%gipaw_ps_vloc(upf%mesh) )
     paw_recon(is)%gipaw_ae_vloc(:upf%mesh) = upf%gipaw_vlocal_ae(:upf%mesh)
     paw_recon(is)%gipaw_ps_vloc(:upf%mesh) = upf%gipaw_vlocal_ps(:upf%mesh)
     
     allocate ( paw_recon(is)%aephi(upf%gipaw_wfs_nchannels) )
     allocate ( paw_recon(is)%psphi(upf%gipaw_wfs_nchannels) )
     
     do nb = 1, upf%gipaw_wfs_nchannels
        allocate ( paw_recon(is)%aephi(nb)%psi(rgrid(is)%mesh) )
        paw_recon(is)%aephi(nb)%label%nt = is
        paw_recon(is)%aephi(nb)%label%n = nb
        paw_recon(is)%aephi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%aephi(nb)%label%m = 
        paw_recon(is)%aephi(nb)%label%nrc = upf%mesh
        paw_recon(is)%aephi(nb)%kkpsi = upf%mesh
        paw_recon(is)%aephi(nb)%label%rc = upf%gipaw_wfs_rcut(nb)
        paw_recon(is)%aephi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ae(:upf%mesh,nb)
        
        allocate ( paw_recon(is)%psphi(nb)%psi(rgrid(is)%mesh) )
        paw_recon(is)%psphi(nb)%label%nt = is
        paw_recon(is)%psphi(nb)%label%n = nb
        paw_recon(is)%psphi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%psphi(nb)%label%m = 
        paw_recon(is)%psphi(nb)%label%nrc = upf%mesh
        paw_recon(is)%psphi(nb)%kkpsi = upf%mesh
        paw_recon(is)%psphi(nb)%label%rc = upf%gipaw_wfs_rcutus(nb)
        paw_recon(is)%psphi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ps(:upf%mesh,nb)
     end do
  end if

  END SUBROUTINE set_paw_upf
  
END MODULE paw_gipaw




! this two routines are for backward compatibility and are deactivated for
! the time being
#if 0 
    !==================================================================
    SUBROUTINE read_recon(filerec_sp, jtyp, paw_recon_sp)
    !==================================================================
    !
    ! Read all-electron and pseudo atomic wavefunctions 
    ! needed for PAW-GIPAW reconstruction
    !
    USE upf_module,         ONLY : scan_begin, scan_end
    USE atom,               ONLY : rgrid
    USE io_global,          ONLY : stdout
    implicit none
    character(len=256), intent(in) :: filerec_sp
    integer, intent(in) :: jtyp
    type(paw_recon_type), intent(inout) :: paw_recon_sp
    integer :: j, i, kkphi
    
    ! if the data has already been read from a UPF file
    if ( paw_recon_sp%gipaw_data_in_upf_file ) return
    
    OPEN ( 14, FILE = filerec_sp )
    CALL scan_begin ( 14, 'PAW', .true. )
    READ(14,*) paw_recon_sp%paw_nbeta
    CALL scan_end ( 14, 'PAW' )
    CLOSE ( 14 )
    
    ALLOCATE ( paw_recon_sp%psphi(paw_recon_sp%paw_nbeta) )
    ALLOCATE ( paw_recon_sp%aephi(paw_recon_sp%paw_nbeta) )
    
    CALL paw_wfc_init ( paw_recon_sp%psphi )
    CALL paw_wfc_init ( paw_recon_sp%aephi )
    
    OPEN ( 14, FILE = filerec_sp )
    
    WRITE (stdout,*) "N_AEwfc atom",jtyp,":", paw_recon_sp%paw_nbeta
    
    recphi_loop: DO i = 1, paw_recon_sp%paw_nbeta
       ALLOCATE ( paw_recon_sp%aephi(i)%psi(rgrid(jtyp)%mesh) )
       paw_recon_sp%aephi(i)%label%nt=jtyp
       paw_recon_sp%aephi(i)%label%n=i
       CALL scan_begin(14,'REC',.false.)
       CALL scan_begin(14,'kkbeta',.false.)
       read(14,*)  kkphi
       CALL scan_end(14,'kkbeta')
       paw_recon_sp%aephi(i)%kkpsi=kkphi
       CALL scan_begin(14,'L',.false.)        
       read(14,*) paw_recon_sp%aephi(i)%label%l
       CALL scan_end(14,'L')
       CALL scan_begin(14,'REC_AE',.false.)
       read(14,*) ( paw_recon_sp%aephi(i)%psi(j),j=1,kkphi)
       CALL scan_end(14,'REC_AE')
       ALLOCATE ( paw_recon_sp%psphi(i)%psi(rgrid(jtyp)%mesh) )
       paw_recon_sp%psphi(i)%label%nt = jtyp
       paw_recon_sp%psphi(i)%label%n = i
       paw_recon_sp%psphi(i)%label%l = paw_recon_sp%aephi(i)%label%l
       paw_recon_sp%psphi(i)%kkpsi = kkphi
       CALL scan_begin(14,'REC_PS',.false.)
       READ(14,*) ( paw_recon_sp%psphi(i)%psi(j),j=1,kkphi)
       CALL scan_end(14,'REC_PS')
       CALL scan_end(14,'REC')
    END DO recphi_loop
    CLOSE(14)
    
  END SUBROUTINE read_recon
  
   subroutine read_recon_paratec ( filerec_sp, jtyp, paw_recon_sp, vloc_set )
    
    !
    ! Read all-electron and pseudo atomic wavefunctions 
    !  needed for PAW reconstruction
    !
    
    use upf_module, only: scan_begin, scan_end
    USE ions_base,          ONLY : ntyp => nsp
    use atom, only: msh, rgrid
    use kinds, only: DP
    use parameters, only : ntypx
    USE io_global,  ONLY : stdout
    use splinelib
    
    IMPLICIT NONE
    
    ! Arguments
    CHARACTER ( LEN = 256 ), INTENT ( IN ) :: filerec_sp
    INTEGER, INTENT ( IN ) :: jtyp
    TYPE ( paw_recon_type ), INTENT ( INOUT ) :: paw_recon_sp
    LOGICAL, INTENT ( OUT ) :: vloc_set
    
    ! Local
    character (len=256) :: readline
    logical :: file_exists, new_file_format
    integer :: iostatus, nlines
    integer :: j,i, kkpsi
    real(dp) :: d1
    INTEGER :: local_component
    
    type at_wfc_r
       type(wfc_label)          :: label
       integer                  :: kkpsi
       real(DP)  , pointer :: r(:)
    end type at_wfc_r
    
    real(dp), allocatable :: xdata(:), tab(:), tab_d2y(:)
    real(dp), allocatable :: gipaw_ae_vloc2 ( :, : )
    real(dp), allocatable :: gipaw_ps_vloc2 ( :, : )
    type(at_wfc_r), allocatable :: aephi2_r ( : ), psphi2_r ( : )
    type(at_wfc),pointer :: aephi2(:), psphi2(:) ! Atom
    
    !
    
    IF ( paw_recon_sp%gipaw_data_in_upf_file ) THEN
       RETURN
    END IF
    
    IF ( TRIM ( filerec_sp ) == "" ) THEN
       ! No reconstruction asked for this species
       paw_recon_sp%paw_nbeta = 0
       nlines = 0
       RETURN
    END IF
    
    inquire ( file = filerec_sp, exist = file_exists )
    if ( .not. file_exists ) then
       call errore ( "reconstruction file does not exist", &
            TRIM(filerec_sp), 1 )
       stop
    end if
    
    open ( 14, file = filerec_sp )
    
    paw_recon_sp%paw_nbeta = 0
    nlines = 0
    do
       READ ( UNIT = 14, FMT = '( 256A )', IOSTAT = iostatus ) readline
       IF ( iostatus /= 0 ) THEN
          EXIT
       END IF
       IF ( INDEX ( readline, "#core wavefunctions" ) /= 0 ) THEN
          EXIT
       END IF
       IF ( INDEX ( readline, "#" ) /= 0 ) THEN
          nlines = 0
       END IF
       IF ( INDEX ( readline, "&" ) /= 0 ) THEN
          paw_recon_sp%paw_nbeta = paw_recon_sp%paw_nbeta + 1
       END IF
       IF ( INDEX ( readline, "&" ) == 0 .AND. &
            INDEX ( readline, "#" ) == 0 ) THEN
          nlines = nlines + 1
       END IF
    end do
    close ( unit = 14 )
    
    allocate( psphi2(paw_recon_sp%paw_nbeta) )
    allocate( aephi2(paw_recon_sp%paw_nbeta) )
    allocate( psphi2_r(paw_recon_sp%paw_nbeta) )
    allocate( aephi2_r(paw_recon_sp%paw_nbeta) )
    
    allocate ( gipaw_ae_vloc2(nlines,paw_recon_sp%paw_nbeta) )
    allocate ( gipaw_ps_vloc2(nlines,paw_recon_sp%paw_nbeta) )
    
    call paw_wfc_init ( psphi2 )
    call paw_wfc_init ( aephi2 )
    
    paw_recon_sp%vloc_present = .FALSE.
    
    open ( 14, file = filerec_sp )
    rewind ( unit = 14 )
    
    recphi_loop: do i = 1, paw_recon_sp%paw_nbeta
       aephi2(i)%label%nt = jtyp
       aephi2(i)%label%n = i
       kkpsi = nlines
       aephi2(i)%kkpsi = kkpsi
       allocate ( aephi2(i)%psi(kkpsi), aephi2_r(i)%r(kkpsi) )
       allocate ( psphi2(i)%psi(kkpsi), psphi2_r(i)%r(kkpsi) )
       read(14,*) readline
       IF ( i == 1 ) new_file_format = .FALSE.
       read(14, FMT = '( 256A )' ) readline
       IF ( readline(1:6) == "#local" ) THEN
          new_file_format = .TRUE.
          paw_recon_sp%vloc_present = .TRUE.
          READ ( readline(8:8), * ) local_component
          read(14, FMT = '( 256A )' ) readline
       END IF
       IF ( readline(1:3) /= "#l=" ) THEN
          WRITE ( UNIT = stdout, FMT = '( 2A, ":" )' ) &
               "Wrong control string in file ", TRIM ( filerec_sp )
          WRITE ( UNIT = stdout, FMT = '( A )' ) TRIM ( readline )
          CALL errore ( "read_recon_paratec", "wrong control string", 1 )
          stop
       END IF
       read(readline(4:4), FMT = * ) aephi2(i)%label%l
       read(readline(12:), FMT = * ) aephi2(i)%label%rc
       
       IF ( new_file_format ) THEN
          DO j = 1, kkpsi
             read(14,*) aephi2_r(i)%r(j), aephi2(i)%psi(j), &
                  psphi2(i)%psi(j), &
                  gipaw_ae_vloc2(j,i), gipaw_ps_vloc2(j,i)
          END DO
       ELSE
          DO j = 1, kkpsi
             read(14,*) aephi2_r(i)%r(j), aephi2(i)%psi(j), &
                  psphi2(i)%psi(j)
          END DO
       END IF
          
       psphi2(i)%label%nt = jtyp
       psphi2(i)%label%n = i
       psphi2(i)%label%l = aephi2(i)%label%l
       psphi2(i)%label%rc = aephi2(i)%label%rc
       psphi2(i)%kkpsi = kkpsi
       psphi2_r(i)%r(:) = aephi2_r(i)%r(:)
    end do recphi_loop
    close ( 14 )
    
    !
    
    allocate ( paw_recon_sp%psphi(paw_recon_sp%paw_nbeta) )
    allocate ( paw_recon_sp%aephi(paw_recon_sp%paw_nbeta) )
    
    call paw_wfc_init ( paw_recon_sp%psphi )
    call paw_wfc_init ( paw_recon_sp%aephi )
    
    IF ( paw_recon_sp%vloc_present ) THEN
       ALLOCATE ( paw_recon_sp%gipaw_ae_vloc(msh(jtyp)) )
       ALLOCATE ( paw_recon_sp%gipaw_ps_vloc(msh(jtyp)) )
    END IF
    
    vloc_set = .FALSE.
    
    do i = 1, paw_recon_sp%paw_nbeta
       
       ! AE
       kkpsi = aephi2(i)%kkpsi
       allocate( xdata(kkpsi), tab(kkpsi), tab_d2y(kkpsi) )
       xdata(:) = aephi2_r(i)%r(:)
       tab(:) = aephi2(i)%psi(:)
       
       ! initialize spline interpolation
       d1 = ( tab(2) - tab(1) ) / ( xdata(2) - xdata(1) )
       call spline ( xdata, tab, 0.0_dp, d1, tab_d2y )
       
       ! use interpolation
       allocate ( paw_recon_sp%aephi(i)%psi(msh(jtyp)) )
       paw_recon_sp%aephi(i)%label%nt = jtyp
       paw_recon_sp%aephi(i)%label%n = i
       paw_recon_sp%aephi(i)%label%l = aephi2(i)%label%l
       paw_recon_sp%aephi(i)%label%rc = aephi2(i)%label%rc
       paw_recon_sp%aephi(i)%kkpsi = msh(jtyp)
       do j = 1, msh(jtyp)
          paw_recon_sp%aephi(i)%psi(j) &
               = splint ( xdata, tab, tab_d2y, rgrid(jtyp)%r(j) )
       end do
       
       ! PS        
       allocate ( paw_recon_sp%psphi(i)%psi(msh(jtyp)) )
       xdata(:) = psphi2_r(i)%r(:)
       tab(:) = psphi2(i)%psi(:)
       
       ! initialize spline interpolation
       d1 = ( tab(2) - tab(1) ) / ( xdata(2) - xdata(1) )
       call spline ( xdata, tab, 0.0_dp, d1, tab_d2y )
       
       ! use interpolation
       allocate ( paw_recon_sp%psphi(i)%psi(msh(jtyp)) )
       paw_recon_sp%psphi(i)%label%nt = jtyp
       paw_recon_sp%psphi(i)%label%n = i
       paw_recon_sp%psphi(i)%label%l = psphi2(i)%label%l
       paw_recon_sp%psphi(i)%label%rc = psphi2(i)%label%rc
       paw_recon_sp%psphi(i)%kkpsi = msh(jtyp)
       !paw_recon_sp%psphi(i)%r(1:msh(jtyp)) = r(1:msh(jtyp),jtyp)
       do j = 1, msh(jtyp)
          paw_recon_sp%psphi(i)%psi(j) = splint(xdata, tab, tab_d2y, rgrid(jtyp)%r(j))
       end do
       
       ! for the local potential; it is the same for all components, choose 1
       IF ( paw_recon_sp%vloc_present .AND. i == 1 ) THEN
          tab(:) = gipaw_ae_vloc2(:kkpsi,i)
          d1 = ( gipaw_ae_vloc2(2,i) - gipaw_ae_vloc2(1,i) ) &
               / ( xdata(2) - xdata(1) )
          call spline(xdata, tab, 0.0_dp, d1, tab_d2y)
          DO j = 1, msh(jtyp)
             paw_recon_sp%gipaw_ae_vloc(j) &
                  = splint ( xdata, tab, tab_d2y, rgrid(jtyp)%r(j) )
          END DO
       END IF
       
       ! This assumes that the FIRST projector on a given channel is the
       !    local component
       IF ( paw_recon_sp%vloc_present &
            .AND. paw_recon_sp%aephi(i)%label%l == local_component &
            .AND. .NOT. vloc_set ) THEN
          tab(:) = gipaw_ps_vloc2(:kkpsi,i)
          d1 = ( gipaw_ps_vloc2(2,i ) - gipaw_ps_vloc2(1,i) ) &
               / ( xdata(2) - xdata(1) )
          call spline ( xdata, tab, 0.0_dp, d1, tab_d2y )
          DO j = 1, msh(jtyp)
             paw_recon_sp%gipaw_ps_vloc(j) &
                  = splint ( xdata, tab, tab_d2y, rgrid(jtyp)%r(j))
          END DO
          
          vloc_set = .TRUE.
       END IF
       
       deallocate( xdata, tab, tab_d2y)
    end do

    IF ( .NOT. vloc_set ) write(stdout,*) &
                  "WARNING: no local potential set in read_recon_paratec." 

    
    deallocate( psphi2 )
    deallocate( aephi2 )
    deallocate( psphi2_r )
    deallocate( aephi2_r )
    
    DEALLOCATE ( gipaw_ae_vloc2 )
    DEALLOCATE ( gipaw_ps_vloc2 )
    
  end subroutine read_recon_paratec
#endif

