!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE write_tensor_field(name, ispin, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the tensor field in VTK format
  !
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout, ionode
  USE cell_base,                   ONLY : at, alat
  USE ions_base,                   ONLY : nat, tau, atm, ityp
  USE fft_base,                    ONLY : dfftp, grid_gather
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  implicit none
  character*(*) name
  integer :: ispin
  real(dp) :: field(dfftp%nnr,3,3)
  !--------------------------------------------------------------------
  integer, parameter :: ounit = 48
  character*80 :: fname
  integer :: ios, ipol
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: aux(:,:,:)
#ifdef __MPI
  integer :: i, j
#endif

 allocate( aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3,3) )

#ifdef __MPI
 ! gather the data 
  do i = 1, 3
    do j = 1, 3
      call grid_gather(field(:,i,j), aux(:,i,j))
    enddo
  enddo
#else
  aux = field
#endif

  if (ionode) then
    do ipol = 1, 3
      ! form the name of the output file
      if (nspin == 1) then
        fname = trim(name)//'_'
      elseif (nspin == 2 .and. ispin == 1) then
        fname = trim(name)//'_UP_'
      elseif (ispin == 2 .and. ispin == 2) then
        fname = trim(name)//'_DW_'
      endif

      if (ipol == 1) fname = trim(fname)//'X.vtk'
      if (ipol == 2) fname = trim(fname)//'Y.vtk'
      if (ipol == 3) fname = trim(fname)//'Z.vtk'
      write(stdout, '(5X,''write_tensor_field: '',A40)') fname

      open(unit=ounit, file=fname, iostat=ios, form='formatted', status='unknown')
      if (ios /= 0) call errore('write_tensor_field', 'error opening '//fname, ounit)

      call vtk_vector_3d(aux(:,:,ipol), dfftp%nr1, dfftp%nr2, dfftp%nr3, at, alat, ounit)

      close(unit=ounit)
    enddo
  endif

  deallocate(aux)

end subroutine write_tensor_field


!-------------------------------------------------------------------
! this routine writes a 3D vector field in VTK format
!-------------------------------------------------------------------
subroutine vtk_vector_3d(vin, nr1, nr2, nr3, at, alat, ounit)
  USE kinds,           only : dp
  USE constants,       only : bohr_radius_angs
  USE fft_base,        only : dfftp

  implicit none
  integer, intent(in) ::  nr1, nr2, nr3, ounit
  real(dp), intent(in) :: at(3,3), alat
  real(dp), intent(in) :: vin(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3)  

  integer :: i1, i2, i3, bla
  real(dp) :: x(3)

  ! header
  write(ounit,'(A)') '# vtk DataFile Version 2.0'
  write(ounit,'(A)') 'created by qe-gipaw'
  write(ounit,'(A)') 'ASCII'
  write(ounit,'(A)') 'DATASET STRUCTURED_GRID'
  write(ounit,'(A,3I5)') 'DIMENSIONS', nr1, nr2, nr3
  write(ounit,'(A,I10,4X,A)') 'POINTS', nr1*nr2*nr3, 'float'

  ! point coordinates
  do i3 = 1, nr3
    do i2 = 1, nr2
      do i1 = 1, nr1
        ! coordinate in angstrom
        x(1) = dble(i1-1)/dble(nr1)     
        x(2) = dble(i2-1)/dble(nr2)     
        x(3) = dble(i3-1)/dble(nr3)
        ! crystal to cartesian
        call cryst_to_cart (1, x, at, 1)
        x = x * alat * BOHR_RADIUS_ANGS
        write(ounit,'(3F15.8)') x
      enddo
    enddo
  enddo

  ! vectors
  write(ounit,'(A,I10)') 'POINT_DATA', nr1*nr2*nr3
  write(ounit,'(A)') 'VECTORS vectors float'
  do i3 = 1, nr3
    do i2 = 1, nr2
      do i1 = 1, nr1
        bla = i1 + (i2-1)*dfftp%nr1x + (i3-1)*dfftp%nr1x*dfftp%nr2x
        write(ounit,'(3E15.8)') vin(bla,1:3)
      enddo
    enddo
  enddo

end subroutine vtk_vector_3d


!-----------------------------------------------------------------------
SUBROUTINE write_nics(filename, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the NICS in PP postproc format
  !
  USE kinds,           ONLY : dp
  USE io_global,       ONLY : stdout, ionode
  USE fft_base,        ONLY : dfftp, grid_gather
  USE gvect,           ONLY : gcutm
  USE cell_base,       ONLY : at, alat, tpiba2, omega, ibrav, celldm
  USE ions_base,       ONLY : zv, ntyp => nsp, nat, ityp, atm, tau
  USE wvfct,           ONLY : ecutwfc
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  implicit none
  character*(*) filename
  real(dp) :: field(dfftp%nnr,3,3,nspin)
  !-- local variables ----------------------------------------------------
  character(75), parameter :: title = 'NICS'
  real(dp), allocatable :: nics(:), aux(:)
  integer :: ispin

  allocate(nics(dfftp%nnr))
  nics = 0.d0
  do ispin = 1,nspin
    nics(:) = nics(:) + (field(:,1,1,ispin) + field(:,2,2,ispin) + &
                         field(:,3,3,ispin))/3.d0
  enddo
  nics = nics * 1d6

  ! gather the data 
  allocate(aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
#ifdef __MPI
  call grid_gather(nics, aux)
#else
  aux = nics
#endif

  if (ionode) then
      write(stdout, '(5X,''writings NICS to: '',A40)') trim(filename)
      call plot_io(filename, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
                   celldm, at, gcutm, dual, ecutwfc, 100, atm, ityp, zv, &
                   tau, aux, 1)
  endif

  deallocate(aux, nics)

end subroutine write_nics


!  NEW_part  ---- highly experimental ------------------

!-----------------------------------------------------------------------
SUBROUTINE write_tensor_field_xsf(name, ipol, ispin, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the tensor field in xcrysden format
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout, ionode
  USE mp_global,                   ONLY : me_pool
  USE cell_base,                   ONLY : at, bg, alat
  USE ions_base,                   ONLY : nat, tau, atm, ityp
  !USE something                   ONLY : nr1x, nr2x, nr3x, nr1, nr2, nr3
       USE fft_base, only : dffts, dfftp, grid_gather
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  character*(*) name
  integer :: ispin
!  real(dp) :: field(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x,3)
  real(dp) :: field(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3)
  !--------------------------------------------------------------------
  integer, parameter :: ounit = 48
  character*80 :: fname
  integer :: ios, ipol

  real(dp) :: vector_field(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3)                                                                                              
#ifdef __MPI
  ! gather the data if necessary
  do ipol = 1, 3
     call grid_gather(field(:,ipol), vector_field(:,ipol))
  enddo
#else
  vector_field = field
#endif

!  vector_field = field

    ! TEST !
  nr1 = dfftp%nr1
  nr2 = dfftp%nr2
  nr3 = dfftp%nr3

  nr1x = dfftp%nr1x
  nr2x = dfftp%nr2x
  nr3x = dfftp%nr3x


  if (me_pool /= 0) return


  if (ionode) then
!  do ipol = 1, 3
    ! form the name (and info) of the output file
    if (ispin == 0) fname = trim(name)//''
    if (ispin == 1) fname = trim(name)//'_UP'
    if (ispin == 2) fname = trim(name)//'_DW'

    if (ipol == 0) fname = trim(fname)//'.xsf'
    if (ipol == 1) fname = trim(fname)//'_X.xsf'
    if (ipol == 2) fname = trim(fname)//'_Y.xsf'
    if (ipol == 3) fname = trim(fname)//'_Z.xsf'
    write(stdout, '(5X,"write_tensor_field: ",A40)') fname

    open(unit=ounit, file=fname, iostat=ios, form='formatted', &
         status='unknown')
    if (ios /= 0) &
      call errore('write_tensor_field', 'error opening '//fname, ounit)

    call xsf_struct (alat, at, nat, tau, atm, ityp, nr1*nr2*nr3, ounit)
    call xsf_vector_3d(vector_field, &
                  nr1, nr2, nr3, nr1x, nr2x, nr3x, at, bg, alat, ounit)
    close(unit=ounit)
!  enddo
  endif
end subroutine write_tensor_field_xsf

!
! Copyright (C) 2003 Tone Kokalj
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
!
! -------------------------------------------------------------------
!   this routine writes the crystal structure in XSF format
! -------------------------------------------------------------------
subroutine xsf_struct (alat, at, nat, tau, atm, ityp, nr, ounit)
  USE kinds, only : DP
  implicit none
  integer          :: nat, ityp (nat), nr, ounit
  character(len=3) :: atm(*)
  real(DP)    :: alat, tau (3, nat), at (3, 3)
  ! --
  integer          :: i, j, n
  real(DP)    :: at1 (3, 3)
  ! convert lattice vectors to ANGSTROM units ...
  do i=1,3
     do j=1,3
        at1(j,i) = at(j,i)*alat*0.529177d0
     enddo
  enddo

  write(ounit,*) 'CRYSTAL'
  write(ounit,*) 'PRIMVEC'
  write(ounit,'(2(3F15.9/),3f15.9)') at1
  write(ounit,*) 'PRIMCOORD'
  write(ounit,*) nat + nr, 1

  do n=1,nat
     ! positions are in Angstroms
     write(ounit,'(a3,3x,3f15.9)') atm(ityp(n)), &
          tau(1,n)*alat*0.529177d0, &
          tau(2,n)*alat*0.529177d0, &
          tau(3,n)*alat*0.529177d0
  enddo
  return
end subroutine xsf_struct




! -------------------------------------------------------------------
!   this routine writes a 3D vector field
! -------------------------------------------------------------------
subroutine xsf_vector_3d(vfield, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                         at, bg, alat, ounit)
  USE kinds, only : DP
  implicit none
  integer  :: nr1x, nr2x, nr3x, nr1, nr2, nr3, ounit
!  real(DP) :: vfield(nr1x*nr2x*nr3x, 3)
  real(DP) :: vfield(nr1x,nr2x,nr3x, 3)
  real(DP) :: x(3), at(3,3), bg(3,3), alat
  integer  :: i1, i2, i3, i_point

  do i1 = 1, nr1
    do i2 = 1, nr2
      do i3 = 1, nr3
        !i_point = i1 + (i2-1)*nr1x + (i3-1)*nr1x*nr2x
        ! coordinate in angstrom
        x(1) = dble(i1-1)/dble(nr1)     
        x(2) = dble(i2-1)/dble(nr2)     
        x(3) = dble(i3-1)/dble(nr3)
        ! crystal to cartesian
        call trnvect (x, at, bg, 1)
        x = x * alat * 0.529177d0
        !write(ounit,'(''X  '',3x,3f15.9,2x,3e12.4)') x, vfield(i_point, 1:3)
        write(ounit,'(''X  '',3x,3f15.9,2x,3e12.4)') x, vfield(i1,i2,i3, 1:3)
      enddo
    enddo
  enddo
  return
end subroutine xsf_vector_3d



!-----------------------------------------------------------------------
subroutine trnvect (vect, at, bg, iflag)
!-----------------------------------------------------------------------
  !
  !  This routine transforms a vector (like forces which in the
  !  crystal axis is represented on the basis of the reciprocal lattice
  !  vectors) from crystal to cartesian axis (iflag.gt.0)
  !  and viceversa (iflag.le.0)
  !
  USE kinds
  implicit none
  integer :: iflag
  ! input: gives the versus of the transformati

  real(DP) :: vect (3), at (3, 3), bg (3, 3)
  ! inp/out: the vector to transform
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  real(DP) :: work (3)
  ! a working array

  integer :: ipol, ialpha
  ! counter on crystal coordinates
  ! counter on cartesian coordinates
  if (iflag.gt.0) then
     !
     !     forward transformation, from crystal to cartesian axis
     !
     do ipol = 1, 3
        work (ipol) = vect (ipol)
     enddo
     do ialpha = 1, 3
        vect (ialpha) = 0.d0
        do ipol = 1, 3
           vect (ialpha) = vect (ialpha) + work (ipol) * bg (ialpha, ipol)
        enddo
     enddo
  else
     !
     !    backward transformation, from cartesian to crystal axis
     !
     do ipol = 1, 3
        work (ipol) = 0.d0
        do ialpha = 1, 3
           work (ipol) = work (ipol) + vect (ialpha) * at (ialpha, ipol)
        enddo
     enddo
     do ipol = 1, 3
        vect (ipol) = work (ipol)
     enddo
  endif
  return
end subroutine trnvect

