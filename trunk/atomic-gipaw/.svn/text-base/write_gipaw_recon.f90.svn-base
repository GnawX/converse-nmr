!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
subroutine write_gipaw_recon()
!--------------------------------------------------------------
      use ld1inc
      implicit none
      integer :: i, j, n, m, l, ios, iae, ncore, thesign

      ! set number of core orbitals
      ncore = nwf - nwfts
      print*, 'nwf  =', nwf
      print*, 'nwfts=', nwfts
      print*, 'ncore=', ncore
      do i = 1, nwfts
        write(*,'(I2,2X,A2,2X,''n='',I1,2X,''l='',I1)') &
              i, el(i+ncore), nn(i+ncore), ll(i+ncore)
        write(*,'(I2,2X,A2,2X,''n='',I1,2X,''l='',I1,4X,''r_c='',F10.4)') &
              i, elts(i), nnts(i), llts(i), rcutts(i)
      enddo

      ! core orbitals
      open(unit=51,file=trim(prefix)//'.gipaw_core',status='unknown',& 
           err=100,iostat=ios,form='formatted')
100   call errore('write_result','opening 51',abs(ios))
      write(51,'(A)') '<PP_GIPAW_CORE_ORBITALS>'
      write(51,*) ncore
      do i = 1, nwf - nwfts
        !if (core_state(i) == .false.) cycle
        write(51,'(A)') '<PP_GIPAW_CORE_ORBITAL>'
        write(51,*) nn(i), ll(i), 0, 0, oc(i)
        call check_sign(mesh, psi(:,1,i), thesign)
        write(51,'(1p4e19.11)') ( thesign*psi(n,1,i), n=1,mesh )
        write(51,'(A)') '</PP_GIPAW_CORE_ORBITAL>'
      enddo
      write(51,'(A)') '</PP_GIPAW_CORE_ORBITALS>'
      close(51)

      ! local potential
      open(unit=51,file=trim(prefix)//'.gipaw_vloc',status='unknown',& 
           err=101,iostat=ios,form='formatted')
101   call errore('write_result','opening 51',abs(ios))
      write(51,'(A)') '<PP_GIPAW_LOCAL_DATA>'
      write(51,'(A)') '<PP_GIPAW_VLOCAL_AE>'
      write(51,'(1p4e19.11)') ( r(n)*vpot(n,1), n=1,mesh )
      write(51,'(A)') '</PP_GIPAW_VLOCAL_AE>'
      write(51,'(A)') '<PP_GIPAW_VLOCAL_PS>'
      write(51,'(1p4e19.11)') ( r(n)*vpstot(n,1), n=1,mesh )
      write(51,'(A)') '</PP_GIPAW_VLOCAL_PS>'
      write(51,'(A)') '</PP_GIPAW_LOCAL_DATA>'
      do n = 1, mesh
        write(70,'(5(E12.6,2X))') r(n), r(n)*vpot(n,1), r(n)*vpstot(n,1)
      enddo
      close(51)


      ! GIPAW orbitals
      do i = 1, nwfts
        iae = ncore + i
        open(unit=51,file=trim(prefix)//'.gipaw_'//el(iae),status='unknown',& 
             err=102,iostat=ios,form='formatted')
102     call errore('write_result','opening 51',abs(ios))
        write(51,'(A)') '<PP_GIPAW_AE_ORBITAL>'
        write(51,*) el(iae), ll(iae)
        call check_sign(mesh, psi(:,1,iae), thesign)
        write(51,'(1p4e19.11)') ( thesign*psi(n,1,iae), n=1,mesh )
        write(80+i,*) '#', el(iae)
        do l = 1, mesh
          write(80+i,*) r(l), thesign*psi(l,1,iae)
        enddo
        write(51,'(A)') '</PP_GIPAW_AE_ORBITAL>'

        write(51,'(A)') '<PP_GIPAW_PS_ORBITAL>'
        write(51,*) rcutts(i), rcutusts(i)
        call check_sign(mesh, phits(:,i), thesign)
        write(51,'(1p4e19.11)') ( thesign*phits(n,i), n=1,mesh )
        write(90+i,*) '#', elts(i)
        do l = 1, mesh
          write(90+i,*) r(l), thesign*phits(l,i)
        enddo
        write(51,'(A)') '</PP_GIPAW_PS_ORBITAL>'
        close(51)
      enddo
      return

contains
   
      ! Check sign of wfct (max should be >0)
      subroutine check_sign(n, psi, thesign)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: psi(ndm)
      integer, intent(out) :: thesign
      integer :: i
      real(dp) :: wmax
           
      wmax = 0.0_dp
      do i = 1, n
        if ( abs(psi(i)) > wmax .and. r(i) < 4.0_dp) then
          wmax = abs(psi(i))
          if ( psi(i) < 0.0_dp) then
            thesign = -1
          else
            thesign = +1
          endif
        endif
      enddo
      end subroutine check_sign

end subroutine write_gipaw_recon
