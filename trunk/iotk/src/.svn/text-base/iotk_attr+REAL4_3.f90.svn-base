# 48 "iotk_attr.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_REAL4
#if 3 <= __IOTK_MAXRANK

# 158 "iotk_attr.spp"

# 233 "iotk_attr.spp"

# 236 "iotk_attr.spp"
subroutine iotk_write_attr_REAL4_3(attr,name,val,dummy,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  REAL(kind=iotk_REAL4), intent(in)  :: val (:,:,:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 257 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 265 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 282 "iotk_attr.spp"
  delim = '"'
# 286 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 288 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 289 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 293 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_REAL4_3

# 308 "iotk_attr.spp"
subroutine iotk_scan_attr_REAL4_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=iotk_REAL4)                        :: val (:,:,:)
#else
  REAL(kind=iotk_REAL4), intent(out)           :: val (:,:,:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  REAL(kind=iotk_REAL4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 339 "iotk_attr.spp"
  integer :: index
  REAL(kind=iotk_REAL4), allocatable :: tmpval (:)
# 342 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 353 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 360 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 366 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 371 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 380 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
  else
    goto 1
  end if
# 401 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 406 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 412 "iotk_attr.spp"
  if(index/=size(val)) then
# 414 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 420 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 422 "iotk_attr.spp"
  deallocate(tmpval)
# 424 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 428 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 441 "iotk_attr.spp"
    val = default
# 443 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_REAL4_3
# 451 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_REAL4_3
  write(0,*)
end subroutine iotk_attr_dummy_REAL4_3

# 45 "iotk_attr.spp"

# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_REAL4
#if 4 <= __IOTK_MAXRANK

# 158 "iotk_attr.spp"

# 233 "iotk_attr.spp"

# 236 "iotk_attr.spp"
subroutine iotk_write_attr_REAL4_4(attr,name,val,dummy,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  REAL(kind=iotk_REAL4), intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 257 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 265 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 282 "iotk_attr.spp"
  delim = '"'
# 286 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 288 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 289 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 293 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_REAL4_4

# 308 "iotk_attr.spp"
subroutine iotk_scan_attr_REAL4_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=iotk_REAL4)                        :: val (:,:,:,:)
#else
  REAL(kind=iotk_REAL4), intent(out)           :: val (:,:,:,:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  REAL(kind=iotk_REAL4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 339 "iotk_attr.spp"
  integer :: index
  REAL(kind=iotk_REAL4), allocatable :: tmpval (:)
# 342 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 353 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 360 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 366 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 371 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 380 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
  else
    goto 1
  end if
# 401 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 406 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 412 "iotk_attr.spp"
  if(index/=size(val)) then
# 414 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 420 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 422 "iotk_attr.spp"
  deallocate(tmpval)
# 424 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 428 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 441 "iotk_attr.spp"
    val = default
# 443 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_REAL4_4
# 451 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_REAL4_4
  write(0,*)
end subroutine iotk_attr_dummy_REAL4_4

# 45 "iotk_attr.spp"

# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_REAL4
#if 5 <= __IOTK_MAXRANK

# 158 "iotk_attr.spp"

# 233 "iotk_attr.spp"

# 236 "iotk_attr.spp"
subroutine iotk_write_attr_REAL4_5(attr,name,val,dummy,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  REAL(kind=iotk_REAL4), intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 257 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 265 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 265 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 282 "iotk_attr.spp"
  delim = '"'
# 286 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 288 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 289 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 293 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 295 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_REAL4_5

# 308 "iotk_attr.spp"
subroutine iotk_scan_attr_REAL4_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=iotk_REAL4)                        :: val (:,:,:,:,:)
#else
  REAL(kind=iotk_REAL4), intent(out)           :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  REAL(kind=iotk_REAL4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 339 "iotk_attr.spp"
  integer :: index
  REAL(kind=iotk_REAL4), allocatable :: tmpval (:)
# 342 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 353 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 353 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 360 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 366 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 371 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 380 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
      goto 1
    end if
  else
    goto 1
  end if
# 401 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 406 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
    goto 1
  end if
# 412 "iotk_attr.spp"
  if(index/=size(val)) then
# 414 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 414 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 414 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 420 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 422 "iotk_attr.spp"
  deallocate(tmpval)
# 424 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.4 ")
# 428 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 428 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 441 "iotk_attr.spp"
    val = default
# 443 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_REAL4_5
# 451 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_REAL4_5
  write(0,*)
end subroutine iotk_attr_dummy_REAL4_5

# 45 "iotk_attr.spp"

