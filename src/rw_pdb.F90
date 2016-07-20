!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2016 by Paolo Raiteri                               !
!                                                                   !
! p.raiteri@curtin.edu.au                                           !
!                                                                   !
! This program is free software; you can redistribute it and/or     !
! modify it under the terms of the GNU General Public License       !
! as published by the Free Software Foundation; either version 2    !
! of the License, or (at your option) any later version.            !
!                                                                   !
! This program is distributed in the hope that it will be useful,   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of    !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !
! GNU General Public License for more details.                      !
!                                                                   !
! For the full text of the GNU General Public License,              !
! write to: Free Software Foundation, Inc.,                         !
!           675 Mass Ave, Cambridge, MA 02139, USA.                 !
!                                                                   !
! The GNU GPL can also be found at http://www.gnu.org               !
!                                                                   !
! No claim is made that this program is free from errors and        !
! no liability will be accepted for any loss or damage that         !
! may result. The user is responsible for checking the validity     !
! of their results.                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_natom_pdb(uinp,natoms)
  use variables, only : ip
  implicit none
  integer(ip), intent(in)  :: uinp
  integer(ip), intent(out) :: natoms
  character(len=80)  :: line
  character(len=10)  :: str

#ifdef DEBUG
  write(0,*)"Entering subroutine :: get_natom_pdb"
#endif

  natoms=0
  do
    read(uinp,'(a80)',end=546)line
    read(line,*)str
    if(str=="ATOM" .or. str=="HETATM")then
      natoms=natoms+1
    elseif(str=="END")then
      exit
    else
      cycle
    endif
  enddo
546 continue
#ifdef DEBUG
  write(0,*)" Leaving subroutine :: get_natom_pdb"
#endif
  return
end subroutine get_natom_pdb

subroutine read_pdb(uinp,natoms,pos,label,chg,cell,go)
  use variables, only : ip, dp, cp, pi
  implicit none
  integer(ip), intent(in) :: uinp
  integer(ip), intent(in) :: natoms
  real(dp), dimension(3,natoms), intent(out) :: pos
  character(cp), dimension(natoms), intent(out) :: label
  real(dp), dimension(natoms), intent(out) :: chg
  real(dp), dimension(6), intent(out) :: cell
  logical, intent(out) :: go

  integer(ip) :: iatom
  real(dp) :: rtmp1, rtmp2
  character(len=100)  :: pdb_line
  character(len=10)  :: str
  logical, save :: first_time_in=.true.
  integer :: ios

#ifdef DEBUG
  write(0,*)"Entering subroutine :: read_pdb"
#endif

  iatom=0
  go=.true.
  cell=0.0_dp
  do
    read(uinp,'(a100)',end=222,err=222)pdb_line
    if (len_trim(pdb_line)==0) goto 222
    read(pdb_line,*)str
    if(str=="ATOM" .or. str=="HETATM")then
      iatom=iatom+1
! Positions & labels
      read(pdb_line(13:16),*)label(iatom) 
      read(pdb_line(31:54),'(3f8.3)')pos(:,iatom)
      chg(iatom)=0.0d0
      read(pdb_line(81:99),'(3f8.3)',end=100,err=100)chg(iatom)
100   continue
! Cell
    elseif(str == "CRYST1") then
      read(pdb_line,*)str,cell
      cell(4) = cell(4) * pi / 180.0_dp
      cell(5) = cell(5) * pi / 180.0_dp
      cell(6) = cell(6) * pi / 180.0_dp
    elseif(str=="END")then
      exit
    else
      cycle
    endif
  enddo
222 continue
  if(iatom/=natoms)go=.false.

#ifdef DEBUG
  write(0,*)" Leaving subroutine :: read_pdb"
#endif
  return
end subroutine read_pdb

subroutine write_pdb(uout,natoms,pos,label,chg,hmat)
  use variables, only : ip, dp, cp, pi
  implicit none
  integer(ip), intent(in) :: uout
  integer(ip), intent(in) :: natoms
  real(dp), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms), intent(in) :: label
  real(dp), dimension(natoms), intent(in) :: chg
  real(dp), intent(in) :: hmat(3,3)

  integer(ip) :: i, j
  character(len=100)  :: line
  real(dp) :: cell(6)
  real(dp)  :: occ, beta, tmp

#ifdef DEBUG
  write(0,*)"Entering subroutine :: write_pdb"
#endif

  occ=1.0_dp
  beta=0.0_dp
  write(uout,'(a6)')"REMARK"
  call get_cell(hmat,cell,"DEG")
  write(uout,'(a6,3f9.3,3f7.2)')"CRYST1",cell

  do i=1,natoms
    line=" "

    write(line( 1: 6),'(a6   )')"ATOM  "          ! a literal "ATOM  " (note two trailing spaces).
    if (natoms<=99999) then
      write(line( 7:11),'(i5   )')i               ! atom serial number, e.g. "   86". 
    else
      write(line( 7:11),'(z5   )')i               ! NON-STANDARD atom serial number, e.g. "   86". 
    endif
    write(line(13:16),'(a4   )')adjustr(label(i)) ! Atom role name, e.g. " CG1;". 

    write(line(18:22),'(a4,1x)')"UNK "            ! amino acid abbreviation, e.g. "ARG". 
    line(23:27) = line( 7:11)                     ! residue sequence number (I4) and insertion code (A1), e.g. "  11 " or " 256C". 

    write(line(31:54),'(3f8.3)')pos(:,i)          ! atom coordinates (X, Y, Z)
    write(line(55:60),'(f6.2 )')occ               ! atom occupancy, usually "  1.00".
    write(line(61:66),'(f6.2 )')beta              ! B value or temperature factor.
    write(line(81:100),'(f12.8)')chg(i)
    write(uout,'(a100)')line
  enddo
  write(uout,'("END")')

#ifdef DEBUG
  write(0,*)" Leaving subroutine :: write_pdb"
#endif

  return
end subroutine write_pdb
