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

program dislocation
  use variables
  implicit none

  real(dp),  dimension(3) :: bvec
  real(dp), dimension(3) :: d0, dvec, nvec
  real(dp), dimension(3) :: p0, p1, rtmp
  real(dp), dimension(3) :: dij, sij

  integer(ip) :: iatm, nat_per_mol, nmols, imol, idx
  real(dp), allocatable, dimension(:,:) :: pfrac
  real(dp), dimension(3) :: xcom
  real(dp) :: dotp
  real(dp) :: dispx, dispy, disp
  real(dp) :: d1, d2, dl, theta
  real(dp) :: x1, x2

! Tapering variables
  character(len=20) :: dtype
  real(dp) :: r1, r2
  integer(ip) :: n1, m1, n2, m2
  
  real(dp) :: cell(6), volume
  logical :: go
  integer(ip) :: inpfile=100
  character(len=100) :: inpfilename
  integer(ip) :: cunit=200
  character(len=100) :: coordfile
  integer(ip) :: outfile=100
  character(len=100) :: outfilename

  call getarg(1,inpfilename)
  
  open (inpfile,file=inpfilename,status='old')
  
  read(inpfile,*)dtype
  read(inpfile,*)bvec
  read(inpfile,*)d0
  read(inpfile,*)dvec
  read(inpfile,*)r1,n1,m1
  read(inpfile,*)coordfile
  read(inpfile,*)nat_per_mol
  read(inpfile,*)outfilename

#ifdef DEBUG
  write(0,*)dtype
  write(0,*)bvec
  write(0,*)d0
  write(0,*)dvec
  write(0,*)r1,n1,m1
  write(0,*)coordfile
#endif

  open (cunit,file=coordfile,status='old',form='formatted')
  open (outfile,file=outfilename,status='unknown',form='formatted')
  call get_natom_pdb(cunit,natoms)
  rewind(cunit)

#ifdef DEBUG
  write(0,*)natoms
#endif

  allocate(pos(3,natoms))
  allocate(lab(natoms))
  call read_pdb(cunit,natoms,pos,lab,cell,go)
  call get_hmat(cell,hmat)
  call get_3x3_inv(hmat,hinv,volume)

#ifdef DEBUG
  write(0,*)go
  write(0,'(6f12.5)')cell
  write(0,'(3f12.5)')hmat
  write(0,'(a4,3f12.5)')lab(1),pos(1:3,1)
  write(0,'(a4,3f12.5)')lab(2),pos(1:3,2)
#endif

! Dislocation line
  dvec = dvec / sqrt(sum(dvec*dvec))

! Unit vector normal to the dislocation line in fractional coordinated
  call vec_prod(bvec,dvec,nvec)

#ifdef DEBUG
  write(0,'(3f5.2)') bvec
  write(0,'(3f5.2)') dvec
  write(0,'(3f5.2)') nvec
#endif

! Transform the dislocation mid point in factional coordinates
  rtmp=d0
  d0(1) = hinv(1,1)*rtmp(1) + hinv(1,2)*rtmp(2) + hinv(1,3)*rtmp(3)
  d0(2) = hinv(2,1)*rtmp(1) + hinv(2,2)*rtmp(2) + hinv(2,3)*rtmp(3)
  d0(3) = hinv(3,1)*rtmp(1) + hinv(3,2)*rtmp(2) + hinv(3,3)*rtmp(3)

! Atomic position in fractional coordinates
  allocate(pfrac(3,natoms))
  do iatm=1,natoms
    pfrac(1,iatm) = hinv(1,1)*pos(1,iatm) + hinv(1,2)*pos(2,iatm) + hinv(1,3)*pos(3,iatm)
    pfrac(2,iatm) = hinv(2,1)*pos(1,iatm) + hinv(2,2)*pos(2,iatm) + hinv(2,3)*pos(3,iatm)
    pfrac(3,iatm) = hinv(3,1)*pos(1,iatm) + hinv(3,2)*pos(2,iatm) + hinv(3,3)*pos(3,iatm)
  enddo

  nmols=natoms/nat_per_mol
  
  iatm=0
  do imol=1,nmols

!! Molecules' centre of mass
    xcom=0.0_dp
    do idx=1,nat_per_mol
      xcom(1:3) = xcom(1:3) + pfrac(1:3,iatm+idx)
    enddo
    xcom(1:3) = xcom(1:3) / nat_per_mol

!!#ifdef DEBUG
!!    write(0,'(i6,3f10.5)')imol,xcom
!!#endif

! Vector distance between the molecule's c.o.m. and the mid point
    dij(1:3) = xcom-d0

! Angle with the dislocation line...
!      d1 = dot_product(dij,dvec) ! projection

! .. and its normal
!      d2 = dot_product(dij,nvec) ! projection
     
    sij(1:3) = dij(1:3) - dot_product(dij,bvec)*bvec(1:3)/sqrt(sum(bvec*bvec))
    rtmp=sij-dot_product(sij,dvec)/sqrt(sum(dvec*dvec))*dvec
    d2=sqrt(sum(rtmp*rtmp))
    dl=dot_product(rtmp,nvec)

! Back to cartesian distances to get the vertical shift
    dispx=sign(1.0_dp,dl)

    sij(1) = hmat(1,1)*dij(1) + hmat(1,2)*dij(2) + hmat(1,3)*dij(3)
    sij(2) = hmat(2,1)*dij(1) + hmat(2,2)*dij(2) + hmat(2,3)*dij(3)
    sij(3) = hmat(3,1)*dij(1) + hmat(3,2)*dij(2) + hmat(3,3)*dij(3)
    x1=sqrt(sum(sij**2))
    disp = dispx * (1. - (x1/r1)**n1) / (1. - (x1/r1)**m1)

!#ifdef DEBUG
!  write(123,'(10f14.5)')xcom(1:3),disp,disp*bvec(1:3)
!#endif

    do idx=1,nat_per_mol
      pfrac(1:3,iatm+idx) = pfrac(1:3,iatm+idx) + disp*bvec(1:3)
    enddo

    iatm=iatm+nat_per_mol
  enddo

  do iatm=1,natoms
    pos(1,iatm) = hmat(1,1)*pfrac(1,iatm) + hmat(1,2)*pfrac(2,iatm) + hmat(1,3)*pfrac(3,iatm)
    pos(2,iatm) = hmat(2,1)*pfrac(1,iatm) + hmat(2,2)*pfrac(2,iatm) + hmat(2,3)*pfrac(3,iatm)
    pos(3,iatm) = hmat(3,1)*pfrac(1,iatm) + hmat(3,2)*pfrac(2,iatm) + hmat(3,3)*pfrac(3,iatm)
  enddo

  call write_pdb(outfile,natoms,pos,lab,hmat)
  deallocate(pfrac)

  stop
end program dislocation
