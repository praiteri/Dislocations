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

module variables
  implicit none
!
! Precision
! cp = character length
! ip = integer precision
! dp = floating point precision
!
  integer, parameter :: cp = 4
  integer, parameter :: ip = 4 
  integer, parameter :: dp = 8 
  real(dp), parameter:: pi = 3.1415926535898_dp
!
! I/O
!
  integer(ip) :: io = 6
!
! atoms
!
  integer(ip) :: natoms
  real(dp), allocatable, dimension(:,:) :: pos
  character(cp), allocatable, dimension(:) :: lab
  real(dp), allocatable, dimension(:) :: chg
  real(dp), dimension(3,3) :: hmat, hinv
!
! Neighbours' list
!
  integer(ip), allocatable, dimension(:) :: nlist
  integer(ip), allocatable, dimension(:,:) :: llist

end module
