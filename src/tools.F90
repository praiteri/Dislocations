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

subroutine get_cell(hmat,cell,angle)
  use variables, only : dp, pi
  implicit none
  real(dp) :: hmat(3,3)
  real(dp) :: cell(6)
  character(len=3)  :: angle

!!p a, b, c
  cell(1)= sqrt(dot_product(hmat(:,1),hmat(:,1)))
  cell(2)= sqrt(dot_product(hmat(:,2),hmat(:,2)))
  cell(3)= sqrt(dot_product(hmat(:,3),hmat(:,3)))
!!p alpha, beta, gamma in degrees
  if(angle=="DEG")then
    cell(4)=acos(dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3))/pi*180_dp
    cell(5)=acos(dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3))/pi*180_dp
    cell(6)=acos(dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2))/pi*180_dp
!!p alpha, beta, gamma in radians
  elseif(angle=="RAD")then
    cell(4)=acos(dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3))
    cell(5)=acos(dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3))
    cell(6)=acos(dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2))
!!p cosines of alpha, beta, gamma
  elseif(angle=="COS")then
    cell(4)=dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3)
    cell(5)=dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3)
    cell(6)=dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2)
  endif
  return
end subroutine

subroutine get_hmat(cell,h)
  use variables, only : dp
  implicit none
  real(dp) :: cell(6)
  real(dp) :: h(3,3)

!!p cell  :: a, b, c, alpha, beta, gamma in radians
  h = 0_dp
  h(1,1) = cell(1)
  h(1,2) = cell(2) * cos(cell(6))
  h(2,2) = cell(2) * sin(cell(6))
  h(1,3) = cell(3) * cos(cell(5))
  h(2,3) = (cell(3) * cell(2) * cos(cell(4)) - h(1,2) * h(1,3)) / h(2,2)
  h(3,3) = sqrt(cell(3)*cell(3) - h(1,3)*h(1,3) - h(2,3)*h(2,3))
!
  if ( abs(h(1,2)) < 1.e-4_dp ) h(1,2) = 0.0_dp
  if ( abs(h(1,3)) < 1.e-4_dp ) h(1,3) = 0.0_dp
  if ( abs(h(2,3)) < 1.e-4_dp ) h(2,3) = 0.0_dp
!
  return
end subroutine get_hmat

subroutine vec_prod(a,b,c)
  use variables, only : dp
  implicit none
  real(dp), dimension (3), intent(in)  :: a
  real(dp), dimension (3), intent(in)  :: b
  real(dp), dimension (3), intent(out) :: c

  c(1) =  a(2)*b(3) - a(3)*b(2)
  c(2) = -a(1)*b(3) + a(3)*b(1)
  c(3) =  a(1)*b(2) - a(2)*b(1)

  return
end subroutine vec_prod

subroutine get_3x3_inv (hmat,hmati,deth)
  use variables, only : dp
  implicit none
  real(dp), dimension(3,3)  :: hmat, hmati
  real(dp) :: deth
  real(dp) :: odet

  deth = &
       hmat(1,1) * ( hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2) ) + &
       hmat(1,2) * ( hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3) ) + &
       hmat(1,3) * ( hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1) )
  if ( deth < 1.e-4_dp ) return
  odet = 1_dp / deth
  hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
  hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
  hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
  hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
  hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
  hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
  hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
  hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
  hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet
  return
end subroutine get_3x3_inv

