#include "constants.h"

module parameters
  implicit none 
  
  integer, parameter :: dimX = 32
  integer, parameter :: dimV = 64
  
  integer, parameter :: sizex = 2*dimX
  integer, parameter :: sizev = 2*dimV

!~   integer, parameter, dimension(3) :: sizex = [2*(dimX(1)), 1, 1]
!~   integer, parameter, dimension(3) :: sizev = [2*(dimV(1)), 1, 1]

  real(kind=DTYPE), parameter :: PI = 4.D0*DATAN(1.D0)
  real(kind=DTYPE), parameter :: xe = 2*PI
  real(kind=DTYPE), parameter :: xb = -2*PI
  
!~   real(kind=DTYPE), parameter, dimension(3) :: xe = [2, 2, 2]
!~   real(kind=DTYPE), parameter, dimension(3) :: xb = [-2, -2, -2]
  
  real(kind=DTYPE), parameter :: ve = 4.5
  real(kind=DTYPE), parameter :: vb = -4.5
  
  real(kind=DTYPE), parameter :: dx = (xe-xb)/dimX
  real(kind=DTYPE), parameter :: dv = (ve-vb)/dimV
  
  real(kind=DTYPE), parameter :: dt = 0.01
  real(kind=DTYPE), parameter :: tmax = 45.
  
  
  real(kind=DTYPE), parameter :: me = 1.
  real(kind=DTYPE), parameter :: mi = 1.
  real(kind=DTYPE), parameter :: qe = -1.
  real(kind=DTYPE), parameter :: qi = 1.
  real(kind=DTYPE), parameter :: Te = 1.
  real(kind=DTYPE), parameter :: Ti = 1.
 
end module parameters
