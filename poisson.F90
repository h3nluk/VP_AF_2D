#include "constants.h"

module poisson
	
implicit none
	
contains
		
!~ real(kind=DTYPE) function reconstruction(L,R,avg,xi)
		
!~ 	implicit none
!~ 	real(kind=DTYPE) :: L,R,avg,xi
			
!~ 	reconstruction = L*(3.*xi**2 -4.*xi+1.)+avg*(6.*xi-6.*xi**2)+R*(3.*xi**2 -2.*xi)
			
!~ end function reconstruction

subroutine boundary_fluid(E,sizex)
	implicit none 
	integer :: sizex
	real(kind=DTYPE) :: E(-2*B:sizex+2*B)
			
	integer :: i
			
	do i=0, 2*B
		E(-i) = E(sizex-i)
		E(sizex+i) = E(i)
	end do
	
end subroutine boundary_fluid
		
subroutine get_E_schwarz(fe,fi,rho,phi,E,sizex,sizev,dx,dv,qe,qi)
  
	implicit none 
  
	integer :: sizex, sizev
	real(kind=DTYPE) :: fe(-2*B:sizex+2*B,0:sizev), fi(-2*B:sizex+2*B,0:sizev)
	real(kind=DTYPE) :: rho(-2*B:sizex+2*B), phi(-2*B:sizex+2*B), E(-2*B:sizex+2*B)
	real(kind=DTYPE) :: dx, dv, qe, qi
  
	integer :: ix, iv
  
	call calcRho(fe,fi,rho,sizex,sizev,dv,qe,qi)
			  
	call poisson_schwarz(rho,E,sizex,dx)
			  
end subroutine get_E_schwarz

subroutine get_E_gauss_seidel(fe,fi,rho,phi,E,sizex,sizev,dx,dv,qe,qi)
			  
	implicit none 
			  
	integer :: sizex, sizev
	real(kind=DTYPE) :: fe(-2*B:sizex+2*B,0:sizev), fi(-2*B:sizex+2*B,0:sizev)
	real(kind=DTYPE) :: rho(-2*B:sizex+2*B), phi(-2*B:sizex+2*B), E(-2*B:sizex+2*B)
	real(kind=DTYPE) :: dx, dv, qe, qi
			  
	integer :: iter
	real(kind=DTYPE) :: eps0 = 1.
			  
	call calcRho(fe,fi,rho,sizex,sizev,dv,qe,qi)
	call boundary_fluid(rho,sizex)
			  
	do iter=1,1000
		call poisson_gauss_seidel2(phi,rho,eps0,sizex,dx)
		call boundary_fluid(phi,sizex)
	enddo
			  
	call calc_E_from_phi2(E,phi,sizex,dx)
			  
	!check poisson solver
	E(:) = -1.*E(:)
			  
end subroutine get_E_gauss_seidel

subroutine calcRho(fe,fi,rho,sizex,sizev,dv,qe,qi)
	implicit none

	integer :: sizex,sizev
	real(kind=DTYPE) :: dv, qe, qi
	real(kind=DTYPE) :: rho(-2*B:sizex+2*B)
	real(kind=DTYPE) :: fe(-2*B:sizex+2*B,0:sizev), fi(-2*B:sizex+2*B,0:sizev)

	integer :: ix, iv

	do ix=0,sizex !0, mx
		rho(ix) = 0.
		do iv=1,sizev-1 !1, mvx
			rho(ix) = rho(ix) + qe*fe(ix,iv) + qi*fi(ix,iv)
		enddo
			rho(ix) = -rho(ix)*dv !TODO: Vorzeichen ???
	enddo
end subroutine calcRho

subroutine poisson_schwarz(rho,E,sizex,dx)
	implicit none

	integer :: sizex
	real(kind=DTYPE) :: dx
	real(kind=DTYPE) :: rho(-2*B:sizex+2*B), E(-2*B:sizex+2*B)

	integer :: ix
	real(kind=DTYPE), parameter :: epsilon = 1.e-6
	real(kind=DTYPE) :: vy,vq
	real(kind=DTYPE) :: rhoCopy(-2*B:sizex+2*B), &
				a(-2*B:sizex+2*B), b(-2*B:sizex+2*B), &
				x(-2*B:sizex+2*B), y(-2*B:sizex+2*B)

	rhoCopy(:) = rho(:)*dx**2

	a(0) = 1./(2.*(-2.-epsilon))
	do ix = 1,sizex-1
		a(ix) = 1. / (-2.-epsilon - a(ix-1))
	enddo

	b(0) = rhoCopy(0)/(2.*(-2.-epsilon))
	do ix = 1,sizex-1
		b(ix) = (rhoCopy(ix) - b(ix-1)) / (-2.-epsilon - a(ix-1))
	enddo
	b(sizex) = (rhoCopy(sizex) - b(sizex-1)) / ((-2.-epsilon + 1/(-2.-epsilon)) - a(sizex-1))

	x(sizex) = b(sizex)
	do ix = sizex-1,0,-1
		x(ix) = b(ix) - a(ix)*x(ix+1)
	enddo

	b(0) = -1./2.
	do ix = 1, sizex-1
		b(ix) = (-1.*b(ix-1)) / (-2.-epsilon - a(ix-1))
	enddo
	b(sizex) = (1. - b(sizex-1)) / ((-2.-epsilon + 1./(-2.-epsilon)) - a(sizex-1))

	y(sizex) = b(sizex)
	do ix = sizex-1, 0, -1
		y(ix) = b(ix) - a(ix)*y(ix+1)
	enddo

	vy = 1.*x(0) - 1./(-2.-epsilon)*x(sizex)
	vq = 1.*y(0) - 1./(-2.-epsilon)*y(sizex)

	do ix = 0, sizex
		x(ix) = x(ix) - vy/(1.+vq)*y(ix)
	enddo

	E(:) = 0.
	E(0) = -1.*(x(2) - x(sizex)) / (2.*dx)
	do ix = 1, sizex-1
		E(ix) = -1.*(x(ix+1) - x(ix-1)) / (2.*dx)
	enddo
	E(sizex) = -1.*(x(0) - x(sizex-2)) / (2.*dx)
	
	! E(-1) = E(mx)
	! E(mx+1) = E(0)

	! write(*,*) "Poisson: E(0)", E(0)
	! write(*,*) "Poisson: rho(0)", rho(0)
			
	! Periodic boundary conditions are not correct

end subroutine poisson_schwarz

subroutine calc_E_from_phi2(E,phi,sizex,dx) ! E = -grad(phi)
	
	!"Second order Central Finite Difference"
			  
	implicit none 
			  
	integer :: sizex
	real(kind=DTYPE) :: dx
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: phi
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: E
			  
	integer :: ix
			  
	do ix= 0, sizex
		E(ix) = (-phi(ix+1) + phi(ix-1))/(2.*dx)
	enddo
			  
end subroutine calc_E_from_phi2

subroutine calc_E_from_phi4(E,phi,sizex,dx) ! E = -grad(phi)
	
	!"Fourth order Central Finite Difference"
			  
	implicit none 
			  
	integer :: sizex
	real(kind=DTYPE) :: dx
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: phi
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: E
			  
	integer :: ix
			  
	do ix= 0, sizex
		!E(ix) = (-phi(ix+1) + phi(ix-1))/(2.*dx)
		E(ix) = -1.*(-phi(ix+2) + 8.*phi(ix+1) - 8.*phi(ix-1) + phi(ix-2))/(12.*dx)
	enddo
			  
end subroutine calc_E_from_phi4

subroutine poisson_gauss_seidel2(phi,rho,eps0,sizex,dx)

	!"Second order Finite Difference Poisson solver"

	implicit none
			  
	integer :: sizex
	real(kind=DTYPE) :: dx, eps0
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: phi
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: rho
			  
	integer :: ix
	real(kind=DTYPE) :: invdx2, suminv
			  
	invdx2 = 1./(dx*dx)
	suminv= 2.*invdx2
			  
	do ix= 0,sizex
		phi(ix) = ((phi(ix+1) + phi(ix-1))*invdx2 + rho(ix)/eps0)/suminv 
	enddo 

end subroutine poisson_gauss_seidel2

subroutine poisson_gauss_seidel4(phi,rho,eps0,sizex,dx)

	!"Fourth order Finite Difference Poisson solver"

	implicit none
			  
	integer :: sizex
	real(kind=DTYPE) :: dx, eps0
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: phi
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: rho
			  
	integer :: ix
	real(kind=DTYPE) :: invdx2, suminv
			  
	!invdx2 = 1./(dx*dx)
			  
	do ix= 0,sizex
		phi(ix) = 0.5*(phi(ix-1)+phi(ix+1)) + (dx**2 / 24.)*(rho(ix-1)+rho(ix+1)+10.*rho(ix))
	enddo 

end subroutine poisson_gauss_seidel4
	
end module poisson
