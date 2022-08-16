#include "constants.h"

module poisson
	
implicit none
	
contains
		
real(kind=DTYPE) function reconstruction1D(f,i,ibeg,iend,xi) result (reco)
		
	implicit none
	
	integer :: i, ibeg, iend
	real(kind=DTYPE) :: f(ibeg:iend), xi
	
			
	reco = f(i-1)*(3.*xi**2 -4.*xi+1.)+f(i)*(6.*xi-6.*xi**2)+f(i+1)*(3.*xi**2 -2.*xi)
			
end function reconstruction1D


real(kind=DTYPE) function reconstruction2D(f,i,j,sizex,sizev,xi,eta) result (reco)

	implicit none
	
	integer :: i,j,sizex,sizev
	real(kind=DTYPE) :: f(-2*B:sizex+2*B, -2*B:sizev+2*B), xi, eta
	
	real(kind=DTYPE) ::  edges, nodes, bubble
	
	nodes =   f(i-1,j-1)*0.25*(xi-1.)*(eta-1.)*xi*eta &
			+ f(i+1,j-1)*0.25*(xi+1.)*(eta-1.)*xi*eta &
			+ f(i-1,j+1)*0.25*(xi-1.)*(eta+1.)*xi*eta &
			+ f(i+1,j+1)*0.25*(xi+1.)*(eta+1.)*xi*eta
	
	edges =   f(i+1, j)*0.5*(xi+1.)*(1.-eta**2)*xi &
			+ f(i-1, j)*0.5*(xi-1.)*(1.-eta**2)*xi &
			+ f(i ,j+1)*0.5*(eta+1.)*(1.-xi**2)*eta &
			+ f(i ,j-1)*0.5*(eta-1.)*(1.-xi**2)*eta 
	
	bubble = (1./16.)*(36.*f(i,j)-(f(i-1,j-1)+f(i+1,j-1)+f(i-1,j+1)+f(i+1,j+1))-4.*(f(i+1, j)+f(i-1, j)+f(i ,j+1)+f(i ,j-1)))*(1.-xi**2)*(1.-eta**2)
	
	reco = nodes+edges+bubble

end function reconstruction2D

subroutine boundary_fluid(E,sizex)
	implicit none 
	integer :: sizex
	real(kind=DTYPE) :: E(-2*B:sizex+2*B)
			
	integer :: i
			
	do i=1, 2*B
		E(-i) = E(sizex-i)
		E(sizex+i) = E(i)
	end do
	
end subroutine boundary_fluid
		

subroutine get_E_gauss_seidel(fe,fi,rho,phi,E,sizex,sizev,dx,dv,qe,qi)
			  
	implicit none 
			  
	integer :: sizex, sizev
	real(kind=DTYPE), dimension(-2*B:sizex+2*B,-2*B:sizev+2*B) :: fe, fi
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: rho, phi, E
	real(kind=DTYPE) :: dx, dv, qe, qi
			  
	integer :: iter
	real(kind=DTYPE) :: eps0 = 1.
			  
	call calcRho(fe,fi,rho,sizex,sizev,dv,qe,qi)
	call boundary_fluid(rho,sizex)
			  
	do iter=1,1000
		call poisson_gauss_seidel4(phi,rho,eps0,sizex,dx)
		call boundary_fluid(phi,sizex)
	enddo
			  
	call calc_E_from_phi4(E,phi,sizex,dx)
			  
	!check poisson solver
	E(:) = -1.*E(:)
			  
end subroutine get_E_gauss_seidel

subroutine get_E_gauss_seidel_point(fe,fi,rho,phi,E,sizex,sizev,dx,dv,qe,qi)
			  
	implicit none 
			  
	integer :: sizex, sizev
	real(kind=DTYPE), dimension(-2*B:sizex+2*B,-2*B:sizev+2*B) :: fe, fi, fe_point
	real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: rho, phi, E
	real(kind=DTYPE) :: dx, dv, qe, qi
			  
	integer :: iter
	real(kind=DTYPE) :: eps0 = 1.
	
	!calculate point values with local reconstruction
	call get_point_values(fe, fe_point, sizex, sizev)
	
	call calcRho(fe_point,fi,rho,sizex,sizev,dv,qe,qi)
	call boundary_fluid(rho,sizex)
			  
	do iter=1,1000
		call poisson_gauss_seidel4(phi,rho,eps0,sizex,dx)
		call boundary_fluid(phi,sizex)
	enddo
			  
	call calc_E_from_phi4(E,phi,sizex,dx)
			  
	!check poisson solver
	E(:) = -1.*E(:)
			  
end subroutine get_E_gauss_seidel_point

subroutine get_point_values(f,fpoint,sizex,sizev)

	implicit none
	
	integer :: sizex, sizev
	real(kind=DTYPE), dimension(-2*B:sizex+2*B,-2*B:sizev+2*B) :: f, fpoint
	
	integer :: i, j
	!real(kind=DTYPE) :: xi, eta
	
	!use reconstruction to get full grid of point values
	
	do i=0, sizex
		
		if (mod(i,2) .eq. 0) then 
			
			do j=0, sizev
				
				if (mod(j,2) .eq. 0) then 
					
					fpoint(i,j) = f(i,j) !nodes
					
				else if (mod(j,2) .eq. 1) then 
					
					!fpoint(i,j) = reconstruction1D(f(i,j-1), f(i,j+1), f(i,j), 0.5) !vertical edge
					fpoint(i,j) = reconstruction1D(f(i,:), j, 0, sizev, 0.5_8) !vertical edge
					
				endif
				
			enddo
			
		else if (mod(i,2) .eq. 1) then
			
			do j=0, sizev
				
				if (mod(j,2) .eq. 0) then 
					
					!fpoint(i,j) = reconstruction1D(f(i-1,j), f(i+1,j), f(i,j), 0.5) !horizontal edge
					fpoint(i,j) = reconstruction1D(f(:,j), i, -2*B, sizex+2*B, 0.5_8) !horizontal edge
					
				else if (mod(j,2) .eq. 1) then 
					
					fpoint(i,j) = reconstruction2D(f, i, j, sizex, sizev, 0._8, 0._8)
					
				endif
				
			enddo
			
		endif
		
	enddo
	
end subroutine get_point_values

subroutine calcRho(fe,fi,rho,sizex,sizev,dv,qe,qi)
	implicit none

	integer :: sizex,sizev
	real(kind=DTYPE) :: dv, qe, qi
	real(kind=DTYPE) :: rho(-2*B:sizex+2*B)
	real(kind=DTYPE), dimension(-2*B:sizex+2*B,-2*B:sizev+2*B) :: fe, fi

	integer :: ix, iv

	do ix=0,sizex !0, mx
		rho(ix) = 0.
		do iv=0,sizev-1 !1, mvx
			rho(ix) = rho(ix) + qe*fe(ix,iv) + qi*fi(ix,iv)
		enddo
			rho(ix) = -rho(ix)*dv !TODO: Vorzeichen ???
	enddo
end subroutine calcRho

!~ subroutine calcRho_reconstruction(fe,fi,rho,sizex,sizev,dv,qe,qi)
!~ 	implicit none

!~ 	integer :: sizex,sizev
!~ 	real(kind=DTYPE) :: dv, qe, qi
!~ 	real(kind=DTYPE) :: rho(-2*B:sizex+2*B)
!~ 	real(kind=DTYPE) :: fe(-2*B:sizex+2*B,0:sizev), fi(-2*B:sizex+2*B,0:sizev)

!~ 	integer :: ix, iv
	
!~ 	!reconstruction
	
!~ 	!nodes
!~ 	do ix=0, sizex, 2
		
!~ 		rho(ix) = 0.
		
!~ 		do iv=0, sizev-1
			
!~ 			rho(ix) = rho(ix) + qe*fe(ix,iv) + qi*fi(ix,iv)
			
!~ 		enddo
		
!~ 		rho(ix) = -rho(ix)*dv
		
!~ 	enddo
	
!~ 	!centers
!~ 	do ix=1, sizex-1, 2
	
!~ 		rho(ix) = 0.
		
!~ 		do iv=0, sizev-1
			
!~ 			if (mod(iv,2) .eq. 0) then !horizontal edges
				
!~ 				center_recon = reconstruction1D(fe(ix-1,iv),fe(ix+1,iv),fe(ix,iv),0.5)
				
!~ 			else if (mod(iv,2) .eq. 1) then !cell centers
				
!~ 				center_recon = &
!~ 				reconstruction2D(fe(ix,iv+1),fe(ix,iv-1),fe(ix-1,iv),fe(ix+1,iv),fe(ix-1,iv+1),fe(ix+1,iv+1),fe(ix-1,iv-1),fe(ix+1,iv-1),fe(ix,iv),0.,0.)
				
!~ 			endif
			
!~ 			rho(ix) = rho(ix) + qe*center_recon + qi*fi(ix,iv)
			
!~ 		enddo
		
!~ 		rho(ix) = -rho(ix)*dv
	
!~ 	enddo
	
!~ end subroutine calcRho_reconstruction


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
