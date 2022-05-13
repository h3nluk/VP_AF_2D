#include "constants.h"

subroutine boundary(u,sizex,sizev)
  !periodic boundary condition
  implicit none

  integer :: sizex, sizev
  real(kind=DTYPE) :: u(-2*B:sizex+2*B,0:sizev)

  integer :: ix,iy

  !x-direction
  do ix=1,2*B
    u(-ix,:) = u(sizex-ix,:)
    u(sizex+ix,:) = u(ix,:)
  enddo
  
  !y-direction
!~   do iy=1,2*B
!~     u(:,-iy) = u(:,sizev-iy)
!~     u(:,sizev+iy) = u(:,iy)
!~   enddo
  
end subroutine boundary

subroutine Evolution(u,sizex,sizev,ax,ay,dx,dy,dt)
  
  implicit none 
  
  integer :: sizex, sizev
  real(kind=DTYPE) :: dx, dy, dt
  real(kind=DTYPE) :: ax(0:sizev) !v
  real(kind=DTYPE) :: ay(-2*B:sizex+2*B) !Ex
  real(kind=DTYPE) :: u(-2*B:sizex+2*B, 0:sizev) !f(x,v)
  
  integer :: i,j,k,a,b
  integer :: Nx,Ny,Sx,Sy,Wx,Wy,Ex,Ey,Cx,Cy
  integer :: x_shift(3), y_shift(3)
  real(kind=DTYPE) :: xi,eta,xi_eta,xisq,etasq,nodes,edges,bubble
  real(kind=DTYPE) :: xi_arr(3), eta_arr(3)
  real(kind=DTYPE) :: newinterfaces(-2*B:sizex+2*B, 0:sizev)
  
  !tracing
  
  xi_arr(1) = -1. !SW
  xi_arr(2) = -1. !W
  xi_arr(3) = 0.  !S
  eta_arr(1) = -1.
  eta_arr(2) = 0.
  eta_arr(3) = -1.
  
!~   x_shift(1) = -1
!~   x_shift(2) = -1
!~   x_shift(3) = 0
!~   y_shift(1) = -1
!~   y_shift(2) = 0
!~   y_shift(3) = -1
  
  x_shift(:) = int(xi_arr(:))
  y_shift(:) = int(eta_arr(:))
  
  do i=-1, sizex+1,2
  do j=3, sizev-3,2
  
    do k=1,3
    
      a = 0
      b = 0
      !characteristic origin
      xi = xi_arr(k) - dt*(2.*ax(j+y_shift(k))/dx)
      eta = eta_arr(k) - dt*(2.*ay(i+x_shift(k))/dy)
      
      !cell shift
      !x-direction
      if (xi < -1.) then 
      a = -2
      xi = xi + 2.
      else if(xi > 1.) then 
      a = 2
      xi = xi - 2.
      endif
      
      !y-direction
      if(eta < -1.) then
      b = -2
      eta = eta + 2.
      else if(eta > 1.) then 
      b = 2
      eta = eta - 2.
      endif  
      
      xi_eta = xi*eta
      xisq = xi**2
      etasq = eta**2
      
      Cx = i+a; Cy = j+b
      Nx = Cx; Ny = j+1+b
      Sx = Cx; Sy = j-1+b
      Wx = i-1+a; Wy = Cy
      Ex = i+1+a; Ey = Cy
      

      !reconstruction
      nodes = 0.25*xi*eta*(u(Wx,Sy)*(xi-1.)*(eta-1.) & 
            + u(Wx,Ny)*(xi-1.)*(eta+1.) & 
            + u(Ex,Ny)*(xi+1.)*(eta+1.) & 
            + u(Ex,Sy)*(xi+1.)*(eta-1.))
      
      edges = 0.5*(u(Cx,Sy)*(1.-xisq)*(eta-1.)*eta &
            + u(Ex,Cy)*(xi+1.)*(1.-etasq)*xi &
            + u(Cx,Ny)*(1.-xisq)*(eta+1.)*eta &
            + u(Wx,Cy)*(xi-1.)*(1.-etasq)*xi)
      
      bubble = (1./16.)*(36.*u(Cx,Cy)-(u(Wx,Sy)+u(Wx,Ny) &
             + u(Ex,Ny)+u(Ex,Sy)) &
             - 4.*(u(Wx,Cy)+u(Ex,Cy)+u(Cx,Ny)+u(Cx,Sy))) &
             * (1-xisq)*(1-etasq)

      newinterfaces(i+x_shift(k),j+y_shift(k)) = nodes + edges + bubble
!~       newinterfaces(i+x_shift(k),j+y_shift(k)) = k
      
    enddo
    
  enddo
  enddo
  
  !overwrite
  do i=1,sizex+1,2
  do j=1,sizev-1,2
  do k=1,3
    u(i+x_shift(k),j+y_shift(k)) = newinterfaces(i+x_shift(k),j+y_shift(k))
  enddo
  enddo
  enddo
  
  !test
  !u(:sizev) = 0.
  !u(:,sizev) = u(:,sizev-2)
  
end subroutine Evolution


subroutine Conservation(unew,uhalf,uold,ax,aynew,ayhalf,ayold,dt,dx,dy,sizex,sizev)
  
  implicit none 
  
  integer :: sizex, sizev
  real(kind=DTYPE) :: dx, dy, dt
  real(kind=DTYPE) :: ax(0:sizev) !v
  real(kind=DTYPE) :: aynew(-2*B:sizex+2*B)  !E^{n+1}(x)
  real(kind=DTYPE) :: ayhalf(-2*B:sizex+2*B) !E^{n+1/2}(x)
  real(kind=DTYPE) :: ayold(-2*B:sizex+2*B)  !E^{n}(x)
  real(kind=DTYPE) :: unew(-2*B:sizex+2*B, 0:sizev)  !f^{n+1}(x,v)
  real(kind=DTYPE) :: uhalf(-2*B:sizex+2*B, 0:sizev) !f^{n+1/2}(x,v)
  real(kind=DTYPE) :: uold(-2*B:sizex+2*B, 0:sizev)  !f^{n}(x,v)
  
  integer :: i, j
  real(kind=DTYPE) :: flux_left, flux_right, flux_top, flux_bottom
  real(kind=DTYPE) :: dt_dxdy
  real(kind=DTYPE) :: newavg(-2*B:sizex+2*B, 0:sizev)
  
  !update average
  dt_dxdy = dt/(dx*dy)
  
  do i=1,sizex-1,2
  do j=1,sizev-1,2
  
    !composite Simpsons rule
    flux_left = (1./36.)*((ax(j-1)*uold(i-1,j-1)+ax(j+1)*uold(i-1,j+1)+ax(j-1)*unew(i-1,j-1)+ax(j+1)*unew(i-1,j+1)) & 
    + 4.* (ax(j)*(uold(i-1,j)+unew(i-1,j))+ax(j-1)*uhalf(i-1,j-1)+ax(j+1)*uhalf(i-1,j+1)) & 
    + 16.*ax(j)*uhalf(i-1,j))
    
    flux_right = (1./36.)*((ax(j-1)*uold(i+1,j-1)+ax(j+1)*uold(i+1,j+1)+ax(j-1)*unew(i+1,j-1)+ax(j+1)*unew(i+1,j+1)) & 
    + 4.* (ax(j)*(uold(i+1,j)+unew(i+1,j))+ax(j-1)*uhalf(i+1,j-1)+ax(j+1)*uhalf(i+1,j+1)) & 
    + 16.*ax(j)*uhalf(i+1,j))
    
    flux_top = (1./36.)*((ayold(i-1)*uold(i-1,j+1)+ayold(i+1)*uold(i+1,j+1)+aynew(i-1)*unew(i-1,j+1)+aynew(i+1)*unew(i+1,j+1)) & 
    + 4.* (ayold(i)*uold(i,j+1)+aynew(i)*unew(i,j+1)+ayhalf(i-1)*uhalf(i-1,j+1)+ayhalf(i+1)*uhalf(i+1,j+1)) & 
    + 16.*ayhalf(i)*uhalf(i,j+1))
    
    flux_bottom = (1./36.)*((ayold(i-1)*uold(i-1,j-1)+ayold(i+1)*uold(i+1,j-1)+aynew(i-1)*unew(i-1,j-1)+aynew(i+1)*unew(i+1,j-1)) & 
    + 4.* (ayold(i)*uold(i,j-1)+aynew(i)*unew(i,j-1)+ayhalf(i-1)*uhalf(i-1,j-1)+ayhalf(i+1)*uhalf(i+1,j-1)) & 
    + 16.*ayhalf(i)*uhalf(i,j-1))
  
    !Finite Volume update
    newavg(i,j) = uold(i,j) - dt_dxdy &
    * (flux_top*dx &
    + flux_left*(-dy)&
    + flux_right*dy &
    + flux_bottom*(-dx))
  
    !unew(i,j) = newavg
  
  enddo
  enddo
  
  !overwrite
  do i=1,sizex-1,2
  do j=1,sizev-1,2
    unew(i,j) = newavg(i,j)
  enddo
  enddo
  
  !test
  !unew(:,sizev-1) = 0.
  
end subroutine Conservation
