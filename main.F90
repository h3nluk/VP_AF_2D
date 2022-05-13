#include "constants.h"

!################################INITIALIZATION#########################

module parameters
  implicit none 
  
  integer, parameter :: dimX = 32
  integer, parameter :: dimV = 32
  
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
  real(kind=DTYPE), parameter :: tmax = 0.01
  
  
  real(kind=DTYPE), parameter :: me = 1.
  real(kind=DTYPE), parameter :: mi = 1.
  real(kind=DTYPE), parameter :: qe = -1.
  real(kind=DTYPE), parameter :: qi = 1.
  real(kind=DTYPE), parameter :: Te = 1.
  real(kind=DTYPE), parameter :: Ti = 1.
 
end module parameters

subroutine init(fe,fi,v, xb, xe, vb, ve, sizex,sizev)
  implicit none
  
  integer :: sizex, sizev
  real(kind=DTYPE) :: xe,xb,ve,vb
  !real(kind=DTYPE) :: fe_val(0:sizex(1),0:sizex(2),0:sizex(3),0:sizev(1),0:sizev(2),0:sizev(3))
  !real(kind=DTYPE) :: fi_val(0:sizex(1),0:sizex(2),0:sizex(3),0:sizev(1),0:sizev(2),0:sizev(3))
  
  real(kind=DTYPE) :: fe(-2*B:sizex+2*B, 0:sizev)
  real(kind=DTYPE) :: fi(-2*B:sizex+2*B, 0:sizev)
  
  real(kind=DTYPE) :: v(0:sizev)
  
  integer :: ix,iv
  real(kind=DTYPE) :: x,vx,v2,v0
  real(kind=DTYPE) :: dx_val,dv_val
  real(kind=DTYPE) :: m_i, m_e, T_e, T_i
  
  m_i = 1.
  m_e = 1.
  T_e = 1.
  T_i = 1.
  
  dx_val = (xe-xb)/sizex
  dv_val= (ve-vb)/sizev
  
  do ix=-2*B,sizex+2*B

    x = xb + ix*dx_val
    
    do iv=0,sizev
    
      vx = vb + iv*dv_val
      v(iv) = vx
      
      v2 = vx**2
      
      !linear LD
      fe(ix,iv) = sqrt(1/(2*M_PI))*exp(-0.5*v2)*(1+0.01*cos(0.5*x))
      fi(ix,iv) = sqrt(1/(2*M_PI))*exp(-0.5*v2)*(1+0.00*cos(0.5*x))
      
    
    enddo

  enddo

end subroutine init

!#########################MAIN##########################################

program AF_Vlasov
  
  use parameters
  use advection
  use poisson
  use output
  
  implicit none
  
  real(kind=DTYPE), dimension(:,:), allocatable :: fe, fehalf, feold
  real(kind=DTYPE), dimension(:,:), allocatable :: fi
  real(kind=DTYPE), dimension(:), allocatable :: v
  real(kind=DTYPE), dimension(:), allocatable :: E, Ehalf, Eold
  
  real(kind=DTYPE), dimension(:), allocatable :: Phi
  real(kind=DTYPE), dimension(:), allocatable :: Rho

  real(kind=DTYPE) :: t = 0.
    
  !integer :: nproc(3)
  !character(len=256) :: outputDirectory
  !nproc(:) = 1
  !outputDirectory = 'AF_Vlasov_Vtki_Test'
  
  !Initialization
  write(*,*) '-------------------------------'
  write (*,*) 'AF_Vlasov'
  write(*,*) '-------------------------------'
  write (*,*) 'Init'
  
  allocate(fe(-2*B:sizex+2*B,0:sizev))
  allocate(fehalf(-2*B:sizex+2*B,0:sizev))
  allocate(feold(-2*B:sizex+2*B,0:sizev))
  allocate(fi(-2*B:sizex+2*B,0:sizev))
  allocate(v(0:sizev))
  allocate(E(-2*B:sizex+2*B))
  allocate(Ehalf(-2*B:sizex+2*B))
  allocate(Eold(-2*B:sizex+2*B))
  allocate(Phi(-2*B:sizex+2*B))
  allocate(Rho(-2*B:sizex+2*B))
  
  call init(fe,fi,v,xb,xe,vb,ve,sizex,sizev)
  
  !call get_E_gauss_seidel(fe,fi,rho,phi,E,sizex,sizev,0.5*dx,0.5*dv,qe,qi)
  
  call boundary_kinetic(fe,sizex,sizev)
  call boundary_kinetic(fi,sizex,sizev)
  call boundary_fluid(E,sizex)
  
  call output_kinetic(fe,sizex,sizev,0)
  
  write(*,*) 'Main Loop'
  write(*,*) '-------------------------------'
  
  do while (t < tmax)
  
    feold(:,:)  = fe(:,:)
    fehalf(:,:) = fe(:,:)
    
    !call get_E_gauss_seidel(fe,fi,rho,phi,E,sizex,sizev,0.5*dx,0.5*dv,qe,qi)
    call boundary_fluid(E,sizex)
    
    !half timestep interface update
    call Evolution(fehalf,sizex,sizev,v,E,dx,dv,0.5*dt)
    call boundary_kinetic(fehalf,sizex,sizev)
    
    !full timestep interface update
    call Evolution(fe,sizex,sizev,v,E,dx,dv,dt)
    call boundary_kinetic(fe,sizex,sizev)

    !update cell average
    call Conservation(fe,fehalf,feold,v,E,E,E,dt,dx,dv,sizex,sizev)
    call boundary_kinetic(fe,sizex,sizev)
    
    !Output
    call outputDat(E,t,dimX,sizex)
    call output_max(E,t,sizex)
    !call outputVtkDistributionFunction(fe,fi,dimX,dimV,nproc,me,outputDirectory)
    
    t = t + dt
    write (*,*) "t = ", t
    write (*,*) '--------------------------------------'
    
    write (*,*) 'max_E', maxval(E)
    write (*,*) 'min_E', minval(E)
    write (*,*) 'max_fe', maxval(fe)
    write (*,*) 'min_fe', minval(fe)
    
  write(*,*) '-----------------------------------------'
  write(*,*) ' '
    
  enddo
  
  write (*,*) 'Complete'
!~   write (*,*) '--------------------------------------'
!~   write (*,*) 'fe_shape:', shape(fe)
!~   write (*,*) 'fi_shape:', shape(fi)
!~   write (*,*) 'E_shape:', shape(E)
!~   write (*,*) 'v_shape:', shape(v)
  write (*,*) '--------------------------------------'
  write (*,*) 'max_v', maxval(v)
  write (*,*) 'min_v', minval(v)
  write (*,*) 'max_E', maxval(E)
  write (*,*) 'min_E', minval(E)
  write (*,*) 'max_fe', maxval(fe)
  write (*,*) 'min_fe', minval(fe)
  write (*,*) 'max_fi', maxval(fi)
  write (*,*) 'min_fi', minval(fi)
  
  call output_kinetic(fe,sizex,sizev,1)
  call output_fluid(E,sizex)
  call output_scalar(Phi,sizex,2,"Phi.dat")
  call output_scalar(Rho,sizex,3,"Rho.dat")
  
  deallocate(fe)
  deallocate(fehalf)
  deallocate(feold)
  deallocate(fi)
  deallocate(v)
  deallocate(E)
  deallocate(Ehalf)
  deallocate(Eold)
  deallocate(Phi)
  deallocate(Rho)
  
end program 
