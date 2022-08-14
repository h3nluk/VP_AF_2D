#include "constants.h"

module output

implicit none
	
contains
		
subroutine output_max(E,time,sizex)

  implicit none 
  integer :: sizex
  real(kind=DTYPE) :: E(-2*B:sizex+2*B), time
  
  open(10, file="Emax.dat", status="unknown", position="append")
  write(10,*) time, maxval(E)
  close(10)

end subroutine output_max

subroutine output_kinetic(fs,sizex,sizev,a)
  
  implicit none 
  integer :: sizex, sizev, a
  real(kind=DTYPE) :: fs(-2*B:sizex+2*B,-2*B:sizev+2*B)
  
  integer :: i, j
  
  if (a == 0) then
    do i=-2*B, sizex+2*B
      open(11, file="f0.dat", status ="unknown", position="append")
      write(11,*) fs(i,:)
      close(11)
    end do
  
  else if (a == 1) then
    do i=-2*B, sizex+2*B 
      open(10, file="f.dat", status ="unknown", position="append")
      write(10,*) fs(i,:)
      close(10)
    end do
  
  endif
  
end subroutine output_kinetic

subroutine output_fluid(E,sizex)
  
  implicit none 
  integer :: sizex
  real(kind=DTYPE) :: E(-2*B:sizex+2*B)
  
  integer :: i
  
  do i=-2*B,sizex+2*B
    !output 1D electric field Ex
    open(10, file="E.dat", status="unknown", position="append")
    write(10,*) E(i)
    close(10)
  end do
  
end subroutine output_fluid

subroutine output_scalar(phi, sizex, fileid, filename)

  implicit none 
  integer :: sizex
  real(kind=DTYPE) :: phi(-2*B:sizex+2*B)
  !character(len=10) :: fileid
  integer :: fileid
  character(len=25) :: filename
  
  integer :: i
  
  do i=-2*B, sizex+2*B
    open(fileid,file=trim(filename), status="unknown", position="append")
    write(fileid,*) phi(i)
    close(fileid)
  end do

end subroutine output_scalar

!~ subroutine output_fluid(E,sizex)

!~   implicit none 
!~   integer :: sizex(3)
!~   real(kind=DTYPE) :: E(-2*B:sizex(1)+2*B, -2*B:sizex(2)+2*B, -2*B:sizex(3)+2*B,3)
  
!~   call write_vti_cells()

!~ end subroutine output_fluid

subroutine outputDat(E,time,dimX,sizex)
    ! for Landau damping E output
    implicit none
    integer :: dimX,sizex
    real(kind=DTYPE)  :: E(-2*B:sizex+2*B), time

    !character(len=*), intent(in) :: output_file
    !integer       :: ierror
    real(kind=DTYPE) :: Esq_sum
    real(kind=DTYPE) :: Esq(-2*B:sizex+2*B)
    !real(kind=DTYPE) :: Esq(-B:dimX-1+B)

!~     Esq(:,:,:) = E(:,:,:,1)*E(:,:,:,1) + &
!~                  E(:,:,:,2)*E(:,:,:,2) + &
!~                  E(:,:,:,3)*E(:,:,:,3)
                 
    !Esq(0:sizex:2) = E(0:sizex:2)*E(0:sizex:2)
    Esq(:) = sqrt(E(:)*E(:))
    !Esq(0:dimX-1) = sqrt(E(1:sizex-1:2)*E(1:sizex-1:2))
    
    Esq_sum = sum(Esq)
!~     call sum_parallel(Esq_sum)
!~     Esq_sum = Esq_sum/(nproc(1)*(dimX(1)+1))/(nproc(2)*(dimX(2)+1))/(nproc(3)*(dimX(3)+1))

!~     write (*,*) 'Esq_sum=', Esq_sum
!~     write (*,*) 'E_max=', maxval(E)

    open(10, file= "Esq.dat", status="unknown", position="append")
    write(10,*) sngl(time), sngl(Esq_sum)
    close(10)
    
end subroutine outputDat

end module output
