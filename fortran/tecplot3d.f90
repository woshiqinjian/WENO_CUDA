program main
implicit none
integer, parameter::  nx=32
integer, parameter::  ny=32
integer, parameter::  nz=32
double precision u(nx,ny,nz)
double precision v(nx,ny,nz)
double precision w(nx,ny,nz)

integer i ,j ,k, ii, jj, kk
double precision xd, yd, zd, ut, vt, wt

xd = 240
yd = 18
zd = 12

open(10,file='../build/HIT_0032_4.0')
do k=1, nz, 1
  do j=1, ny, 1
    do i=1,nx, 1
      read(10,*) u(i,j,k), v(i,j,k), w(i,j,k)
    enddo 
  enddo
enddo
close(10)

open(10,file='tec_uvw.dat')
write(10,*) ' TITLE = "Example: Simple 3D Volume Data" '
write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
write(10,*) ' ZONE I=',nx, ', J=',ny,', K=',nz,', DATAPACKING=POINT '
!write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
!write(10,*) ' ZONE I=',nx, ', J=',ny,', DATAPACKING=POINT '
do k=1, nz
  do j=1, ny, 1
    do i=1, nx, 1
      write(10,*) i,j,k, u(i,j,k), v(i,j,k), w(i,j,k)
    enddo
  enddo
enddo
close(10)


stop
end
