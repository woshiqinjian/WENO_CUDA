program main
implicit none
integer, parameter::  nx=192
integer, parameter::  ny=96
integer, parameter::  nz=32
double precision u(nx,ny,nz)
double precision v(nx,ny,nz)
double precision w(nx,ny,nz)

integer i ,j ,k, ii, jj, kk
double precision xd, yd, zd, ut, vt, wt

xd = 240
yd = 18
zd = 12

open(10,file='u_dump_it12000.dat')
read(10,*)
do i=1, nx, 1
  do j=1, ny, 1
    do k=1,nz, 1
      read(10,*) u(i,j,k)
    enddo 
  enddo
enddo
close(10)
open(10,file='v_dump_it12000.dat')
read(10,*)
do i=1, nx, 1
  do j=1, ny-1, 1
    do k=1,nz, 1
      read(10,*) v(i,j,k)
    enddo
  enddo
enddo
v(:,ny,:)= v(:,ny-1,:)
close(10)
open(10,file='w_dump_it12000.dat')
read(10,*)
do i=1, nx, 1
  do j=1, ny, 1
    do k=1,nz, 1
      read(10,*) w(i,j,k)
    enddo
  enddo
enddo
close(10)
open(10,file='tec_uvw.dat')
write(10,*) ' TITLE = "Example: Simple 3D Volume Data" '
write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
write(10,*) ' ZONE I=',nx-1, ', J=',1,', K=',nz-1,', DATAPACKING=POINT '
!write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
!write(10,*) ' ZONE I=',nx, ', J=',ny,', DATAPACKING=POINT '
do k=1, nz-1, 1
!   k = 32
!  do j=1, ny, 1
      j = 1
    do i=1,nx-1, 1
	  ut = ( u(i,j,k)+u(i+1,j,k) )/2.d0
	  vt =  v(i,j,k) / 2.d0
	  wt =  w(i,j,k) 
      write(10,*) xd/nx*i, yd/ny*j, zd/nz*k, ut, vt, wt
    enddo
  enddo
!enddo
close(10)

open(10,file='separation_tec_uvw.dat')
write(10,*) ' TITLE = "Example: Simple 3D Volume Data" '
write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
write(10,*) ' ZONE I=',16, ', J=',1,', K=',nz-1,', DATAPACKING=POINT '
!write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
!write(10,*) ' ZONE I=',nx, ', J=',ny,', DATAPACKING=POINT '
do k=1, nz-1, 1
!   k = 32
!  do j=1, ny, 1
      j = 1
    do i=45,60, 1
	  ut = ( u(i,j,k)+u(i+1,j,k) )/2.d0
	  vt =  v(i,j,k) / 2.d0
	  wt =  w(i,j,k) 
      write(10,*) xd/nx*i, yd/ny*j, zd/nz*k, ut, vt, wt
    enddo
  enddo
!enddo
close(10)

open(10,file='separation2_tec_uvw.dat')
write(10,*) ' TITLE = "Example: Simple 3D Volume Data" '
write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
write(10,*) ' ZONE I=',5, ', J=',1,', K=',nz-1,', DATAPACKING=POINT '
!write(10,*) ' VARIABLES = "X", "Y", "Z", "u", "v", "w" '
!write(10,*) ' ZONE I=',nx, ', J=',ny,', DATAPACKING=POINT '
do k=1, nz-1, 1
!   k = 32
!  do j=1, ny, 1
      j = 1
    do i=56,60, 1
	  ut = ( u(i,j,k)+u(i+1,j,k) )/2.d0
	  vt =  v(i,j,k) / 2.d0
	  wt =  w(i,j,k) 
      write(10,*) xd/nx*i, yd/ny*j, zd/nz*k, ut, vt, wt
    enddo
  enddo
!enddo
close(10)

stop
end
