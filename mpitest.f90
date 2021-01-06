program mpitest

! before compiling
! module load netCDF-Fortran/4.5.3-gmpich-2021.01

! Compile line
! mpifort -g -Wall -fcheck=all -o paralleltest parametersmod.f90 coordsmod.f90 outputmod.f90 averagemod.f90 mpitest.f90 -I/home/public/easybuild/software/netCDF-Fortran/4.5.3-gmpich-2021.01/include -L/home/public/easybuild/software/netCDF-Fortran/4.5.3-gmpich-2021.01/lib -lnetcdff

! Run line
! mpirun -np 2 ./paralleltest 1871/140 /home/terraces/datasets/dgvm_input/climate/transient1871-2010_list-formatted.nc test.nc

use parametersmod,  only : i1,i2,i4,sp,dp
use coordsmod,      only : index,parsecoords
use outputmod,      only : getoutfile,handle_err
use averagemod,     only : average,infompi
use netcdf
use mpi

implicit none

!--------------------
! Parallel program to read in list-formatted ncfile an calculate the average of each validcell across specified timeframe

integer :: rank
integer :: numtasks
integer :: ierr

integer :: infosize
integer :: sendsize

type(infompi) :: info

integer(i4), allocatable, dimension(:) :: srt
integer(i4), allocatable, dimension(:) :: cnt

integer(i1), allocatable, dimension(:) :: ob

integer(i4), dimension(2) :: job

real :: start_time
real :: end_time

!--------------------

call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

!--------------------

call CPU_TIME(start_time)

! ---

allocate(srt(numtasks))
allocate(cnt(numtasks))

srt = 0
cnt = 0

info%nproc = numtasks

infosize = sizeof(info)

sendsize = infosize + (8*numtasks)

allocate(ob(sendsize))

ob = 0

!--------------------

if (rank == 0) then

  call startmpi(info, srt, cnt)

  info%validcell = sum(cnt)

  call getoutfile(info%outfile, info%validcell)

  call infotobyte(info, srt, cnt, ob)

!   do i = 1, numtasks-1
!
!     call MPI_SEND(ob, sendsize, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
!
!   end do
!
!   call bytetoinfo(info, srt, cnt, ob)
!
!   write(0,*) 'srt from rank 0', srt
!
! else
!
!   call MPI_RECV(ob, sendsize, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)
!
!   write(0,*) 'timestring:', info%timestring, 'sendsize:', sendsize
!
!   call bytetoinfo(info, srt, cnt, ob)

end if

! ---
! Broadcast all info to all processes

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_BCAST(ob, (sendsize/4), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) ! For some reason, count is calculated in byte????

! ---

if (rank == 0) write(0,*) 'Broadcast complete'

if (rank /= 0) then

  call bytetoinfo(info, srt, cnt, ob)

end if

! ---

job = [srt(rank+1), cnt(rank+1)]

write(0,*) 'Rank', rank, 'recieved srt and cnt: ', job

!--------------------

call average(info, job, rank)

!---

call CPU_TIME(end_time)

if (rank == 0) write(0,*) "Total CPU time used by Rank 0:", end_time - start_time

!---

call MPI_FINALIZE(ierr)



!-------------------------------------------------------

contains

subroutine startmpi(info,srt,cnt)

type(infompi), target              , intent(inout) :: info
integer(i4),           dimension(:), intent(inout) :: srt
integer(i4),           dimension(:), intent(inout) :: cnt

integer :: i
integer :: n

integer :: status
integer :: ifid
integer :: dimid
integer :: ilen
integer :: tlen

!--- File variables
character(100), pointer :: infile
character(100), pointer :: outfile
character(100), pointer :: timestring
integer(i4)   , pointer :: nproc
integer(i4)   , pointer :: t0
integer(i4)   , pointer :: nt

integer :: len
integer :: remainder

type(index) :: timevals

integer(i4) :: baseyr
integer(i4) :: startyr
integer(i4) :: calcyrs
integer(i4) :: t1


infile     => info%infile
outfile    => info%outfile
timestring => info%timestring
nproc      => info%nproc
t0         => info%t0
nt         => info%nt


call getarg(2,infile)

status = nf90_open(infile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'index',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ilen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)


!--------------------

baseyr = 1871

call getarg(1,timestring)

call parsecoords(timestring,timevals)

startyr = timevals%minlon
calcyrs = timevals%minlat

t0 = 1 + 12 * (startyr - baseyr)
t1 = t0 + 12 * calcyrs - 1

nt = t1 - t0 + 1

!--------------------

len = floor(real(ilen) / nproc)

remainder = ilen - (len * nproc)

n = 1
do i = 1, nproc

  srt(n)   = (i-1) * len + 1
  cnt(n) = len

  if (i == nproc) then

    cnt(n) = len + remainder

  end if

  n = n + 1

end do

call getarg(3, outfile)


end subroutine startmpi

!-------------------------------------------------------

subroutine infotobyte(info,srt,cnt,ob)

type(infompi)              , intent(in)    :: info
integer(i4)  , dimension(:), intent(in)    :: srt
integer(i4)  , dimension(:), intent(in)    :: cnt
integer(i1)  , dimension(:), intent(inout) :: ob

integer :: len
integer :: ilen
integer :: alen
integer :: n
integer :: m

ilen = sizeof(info)

alen = sizeof(srt)

len = sizeof(ob)

n = ilen

ob(1:n) = transfer(info, ob(1:n))

m = n + alen

ob((n+1):m) = transfer(srt, ob((n+1):m))

n = m + alen

ob((m+1):n) = transfer(cnt, ob((m+1):n))

end subroutine infotobyte

!-------------------------------------------------------

subroutine bytetoinfo(info,srt,cnt,ob)

type(infompi)              , intent(inout)    :: info
integer(i4)  , dimension(:), intent(inout)    :: srt
integer(i4)  , dimension(:), intent(inout)    :: cnt
integer(i1)  , dimension(:), intent(in)       :: ob

integer :: len
integer :: ilen
integer :: alen
integer :: n
integer :: m

ilen = sizeof(info)

alen = sizeof(srt)

len = sizeof(ob)

n = ilen

info = transfer(ob(1:n), info)

m = n + alen

srt = transfer(ob((n+1):m), srt)

n = m + alen

cnt = transfer(ob((m+1):n), cnt)

end subroutine bytetoinfo

!-------------------------------------------------------

end program mpitest
