module averagemod

use parametersmod, only : i2,i4,sp,dp
use outputmod,     only : getoutfile,handle_err
use netcdf
use mpi

implicit none

!--------------------

type infompi
  character(100) :: infile
  character(100) :: outfile
  character(100) :: timestring
  integer(i4)    :: ofid
  integer(i2)    :: nproc
  integer(i4)    :: t0
  integer(i4)    :: nt
end type

!--------------------

contains

!-------------------------------------------------------

subroutine average(info,job)

type(infompi), intent(in) :: info
integer(i4)  , dimension(2), intent(in) :: job

integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt

character(100) :: infile
character(100) :: outfile

integer :: status
integer :: ifid
integer :: ofid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen
integer :: tlen

integer :: rank
integer :: ierr

integer(i2) :: missing_value
real(sp) :: scale_factor
real(sp) :: add_offset

real(sp), allocatable, dimension(:,:) :: var_in
integer(i4), allocatable, dimension(:) :: var_out

integer :: i

! --------------------

infile = info%infile
outfile = info%outfile
ofid = info%ofid

srt = [job(1), info%t0]

cnt = [job(2), info%nt]

allocate(var_in(cnt(1),cnt(2)))

! --------------------

status = nf90_open(infile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,var_in,start=srt,count=cnt)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ifid,varid,"missing_value",missing_value)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ifid,varid,"scale_factor",scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ifid,varid,"add_offset",add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

where (var_in /= missing_value)
  var_in = real(var_in) * scale_factor + add_offset
elsewhere
  var_in = -9999.
end where

! --------------------

allocate(var_out(size(var_in,1)))

do i = 1, size(var_out)

  var_out(i) = nint(sum(var_in(i,:)) / real(size(var_in,2)) * 10)

end do

! --------------------

! call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
!
! if(rank == 0) then
!
!   call getoutfile(outfile, ofid)
status = nf90_open_par(outfile,nf90_write,MPI_COMM_WORLD,MPI_INFO_NULL,ofid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'ave_tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,var_out)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

! end if

! --------------------



end subroutine average



end module averagemod
