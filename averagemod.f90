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
  integer(i4)    :: nproc
  integer(i4)    :: validcell
  integer(i4)    :: t0
  integer(i4)    :: nt
end type

!--------------------

contains

!-------------------------------------------------------

subroutine average(info,job,rank)

type(infompi), intent(in) :: info
integer(i4)  , dimension(2), intent(in) :: job


integer(i4), dimension(2) :: srt
integer(i4), dimension(2) :: cnt

character(100) :: infile
character(100) :: outfile

integer :: status
integer :: ifid
integer :: ofid
integer :: varid

integer :: rank
integer :: ierr

integer(i4) :: mpi_info

integer(i2) :: missing_value
real(sp) :: scale_factor
real(sp) :: add_offset

real(sp)   , allocatable, dimension(:,:) :: var_in
integer(i2), allocatable, dimension(:)   :: var_out

integer :: i

! --------------------

infile = info%infile
outfile = info%outfile

srt = [job(1), info%t0]

cnt = [job(2), info%nt]

allocate(var_in(cnt(1),cnt(2)))

! --------------------
! Read temperature variable from infile

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
! Calculate average temperature of each validcell over entire time period

allocate(var_out(size(var_in,1)))

do i = 1, size(var_out)

  var_out(i) = nint(sum(var_in(i,:)) / real(size(var_in,2)) * 10)

end do

write(0,*) 'Rank', rank, 'var_in min/max check:',minval(var_in), maxval(var_in)

! --------------------
! Write average temperature (var_out) into the outfile

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

status = nf90_open(outfile,nf90_write,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'ave_tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,var_out,start=srt,count=cnt)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Data written successfully by rank', rank

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

! --------------------



end subroutine average



end module averagemod
