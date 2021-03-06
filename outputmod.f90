module outputmod

use parametersmod, only : i2,i4,sp,dp
use netcdf
use mpi

implicit none

real(sp), parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

contains

!-------------------------------------------------------

subroutine getoutfile(outfile,validcell)

character(*), intent(in) :: outfile
integer(i4) , intent(in) :: validcell

integer :: ofid
integer :: status
integer :: varid
integer :: dimid

integer :: i

character(8) :: today
character(10) :: now

write(0,*) 'Creating outfile'

status = nf90_create(outfile,nf90_netcdf4,ofid)
if (status/=nf90_noerr) call handle_err(status)

! status = nf90_create_par(outfile,ior(nf90_netcdf4,nf90_mpiio),MPI_COMM_WORLD,MPI_INFO_NULL,ofid)
! if (status/=nf90_noerr) call handle_err(status)

write(0,*) 'Create outfile: success'

status = nf90_put_att(ofid,nf90_global,'title','paralleltest output file')
if (status/=nf90_noerr) call handle_err(status)

call date_and_time(today,now)

status = nf90_put_att(ofid,nf90_global,'created',today//' '//now(1:4))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (status/=nf90_noerr) call handle_err(status)


!----
! index

status = nf90_def_dim(ofid,'index',validcell,dimid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ofid,'index',nf90_int,dimid,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,(/(i,i=1,validcell,1)/))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'long_name','index of lon and lat')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'units','1 to length')
if (status/=nf90_noerr) call handle_err(status)


!----
!average temperature

status = nf90_def_var(ofid,'ave_tmp',nf90_short,dimid,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'long_name','test average temperature')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'units','degC')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (status/=nf90_noerr) call handle_err(status)

!----

status = nf90_enddef(ofid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Finished create outfile'



end subroutine getoutfile


!-------------------------------------------------------

subroutine handle_err(status)

!   Internal subroutine - checks error status after each netcdf call,
!   prints out text message each time an error code is returned.

integer, intent (in) :: status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

!-------------------------------------------------------

end module outputmod
