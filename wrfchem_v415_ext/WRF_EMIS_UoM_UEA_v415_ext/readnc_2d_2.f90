subroutine handle_err(status)
integer:: status
  return
end subroutine handle_err

subroutine readnc_2d_2(   &
     twoDfield     &  !O A three dimensional field
     ,fieldname      &  !I The name of the field
     ,longi           &
     ,lati            &
     ,times           &
     ,filename       &  !I The name of the file to read from
     ,IM             &  !I Number of longitudes
     ,JM             &  !I Number of latitudes
     ,ntimestep      &  !I The timestep to read
     )

  !Purpose: Open a netCDF file and get 3D field from the file
  !The field you want to get has the name "fieldname", and the file
  !has the name "filename". 
  !The timestep you want to read is the input parameter ntimestep

  !Remember to include -I${NETCDF_INC} -I${NETCDF_LIB} for compiling
  !And for linking -L${NETCDF_LIB} and -lnetcdf in your makefile
  !The values of these two should be set in your .bashrc file like this: 
  !export NETCDF_INC = /mn/hox/u1/jsundet/include/  
  !export NETCDFid(_LIB = /mn/hox/u1/jsundet/lib/
  !Or if you are using C-shell put the following in your .cshrc file
  !setenv NETCDF_INC /mn/hox/u1/jsundet/include/
  !setenv NETCDF_LIB /mn/hox/u1/jsundet/lib/
  !If this doesn't work, the explicitly link $NETCDF_INC/netcdf.mod to your run-dir

  !Author: Alf Grini, alf.grini@geofysikk.uio.no

  use netcdf

  implicit none
  
  !INPUT
  character*(*), intent(in)           :: fieldname  !I Name of field
  character*(*), intent(in)           :: filename   !I Name of netCDFfile
  character*(*), intent(in)           :: longi  !I Name of field
  character*(*), intent(in)           :: lati   !I Name of netCDFfile
  character*(*), intent(in)           :: times   !I Name of netCDFfile
  integer,intent(in)                  :: IM         !I Number of longitudes
  integer,intent(in)                  :: JM         !I Number of latitudes
  integer,intent(in)                  :: ntimestep  !I The timestep to be read

  !OUTPUT
  real*8, intent(out)        :: twoDfield(IM,JM)  !Three dimensional field

  !LOCAL
  !LOCAL NETCDF DIMENSION IDs ETCETERA
  integer                  :: lon_dim_id      !Id for longitude dimension
  integer                  :: lon_id          !Id for variable longitude
  real                     :: lon(IM)         !variable lon (in file)
  integer                  :: lat_dim_id      !Id for latitude dimension
  integer                  :: lat_id          !Id for latitude
  real                     :: lat(JM)         !Variable for latitude
  integer                  :: time_dim_id     !Id for time dimension
  integer                  :: time_id         !Id for time
  integer                  :: field_dim_id(3) !Dimension id for field
  integer                  :: field_id        !Variable id for field
  integer                  :: srt_lon_lat_time(3) !Start point 
  integer                  :: cnt_lon_lat_time(3) !Count indexes
  integer                  :: nlons                   !Longitudes in file
  integer                  :: nlats                   !Latitudes in file
  integer                  :: nsteps                  !Timesteps avaiable in file
  integer                  :: status                  !status of process (0=OK)
  integer                  :: ncid                    !file id 

  !Array which tells you where to start picking your 3D field
  srt_lon_lat_time= (/ 1 , 1 , ntimestep /)    !Start array
  !Array which tells you how far to count when picking it
  cnt_lon_lat_time= (/ IM , JM , 1 /)         !Count array

  !**********START CODE************************''

  status=nf90_noerr  !Status is 0 and should be kept that way !!




  !Open the existing file
  status=nf90_open(filename, nf90_nowrite, ncid)
  if(status/=nf90_noerr)call handle_err(status)

  !Inquire dimension ids

  status = nf90_inq_dimid(ncid,lati,lat_dim_id)
  if(status/=nf90_noerr)call handle_err(status)

  status = nf90_inq_dimid(ncid,longi,lon_dim_id)
  if(status/=nf90_noerr)call handle_err(status)
  status = nf90_inq_dimid(ncid,times,time_dim_id)
  if(status/=nf90_noerr)call handle_err(status)

  !Dimension id for 3D field /lon/lat/lev/time

  field_dim_id(1)=lon_dim_id
  field_dim_id(2)=lat_dim_id
  field_dim_id(3)=time_dim_id


  !Inquire dimensions
  status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
  if(status/=nf90_noerr)call handle_err(status)
  if(nlats/=JM)then
     write(6,*)'file'//filename//' reports JM = ',nlats
     write(6,*)'your array has dimension ',JM
     stop
  endif



  status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
  if(status/=nf90_noerr)call handle_err(status)
  if(nlons/=IM)then
     write(6,*)'file'//filename//'file reports IM = ',nlons
     write(6,*)'your array has dimension',IM
     stop
  endif
 
  status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
  if(status/=nf90_noerr)call handle_err(status)
  if(ntimestep.gt.nsteps.or.nsteps.le.0)then
     write(6,*)'file'//filename//'file reports nsteps = ',nsteps
     write(6,*)'you try to read timestep',ntimestep
     stop
  endif

  
  !Get variable ID
  status=nf90_inq_varid(ncid,fieldname,field_id)
  if(status/=nf90_noerr)call handle_err(status)

  !Finally after all this, you can get the variable you want !!
  !and put it in the threeDfield array
  status=nf90_get_var(ncid,field_id,twoDfield &
       ,start=srt_lon_lat_time   &
       ,count=cnt_lon_lat_time )
  if(status/=nf90_noerr)call handle_err(status)

  write(6,*)'got variable ',fieldname

  !Closing file
  status=nf90_close(ncid)
  if(status/=nf90_noerr)call handle_err(status)
  
  return
  end subroutine readnc_2d_2
