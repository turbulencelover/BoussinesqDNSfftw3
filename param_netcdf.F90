module param_netcdf
  implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Netcdf i/o Parameters          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -----------------------------------------------------------------
! Netcdf Stuff for Real Space Output

  integer, parameter :: iputu=1,iputv=1,iputw=1,iputth=1,iputzx=1,iputzy=1,iputzz=1  

! NETCDF IDs: every real space variable at a timestep has its own file 
! indnetcdf is a (7 x 6 x nrsp) matrix
! for storing values (U,V,W,TH,ZX,ZY,ZZ) x (ncid, varid, xid, yid, zid,vartid) x (nrsp time points) 
  integer, dimension(:,:,:), save, allocatable :: idnetcdf

! start and count arrays
  integer, save      :: ncstart(3),nccount(3)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Netcdf stuff for slice outputs
  integer, parameter :: islputxu=1,islputxv=1,islputxw=1,islputxth=1,islputxzx=1,islputxzy=1,islputxzz=1

  integer, parameter :: islputyu=1,islputyv=1,islputyw=1,islputyth=1,islputyzx=1,islputyzy=1,islputyzz=1

  integer, parameter :: islputzu=1,islputzv=1,islputzw=1,islputzth=1,islputzzx=1,islputzzy=1,islputzzz=1

! idncsl is a (7 x 6 x len(xslice)) matrix
! for storing values (U,V,W,TH,ZX,ZY,ZZ) x (ncid,varid,yid,zid,tid,vartid) x (yslice)
  integer, dimension(:,:,:), save, allocatable :: idncslx,idncsly,idncslz

! start and count arrays
  integer, save      :: ncslstart(3),ncslcount(3)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Netcdf Stuff for Restart File:  
! every variable at a timestep has its own file
! idncrst is a (4 x 7 x nrst) matrix
! for storing values (ZXK,ZYK,ZZK,ZZ) x (ncid,varid,kxid,kyid,kzid,kri,vartid) x (nrst time points)
  integer, dimension(:,:,:), save, allocatable :: idncrst

! start and count arrays
  integer, save :: ncstartk(4),nccountk(4)
! -----------------------------------------------------------------
end module param_netcdf
