module param_fftw
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03' 
  include 'mpif.h'

! -----------------------------------------------------------------
! FFTW 3 Plans
  type(C_PTR), save :: plan3_ur_uk,plan3_vr_vk,plan3_wr_wk
  type(C_PTR), save :: plan3_zxnr_zxnk,plan3_zynr_zynk,plan3_zznr_zznk,plan3_ttnr_ttnk  
  type(C_PTR), save :: plan3_nzxr_nzxk,plan3_nzyr_nzyk,plan3_nzzr_nzzk

  type(C_PTR), save :: plan3_uk_ur,plan3_vk_vr,plan3_wk_wr
  type(C_PTR), save :: plan3_zxnk_zxnr,plan3_zynk_zynr,plan3_zznk_zznr,plan3_ttnk_ttnr  
! -----------------------------------------------------------------

end module param_fftw
