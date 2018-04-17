! Triply Periodic Boussinesq
! Original Code by Peter Bartello
! Rewritten by Michael Waite
! f90, mpi: Winter 2008-2009
! Cleaned up Winter 2012
! Summer 2012: Upgraded netcdf calls to nf90; upgraded to module param
! Modified by Karthik Velakur, 2013-2014: updated to FFTW3, output slices
! To compile: have to link to FFTW and NETCDF libraries


!     NOTES: 
!     kx,ky,kz are wavenumbers. 
!     ikx,y,z  are corresponding array indices.
!     KTX,Y,Z  are truncation wavenumbers.  
!     IKTX,Y,Z are corresponding indices.

!     Code solves Boussinesq equations in vorticity form.
!     zx,y,z is vorticity. u,v,w is velocity.  t is potential temp.
!     nzx,y,z is nonlinear term in vorticity equation.

!     zxnk,zynk,zznk,ttnk are fields at current time level (n)
!     zxok,zyok,zzok,ttok are fields at time level n-1. **
!     zx1k,zy1k,zz1k,tt1k are fields at time level n-2. ** (AB3 only)
!     ** NOTE: for AB3, store right-hand-side, not basic fields, at times n-1 and n-2.
!     FFTW v 3.3.3 Manual can be found at http://www.fftw.org


PROGRAM MAIN
  use param
  use param_fftw 
  implicit none 
  

! -----------------------------------------------------------------
! Additional Parameters
! -----------------------------------------------------------------

  ! DAMPING
  real, parameter    :: ek      =   1./(1000.)*delt      ! kh=0 damping time scale

  ! INITIAL CONDITIONS
  integer, parameter :: irest  =   0                     ! flag for starting from restart file
  integer, parameter :: myrest =   irest                 ! used in random seed
  real,    parameter :: ki     =   3.                    ! wavenumber for initial conditions 
  real,    parameter :: etot   =   1.e-5                 ! initial energy                    
  real,    parameter :: pe     =   etot/3.               ! initial potential energy
  real,    parameter :: ke     =   etot - pe             ! initial kinetic energy
  integer, parameter :: linear  =   0                    ! flag for omitting nonlinear terms

  ! FORCING
  real,    parameter :: kf     =   4.                    ! wavenumber for forcing
  real,    parameter :: tau    =   10.*delt              ! time scale for forcing
  real,    parameter :: alpha  = exp(-delt/tau)          ! memory parameter
  real,    parameter :: eps    =   6.e-9                 ! approx target dissipation rate for forcing amplitude
  real,    parameter :: ampv   =   sqrt(eps/tau)         ! random forcing amplitude
  real,    parameter :: ampw   =   0.                    ! random wave forcing amplitude (when used)
  integer, parameter :: nfmax  = 1200                    ! max number of forced modes

  ! MISC
  real, parameter    :: d2=2*delt,d12=delt/12     
  real, parameter    :: k2h=visch*delt,k2z=viscz*delt
  real, parameter    :: v2h=visch*delt,v2z=viscz*delt


! -----------------------------------------------------------------
! Declaration of variables
! -----------------------------------------------------------------

  complex, dimension(:,:,:), pointer :: zxok,zyok,zzok,ttok,zxnk,zynk,zznk,ttnk
  real,    dimension(:,:,:), pointer :: zxor,zyor,zzor,ttor,zxnr,zynr,zznr,ttnr
  type(C_PTR) :: zxo_p,zyo_p,zzo_p,tto_p,zxn_p,zyn_p,zzn_p,ttn_p

  complex, dimension(:,:,:), pointer :: nzxk,nzyk,nzzk,nttk,uk,vk,wk
  real,    dimension(:,:,:), pointer :: nzxr,nzyr,nzzr,nttr,ur,vr,wr
  type(C_PTR) :: nzx_p,nzy_p,nzz_p,ntt_p,u_p,v_p,w_p  

  complex, dimension(:,:,:), allocatable :: rhzx,rhzy,rhzz,rhtt,geok,gw1k,gw2k
  complex, dimension(:,:,:), allocatable :: zx1k,zy1k,zz1k,tt1k

  complex :: termzx,termzy,termzz,termtt,tzx,tzy,tzz,ttt
  complex :: u,v,w,c1,c2,c3,gtau(nfmax,4)
  real :: dmz,dpz,r1,r2,kx,ky,kz,kh,khn,wk2,k
  real :: ranno,fseed,ts=0,etime
  real*4 :: time1(2),time2(2),time3
      
  integer :: ikx,iky,ikz,ikza
  integer :: iseed,nt=0,nt0,ntdump

  external :: init,out,convol,spec,constr,velo,transf,ranno,forset,force

! -----------------------------------------------------------------
! Initialize MPI
! -----------------------------------------------------------------

  time3 = etime(time1)

  call mpi_init(istatus)
  call mpi_comm_size(mpi_comm_world,npe,istatus)
  call mpi_comm_rank(mpi_comm_world,mype,istatus)
  call fftwf_mpi_init()
  
  ! checks on the number of processors
  if (mype.eq.0) print*,'npe = ',npe
  if(mod(n2,npe).ne.0 .or. mod(iktz,npe).ne.0) then
    print*,'npe must divide n2 and iktz'
    stop
  endif

  ! array sizes based on npe
  iktzp=iktz/npe
  n2dp = n2d/npe; n2p = n2/npe;

! -----------------------------------------------------------------
! Allocate arrays not used in FFTW
! -----------------------------------------------------------------

  allocate(L(iktx,ikty,iktzp))
  allocate(kxa(iktx)) 
  allocate(kya(ikty)) 
  allocate(kza(iktz))
  allocate(rhzx(iktx,ikty,iktzp))
  allocate(rhzy(iktx,ikty,iktzp))
  allocate(rhzz(iktx,ikty,iktzp))
  allocate(rhtt(iktx,ikty,iktzp))
  allocate(geok(iktx,ikty,iktzp))
  allocate(gw1k(iktx,ikty,iktzp))
  allocate(gw2k(iktx,ikty,iktzp))
  allocate(zx1k(iktx,ikty,iktzp))
  allocate(zy1k(iktx,ikty,iktzp))
  allocate(zz1k(iktx,ikty,iktzp))
  allocate(tt1k(iktx,ikty,iktzp))

! -----------------------------------------------------------------
! Allocate arrays used in FFTW
! -----------------------------------------------------------------

! Get local data size (note dimension reversal), and allocate. 
! Use in-place transforms.
! In Fortran order we see real space to Fourier space dimensions in order (kx,ky,kz) <-> (x,z,y)
! For MPI, in C order we have to pass dimensions in reverse order (kz,ky,kx) <-> (y,z,x). 
! locz is local block size of z, and loczstart is local starting value of z. 
! In our case, since npe|n3, we will have locz=iktzp over all processors.
! Also rows are always divided in kz by rank order. So proc of rank 0 will have first
! block, proc 1 will have second, and so on. 
! cf. Sec 6.4, 6.5, 6.13 and 7.2 of FFTW manual

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zxo_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zxo_p, zxok, [iktx,ikty,iktzp])
  call c_f_pointer(zxo_p, zxor, [n1d,n3d,n2p])

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zyo_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zyo_p, zyok, [iktx,ikty,iktzp])
  call c_f_pointer(zyo_p, zyor, [n1d,n3d,n2p])

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zzo_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zzo_p, zzok, [iktx,ikty,iktzp])
  call c_f_pointer(zzo_p, zzor, [n1d,n3d,n2p])

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  tto_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(tto_p, ttok, [iktx,ikty,iktzp])
  call c_f_pointer(tto_p, ttor, [n1d,n3d,n2p])

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zxn_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zxn_p, zxnk, [iktx,ikty,iktzp])
  call c_f_pointer(zxn_p, zxnr, [n1d,n3d,n2p])

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zyn_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zyn_p, zynk, [iktx,ikty,iktzp])
  call c_f_pointer(zyn_p, zynr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  zzn_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(zzn_p, zznk, [iktx,ikty,iktzp])
  call c_f_pointer(zzn_p, zznr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  ttn_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(ttn_p, ttnk, [iktx,ikty,iktzp])
  call c_f_pointer(ttn_p, ttnr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  nzx_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(nzx_p, nzxk, [iktx,ikty,iktzp])
  call c_f_pointer(nzx_p, nzxr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  nzy_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(nzy_p, nzyk, [iktx,ikty,iktzp])
  call c_f_pointer(nzy_p, nzyr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  nzz_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(nzz_p, nzzk, [iktx,ikty,iktzp])
  call c_f_pointer(nzz_p, nzzr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  ntt_p       = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(ntt_p, nttk, [iktx,ikty,iktzp])
  call c_f_pointer(ntt_p, nttr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  u_p         = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(u_p, uk, [iktx,ikty,iktzp])
  call c_f_pointer(u_p, ur, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  v_p         = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(v_p, vk, [iktx,ikty,iktzp])
  call c_f_pointer(v_p, vr, [n1d,n3d,n2p])  

  alloc_local = fftwf_mpi_local_size_many(3,(/iktz,ikty,iktx/),1_8,iktzp,mpi_comm_world,locz,loczstart)
  w_p   = fftwf_alloc_complex(alloc_local)
  call c_f_pointer(w_p, wk, [iktx,ikty,iktzp])
  call c_f_pointer(w_p, wr, [n1d,n3d,n2p])

! -----------------------------------------------------------------
! INITIALIZE FFTW. Dimesension are passed in reverse to MPI C routines. 
! -----------------------------------------------------------------

! FFTW_MPI_TRANSPOSED_OUT flag to compute forward  DFT: (y,z,x) -> (kz,ky,kx). Fortran sees this as (kx,ky,kz).
! FFTW_MPI_TRANSPOSED_IN  flag to compute backward DFT: (kz,ky,kx), viewed as (ky,kz,kx) -> (y,z,x). Fortran sees this as (x,z,y).

  if (mype.eq.0) print*,'Initializing FFTW'
  plan3_ur_uk     = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,ur,  uk,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_vr_vk     = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,vr,  vk,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_wr_wk     = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,wr,  wk,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))  
  plan3_zxnr_zxnk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,zxnr,zxnk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_zynr_zynk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,zynr,zynk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_zznr_zznk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,zznr,zznk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_ttnr_ttnk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,ttnr,ttnk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_nzxr_nzxk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,nzxr,nzxk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_nzyr_nzyk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,nzyr,nzyk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
  plan3_nzzr_nzzk = fftwf_mpi_plan_dft_r2c_3d(n2,n3,n1,nzzr,nzzk,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))

  plan3_uk_ur     = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,uk,  ur,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))
  plan3_vk_vr     = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,vk,  vr,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))
  plan3_wk_wr     = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,wk,  wr,  mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))  
  plan3_zxnk_zxnr = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,zxnk,zxnr,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))
  plan3_zynk_zynr = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,zynk,zynr,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))
  plan3_zznk_zznr = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,zznk,zznr,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))
  plan3_ttnk_ttnr = fftwf_mpi_plan_dft_c2r_3d(n2,n3,n1,ttnk,ttnr,mpi_comm_world,IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN))

! INITIALIZE I/O  
  call io_prep()

! Open files for diagnostics
  if (mype.eq.0) then
    open (51,file='spcz.dat', form='formatted')
    open (52,file='spch.dat', form='formatted')
    open (53,file='spc.dat',  form='formatted')
    open (46,file='eng.dat',  form='formatted')
    open (48,file='trnh.dat', form='formatted')
    open (49,file='trnz.dat', form='formatted')
    open (50,file='trn.dat',  form='formatted')
    open (79,file='eps.dat',  form='formatted')
  endif
      
! Set random seed
  fseed = ranno(1946 + myrest + mype)
  fseed = ranno(0)/6.29 + 0.5
  iseed = int(10000.*fseed) + 1
  fseed = ranno(iseed)

! -----------------------------------------------------------------
! Initialize fields, L, etc.
! -----------------------------------------------------------------

  call init(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,zxor,zyor,zzor,ur,vr,wr,ttor,ke,pe,irest,ts,ki)

  ! initialize a few last things
  time = ts + nt*delt
  gtau = cmplx(0.,0.)
  zx1k=cmplx(0.,0.)
  zy1k=cmplx(0.,0.)
  zz1k=cmplx(0.,0.)
  tt1k=cmplx(0.,0.)


! -----------------------------------------------------------------
! Print out all Parameters
! -----------------------------------------------------------------
  if (mype.eq.0) then
    print*,'                '
    print*,'                '
    print*,'Parameters: --------------------------------'
    print*,'              N1,N2,N3 =  ', n1,n2,n3
    print*,'              L1,L2,L3 =  ', L1,L2,L3
    print*,'              dx,dy,dz =  ', L1/n1,L2/n2,L3/n3
    print*,'              IKTX,Y,Z =  ', iktx,ikty,iktz
    print*,'                     KT = ', ktrunc_x,ktrunc_y,ktrunc_z
    print*,'Order of Laplacian diss.=  ', ilap 
    print*,'                  VISCH = ', visch
    print*,'                  VISCZ = ', viscz
    print*,'              tau_VISCH = ', (visch*ktrunc_x**ilap2)**(-1)
    print*,'              tau_VISCZ = ', (viscz*ktrunc_z**ilap2)**(-1)
    if (linear.eq.1) print*,'Nonlinear terms switched off.' 
    print*,'               Timestep = ', delt
    print*,'     Integration length = ', nstop*delt,' = ', nstop,' DT.'
    print*,'       Output frequency = ', nout*delt, ' = ', nout,' dt.'
    print*,'  Rspace dump frequency = ', nbig*delt, ' = ', nbig,' dt.'
    print*,'  '
    if (irest.eq.0) then
      print*,' Starting from ICs '
      print*,' Initial kinetic energy = ', ke
      print*,'   "    potential  "    = ', pe
      print*,'   "      total    "    = ', ke + pe
    else
      print*,'  Starting from restart file'
      print*,'  Restart from rec.     = ', irest
      print*,'  Restart from time     = ', ts
    endif
    print*,'    '
    print*,'    Thermal expansivity = ', aj
    print*,'    Vertical T gradient = ', bj
    print*,'    Brunt-Vaisala freq. = ', sqrt(bf2)
    print*,'     Coriolis parameter = ', cor
    print*,'   '
    if (max(ampv,ampw).ne.0.) then
      print*,' Random forcing: ' 
      print*,'     Amplitude (vortical) = ', ampv
      print*,'     Amplitude (wave)     = ', ampw
      print*,'               Wavenumber = ', kf
      print*,'               Memory = ', tau,' = ', int(tau/delt),' DT.'
    else
      print*,'      No forcing.' 
    endif
    print*,'                '
  endif

! -----------------------------------------------------------------
! Diagnostics on ICs
! -----------------------------------------------------------------


  if (irest.ge.0) then
    call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk)
    call proj(zxok,zyok,zzok)  
    call out(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k) 
    call spec(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,2,52)
    call spec(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,3,51)
    call spec(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,1,53)
    call cfl(zxok,zyok,zzok,uk,vk,wk,ur,vr,wr)
    if (rsflag.eq.1) call dumpreal(zxok,zyok,zzok,ttok,zxor,zyor,zzor,ttor,uk,vk,wk,ur,vr,wr,1,FULLFIELD_MODE)
    if (slxflag.eq.1 .or. slyflag.eq.1 .or. slzflag.eq.1) call dumpreal(zxok,zyok,zzok,ttok,zxor,zyor,zzor,ttor,uk,vk,wk,ur,vr,wr,1,SLICE_MODE)
  endif ! irest

  time3 = etime(time2)
  time3 = time2(1) - time1(1)
  if (mype.eq.0) then
    print*,' Initialization time:   '
    write(6,5001) time3, time3/60., time3/3600., time3/86400.
    print*,'    '
  endif


! -----------------------------------------------------------------
! Beginning of the first timestep.  Use explicit trapezoidal.
! -----------------------------------------------------------------
  nt   = 1
  time = ts + nt*delt

  if (linear.ne.1) then
    call constr(zxok,zyok,zzok,ttok,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr,zxor,zyor,zzor,ttor,nzxr,nzyr,nzzr,nttr)
  else
    nzxk = cmplx(0.,0.)
    nzyk = cmplx(0.,0.)
    nzzk = cmplx(0.,0.)
    nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)
  endif

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        if ( L(ikx,iky,ikz).ne.1 ) then
          zxnk(ikx,iky,ikz) = cmplx(0.,0.)
          zynk(ikx,iky,ikz) = cmplx(0.,0.)
          zznk(ikx,iky,ikz) = cmplx(0.,0.)
          ttnk(ikx,iky,ikz) = cmplx(0.,0.)
        else
          kx  = kxa(ikx)
          kh  = sqrt( kx*kx + ky*ky )
          khn = kh * L1/twopi
          wk2 = kx*kx + ky*ky + kz*kz
          k  = sqrt(wk2)

          c1 = +  ky * zzok(ikx,iky,ikz) -  kz * zyok(ikx,iky,ikz)
          c2 = +  kz * zxok(ikx,iky,ikz) -  kx * zzok(ikx,iky,ikz)
          c3 = +  kx * zyok(ikx,iky,ikz) -  ky * zxok(ikx,iky,ikz)
          u = zi * c1 / wk2
          v = zi * c2 / wk2
          w = zi * c3 / wk2

          ! spherical dissipation
          r1 = v2h/delt*wk2**ilap
          r2 = k2h/delt*wk2**ilap

          ! cylindrical dissipation
!          r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
!          r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2


          if (khn.lt.0.1) then
            r1 = r1 + ek/delt
            r2 = r2 + ek/delt
          endif

          termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttok(ikx,iky,ikz) + cor*zi*kz*u - zxok(ikx,iky,ikz)*r1
          termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttok(ikx,iky,ikz) + cor*zi*kz*v - zyok(ikx,iky,ikz)*r1
          termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zzok(ikx,iky,ikz)*r1
          termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

          rhzx(ikx,iky,ikz) = termzx
          zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) + delt*termzx
          rhzy(ikx,iky,ikz) = termzy
          zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) + delt*termzy
          rhzz(ikx,iky,ikz) = termzz
          zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) + delt*termzz
          rhtt(ikx,iky,ikz) = termtt
          ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*termtt

        endif ! L
      enddo
    enddo
  enddo

  if (linear.ne.1) then
    call constr(zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr)
  else
    nzxk = cmplx(0.,0.)
    nzyk = cmplx(0.,0.)
    nzzk = cmplx(0.,0.)
    nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)
  endif

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        if ( L(ikx,iky,ikz).ne.1 ) then
          zxnk(ikx,iky,ikz) = cmplx(0.,0.)
          zynk(ikx,iky,ikz) = cmplx(0.,0.)
          zznk(ikx,iky,ikz) = cmplx(0.,0.)
          ttnk(ikx,iky,ikz) = cmplx(0.,0.)
        else
          kx  = kxa(ikx)             
          kh  = sqrt( kx*kx + ky*ky )
          khn = kh * L1/twopi
          wk2 = kx*kx + ky*ky + kz*kz
          k   = sqrt( wk2 )

          c1 = +  ky*zznk(ikx,iky,ikz) -  kz*zynk(ikx,iky,ikz)
          c2 = +  kz*zxnk(ikx,iky,ikz) -  kx*zznk(ikx,iky,ikz)
          c3 = +  kx*zynk(ikx,iky,ikz) -  ky*zxnk(ikx,iky,ikz)
          u  = zi * c1 / wk2
          v  = zi * c2 / wk2
          w  = zi * c3 / wk2

          ! spherical dissipation
          r1 = v2h/delt*wk2**ilap
          r2 = k2h/delt*wk2**ilap

          ! cylindrical dissipation
!          r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
!          r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2


          if (khn.lt.0.1) then
            r1 = r1 + ek/delt
            r2 = r2 + ek/delt
          endif

          termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttnk(ikx,iky,ikz) + cor*zi*kz*u - zxnk(ikx,iky,ikz)*r1
          termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttnk(ikx,iky,ikz) + cor*zi*kz*v - zynk(ikx,iky,ikz)*r1
          termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zznk(ikx,iky,ikz)*r1
          termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttnk(ikx,iky,ikz)*r2

          zx1k(ikx,iky,ikz) = rhzx(ikx,iky,ikz) + zxok(IKX,IKY,IKZ)*r1
          zy1k(ikx,iky,ikz) = rhzy(ikx,iky,ikz) + zyok(IKX,IKY,IKZ)*r1
          zz1k(ikx,iky,ikz) = rhzz(ikx,iky,ikz) + zzok(IKX,IKY,IKZ)*r1
          tt1k(ikx,iky,ikz) = rhtt(ikx,iky,ikz) + ttok(IKX,IKY,IKZ)*r2

          zxok(ikx,iky,ikz) = zxok(ikx,iky,ikz)  +   delt*(termzx + rhzx(ikx,iky,ikz))/2.
          zyok(ikx,iky,ikz) = zyok(ikx,iky,ikz)  +   delt*(termzy + rhzy(ikx,iky,ikz))/2.
          zzok(ikx,iky,ikz) = zzok(ikx,iky,ikz)  +   delt*(termzz + rhzz(ikx,iky,ikz))/2.
          ttok(ikx,iky,ikz) = ttok(ikx,iky,ikz)  +   delt*(termtt + rhtt(ikx,iky,ikz))/2.

        endif ! L
      enddo
    enddo
  enddo


! -----------------------------------------------------------------
! Beginning of the second timestep (only necessary for AB3).
! -----------------------------------------------------------------
  nt   = 2
  time = ts + nt*delt

  if (LINEAR.ne.1) then
    call constr (zxok,zyok,zzok,ttok,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr,zxor,zyor,zzor,ttor,nzxr,nzyr,nzzr,nttr)
  else
    nzxk = cmplx(0.,0.)
    nzyk = cmplx(0.,0.)
    nzzk = cmplx(0.,0.)
    nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)
  endif

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx

        if ( L(ikx,iky,ikz).ne.1 ) then
          zxnk(ikx,iky,ikz) = cmplx(0.,0.)
          zynk(ikx,iky,ikz) = cmplx(0.,0.)
          zznk(ikx,iky,ikz) = cmplx(0.,0.)
          ttnk(ikx,iky,ikz) = cmplx(0.,0.)
        else
          kx  = kxa(ikx)
          kh  = sqrt( kx*kx + ky*ky )
          khn = kh * L1/twopi
          wk2 = kx*kx + ky*ky + kz*kz
          k  = sqrt(wk2)

          c1 = +  ky * zzok(ikx,iky,ikz) -  kz * zyok(ikx,iky,ikz)
          c2 = +  kz * zxok(ikx,iky,ikz) -  kx * zzok(ikx,iky,ikz)
          c3 = +  kx * zyok(ikx,iky,ikz) -  ky * zxok(ikx,iky,ikz)
          u = ZI * c1 / wk2
          v = ZI * c2 / wk2
          w = ZI * c3 / wk2

          ! spherical dissipation
          r1 = V2H/delt*wk2**ILAP
          r2 = K2H/delt*wk2**ILAP

          ! cylindrical dissipation
!          r1 = V2H/delt*kh**ILAP2 + V2Z/delt*kz**ILAP2
!          r2 = K2H/delt*kh**ILAP2 + K2Z/delt*kz**ILAP2


          if (khn.lt.0.1) then
            r1 = r1 + ek/delt
            r2 = r2 + ek/delt
          endif

          termzx = nzxk(ikx,iky,ikz) + aj*ZI*ky*ttok(ikx,iky,ikz) + cor*ZI*kz*u - zxok(IKX,IKY,IKZ)*r1
          termzy = nzyk(ikx,iky,ikz) - aj*ZI*kx*ttok(ikx,iky,ikz) + cor*ZI*kz*v - zyok(IKX,IKY,IKZ)*r1
          termzz = nzzk(ikx,iky,ikz)                              + cor*ZI*kz*w - zzok(ikx,iky,ikz)*r1
          termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

          rhzx(ikx,iky,ikz) = termzx
          zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) + delt*termzx
          rhzy(ikx,iky,ikz) = termzy
          zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) + delt*termzy
          rhzz(ikx,iky,ikz) = termzz
          zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) + delt*termzz
          rhtt(ikx,iky,ikz) = termtt
          ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*termtt
        endif ! L
      enddo
    enddo
  enddo

  if (linear.ne.1) then
    call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr)
  else
    nzxk = cmplx(0.,0.)
    nzyk = cmplx(0.,0.)
    nzzk = cmplx(0.,0.)
    nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)
  endif

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        if ( L(ikx,iky,ikz).ne.1 ) then
          zxnk(ikx,iky,ikz) = cmplx(0.,0.)
          zynk(ikx,iky,ikz) = cmplx(0.,0.)
          zznk(ikx,iky,ikz) = cmplx(0.,0.)
          ttnk(ikx,iky,ikz) = cmplx(0.,0.)
        else
          kx  = kxa(ikx)             
          kh  = sqrt( kx*kx + ky*ky )
          khn = kh * L1/twopi
          wk2 = kx*kx + ky*ky + kz*kz
          k   = sqrt( wk2 )

          c1 = +  ky*zznk(ikx,iky,ikz) -  kz*zynk(ikx,iky,ikz)
          c2 = +  kz*zxnk(ikx,iky,ikz) -  kx*zznk(ikx,iky,ikz)
          c3 = +  kx*zynk(ikx,iky,ikz) -  ky*zxnk(ikx,iky,ikz)
          u  = zi * c1 / wk2
          v  = zi * c2 / wk2
          w  = zi * c3 / wk2

          ! spherical dissipation
          r1 = v2h/delt*wk2**ilap
          r2 = k2h/delt*wk2**ilap

          ! cylindrical dissipation
!          r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
!          r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2


          if (khn.lt.0.1) then
            r1 = r1 + ek/delt
            r2 = r2 + ek/delt
          endif
          
          termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttnk(ikx,iky,ikz) + cor*zi*kz*u - zxnk(ikx,iky,ikz)*r1
          termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttnk(ikx,iky,ikz) + cor*zi*kz*v - zynk(ikx,iky,ikz)*r1
          termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zznk(ikx,iky,ikz)*r1
          termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttnk(ikx,iky,ikz)*r2

          zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz)  +   delt*(termzx + rhzx(ikx,iky,ikz))/2.
          zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz)  +   delt*(termzy + rhzy(ikx,iky,ikz))/2.
          zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz)  +   delt*(termzz + rhzz(ikx,iky,ikz))/2.
          ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz)  +   delt*(termtt + rhtt(ikx,iky,ikz))/2.

          zxok(ikx,iky,ikz) = rhzx(ikx,iky,ikz)  + zxok(IKX,IKY,IKZ)*r1
          zyok(ikx,iky,ikz) = rhzy(ikx,iky,ikz)  + zyok(IKX,IKY,IKZ)*r1
          zzok(ikx,iky,ikz) = rhzz(ikx,iky,ikz)  + zzok(IKX,IKY,IKZ)*r1
          ttok(ikx,iky,ikz) = rhtt(ikx,iky,ikz)  + ttok(IKX,IKY,IKZ)*r2
                  
        endif ! L
      enddo
    enddo
  enddo

  time3 = etime(time1)

! -----------------------------------------------------------------
!                        Subsequent Timesteps
! -----------------------------------------------------------------

  nt0=3   ! AB3: start at nt=3

  do nt = nt0,NSTOP
    time = ts + nt*delt

    if (LINEAR.ne.1) then
      call convol (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr,rhzx,rhzy,rhzz,rhtt,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr)
    else
      nzxk = cmplx(0.,0.)
      nzyk = cmplx(0.,0.)
      nzzk = cmplx(0.,0.)
      nttk = cmplx(0.,0.)
    endif

    if (max(ampv,ampw).ne.0.) then
      call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)
    endif

    call velo (zxnk,zynk,zznk,rhzx,rhzy,rhzz)

    do ikz = 1,iktzp
      ikza = mype*iktzp+ikz
      kz = kza(ikza)
      do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
          kx  = kxa(ikx)
          kh  = sqrt( kx*kx + ky*ky )
          khn = kh * L1/twopi
          wk2 = kx*kx + ky*ky + kz*kz
          k   = sqrt( wk2 )
          ! Equations:
          termzx = nzxk(ikx,iky,ikz) +  aj*ZI*ky*ttnk(ikx,iky,ikz) + cor*ZI*kz*rhzx(ikx,iky,ikz)
          termzy = nzyk(ikx,iky,ikz) -  aj*ZI*kx*ttnk(ikx,iky,ikz) + cor*ZI*kz*rhzy(ikx,iky,ikz)
          termzz = nzzk(ikx,iky,ikz)                               + cor*ZI*kz*rhzz(ikx,iky,ikz)
          termtt = nttk(ikx,iky,ikz)                               -        bj*rhzz(ikx,iky,ikz)

          ! spherical dissipation
          r1 = v2h*wk2**ilap 
          r2 = k2h*wk2**ilap

          ! cylindrical dissipation
!          r1 = v2h*kh**ilap2 + v2z*kz**ilap2 
!          r2 = k2h*kh**ilap2 + k2z*kz**ilap2


          if (khn.lt.0.1) then
            r1 = r1 + ek
            r2 = r2 + ek
          endif

          dpz = 1. + r1/2.
          dmz = 1. - r1/2.
          tzx = ( zxnk(ikx,iky,ikz)*dmz + D12*( 23.*termzx - 16.*zxok(ikx,iky,ikz) & 
              + 5.*zx1k(ikx,iky,ikz) ) )/dpz
          tzy = ( zynk(ikx,iky,ikz)*dmz + D12*( 23.*termzy - 16.*zyok(ikx,iky,ikz) & 
              + 5.*zy1k(ikx,iky,ikz) ) )/dpz
          tzz = ( zznk(ikx,iky,ikz)*dmz + D12*( 23.*termzz - 16.*zzok(ikx,iky,ikz) & 
              + 5.*zz1k(ikx,iky,ikz) ) )/dpz

          dpz = 1. + r2/2.
          dmz = 1. - r2/2.
          ttt = ( ttnk(ikx,iky,ikz)*dmz + D12*( 23.*termtt - 16.*ttok(ikx,iky,ikz) & 
              + 5.*tt1k(ikx,iky,ikz) ) )/dpz

          zx1k(ikx,iky,ikz) = zxok(ikx,iky,ikz)*L(ikx,iky,ikz)
          zy1k(ikx,iky,ikz) = zyok(ikx,iky,ikz)*L(ikx,iky,ikz)
          zz1k(ikx,iky,ikz) = zzok(ikx,iky,ikz)*L(ikx,iky,ikz)
          tt1k(ikx,iky,ikz) = ttok(ikx,iky,ikz)*L(ikx,iky,ikz)

          zxok(ikx,iky,ikz) = termzx*L(ikx,iky,ikz)
          zyok(ikx,iky,ikz) = termzy*L(ikx,iky,ikz)
          zzok(ikx,iky,ikz) = termzz*L(ikx,iky,ikz)
          ttok(ikx,iky,ikz) = termtt*L(ikx,iky,ikz)

          zxnk(ikx,iky,ikz) =    tzx*L(ikx,iky,ikz)
          zynk(ikx,iky,ikz) =    tzy*L(ikx,iky,ikz)
          zznk(ikx,iky,ikz) =    tzz*L(ikx,iky,ikz)
          ttnk(ikx,iky,ikz) =    ttt*L(ikx,iky,ikz)

        enddo
      enddo
    enddo


! -----------------------------------------------------------------
! Diagnostics and i/o
! -----------------------------------------------------------------

    ! Write to restart file
    if ( mod(nt,NDUMP).eq. 0 ) then
      ntdump=nt/ndump
      if (rstflag.eq.1) call ncdumprst(zxnk,zynk,zznk,ttnk,ntdump)
    endif ! ndump

    ! Compute diagnostics
    if ( mod(nt,nout).eq. 0 ) then
      call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk)
      call proj(zxnk,zynk,zznk)  
      call out(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k) 
      call spec(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,2,52)
      call spec(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,3,51)
      call spec(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,1,53)
      call cfl(zxnk,zynk,zznk,uk,vk,wk,ur,vr,wr)
      if (linear.ne.1) then
        ! Compute basic transfer   
        call convol(zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,uk,vk,wk,ur,vr,wr, &
             rhzx,rhzy,rhzz,rhtt,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr)
        call transf(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,2,48)
        call transf(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,3,49)
        call transf(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,1,50)
      endif ! linear
    endif ! nout 

    ! Write real space fields                                                               
    if (rsflag.eq.1 .and. mod(nt,nbig).eq.0) then
      ntdump=nt/nbig+1
      call dumpreal(zxnk,zynk,zznk,ttnk,zxnr,zynr,zznr,ttnr,uk,vk,wk,ur,vr,wr,ntdump,FULLFIELD_MODE)
    endif

    ! output slices
    if ( (slxflag.eq.1 .or. slyflag.eq.1 .or. slzflag.eq.1) .and. mod(nt,nslout).eq.0 ) then
      ntdump=nt/nslout+1
      call dumpreal(zxnk,zynk,zznk,ttnk,zxnr,zynr,zznr,ttnr,uk,vk,wk,ur,vr,wr,ntdump,SLICE_MODE)    
    endif
  enddo ! nt

  time3 = etime(time2)
  time3 = time2(1) - time1(1)

! CLOSE I/O
  call io_close()

  if (mype.eq.0) then
    print*,'    '
    write(6,5000) time3, time3/60., time3/3600., time3/86400.
    print*,'    '
  endif

  ! destroy FFTW plans, free allocated memory.
  call fftwf_destroy_plan(plan3_ur_uk) 
  call fftwf_destroy_plan(plan3_vr_vk) 
  call fftwf_destroy_plan(plan3_wr_wk) 
  call fftwf_destroy_plan(plan3_zxnr_zxnk)
  call fftwf_destroy_plan(plan3_zynr_zynk)
  call fftwf_destroy_plan(plan3_zznr_zznk)
  call fftwf_destroy_plan(plan3_ttnr_ttnk)  
  call fftwf_destroy_plan(plan3_nzxr_nzxk)
  call fftwf_destroy_plan(plan3_nzyr_nzyk)
  call fftwf_destroy_plan(plan3_nzzr_nzzk)
  call fftwf_destroy_plan(plan3_uk_ur)
  call fftwf_destroy_plan(plan3_vk_vr)
  call fftwf_destroy_plan(plan3_wk_wr)
  call fftwf_destroy_plan(plan3_zxnk_zxnr)
  call fftwf_destroy_plan(plan3_zynk_zynr)
  call fftwf_destroy_plan(plan3_zznk_zznr)
  call fftwf_destroy_plan(plan3_ttnk_ttnr)  
  call fftwf_mpi_cleanup()
  call fftwf_free(zxo_p)
  call fftwf_free(zyo_p)
  call fftwf_free(zzo_p)
  call fftwf_free(tto_p)
  call fftwf_free(zxn_p)
  call fftwf_free(zyn_p)
  call fftwf_free(zzn_p)
  call fftwf_free(ttn_p)
  call fftwf_free(nzx_p)
  call fftwf_free(nzy_p)
  call fftwf_free(nzz_p)
  call fftwf_free(ntt_p)
  call fftwf_free(u_p)
  call fftwf_free(v_p)
  call fftwf_free(w_p)  
  deallocate(L)
  deallocate(kxa)
  deallocate(kya) 
  deallocate(kza)
  deallocate(rhzx) 
  deallocate(rhzy) 
  deallocate(rhzz) 
  deallocate(rhtt) 
  deallocate(geok) 
  deallocate(gw1k) 
  deallocate(gw2k)  
  deallocate(zx1k)
  deallocate(zy1k) 
  deallocate(zz1k) 
  deallocate(tt1k)
  call mpi_finalize(istatus)

5001 format(1x,'CPU time - initialization = ',F7.0,' s = ', F7.1,' m = ',F7.2,' h =',F7.3,' d.')
5000 format(1x,'CPU time - main loop = ',F7.0,' s = ', F7.1,' m = ',F7.2,' h =',F7.3,' d.')
end PROGRAM MAIN



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! INITIALIZATION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine init (zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,zxr,zyr,zzr,ur,vr,wr,tr,      & 
                 ampk,ampp,irest,ts,ki)

! Initialises stuff like wavenumbers, indices, spectra, phases, etc.

  use param
  use param_fftw

  implicit none
  integer, intent(in) :: irest
  real,    intent(in) :: ampk,ampp,ki

  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz, tt,uu,vv,ww
  complex, intent(out), dimension(iktx,ikty,iktzp) :: ge,g1,g2
  real,    intent(out), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,tr,ur,vr,wr
  real,    intent(out) :: ts

  integer :: ikx,iky,ikz,ikza,i,j,k,ja,jj,inormke,inormve
  real ::  kx,ky,kz,wk,kh,khn,wkn,kzn,ek
  real ::  phase,kinen,poten,vorten,waven
  real ::  ranno,noiseamp
  real ::  r0,u0,x,y,z,r,dxtan,dytan,tanhx,tanhy
  real ::  bessj0,bessj1,x0,y0,pert,kzpert,xr,yr
  real, parameter :: mu1=3.83170597020751
  external :: ranno,proj,realit

! Initialize wavenumber arrays
  do  ikx = 1,iktx
    kxa(ikx) = float(ikx-1) * twopi/L1
  enddo
  do iky = 1,ikty
    jj = iky - 1
    if (iky.gt.kty)   jj = jj - 2*kty
    if (iky.eq.kty+1) jj = 0
    if (iky.gt.2*kty) jj = 0
    kya(iky) = float(jj) * twopi/L2
  enddo
  do ikz = 1,iktz
    jj = ikz - 1
    if (ikz.gt.ktz)   jj = jj - 2*ktz
    if (ikz.eq.ktz+1) jj = 0
    if (ikz.gt.2*ktz) jj = 0
    kza(ikz) = float(jj) * twopi/L3
    ! if (mype.eq.0) print*,ikz,kza(ikz)
  enddo

  ! L(ikx,iky,ikz) is unity for retained modes and zero beyond the truncation, 
  ! at k=0 and for half the modes on the plane kx=0.

  L  = 1
  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)
  tt = cmplx(0.,0.)
  uu = cmplx(0.,0.)
  vv = cmplx(0.,0.)
  ww = cmplx(0.,0.)
  ge = cmplx(0.,0.)
  g1 = cmplx(0.,0.)
  g2 = cmplx(0.,0.)
  zxr = 0.
  zyr = 0.
  zzr = 0.
  tr  = 0.
  ur  = 0.
  vr  = 0.
  wr  = 0.

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        kx = kxa(ikx)
        kh = sqrt(kx*kx + ky*ky)
        wk = kx*kx + ky*ky + kz*kz
        wk = sqrt( wk )
           
        ! Set L=0 where necessary:
        if (iky.eq.kty+1)                         L(ikx,iky,ikz) = 0
        if (ikza.eq.ktz+1)                        L(ikx,iky,ikz) = 0
        if (kx.lt.0)                              L(ikx,iky,ikz) = 0
        if (kx.eq.0 .and. ky.lt.0)                L(ikx,iky,ikz) = 0
        if (kx.eq.0 .and. ky.eq.0 .and. kz.lt.0)  L(ikx,iky,ikz) = 0
        if (iky.gt.2*kty)                         L(ikx,iky,ikz) = 0
        if (ikza.gt.2*ktz)                        L(ikx,iky,ikz) = 0


        ! spherical truncation
        !if (wk.eq.0. .or. wk.gt.ifix(KTRUNC_X + 0.5) - 0.5)                L(ikx,iky,ikz) = 0
        !if (wk.eq.0. .or. (wk*L1/twopi).gt.ifix(float(n1)/3.+0.5)-0.5)     L(ikx,iky,ikz) = 0

        ! cylindrical truncation
        !if (wk.eq.0. .or. (kh*L1/twopi).gt.ifix(float(n1)/3.+0.5)-0.5)     L(ikx,iky,ikz) = 0         
        !if (abs(kz*L3/twopi).gt.ifix(float(n3)/3. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          

        ! cubic truncation
         if (wk.eq.0.)                                                     L(ikx,iky,ikz) = 0         
!         if (abs(kx*L1/twopi).gt.int(float(n1)/3. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          
!         if (abs(ky*L2/twopi).gt.int(float(n2)/3. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          
!         if (abs(kz*L3/twopi).gt.int(float(n3)/3. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          

!!!    ! cubic truncation with 8/9 instead of 2/3 truncation
         if (abs(kx*L1/twopi).gt.int(float(n1)*4./9. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          
         if (abs(ky*L2/twopi).gt.int(float(n2)*4./9. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          
         if (abs(kz*L3/twopi).gt.int(float(n3)*4./9. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          


!         if (abs(kx*L1/twopi).gt.255.)              L(ikx,iky,ikz) = 0          
!         if (abs(ky*L2/twopi).gt.255.)              L(ikx,iky,ikz) = 0          
!         if (abs(kz*L3/twopi).gt.255.)              L(ikx,iky,ikz) = 0          


      enddo
    enddo
  enddo

  !  Restart from output file?
  if (irest.ne.0) then
    call ncreadrst(zx,zy,zz,tt,ur,vr,ts)

     
  else ! irest = 0, Initialize all the arrays: 

    !!! Initialize in Physical space, standing wave
!    do i=1,n1d
!      x = float(i-1)/float(n1)*L1
!      do j=1,n2dp
!        ja = mype*n2p+j
!        y = float(ja-1)/float(n2)*L2
!        do k=1,n3d
!          z = float(k-1)/float(n3)*L3
!          ur(i,k,j) =  0.
!          vr(i,k,j) =  0.
!          wr(i,k,j) =  0.
!          tr(i,k,j) =  cos(x)*cos(z)
!        enddo
!      enddo
!    enddo

 
!!! Initialize in Fourier space     
     do ikx = 1,iktx
        kx = kxa(ikx)
        do iky = 1,ikty
           ky = kya(iky)
           kh = sqrt(kx*kx + ky*ky)
           do ikz = 1,iktzp
              ikza = mype*iktzp+ikz
              kz = kza(ikza)

              wk = sqrt(kx*kx + ky*ky + kz*kz)
              khn = kh * L1/twopi
              wkn = wk * L1/twopi
              kzn = kz * L3/twopi
              if (L(ikx,iky,ikz).eq.1) then

                EK = 0
                if (wkn.gt.ki-.5 .and. wkn.le.ki+.5) then 
                   EK = 1. 
                endif

                phase = ranno(0)
                zx(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
                phase = ranno(0)
                zy(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
                phase = ranno(0)
                zz(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
                phase = ranno(0)
                tt(ikx,iky,ikz) =    sqrt(EK)*cexp(ZI*phase)

!  Seed each mode with a small, equal amount of KE and PE: 
                 noiseamp = 1.e-3

                 if (wkn.lt.5) then
                    phase = ranno(0)   
                    zx(ikx,iky,ikz) = zx(ikx,iky,ikz) + noiseamp*wk*cexp(ZI*phase)
                    phase = ranno(0)
                    zy(ikx,iky,ikz) = zy(ikx,iky,ikz) + noiseamp*wk*cexp(ZI*phase)
                    phase = ranno(0)
                    zz(ikx,iky,ikz) = zz(ikx,iky,ikz) + noiseamp*wk*cexp(ZI*phase)
                    !phase = ranno(0)
                    !tt(ikx,iky,ikz) = tt(ikx,iky,ikz) + noiseamp*cexp(ZI*phase)
                 endif

                 if (kx.eq.0.) then
                    zx(ikx,iky,ikz) = cmplx(real(zx(ikx,iky,ikz)),0.)
                    zy(ikx,iky,ikz) = cmplx(real(zy(ikx,iky,ikz)),0.)
                    zz(ikx,iky,ikz) = cmplx(real(zz(ikx,iky,ikz)),0.)
                 endif

              endif
           enddo
        enddo
     enddo

     ! KILL wave energy in ICs
     call wtoab(zx,zy,zz,tt,ge,g1,g2,uu,vv,ww)
     g1 = 0.
     g2 = 0.
     call atowb(ge,g1,g2,zx,zy,zz,tt,uu,vv,ww)

     call proj(zx,zy,zz)                          

     call energy_full(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2, &
          kinen,poten,vorten,waven)

    inormke = 1
    inormve = 0

    if (inormke.eq.1) then  
      zx = zx * sqrt(ampk/(kinen))
      zy = zy * sqrt(ampk/(kinen))
      zz = zz * sqrt(ampk/(kinen))
      tt = tt * sqrt(ampp/(poten))
    endif
     
    if (inormve.eq.1) then  
      call wtoab(zx,zy,zz,tt,ge,g1,g2,uu,vv,ww)
      ge = ge * sqrt(ampk/vorten)            
      g1 = g1 * sqrt(ampp/waven)            
      g2 = g2 * sqrt(ampp/waven)            
      call atowb(ge,g1,g2,zx,zy,zz,tt,uu,vv,ww)
    endif ! inormve

  endif ! irest      
end subroutine init


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! FORCING
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,nfmax,alpha,kf)

! Use this after call convol.  Forcing is specified to be divergence-free
! Forcing follows Herring & Metais, 1989, J. Fluid Mech 202, pp. 97-115

  use param

  implicit none
  integer, intent(in) :: nfmax
  real,    intent(in) :: ampv,ampw,alpha,kf
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: nzxk,nzyk,nzzk,nttk
  complex, intent(inout) :: gtau(nfmax,4)

  integer :: ikx,iky,ikz,ikza,ik
  complex :: fu,fv,fw,ft,f1k,f2k,f3k,ftk
  real :: rang,beta,c,g,tf,theta
  real :: kx,ky,kz,wkh,wk,wkh2,sk
  external :: proj,rang

  beta  = sqrt(1.-alpha**2)
  c     = sqrt(bj/(2.*aj))
  ik    = 0

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        if (L(ikx,iky,ikz).eq.1) then 
          kx = kxa(ikx)
          wkh2 = kx**2 + ky**2
          wkh  = sqrt(wkh2) 
          wk   = sqrt(kx**2 + ky**2 + kz**2)
          theta = wkh/wk
          tf    = 1./sqrt2
       

          ! FORCE modes inside spherical shell centered on k=kf, AND exclude VSHF   
          if (abs(wk*L1/twopi-KF).le.1. .and. (wkh*L1/twopi).gt.0.1) then
       
            ik=ik+1                   
            gtau(ik,1) = alpha*gtau(ik,1) + beta*cmplx(rang(0),rang(0))
            gtau(ik,2) = alpha*gtau(ik,2) + beta*cmplx(rang(0),rang(0))
            gtau(ik,3) = alpha*gtau(ik,3) + beta*cmplx(rang(0),rang(0))
            gtau(ik,4) = alpha*gtau(ik,4) + beta*cmplx(rang(0),rang(0))

            !!! specify amplitudes of vort and wave forcing separately
            G = (wk*L1/twopi-(kf-1.))*(kf+1.-wk*L1/twopi) 
            Fv = ampv * G * gtau(ik,1) 
            Fw = ampw * G * ( gtau(ik,2) + gtau(ik,3))
            Ft = ampw * G * (-gtau(ik,2) + gtau(ik,3))
            sk = sqrt(BF2*wkh2 + cor2*kz**2)

            f1k  = - kx*kz*BF/sk*Fv - wk*ky/wkh/sqrt2*Fw + ZI*cor*kx*kz**2/sqrt2/sk/wkh*Ft
            f2k  = - ky*kz*BF/sk*Fv + wk*kx/wkh/sqrt2*Fw + ZI*cor*ky*kz**2/sqrt2/sk/wkh*Ft
            f3k  =   BF*wkh**2/sk*Fv - ZI*cor*wkh*kz/sqrt2/sk*Ft
            ftk  = - C*sqrt2*ZI*cor*kz/sk*Fv + C*BF*wkh/sk*Ft

            nzxk(ikx,iky,ikz) = nzxk(ikx,iky,ikz) + f1k
            nzyk(ikx,iky,ikz) = nzyk(ikx,iky,ikz) + f2k
            nzzk(ikx,iky,ikz) = nzzk(ikx,iky,ikz) + f3k
            nttk(ikx,iky,ikz) = nttk(ikx,iky,ikz) + ftk

            !!! new: force everything isotropically. Use ampv as global amplitude.
            ! G = (wk*L1/twopi-(kf-1.))*(kf+1.-wk*L1/twopi) 
            ! Fu = ampv * G * gtau(ik,1) 
            ! Fv = ampv * G * gtau(ik,2) 
            ! Fw = ampv * G * gtau(ik,3) 
            ! Ft = ampv * G * gtau(ik,4) *sqrt(bj/aj)
            ! nzxk(ikx,iky,ikz) = nzxk(ikx,iky,ikz) + zi*(ky*Fw-kz*Fv)
            ! nzyk(ikx,iky,ikz) = nzyk(ikx,iky,ikz) + zi*(kz*Fu-kx*Fw)
            ! nzzk(ikx,iky,ikz) = nzzk(ikx,iky,ikz) + zi*(kx*Fv-ky*Fu)
            ! nttk(ikx,iky,ikz) = nttk(ikx,iky,ikz) + Ft


            if (ik.eq.nfmax) then
              print*,'Forcing too many modes!'
              stop
            endif

          endif
        endif
      enddo
    enddo
  enddo
        
  return
end subroutine force



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! DIAGNOSTICS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spec (zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,ispec,iu)

! Calculates spectra using: 
! ispec = 1: total wavenumber k
! ispec = 2: horizontal wavenumber kh
! ispec = 3: vertical wavenumber kz

  use param
  implicit none 
  include 'mpif.h'

  integer, intent(in) :: ispec,iu
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz
  
  integer :: ikx,iky,ikz,ikza,j,i,j0
  integer :: nspz,nspz11,nspz13
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2
  real    :: vt,kvisc,vzx,vzy,vzz
  complex :: div
  real, dimension(0:kts,13)       :: spz,spztot
  integer, dimension(0:kts)       :: n,ntot
  external :: velo

  nspz     = kts+1
  nspz11   = (kts+1)*11
  nspz13   = (kts+1)*13

  if (ispec.eq.1) then
    j0 = 0 !1
  elseif (ispec.eq.2) then
    j0 = 0
  elseif (ispec.eq.3) then
    j0 = 0
  else
    print*,"ispec error"
    stop
  endif

  N     = 0
  spz   = 0.

  call velo (zx,zy,zz,ux,uy,uz)

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    kz2 = kz*kz
    kzn = kz*L3/twopi
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        wkh2  = kx*kx+ky*ky
        wkh2n = wkh2 * (L1/twopi)**2
        wkh   = sqrt(wkh2)
        wk    = sqrt(kx*kx+ky*ky+kz*kz)
        kvisc = visch*wkh2**ilap+viscz*kz2**ilap
           
        if (ispec.eq.1) then
          j = int(wk*L3/twopi+0.5)
        elseif (ispec.eq.2) then
          j = int(wkh*L1/twopi+0.5)
        elseif (ispec.eq.3) then
          j = int(abs(kz)*L3/twopi+0.5)
          ! print*,kz,j
        endif

        if (L(ikx,iky,ikz).eq.1) then
          if (j.lt.j0 .or. j.gt.kts) print*,'SPEC: SCREW-UP.',j,'ktx= ',ktx
          wk   = max(wk,  1.e-15)
          wkh2 = max(wkh2,1.e-15)

          ! Kinetic and potential energy.
          vzx      = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
          vzy      = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
          vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
          vt       = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )

          spz(j,1) = spz(j,1) + vzx/wk**2 + vzy/wk**2 + vzz/wk**2
          spz(j,2) = spz(j,2) + vt*aj/bj

          ! KE and PE dissipation
          spz(j,7) = spz(j,7) + kvisc*(vzx/wk**2 + vzy/wk**2 + vzz/wk**2)     
          spz(j,8) = spz(j,8) + kvisc*vt*aj/bj


          ! Geo, ageo decompostition.
          ! k \in R_k
          if (wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
            vzx      = real( ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)) )
            vzy      = real( g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)) )
            vzz      = real( g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)) )

            spz(j,4) = spz(j,4) + vzx/wkh2
            spz(j,5) = spz(j,5) + vzy/wkh2 + vzz/wkh2
            ! spz(j,7) = spz(j,7) + kvisc*vzx/wkh2
            ! spz(j,8) = spz(j,8) + kvisc*vzy/wkh2 + kvisc*vzz/wkh2
          endif

          ! Special cases: i) k_h=0, ii) k_z=0.
          vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
          vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
          vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

          ! k \in V_k
          if(wkh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
            spz(j,4) = spz(j,4) + vzz + vt*aj/bj
            spz(j,5) = spz(j,5) + vzx + vzy
            ! spz(j,7) = spz(j,7) + kvisc*(vzz + vt*aj/bj)
            ! spz(j,8) = spz(j,8) + kvisc*(vzx + vzy)     
          endif

          ! k \in B_k
          if(abs(kzn).lt.1.e-10.and.wkh2n.gt.1.e-10) then
            spz(j,4) = spz(j,4) + vzx + vzy
            spz(j,5) = spz(j,5) + vzz + vt*aj/bj
            ! spz(j,7) = spz(j,7) + kvisc*(vzx + vzy)     
            ! spz(j,8) = spz(j,8) + kvisc*(vzz + vt*aj/bj)
          endif

          ! Buoyancy Flux
          if (aj.gt.0.) then 
            vt = aj*real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
          else
            vt = real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
          endif
          spz(j,6) = spz(j,6) + vt

          ! KE decomposed into u/v/w
          spz(j,9)  = spz(j,9)  + vzx
          spz(j,10) = spz(j,10) + vzy
          spz(j,11) = spz(j,11) + vzz

	  ! Rotational and divergent KE
          div      = kx*ux(ikx,iky,ikz)+ky*uy(ikx,iky,ikz)
          vzx      = real( div*conjg(div) )/wkh2
          vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )/wkh2
	  spz(j,12) = spz(j,12) + vzz
	  spz(j,13) = spz(j,13) + vzx           

          n(j) = n(j) + 2 

        endif
      enddo
    enddo
  enddo

  spz(:,3) = spz(:,1) + spz(:,2)

  call mpi_reduce(spz,spztot,nspz13,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,istatus)
  call mpi_reduce(  n,  ntot,  nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,istatus)

  if (mype.eq.0) then
    do j=j0,kts-1
      if (ntot(j).ne.0) then
        write(iu,5000) float(j),(spztot(j,i),i=1,13),ntot(j)
      endif
    enddo
    write(iu,*) '           '
    write(iu,*) '           '
    call flush(iu)
  endif

  return
5000  format(1X,F4.0,4X,13(E13.6,1x),6X,I10)
end subroutine spec


subroutine transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k,nk1,nk2,nk3,ispec,iu)

! Calculates transfer spectra.
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  use param
  implicit none 
  include 'mpif.h'

  integer, intent(in) :: ispec,iu
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,nzx,nzy,nzz,ntt
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: geok,gw1k,gw2k,ngeok,ngw1k,ngw2k,nk1,nk2,nk3

  integer :: ikx,iky,ikz,ikza,j,j0,nspz,nspz6
  integer, dimension(0:kts) :: n,ntot
  real :: kx,ky,kz,k,k2,kh2,wkh,kh2n,kzn
  real :: vzx,vzy,vzz,vtt
  real, dimension(0:kts,6)  :: spz,spztot
  complex :: u,v,w,c1,c2,c3

  nspz  = kts+1
  nspz6 = (kts+1)*6

  if (ispec.eq.1) then
    j0 = 1
  elseif (ispec.eq.2) then
    j0 = 0
  elseif (ispec.eq.3) then
    j0 = 0
  else
    print*,"ispec error"
    stop
  endif
  
  call wtoab( zx, zy, zz, tt, geok, gw1k, gw2k,nk1,nk2,nk3)
  call wtoab(nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k,nk1,nk2,nk3)
  call velo(nzx,nzy,nzz,nk1,nk2,nk3)

  N   = 0
  spz = 0.

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    kzn = kz*L3/twopi
    do iky = 1,ikty
      ky  = kya(iky)
      do ikx = 1,iktx
        kx   = kxa(ikx)
        kh2  = kx*kx + ky*ky
        kh2n = kh2 * (L1/twopi)**2
        kh2  = max(1.e-15,kh2)
        wkh  = sqrt(kh2)
        k2   = kx*kx + ky*ky + kz*kz
        k    = sqrt(k2)
        k2   = max(1.e-15,k2)
           
        if (ispec.eq.1) then
          j = int(k*L1/twopi+0.5)
        elseif (ispec.eq.2) then
          j = int(wkh*L1/twopi+0.5)
        elseif (ispec.eq.3) then
          j = int(abs(kz)*L3/twopi+0.5)
        endif

        if (L(ikx,iky,ikz).eq.1)  then
          if (j.lt.j0 .or. j.gt.kts) then
            print*,'transf: screw-up.   k= ',j,kx,ky,kz,L(ikx,iky,ikz)
          endif

          c1 = ky*zz(ikx,iky,ikz) - kz*zy(ikx,iky,ikz)
          c2 = kz*zx(ikx,iky,ikz) - kx*zz(ikx,iky,ikz)
          c3 = kx*zy(ikx,iky,ikz) - ky*zx(ikx,iky,ikz)
          u = zi * c1 / k2
          v = zi * c2 / k2
          w = zi * c3 / k2

          ! k \in R_k
          if (kh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
            spz(j,1) = spz(j,1) + real(geok(ikx,iky,ikz)*conjg(ngeok(ikx,iky,ikz)))/kh2
            spz(j,2) = spz(j,2) + real(gw1k(ikx,iky,ikz)*conjg(ngw1k(ikx,iky,ikz)))/kh2
            spz(j,2) = spz(j,2) + real(gw2k(ikx,iky,ikz)*conjg(ngw2k(ikx,iky,ikz)))/kh2
          endif

          ! Special cases: i) k_h=0, ii) k_z=0.
          vzx = real(u              *conjg( nk1(ikx,iky,ikz)))
          vzy = real(v              *conjg( nk2(ikx,iky,ikz)))
          vzz = real(w              *conjg( nk3(ikx,iky,ikz)))
          vtt = real(tt(ikx,iky,ikz)*conjg(ntt(ikx,iky,ikz)))

          !  k \in V_k
          if (kh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
            spz(j,1) = spz(j,1) + vtt
            spz(j,2) = spz(j,2) + vzx + vzy            
          endif

          !  k \in B_k
          if (kh2n.gt.1.e-10.and.abs(kzn).lt.1.e-10) then
            spz(j,1) = spz(j,1) + vzx + vzy 
            spz(j,2) = spz(j,2) + vzz + vtt*aj/bj
          endif

          !  Now, compute KE and PE transfer
          spz(j,4) = spz(j,4) + vzx + vzy + vzz
          spz(j,5) = spz(j,5) + vtt*aj/bj

          n(j)   = n(j) + 2
        endif
      enddo
    enddo
  enddo

  spz(:,3) = spz(:,1) + spz(:,2)
  spz(:,6) = spz(:,4) + spz(:,5)

  call mpi_reduce(spz,spztot,nspz6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,istatus)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,istatus)

  if (mype.eq.0) then  
    do j = j0,kts-1
      if (ntot(j).ne.0) then
        write(iu,5000) float(j),(spztot(j,iky),iky=1,6),ntot(j)
      endif
    enddo
    write(iu,*) '     '
    write(iu,*) '     '
  endif

  return
5000 format(1X,F4.0,4X,6(E11.4,1x),10X,I6)
end subroutine transf



subroutine out(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2)

  use param
  implicit none
  include 'mpif.h'

  complex, intent(in), dimension(iktx,ikty,iktzp) ::  zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz

  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,wk,wkh2,wkh2n,kzn
  real :: sigma
  real :: vh,vz,vzx,vzy,vzz
  real :: vort2,div2,rossby,fr_z,fr_h,ke,pe,e,eg,ea,h,v,epsk,epsp,eps,epskh,epskv,epsph,epspv,tmp
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo
  complex:: der
  external :: proj, velo

  call velo(zx,zy,zz,ux,uy,uz)

  vort2  = 0.
  div2   = 0.
  Rossby = 0.
  Fr_z   = 0.
  Fr_h   = 0.
  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea = 0.
  h  = 0.
  v  = 0.
  epsk = 0.
  epsp = 0.
  epskh = 0.
  epsph = 0.
  epskv = 0.
  epspv = 0.
  eps  = 0.

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    kzn = kz * (L3/twopi)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        kx = kxa(ikx)
        wkh2 = kx*kx + ky*ky
         wk = kx*kx + ky*ky + kz*kz
        if (L(ikx,iky,ikz).eq.1) then
          wkh2n = wkh2 * (L1/twopi)**2
              
          vzx   = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)))
          ke    = ke + vzx/wk
          
          epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzx/wk
          epskh = epskh + (visch*wkh2**ilap)*vzx/wk
          epskv = epskv + (viscz*kz**(2*ilap))*vzx/wk
          v     = v  + vzx
          vzy   = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)))
          ke    = ke + vzy/wk

          epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzy/wk
          epskh = epskh + (visch*wkh2**ilap)*vzy/wk
          epskv = epskv + (viscz*kz**(2*ilap))*vzy/wk
          v     = v + vzy
          vzz   = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)))
          ke    = ke + vzz/wk
         
          epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzz/wk
          epskh = epskh + (visch*wkh2**ilap)*vzz/wk
          epskv = epskv + (viscz*kz**(2*ilap))*vzz/wk
          v     = v  + vzz
          vort2 = vort2 + vzz

          vh    = real( zx(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)))
          h     = h + vh
          vh    = real( zy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)))
          h     = h + vh
          vh    = real( zz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)))
          h     = h + vh
          vh    = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)))
          pe    = pe  + vh
          epsp  = epsp + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vh*aj/bj
          epsph = epsph + (visch*wkh2**ilap)*vh*aj/bj
          epspv = epspv + (viscz*kz**(2*ilap))*vh*aj/bj

          vzx   = real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)))
          vzy   = real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)))
          vzz   = real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)))

          if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
            zero_kh_grv = zero_kh_grv + vzx + vzy
            zero_kh_geo = zero_kh_geo + vh*aj/bj
          endif
              
          if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
            zero_kz_geo = zero_kz_geo + vzx + vzy
            zero_kz_grv = zero_kz_grv + vzz + vh*aj/bj
          endif
              
          if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
            vz    = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
            eg    = eg  + vz/wkh2
            sigma = wkh2*bf2 + kz*kz*cor2
            vz    = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
            ea    = ea  + vz/wkh2
            vz    = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
            ea    = ea  + vz/wkh2
          endif
             
          rossby = rossby + real(zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)))/(cor2)
          der    =          real(zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)))
          der    = (der +   real(zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz))))
          fr_z   = fr_z + der/bf2
          der    = zi*(kx*ux(ikx,iky,ikz) + ky*uy(ikx,iky,ikz))
          div2   = div2 + real(der*conjg(der))
        endif ! L
      enddo
    enddo
  enddo

  v      = 2.*v
  epsk   = 2.*epsk
  epsp   = 2.*epsp
  epskh  = 2.*epskh
  epsph  = 2.*epsph
  epskv  = 2.*epskv
  epspv  = 2.*epspv
  eg     = eg + zero_kz_geo + zero_kh_geo
  ea     = ea + zero_kz_grv + zero_kh_grv
  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj

  call mpi_reduce(ke,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          ke=tmp
  call mpi_reduce(pe,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          pe=tmp
  call mpi_reduce(eg,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          eg=tmp
  call mpi_reduce(ea,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          ea=tmp
  call mpi_reduce(rossby,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);      rossby=tmp
  call mpi_reduce(fr_z,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        fr_z=tmp
  call mpi_reduce(vort2,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       vort2=tmp
  call mpi_reduce(v,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);           v=tmp
  call mpi_reduce(zero_kz_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); zero_kz_grv=tmp
  call mpi_reduce(zero_kz_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); zero_kz_geo=tmp
  call mpi_reduce(zero_kh_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); zero_kh_grv=tmp
  call mpi_reduce(zero_kh_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); zero_kh_geo=tmp
  call mpi_reduce(epsk,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        epsk=tmp
  call mpi_reduce(epsp,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        epsp=tmp
  call mpi_reduce(epskh,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epskh=tmp
  call mpi_reduce(epsph,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epsph=tmp
  call mpi_reduce(epskv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epskv=tmp
  call mpi_reduce(epspv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epspv=tmp


  if (mype.eq.0) then ! prep for output

    if (cor.gt.1.e-8) then
      rossby = sqrt(2.*rossby)
    else
    rossby = - 999.
    endif
    if (bf.gt.1.e-8) then
      fr_z = sqrt(fr_z)
      fr_h = sqrt(2.*vort2)/bf
    else
      fr_z = - 999.
      fr_h = - 999.
    endif
    eps    = epsk + epsp
    if (aj.ne.0. .and. bj.ne.0.) then
      e  = (pe + ke)
    else
      e = 99999999.
    endif
    write(79,5046) time, epsk, epsp, eps, epskh, epskv, epsph, epspv
    write(46,5045) time,ke,pe,e,eg,ea,rossby,fr_z,fr_h,v
    write( 6,5044) time,ke,pe,e,eg,ea,rossby,fr_z,fr_h,v
    call flush(46)
    call flush(79)
  endif
  return

5042  format(1x,a91)
5043  format(7x,' T',9x,'KE',8x,'PE',8x,'E',9x,'GE',8x,'AE',8x,'Ro',8x,&
             'Fr_z',7x,'Fr_h',5x,'V')
5044  format(1x,f12.2,2x,9(f8.3,2x))
5045  format(1x,f12.2,2x,15(e11.4,1x))
5046  format(1x,f12.2,2x,7(e21.14,1x))

end subroutine out


subroutine energy_full(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,ke,pe,eg,ea)
! Computes total energy (KE, PE, VE, WE) 

  use param
  implicit none 
  include 'mpif.h'

  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,ge,g1,g2
  real,    intent(out) :: ke,pe,eg,ea

  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,wk,wkh2,wkh2n,kzn
  real :: vh,vzx,vzy,vzz,vz
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,tmp
  external proj, velo

  call velo(zx,zy,zz,ux,uy,uz)
  call wtoab(zx,zy,zz,tt,ge,g1,g2,ux,uy,uz)

  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea  = 0.

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wkh2n = wkh2*(L1/twopi)**2
           wk = kx*kx + ky*ky + kz*kz

           if (L(ikx,iky,ikz).eq.1) then
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)
              
              vzx = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              ke = ke + vzx/wk
              vzy = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              ke = ke + vzy/wk
              vzz = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              ke = ke + vzz/wk
              vh = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              pe = pe  + vh
            
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
            
              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv+vzx+vzy
                 zero_kh_geo = zero_kh_geo+VH*aj/bj
              endif
            
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo+VZX+VZY
                 zero_kz_grv = zero_kz_grv+VZZ+VH*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
                 eg = eg  + vz/wkh2
                 vz = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
                 vz = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
              endif

           endif
        enddo
     enddo
  enddo

  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj
  eg = eg + zero_kz_geo + zero_kh_geo
  ea = ea + zero_kz_grv + zero_kh_grv

  call mpi_allreduce(ke,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  ke=tmp
  call mpi_allreduce(pe,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  pe=tmp
  call mpi_allreduce(eg,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  eg=tmp
  call mpi_allreduce(ea,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  ea=tmp

end subroutine energy_full


subroutine cfl(zx,zy,zz,ux,uy,uz,ur,vr,wr)
! Computes and prints cfl number based on max velo

  use param
  use param_fftw
  implicit none 

  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz
  real,    intent(inout), dimension(n1d,n3d,n2dp) :: ur,vr,wr
  
  real :: umax,wmax,tmp,dx,dz

  dx = L1/float(n1)*1.5
  dz = L3/float(n3)*1.5

  call velo(zx,zy,zz,ux,uy,uz)
  call fftwkr(plan3_uk_ur,ux,ur)
  call fftwkr(plan3_vk_vr,uy,vr)
  call fftwkr(plan3_wk_wr,uz,wr)

  umax = maxval(abs(ur))
  umax = max(umax,maxval(abs(vr)))
  wmax = maxval(abs(wr))

  call mpi_reduce(umax,tmp,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,istatus)
  umax=tmp
  call mpi_reduce(wmax,tmp,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,istatus)
  wmax=tmp

  umax = umax*delt/dx
  wmax = wmax*delt/dz
  
  if (mype.eq.0) write(6,5555) umax, wmax
     
5555 format('Max CFL (x,z): ',2(f12.6,1x))

end subroutine cfl



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! NORMAL MODE CONVERSION
! See Bartello, 1995, J. Atmos. Sci. 52, pp. 4410-4428 for details
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine wtoab(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w)

! Converts from (vorticity, buoyancy) to (geok,grav.wave_1k,grav.wave_2k)

  use param
                  
  implicit none 
  complex, intent(in),  dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: u,v,w
  complex, intent(out), dimension(iktx,ikty,iktzp) :: geok,gw1k,gw2k

  integer :: ikx,iky,ikz,ikza,ikz0
  real :: omega,norm
  real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
  complex :: div,zk,tk,dk,bk
   
  geok = cmplx(0.,0.)
  gw1k = cmplx(0.,0.)
  gw2k = cmplx(0.,0.)

  call velo(zx,zy,zz,u,v,w)

  if (mype.eq.0) then 
    ikz0 = 2
  else
  ikz0 = 1
  endif

  do iky=1,ikty
    ky = kya(iky)
    do ikx=1,iktx
      kx = kxa(ikx)
      wkh2 = kx*kx + ky*ky
      wkh  = sqrt(wkh2)
      wkhn = wkh*L1/twopi
      if (wkhn.lt.1.e-10) then
        geok(ikx,iky,ikz0:iktzp) = aj*tt(ikx,iky,ikz0:iktzp)
        gw1k(ikx,iky,ikz0:iktzp) = u(ikx,iky,ikz0:iktzp)-zi*v(ikx,iky,ikz0:iktzp)
        gw2k(ikx,iky,ikz0:iktzp) = u(ikx,iky,ikz0:iktzp)+zi*v(ikx,iky,ikz0:iktzp)
      else
        do ikz=ikz0,iktzp
          if (L(ikx,iky,ikz).eq.1) then
            ikza  = mype*iktzp+ikz
            kz    = kza(ikza)
            wk2   = wkh2 + kz*kz
            wk    = sqrt(wk2)                
            wk    = max(wk,1.e-15)
            omega = sqrt(cor2*kz*kz + bf2*wkh2)/wk 
            div   = zi*(kx*u(ikx,iky,ikz)+ky*v(ikx,iky,ikz))
            bk    = aj*tt(ikx,iky,ikz)
            
            zk = zz(ikx,iky,ikz)
            dk = (wk/kz) * div
            tk = (wkh/bf) * bk
            
            norm  = omega*wk
            geok(ikx,iky,ikz) = (bf*wkh*zk + zi*cor*kz*tk)/norm
            
            norm  = sqrt2*omega*wk
            gw1k(ikx,iky,ikz) = (-zi*cor*kz*zk + omega*wk*dk - bf*wkh*tk)/norm
            
            norm  = sqrt2*omega*wk
            gw2k(ikx,iky,ikz) = (+zi*cor*kz*zk + omega*wk*dk + bf*wkh*tk)/norm
          endif
        enddo
      endif
    enddo
  enddo

  if (mype.eq.0) then  ! do kz=0 modes
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        if (L(ikx,iky,1).eq.1) then
          bk=aj*tt(ikx,iky,1)
          geok(ikx,iky,1)=zz(ikx,iky,1)
          gw1k(ikx,iky,1)=w(ikx,iky,1)-zi*bk/bf
          gw2k(ikx,iky,1)=w(ikx,iky,1)+zi*bk/bf
        endif
      enddo
    enddo
  endif
   
  return
end subroutine wtoab
 

subroutine atowb(geok,gw1k,gw2k,zx,zy,zz,tt,u,v,w)

! Converts from (geok,grav.wave_1k,grav.wave_2k) to (zeta,d,t).

  use param

  implicit none 
  complex, intent(in), dimension(iktx,ikty,iktzp) :: geok,gw1k,gw2k
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: u,v,w
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt

  integer :: ikx,iky,ikz,ikz0,ikza
  real :: omega
  real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
  complex :: zk,dk,tk,bk,gn,psi,div
  complex :: zgk,z1k,z2k,dgk,d1k,d2k,tgk,t1k,t2k

  u  = cmplx(0.,0.)
  v  = cmplx(0.,0.)
  w  = cmplx(0.,0.)
  tt = cmplx(0.,0.)
  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)

  if (mype.eq.0) then 
    ikz0 = 2
  else
    ikz0 = 1
  endif
   
  do ikz=ikz0,iktzp
    ikza = mype*iktzp+ikz
    kz   = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        wkh2 = kx*kx + ky*ky
        wkh  = sqrt(wkh2)
        wkhn = wkh*L1/twopi

        if (L(ikx,iky,ikz).eq.1) then 
          wk2 = wkh2 + kz*kz
          wk  = sqrt(wk2)
          wk  = max(wk,1.e-15)
          omega = sqrt(cor2*kz*kz + bf2*wkh2)/wk
          
          gn  = geok(ikx,iky,ikz) / (omega*wk)
          zgk = bf * wkh       * gn
          dgk = cmplx(0.,0.)
          tgk = - zi * cor * kz * gn
          
          gn  = gw1k(ikx,iky,ikz) / (sqrt2*omega*wk)
          z1k = + zi  * cor * kz * gn
          d1k = omega * wk     * gn
          t1k = - bf   * wkh    * gn
          gn  = gw2k(ikx,iky,ikz) / (sqrt2*omega*wk)
          z2k = - zi  * cor * kz * gn
          d2k = omega * wk     * gn
          t2k = + bf   * wkh   * gn
          
          zk = zgk + z1k + z2k
          dk = dgk + d1k + d2k
          tk = tgk + t1k + t2k
               
          if (wkhn.gt.1.e-10) then
            div             = dk*kz/wk
            u(ikx,iky,ikz)  = +zi*(ky*zk-kx*div)/wkh2             
            v(ikx,iky,ikz)  = -zi*(kx*zk+ky*div)/wkh2             
            w(ikx,iky,ikz)  = zi*div/kz
            tt(ikx,iky,ikz) = bf*tk/(aj*wkh)
          else
            u(ikx,iky,ikz) =     0.5*(gw1k(ikx,iky,ikz)+gw2k(ikx,iky,ikz))
            v(ikx,iky,ikz) = -zi*0.5*(gw1k(ikx,iky,ikz)-gw2k(ikx,iky,ikz))
            w(ikx,iky,ikz) = cmplx(0.,0.)
            tt(ikx,iky,ikz) = geok(ikx,iky,ikz)/aj
          endif
        endif
      enddo
    enddo
  enddo

  if (mype.eq.0) then  ! do kz=0 modes
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        if (L(ikx,iky,1).eq.1) then
          psi = -geok(ikx,iky,1)/(kx*kx+ky*ky)
          u(ikx,iky,1) = -zi*ky*psi
          v(ikx,iky,1) = +zi*kx*psi
          w(ikx,iky,1) = 0.5*(gw1k(ikx,iky,1)+gw2k(ikx,iky,1))
          bk = zi*bf*0.5*(gw1k(ikx,iky,1)-gw2k(ikx,iky,1))
          tt(ikx,iky,1) = bk/aj
        endif
      enddo
    enddo
  endif
     
  call vort(u,v,w,zx,zy,zz)
  
  return
end subroutine atowb


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! VELOCITY <--> VORTICITY
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine velo(zx,zy,zz,u,v,w)

 ! Calculates k-space velocity from k-space vorticity.
 ! curl (vorticity) = - laplacian (velocity) if velocity is solenoidal.

  use param

  implicit none
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz
  complex, intent(out), dimension(iktx,ikty,iktzp) :: u,v,w
  
  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,k2
  complex :: c1,c2,c3

  u = cmplx(0.,0.)
  v = cmplx(0.,0.)
  w = cmplx(0.,0.)

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        if (L(ikx,iky,ikz).eq.1) then
          kx = kxa(IKX)
          k2 = max(kx*kx+ky*ky+kz*kz,1.e-15)
          c1 = ky*zz(ikx,iky,ikz) - kz*zy(ikx,iky,ikz)               
          c2 = kz*zx(ikx,iky,ikz) - kx*zz(ikx,iky,ikz)
          c3 = kx*zy(ikx,iky,ikz) - ky*zx(ikx,iky,ikz)               
          u(ikx,iky,ikz) = zi * c1 / k2
          v(ikx,iky,ikz) = zi * c2 / k2
          w(ikx,iky,ikz) = zi * c3 / k2
        endif
      enddo
    enddo
  enddo

  return
end subroutine velo


subroutine vort(u,v,w,zx,zy,zz)

! Calculates k-space vortcity from k-space velocity.

  use param

  implicit none
  complex, intent(in), dimension(iktx,ikty,iktzp) :: u,v,w
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz

  integer :: ikx,iky,ikz,ikza
  real ::  kx,ky,kz
  complex :: c1,c2,c3

  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        if (L(ikx,iky,ikz).eq.1) then
          c1 = ky*w(ikx,iky,ikz) - kz*v(ikx,iky,ikz)               
          c2 = kz*u(ikx,iky,ikz) - kx*w(ikx,iky,ikz)
          c3 = kx*v(ikx,iky,ikz) - ky*u(ikx,iky,ikz)               
          zx(ikx,iky,ikz) = zi * c1
          zy(ikx,iky,ikz) = zi * c2
          zz(ikx,iky,ikz) = zi * c3
        endif
      enddo
    enddo
  enddo
  
  return
end subroutine vort


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c NONLINEAR TERMS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine convol(zxk,zyk,zzk,ttk,nxk,nyk,nzk,ntk,uk,vk,wk,ur,vr,wr, &
                  wxk,wyk,wzk,wtk,zxr,zyr,zzr,ttr,nxr,nyr,nzr,ntr)

! Nonlinear term in vorticity/temperature equation
! Uses scratch arrays wxk,wyk,wzk,wtk to avoid some ffts.
! Calculates convolution sums, calls ffts, etc.

  use param
  use param_fftw

  implicit none 
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,ttr,ur,vr,wr
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: nxr,nyr,nzr,ntr
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,uk,vk,wk,wxk,wyk,wzk,wtk
  complex, intent(out), dimension(iktx,ikty,iktzp)   :: nxk,nyk,nzk,ntk
  integer :: ikx,iky,ikz,ikza
  real    :: kx,ky,kz
  complex :: c1,c2,c3
  external :: velo

  ntk = cmplx(0.,0.)
  wxk = zxk
  wyk = zyk
  wzk = zzk
  wtk = ttk

! Nonlinear term in the temperature equation.

  call velo(zxk,zyk,zzk,uk,vk,wk)
  call fftwkr(plan3_uk_ur,uk,ur)
  call fftwkr(plan3_vk_vr,vk,vr)
  call fftwkr(plan3_wk_wr,wk,wr)
  call fftwkr(plan3_ttnk_ttnr,ttk,ttr)

  nxr = ur*ttr
  nyr = vr*ttr
  nzr = wr*ttr

  call fftwrk(plan3_nzxr_nzxk,nxr,nxk) 
  call fftwrk(plan3_nzyr_nzyk,nyr,nyk) 
  call fftwrk(plan3_nzzr_nzzk,nzr,nzk) 
  
  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        ntk(ikx,iky,ikz) = - zi * ( kx*nxk(ikx,iky,ikz) &
                                  + ky*nyk(ikx,iky,ikz) &
                                  + kz*nzk(ikx,iky,ikz) )*L(ikx,iky,ikz)
      enddo
    enddo
  enddo

! Nonlinear term in the vorticity equation.
  
  call fftwkr(plan3_zxnk_zxnr,zxk,zxr)
  call fftwkr(plan3_zynk_zynr,zyk,zyr)
  call fftwkr(plan3_zznk_zznr,zzk,zzr)

  nxr = wr*zyr - vr*zzr
  nyr = ur*zzr - wr*zxr
  nzr = vr*zxr - ur*zyr

  call fftwrk(plan3_nzxr_nzxk,nxr,nxk) 
  call fftwrk(plan3_nzyr_nzyk,nyr,nyk) 
  call fftwrk(plan3_nzzr_nzzk,nzr,nzk) 

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        c1 = ky*nzk(ikx,iky,ikz) - kz*nyk(ikx,iky,ikz)               
        c2 = kz*nxk(ikx,iky,ikz) - kx*nzk(ikx,iky,ikz)
        c3 = kx*nyk(ikx,iky,ikz) - ky*nxk(ikx,iky,ikz)               
        nxk(ikx,iky,ikz) = - zi*c1*L(ikx,iky,ikz)
        nyk(ikx,iky,ikz) = - zi*c2*L(ikx,iky,ikz)
        nzk(ikx,iky,ikz) = - zi*c3*L(ikx,iky,ikz)
        enddo
    enddo
  enddo

  zxk = wxk
  zyk = wyk
  zzk = wzk
  ttk = wtk
  
  return
end subroutine convol


subroutine constr(zxk,zyk,zzk,ttk,nxk,nyk,nzk,ntk,uk,vk,wk,ur,vr,wr, &
                  zxr,zyr,zzr,ttr,nxr,nyr,nzr,ntr)

! Nonlinear term in vorticity/temperature equation
! Like CONVOL but does not use scratch arrays, so requires a few extra FFTs.
! Calculates convolution sums, calls ffts, etc. without scratch arrays.

  use param
  use param_fftw

  implicit none 
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,ttr,ur,vr,wr
  real,    intent(out), dimension(n1d,n3d,n2dp)      :: nxr,nyr,nzr,ntr
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,uk,vk,wk
  complex, intent(out), dimension(iktx,ikty,iktzp)   :: nxk,nyk,nzk,ntk  

  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz
  complex :: c1,c2,c3

  external :: velo
  
  ntr = cmplx(0.,0.)

! Nonlinear term in temperature equation.

  call velo(zxk,zyk,zzk,uk,vk,wk)  
  call fftwkr(plan3_uk_ur,uk,ur)
  call fftwkr(plan3_vk_vr,vk,vr)
  call fftwkr(plan3_wk_wr,wk,wr)
  call fftwkr(plan3_ttnk_ttnr,ttk,ttr)

  ur = ur*ttr
  vr = vr*ttr
  wr = wr*ttr

! these calls not in convol:  
  call fftwrk(plan3_ur_uk,ur,uk) 
  call fftwrk(plan3_vr_vk,vr,vk) 
  call fftwrk(plan3_wr_wk,wr,wk) 
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk) 

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        ntk(ikx,iky,ikz) = - ZI * ( kx*uk(ikx,iky,ikz)  &
                                  + ky*vk(ikx,iky,ikz)  &
                                  + kz*wk(ikx,iky,ikz)  ) * L(ikx,iky,ikz)
      enddo
    enddo
  enddo


! Nonlinear term in vorticity equation.

  call velo(zxk,zyk,zzk,uk,vk,wk)  
  call fftwkr(plan3_uk_ur,uk,ur)
  call fftwkr(plan3_vk_vr,vk,vr)
  call fftwkr(plan3_wk_wr,wk,wr)
  call fftwkr(plan3_zxnk_zxnr,zxk,zxr)
  call fftwkr(plan3_zynk_zynr,zyk,zyr)
  call fftwkr(plan3_zznk_zznr,zzk,zzr)

  nxr = wr*zyr - vr*zzr
  nyr = ur*zzr - wr*zxr
  nzr = vr*zxr - ur*zyr
  
  call fftwrk(plan3_nzxr_nzxk,nxr,nxk) 
  call fftwrk(plan3_nzyr_nzyk,nyr,nyk) 
  call fftwrk(plan3_nzzr_nzzk,nzr,nzk) 
  call fftwrk(plan3_zxnr_zxnk,zxr,zxk) 
  call fftwrk(plan3_zynr_zynk,zyr,zyk) 
  call fftwrk(plan3_zznr_zznk,zzr,zzk) 

  do ikz=1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
        kx = kxa(ikx)
        c1 = ky*nzk(ikx,iky,ikz) - kz*nyk(ikx,iky,ikz)               
        c2 = kz*nxk(ikx,iky,ikz) - kx*nzk(ikx,iky,ikz)
        c3 = kx*nyk(ikx,iky,ikz) - ky*nxk(ikx,iky,ikz)               
        nxk(ikx,iky,ikz) = - zi * c1 * L(ikx,iky,ikz)
        nyk(ikx,iky,ikz) = - zi * c2 * L(ikx,iky,ikz)
        nzk(ikx,iky,ikz) = - zi * c3 * L(ikx,iky,ikz)
      enddo
    enddo
  enddo

  return
end subroutine constr


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c FFTS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fftwrk(plan,zr,zk)

! wrapper for real -> complex FFT.  Transforms zr -> zk
! Note, transforms are in place.  zr and zk must be equivalenced.
! So transform destroys zr

  use param
  use param_fftw

  implicit none
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: zr
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zk  
  type(C_PTR), intent(in) :: plan
       
  call fftwf_execute_dft_r2c(plan,zr,zk) 

  ! Normalize
  zk=zk/fftnorm
end subroutine fftwrk


subroutine fftwkr(plan,zk,zr)

! wrapper for complex -> real IFFT.  Transforms zk -> zr
! Note, transforms are in place.  zr and zk must be equivalenced.
! So transform destroys zk

  use param
  use param_fftw

  implicit none
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zk
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: zr
  type(C_PTR),intent(in) :: plan

  call realit(zk) !<---- IMPORTANT! Necessary to populate some zeroed-out wavenumbers before transforming.
  call fftwf_execute_dft_c2r(plan,zk,zr)

end subroutine fftwkr


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RANDOM NUMBERS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function rang(i)

  ! If i ne 0 then initializes RANNO with seed = i
  ! If i eq 0 then draws a random GAUSSIAN number with 
  ! mean and std = 1

  implicit none
  real :: rang
  integer :: i 
  real :: v1,v2,R,FAC,twopi,ranno
  external :: ranno

  twopi = 4.*asin(1.)
  
  if (i.ne.0) then
    v1 = ranno(i)
  else
200 v1 = 2.*(ranno(0)+twopi/2.)/twopi -1.
    v2 = 2.*(ranno(0)+twopi/2.)/twopi -1.
    r = v1**2. + v2**2.
    if (r.gt.1.) goto 200
    fac = sqrt( -2.*log(r)/r)
    rang = v1*fac
  endif
  return
end function rang


function ranno (i)

! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  implicit none
  real :: ranno
  integer :: i
  integer :: junk,ihold
  real :: twopi,ran1
  save junk

  twopi = 4.*asin(1.)
  if (i.ne.0) then
    if (i.gt.0) i = - i
    junk  = i
    ranno = (ran1(i)-0.5)*twopi
  else
    junk  = junk - 1
    ihold = junk
    ranno = (ran1(ihold)-0.5)*twopi
  endif
  return
end function ranno


function ran1(idum)
  implicit none
  real :: ran1
  integer :: idum

  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k,iv(32),iy
  save iv,iy
  data iv /ntab*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
    idum=max(-idum,1)
    do j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      if (j.le.ntab) iv(j)=idum
    enddo
    iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum.lt.0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
end function ran1



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c PHYSICAL SPACE DIAGNOSTICS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine dumpreal(zxk,zyk,zzk,ttk,zxr,zyr,zzr,ttr,uk,vk,wk,ur,vr,wr,ntdump,sl_rs)

  use param
  use param_fftw
  
  implicit none
  integer, intent(in) :: ntdump,sl_rs
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: ur,vr,wr,zxr,zyr,zzr,ttr
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,uk,vk,wk
  
  integer :: ikx,iky,ikz,ikza,i,j,k
  real :: kx,ky,kz
  external :: velo

  call velo(zxk,zyk,zzk,uk,vk,wk)

! Move vorticity and velocity to physical space.
  call fftwkr(plan3_zxnk_zxnr,zxk,zxr)
  call fftwkr(plan3_zynk_zynr,zyk,zyr)
  call fftwkr(plan3_zznk_zznr,zzk,zzr)
  call fftwkr(plan3_ttnk_ttnr,ttk,ttr)  
  call fftwkr(plan3_uk_ur,uk,ur)
  call fftwkr(plan3_vk_vr,vk,vr)
  call fftwkr(plan3_wk_wr,wk,wr)

  if(sl_rs.eq.FULLFIELD_MODE) call ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,ntdump)
  if(sl_rs.eq.SLICE_MODE)     call ncdumpslice(ur,vr,wr,zxr,zyr,zzr,ttr,ntdump)

! Put everything back to fourier for future timestepping  
  call fftwrk(plan3_zxnr_zxnk,zxr,zxk) 
  call fftwrk(plan3_zynr_zynk,zyr,zyk) 
  call fftwrk(plan3_zznr_zznk,zzr,zzk) 
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk) 

  uk = cmplx(0.,0.)
  vk = cmplx(0.,0.)
  wk = cmplx(0.,0.)

  return
end subroutine dumpreal


subroutine io_prep()

  use netcdf
  use param_netcdf
  use param
  implicit none
  include 'mpif.h'

  integer :: i,idx,idy,idz,idt,ncdims(3)          ! physical space dimensions
  integer :: idkx,idky,idkz,idkri,idkt,ncdimsk(5) ! Fourier space dimensions
  character(len=20) :: fname

  allocate(idnetcdf(7,6,nrsp))
  if (nrsp.gt.0 .and. rsflag.ne.0) then    
    if (mype.eq.0) print*,'Creating netcdf files for real space output'     
    do i=1,nrsp,1
      if (iputu.eq.1)  call ncf_file_gen( "U",i,idnetcdf(1,1,i),idnetcdf(1,2,i),"X",int(n1),idnetcdf(1,3,i),"Z",int(n3),idnetcdf(1,5,i),"Y",int(n2),idnetcdf(1,4,i),1,idnetcdf(1,6,i))
      if (iputv.eq.1)  call ncf_file_gen( "V",i,idnetcdf(2,1,i),idnetcdf(2,2,i),"X",int(n1),idnetcdf(2,3,i),"Z",int(n3),idnetcdf(2,5,i),"Y",int(n2),idnetcdf(2,4,i),1,idnetcdf(2,6,i))
      if (iputw.eq.1)  call ncf_file_gen( "W",i,idnetcdf(3,1,i),idnetcdf(3,2,i),"X",int(n1),idnetcdf(3,3,i),"Z",int(n3),idnetcdf(3,5,i),"Y",int(n2),idnetcdf(3,4,i),1,idnetcdf(3,6,i))
      if (iputth.eq.1) call ncf_file_gen("TH",i,idnetcdf(4,1,i),idnetcdf(4,2,i),"X",int(n1),idnetcdf(4,3,i),"Z",int(n3),idnetcdf(4,5,i),"Y",int(n2),idnetcdf(4,4,i),1,idnetcdf(4,6,i))
      if (iputzx.eq.1) call ncf_file_gen("ZX",i,idnetcdf(5,1,i),idnetcdf(5,2,i),"X",int(n1),idnetcdf(5,3,i),"Z",int(n3),idnetcdf(5,5,i),"Y",int(n2),idnetcdf(5,4,i),1,idnetcdf(5,6,i))
      if (iputzy.eq.1) call ncf_file_gen("ZY",i,idnetcdf(6,1,i),idnetcdf(6,2,i),"X",int(n1),idnetcdf(6,3,i),"Z",int(n3),idnetcdf(6,5,i),"Y",int(n2),idnetcdf(6,4,i),1,idnetcdf(6,6,i))
      if (iputzz.eq.1) call ncf_file_gen("ZZ",i,idnetcdf(7,1,i),idnetcdf(7,2,i),"X",int(n1),idnetcdf(7,3,i),"Z",int(n3),idnetcdf(7,5,i),"Y",int(n2),idnetcdf(7,4,i),1,idnetcdf(7,6,i))
    enddo ! nrsp loop    
    ncstart(1) = 1
    ncstart(2) = 1      
    nccount = (/n1,n3,n2p/)
  endif ! nrsp 
   
  allocate(idncslx(7,6,size(xslice)))    
  if(size(xslice).gt.0 .and. slxflag.ne.0) then    
    if (mype.eq.0) print*,'Creating netcdf files for slice outputs in X'      
    do i=1,size(xslice),1
      if (islputxu.eq.1)  call ncf_file_gen( 'SLX_U',i,idncslx(1,1,i),idncslx(1,2,i),"Y",int(n3),idncslx(1,4,i),"Z",int(n1),idncslx(1,3,i),"T",nsl,idncslx(1,5,i),nsl,idncslx(1,6,i))
      if (islputxv.eq.1)  call ncf_file_gen( 'SLX_V',i,idncslx(2,1,i),idncslx(2,2,i),"Y",int(n3),idncslx(2,4,i),"Z",int(n1),idncslx(2,3,i),"T",nsl,idncslx(2,5,i),nsl,idncslx(2,6,i))
      if (islputxw.eq.1)  call ncf_file_gen( 'SLX_W',i,idncslx(3,1,i),idncslx(3,2,i),"Y",int(n3),idncslx(3,4,i),"Z",int(n1),idncslx(3,3,i),"T",nsl,idncslx(3,5,i),nsl,idncslx(3,6,i))
      if (islputxth.eq.1) call ncf_file_gen('SLX_TH',i,idncslx(4,1,i),idncslx(4,2,i),"Y",int(n3),idncslx(4,4,i),"Z",int(n1),idncslx(4,3,i),"T",nsl,idncslx(4,5,i),nsl,idncslx(4,6,i))
      if (islputxzx.eq.1) call ncf_file_gen('SLX_ZX',i,idncslx(5,1,i),idncslx(5,2,i),"Y",int(n3),idncslx(5,4,i),"Z",int(n1),idncslx(5,3,i),"T",nsl,idncslx(5,5,i),nsl,idncslx(5,6,i))
      if (islputxzy.eq.1) call ncf_file_gen('SLX_ZY',i,idncslx(6,1,i),idncslx(6,2,i),"Y",int(n3),idncslx(6,4,i),"Z",int(n1),idncslx(6,3,i),"T",nsl,idncslx(6,5,i),nsl,idncslx(6,6,i))
      if (islputxzz.eq.1) call ncf_file_gen('SLX_ZZ',i,idncslx(7,1,i),idncslx(7,2,i),"Y",int(n3),idncslx(7,4,i),"Z",int(n1),idncslx(7,3,i),"T",nsl,idncslx(7,5,i),nsl,idncslx(7,6,i))
    enddo 
  endif ! size(xslice) slxflag

  allocate(idncsly(7,6,size(yslice)))    
  if(size(yslice).gt.0 .and. slyflag.ne.0) then    
    if (mype.eq.0) print*,'Creating netcdf files for slice outputs in Y'      
    do i=1,size(yslice),1
      if (islputyu.eq.1)  call ncf_file_gen( 'SLY_U',i,idncsly(1,1,i),idncsly(1,2,i),"Z",int(n3),idncsly(1,4,i),"X",int(n1),idncsly(1,3,i),"T",nsl,idncsly(1,5,i),nsl,idncsly(1,6,i))
      if (islputyv.eq.1)  call ncf_file_gen( 'SLY_V',i,idncsly(2,1,i),idncsly(2,2,i),"Z",int(n3),idncsly(2,4,i),"X",int(n1),idncsly(2,3,i),"T",nsl,idncsly(2,5,i),nsl,idncsly(2,6,i))
      if (islputyw.eq.1)  call ncf_file_gen( 'SLY_W',i,idncsly(3,1,i),idncsly(3,2,i),"Z",int(n3),idncsly(3,4,i),"X",int(n1),idncsly(3,3,i),"T",nsl,idncsly(3,5,i),nsl,idncsly(3,6,i))
      if (islputyth.eq.1) call ncf_file_gen('SLY_TH',i,idncsly(4,1,i),idncsly(4,2,i),"Z",int(n3),idncsly(4,4,i),"X",int(n1),idncsly(4,3,i),"T",nsl,idncsly(4,5,i),nsl,idncsly(4,6,i))
      if (islputyzx.eq.1) call ncf_file_gen('SLY_ZX',i,idncsly(5,1,i),idncsly(5,2,i),"Z",int(n3),idncsly(5,4,i),"X",int(n1),idncsly(5,3,i),"T",nsl,idncsly(5,5,i),nsl,idncsly(5,6,i))
      if (islputyzy.eq.1) call ncf_file_gen('SLY_ZY',i,idncsly(6,1,i),idncsly(6,2,i),"Z",int(n3),idncsly(6,4,i),"X",int(n1),idncsly(6,3,i),"T",nsl,idncsly(6,5,i),nsl,idncsly(6,6,i))
      if (islputyzz.eq.1) call ncf_file_gen('SLY_ZZ',i,idncsly(7,1,i),idncsly(7,2,i),"Z",int(n3),idncsly(7,4,i),"X",int(n1),idncsly(7,3,i),"T",nsl,idncsly(7,5,i),nsl,idncsly(7,6,i))
    enddo 
  endif ! size(yslice) slyflag

  allocate(idncslz(7,6,size(zslice)))    
  if(size(zslice).gt.0 .and. slzflag.ne.0) then    
    if (mype.eq.0) print*,'Creating netcdf files for slice outputs in Z'
    do i=1,size(zslice),1
      if (islputxu.eq.1)  call ncf_file_gen( 'SLZ_U',i,idncslz(1,1,i),idncslz(1,2,i),"Y",int(n3),idncslz(1,4,i),"X",int(n1),idncslz(1,3,i),"T",nsl,idncslz(1,5,i),nsl,idncslz(1,6,i))
      if (islputxv.eq.1)  call ncf_file_gen( 'SLZ_V',i,idncslz(2,1,i),idncslz(2,2,i),"Y",int(n3),idncslz(2,4,i),"X",int(n1),idncslz(2,3,i),"T",nsl,idncslz(2,5,i),nsl,idncslz(2,6,i))
      if (islputxw.eq.1)  call ncf_file_gen( 'SLZ_W',i,idncslz(3,1,i),idncslz(3,2,i),"Y",int(n3),idncslz(3,4,i),"X",int(n1),idncslz(3,3,i),"T",nsl,idncslz(3,5,i),nsl,idncslz(3,6,i))
      if (islputxth.eq.1) call ncf_file_gen('SLZ_TH',i,idncslz(4,1,i),idncslz(4,2,i),"Y",int(n3),idncslz(4,4,i),"X",int(n1),idncslz(4,3,i),"T",nsl,idncslz(4,5,i),nsl,idncslz(4,6,i))
      if (islputxzx.eq.1) call ncf_file_gen('SLZ_ZX',i,idncslz(5,1,i),idncslz(5,2,i),"Y",int(n3),idncslz(5,4,i),"X",int(n1),idncslz(5,3,i),"T",nsl,idncslz(5,5,i),nsl,idncslz(5,6,i))
      if (islputxzy.eq.1) call ncf_file_gen('SLZ_ZY',i,idncslz(6,1,i),idncslz(6,2,i),"Y",int(n3),idncslz(6,4,i),"X",int(n1),idncslz(6,3,i),"T",nsl,idncslz(6,5,i),nsl,idncslz(6,6,i))
      if (islputxzz.eq.1) call ncf_file_gen('SLZ_ZZ',i,idncslz(7,1,i),idncslz(7,2,i),"Y",int(n3),idncslz(7,4,i),"X",int(n1),idncslz(7,3,i),"T",nsl,idncslz(7,5,i),nsl,idncslz(7,6,i))
    enddo 
  endif ! size(zslice) slzflag
  
  !!! next, prep restart file
  allocate(idncrst(4,7,nrst))  
  if(rstflag.eq.1) then    
    if (mype.eq.0) print*,'Creating netcdf restart files'
    do i=1,nrst,1
      call ncf_rst_file_gen('RST_ZXK',i,idncrst(1,1,i),idncrst(1,2,i),"KX",int(iktx),idncrst(1,3,i),"KY",int(ikty),idncrst(1,4,i),"KZ",int(iktz),idncrst(1,5,i),"RI",2,idncrst(1,6,i),idncrst(1,7,i))
      call ncf_rst_file_gen('RST_ZYK',i,idncrst(2,1,i),idncrst(2,2,i),"KX",int(iktx),idncrst(2,3,i),"KY",int(ikty),idncrst(2,4,i),"KZ",int(iktz),idncrst(2,5,i),"RI",2,idncrst(2,6,i),idncrst(2,7,i))
      call ncf_rst_file_gen('RST_ZZK',i,idncrst(3,1,i),idncrst(3,2,i),"KX",int(iktx),idncrst(3,3,i),"KY",int(ikty),idncrst(3,4,i),"KZ",int(iktz),idncrst(3,5,i),"RI",2,idncrst(3,6,i),idncrst(3,7,i))
      call ncf_rst_file_gen('RST_TTK',i,idncrst(4,1,i),idncrst(4,2,i),"KX",int(iktx),idncrst(4,3,i),"KY",int(ikty),idncrst(4,4,i),"KZ",int(iktz),idncrst(4,5,i),"RI",2,idncrst(4,6,i),idncrst(4,7,i))
    enddo      
    ncstartk = (/1,1,1,1/)
    nccountk = (/int(iktx),int(ikty),int(iktzp),1/)
  end if ! rstflag
end subroutine io_prep


subroutine ncf_file_gen(field_name,index,idnc,idvar,varx,lenx,idx,vary,leny,idy,varz,lenz,idz,lent,idtimes) 

! Takes field_name and index, and creates .ncf file, and returns file... 
! ncf index idnc, variable index idvar, and dimension indicies idx,idy,idz 
! eg. takes ("U",1,...), creates U.01.ncf, and sets indicies 

  use netcdf
  use param
  implicit none
  include 'mpif.h'

  character(len=*),intent(in) :: field_name,varx,vary,varz  
  integer,intent(in) :: index,lenx,leny,lenz,lent
  integer,intent(inout) :: idnc,idvar,idx,idy,idz,idtimes
  integer :: ncdims(3),idt
  character(len=20) :: fname
  call filename_gen(field_name,index,fname)  
  istatus = nf90_create(adjustl(trim(fname)),IOR(NF90_NETCDF4,NF90_MPIIO),idnc,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
  if (istatus.ne.0) print*,'error in nf_create '//fname

  istatus = nf90_def_dim(idnc,'Time',lent,idt)
  if (istatus.ne.0) print*,'error in nf_def_dim T in '//field_name
  istatus = nf90_def_var(idnc,'TIMES',NF90_FLOAT,idt,idtimes)
  if (istatus.ne.0) print*,'error nf_def_var TIMES in '//field_name    
  
  istatus = nf90_def_dim(idnc,varx,lenx,idx)
  if (istatus.ne.0) print*,'error in nf_def_dim x for '//adjustl(trim(fname))
  istatus = nf90_def_dim(idnc,vary,leny,idy)
  if (istatus.ne.0) print*,'error in nf_def_dim y for '//adjustl(trim(fname))
  istatus = nf90_def_dim(idnc,varz,lenz,idz)
  if (istatus.ne.0) print*,'error in nf_def_dim z for '//adjustl(trim(fname))
  ncdims = (/idx,idy,idz/)
  istatus = nf90_def_var(idnc,field_name,NF90_FLOAT,ncdims,idvar)
  if (istatus.ne.0) print*,'error nf_def_var '//field_name
  istatus = nf90_enddef(idnc)
  if (istatus.ne.0) print*,'error enddef '//adjustl(trim(fname))
end subroutine ncf_file_gen


subroutine ncf_rst_file_gen(field_name,index,idnc,idvar,varx,lenx,idx,vary,leny,idy,varz,lenz,idz,varri,lenri,idri,idtimes) 

! Takes field_name and index, and creates .ncf file, and returns file...
! ncf index idnc, variable index idvar, and dimension indicies idkx,idky,idkz,idkr 
! eg. takes ("ZX",1,...), creates ZX.01.ncf, and sets indicies 

  use netcdf
  use param
  implicit none
  include 'mpif.h'

  character(len=*),intent(in) :: field_name,varx,vary,varz,varri  
  integer,intent(in) :: index,lenx,leny,lenz,lenri
  integer,intent(inout) :: idnc,idvar,idx,idy,idz,idri,idtimes
  integer :: ncdims(4),idt
  character(len=20) :: fname
  call filename_gen(field_name,index,fname)  
  istatus = nf90_create(adjustl(trim(fname)),IOR(NF90_NETCDF4,NF90_MPIIO),idnc,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
  if (istatus.ne.0) print*,'error in nf_create '//fname

  istatus = nf90_def_dim(idnc,'T',1,idt)
  if (istatus.ne.0) print*,'error in nf_def_dim T in rst file'
  istatus = nf90_def_var(idnc,'TIMES',NF90_FLOAT,idt,idtimes)
  if (istatus.ne.0) print*,'error nf_def_var TIMES in rst file'  
  
  istatus = nf90_def_dim(idnc,varx,lenx,idx)
  if (istatus.ne.0) print*,'error in nf_def_dim x for '//adjustl(trim(fname))
  istatus = nf90_def_dim(idnc,vary,leny,idy)
  if (istatus.ne.0) print*,'error in nf_def_dim y for '//adjustl(trim(fname))
  istatus = nf90_def_dim(idnc,varz,lenz,idz)
  if (istatus.ne.0) print*,'error in nf_def_dim z for '//adjustl(trim(fname))
  istatus = nf90_def_dim(idnc,varri,lenri,idri)
  if (istatus.ne.0) print*,'error in nf_def_dim ri for '//adjustl(trim(fname))
  ncdims = (/idx,idy,idz,idri/)
  istatus = nf90_def_var(idnc,field_name,NF90_FLOAT,ncdims,idvar)
  if (istatus.ne.0) print*,'error nf_def_var '//field_name

  istatus = nf90_enddef(idnc)
  if (istatus.ne.0) print*,'error enddef '//adjustl(trim(fname))
end subroutine ncf_rst_file_gen


subroutine filename_gen(uncf, tstep, fname)
! takes ("U",1,fname) and sets fname = 'U.01.ncf'

  implicit none
  character(len=20) :: fname, tstep_char
  character(len=*) :: uncf
  integer :: tstep
  write(tstep_char,"(I0.2)") tstep    
  fname = trim(trim(uncf)//"."//trim(tstep_char)//".ncf")    
end subroutine filename_gen


subroutine io_close()
  use netcdf
  use param_netcdf
  use param
  implicit none
  integer :: i
  
  do i=1,nrsp,1 ! close real space files 
    if(rsflag.ne.0) then
      if (iputu.eq.1)  istatus = nf90_close(idnetcdf(1,1,i))
      if (iputv.eq.1)  istatus = nf90_close(idnetcdf(2,1,i))
      if (iputw.eq.1)  istatus = nf90_close(idnetcdf(3,1,i))
      if (iputth.eq.1) istatus = nf90_close(idnetcdf(4,1,i))
      if (iputzx.eq.1) istatus = nf90_close(idnetcdf(5,1,i))
      if (iputzy.eq.1) istatus = nf90_close(idnetcdf(6,1,i))
      if (iputzz.eq.1) istatus = nf90_close(idnetcdf(7,1,i))
    endif
  enddo
  do i=1,size(xslice),1 ! close X slice files
    if(slxflag.ne.0) then
      if (islputxu.eq.1)  istatus = nf90_close(idncslx(1,1,i))
      if (islputxv.eq.1)  istatus = nf90_close(idncslx(2,1,i))
      if (islputxw.eq.1)  istatus = nf90_close(idncslx(3,1,i))
      if (islputxth.eq.1) istatus = nf90_close(idncslx(4,1,i))
      if (islputxzx.eq.1) istatus = nf90_close(idncslx(5,1,i))
      if (islputxzy.eq.1) istatus = nf90_close(idncslx(6,1,i))
      if (islputxzz.eq.1) istatus = nf90_close(idncslx(7,1,i))
    endif
  enddo
  do i=1,size(yslice),1 ! close Y slice files
    if(slyflag.ne.0) then
      if (islputyu.eq.1)  istatus = nf90_close(idncsly(1,1,i))
      if (islputyv.eq.1)  istatus = nf90_close(idncsly(2,1,i))
      if (islputyw.eq.1)  istatus = nf90_close(idncsly(3,1,i))
      if (islputyth.eq.1) istatus = nf90_close(idncsly(4,1,i))
      if (islputyzx.eq.1) istatus = nf90_close(idncsly(5,1,i))
      if (islputyzy.eq.1) istatus = nf90_close(idncsly(6,1,i))
      if (islputyzz.eq.1) istatus = nf90_close(idncsly(7,1,i))
    endif
  enddo
  do i=1,size(zslice),1 ! close Z slice files
    if(slzflag.ne.0) then
      if (islputzu.eq.1)  istatus = nf90_close(idncslz(1,1,i))
      if (islputzv.eq.1)  istatus = nf90_close(idncslz(2,1,i))
      if (islputzw.eq.1)  istatus = nf90_close(idncslz(3,1,i))
      if (islputzth.eq.1) istatus = nf90_close(idncslz(4,1,i))
      if (islputzzx.eq.1) istatus = nf90_close(idncslz(5,1,i))
      if (islputzzy.eq.1) istatus = nf90_close(idncslz(6,1,i))
      if (islputzzz.eq.1) istatus = nf90_close(idncslz(7,1,i))
    endif
  enddo
  do i=1,nrst,1 ! close restart files
    if (rstflag.eq.1) then
      istatus = nf90_close(idncrst(1,1,i))
      istatus = nf90_close(idncrst(2,1,i))
      istatus = nf90_close(idncrst(3,1,i))
      istatus = nf90_close(idncrst(4,1,i))
    end if
  enddo  
  deallocate(idnetcdf)
  deallocate(idncslx)
  deallocate(idncsly)
  deallocate(idncslz)
  deallocate(idncrst)  
end subroutine io_close


subroutine ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,ntdump)
! Dump real space fields into X.*.ncf files

  use netcdf
  use param_netcdf
  use param
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ntdump
  real,intent(in), dimension(n1d,n3d,n2dp) :: ur,vr,wr,zxr,zyr,zzr,ttr
  integer, dimension(1) :: nctimestart = (/1/)

  if (mype.eq.0) then
    print*,'Writing to *.ncf at time ',time   
    if (iputu.eq.1)  istatus = nf90_put_var(idnetcdf(1,1,ntdump),idnetcdf(1,6,ntdump),time,nctimestart)
    if (iputv.eq.1)  istatus = nf90_put_var(idnetcdf(2,1,ntdump),idnetcdf(2,6,ntdump),time,nctimestart)
    if (iputw.eq.1)  istatus = nf90_put_var(idnetcdf(3,1,ntdump),idnetcdf(3,6,ntdump),time,nctimestart)
    if (iputth.eq.1) istatus = nf90_put_var(idnetcdf(4,1,ntdump),idnetcdf(4,6,ntdump),time,nctimestart)
    if (iputzx.eq.1) istatus = nf90_put_var(idnetcdf(5,1,ntdump),idnetcdf(5,6,ntdump),time,nctimestart)
    if (iputzy.eq.1) istatus = nf90_put_var(idnetcdf(6,1,ntdump),idnetcdf(6,6,ntdump),time,nctimestart)
    if (iputzz.eq.1) istatus = nf90_put_var(idnetcdf(7,1,ntdump),idnetcdf(7,6,ntdump),time,nctimestart)
  endif
  ncstart(3) = mype*n2p+1
  if (iputu.eq.1)  istatus = nf90_put_var(idnetcdf(1,1,ntdump),idnetcdf(1,2,ntdump), ur(1:n1,:,:),ncstart,nccount)
  if (iputv.eq.1)  istatus = nf90_put_var(idnetcdf(2,1,ntdump),idnetcdf(2,2,ntdump), vr(1:n1,:,:),ncstart,nccount)
  if (iputw.eq.1)  istatus = nf90_put_var(idnetcdf(3,1,ntdump),idnetcdf(3,2,ntdump), wr(1:n1,:,:),ncstart,nccount)
  if (iputth.eq.1) istatus = nf90_put_var(idnetcdf(4,1,ntdump),idnetcdf(4,2,ntdump),ttr(1:n1,:,:),ncstart,nccount)
  if (iputzx.eq.1) istatus = nf90_put_var(idnetcdf(5,1,ntdump),idnetcdf(5,2,ntdump),zxr(1:n1,:,:),ncstart,nccount)
  if (iputzy.eq.1) istatus = nf90_put_var(idnetcdf(6,1,ntdump),idnetcdf(6,2,ntdump),zyr(1:n1,:,:),ncstart,nccount)
  if (iputzz.eq.1) istatus = nf90_put_var(idnetcdf(7,1,ntdump),idnetcdf(7,2,ntdump),zzr(1:n1,:,:),ncstart,nccount)
end subroutine ncdumprsp


subroutine ncdumpslice(ur,vr,wr,zxr,zyr,zzr,ttr,ntdump)
! Dump real field y-plane slices at planes specified in yslice
  use param
  use param_netcdf
  use netcdf
  implicit none
  include 'mpif.h'

  integer,intent(in) :: ntdump
! real, intent(in),dimension(n1d,n3d,n2dp),target :: ur,vr,wr,zxr,zyr,zzr,ttr    !!! took out "target" to compile on gfortran
  real, intent(in),dimension(n1d,n3d,n2dp) :: ur,vr,wr,zxr,zyr,zzr,ttr  
  integer,dimension(1) :: nctimestart
  integer :: i,j

  nctimestart(1) = ntdump  
  ncslstart = (/1,int(mype*n2dp+1),ntdump/)
  ncslcount = (/int(n3),int(n2dp),1/)
  do i=1,size(xslice),1    
    if (slxflag.eq.1) then 
      j=xslice(i)
      if (islputxu.eq.1)  istatus = nf90_put_var(idncslx(1,1,i),idncslx(1,2,i), ur(j,:,:),ncslstart,ncslcount)      
      if (islputxv.eq.1)  istatus = nf90_put_var(idncslx(2,1,i),idncslx(2,2,i), vr(j,:,:),ncslstart,ncslcount)      
      if (islputxw.eq.1)  istatus = nf90_put_var(idncslx(3,1,i),idncslx(3,2,i), wr(j,:,:),ncslstart,ncslcount)      
      if (islputxth.eq.1) istatus = nf90_put_var(idncslx(4,1,i),idncslx(4,2,i),ttr(j,:,:),ncslstart,ncslcount)      
      if (islputxzx.eq.1) istatus = nf90_put_var(idncslx(5,1,i),idncslx(5,2,i),zxr(j,:,:),ncslstart,ncslcount)      
      if (islputxzy.eq.1) istatus = nf90_put_var(idncslx(6,1,i),idncslx(6,2,i),zyr(j,:,:),ncslstart,ncslcount)      
      if (islputxzz.eq.1) istatus = nf90_put_var(idncslx(7,1,i),idncslx(7,2,i),zzr(j,:,:),ncslstart,ncslcount)      
      if (mype.eq.0) then 
        if(i.eq.1) print*,'Writing X slice *.ncf at time ',time        
        if (islputxu.eq.1)  istatus = nf90_put_var(idncslx(1,1,i),idncslx(1,6,i),time,nctimestart)
        if (islputxv.eq.1)  istatus = nf90_put_var(idncslx(2,1,i),idncslx(2,6,i),time,nctimestart)
        if (islputxw.eq.1)  istatus = nf90_put_var(idncslx(3,1,i),idncslx(3,6,i),time,nctimestart)
        if (islputxth.eq.1) istatus = nf90_put_var(idncslx(4,1,i),idncslx(4,6,i),time,nctimestart)
        if (islputxzx.eq.1) istatus = nf90_put_var(idncslx(5,1,i),idncslx(5,6,i),time,nctimestart)
        if (islputxzy.eq.1) istatus = nf90_put_var(idncslx(6,1,i),idncslx(6,6,i),time,nctimestart)
        if (islputxzz.eq.1) istatus = nf90_put_var(idncslx(7,1,i),idncslx(7,6,i),time,nctimestart)
      endif
    endif
  enddo
  
  ncslstart = (/1,1,ntdump/)
  ncslcount = (/int(n1),int(n3),1/)
  do i=1,size(yslice),1
    if (slyflag.eq.1 .and. int((yslice(i)-1)/n2dp).eq.mype) then
      j=yslice(i)-(mype*n2dp)
      if (islputyu.eq.1)  istatus = nf90_put_var(idncsly(1,1,i),idncsly(1,2,i), ur(1:n1,:,j),ncslstart,ncslcount) 
      if (islputyv.eq.1)  istatus = nf90_put_var(idncsly(2,1,i),idncsly(2,2,i), vr(1:n1,:,j),ncslstart,ncslcount)
      if (islputyw.eq.1)  istatus = nf90_put_var(idncsly(3,1,i),idncsly(3,2,i), wr(1:n1,:,j),ncslstart,ncslcount) 
      if (islputyth.eq.1) istatus = nf90_put_var(idncsly(4,1,i),idncsly(4,2,i),ttr(1:n1,:,j),ncslstart,ncslcount)
      if (islputyzx.eq.1) istatus = nf90_put_var(idncsly(5,1,i),idncsly(5,2,i),zxr(1:n1,:,j),ncslstart,ncslcount) 
      if (islputyzy.eq.1) istatus = nf90_put_var(idncsly(6,1,i),idncsly(6,2,i),zyr(1:n1,:,j),ncslstart,ncslcount) 
      if (islputyzz.eq.1) istatus = nf90_put_var(idncsly(7,1,i),idncsly(7,2,i),zzr(1:n1,:,j),ncslstart,ncslcount)       
      if (i.eq.1) print*,'Writing Y slice *.ncf at time ',time        
      if (islputyu.eq.1)  istatus = nf90_put_var(idncsly(1,1,i),idncsly(1,6,i),time,nctimestart)
      if (islputyv.eq.1)  istatus = nf90_put_var(idncsly(2,1,i),idncsly(2,6,i),time,nctimestart)
      if (islputyw.eq.1)  istatus = nf90_put_var(idncsly(3,1,i),idncsly(3,6,i),time,nctimestart)
      if (islputyth.eq.1) istatus = nf90_put_var(idncsly(4,1,i),idncsly(4,6,i),time,nctimestart)
      if (islputyzx.eq.1) istatus = nf90_put_var(idncsly(5,1,i),idncsly(5,6,i),time,nctimestart)
      if (islputyzy.eq.1) istatus = nf90_put_var(idncsly(6,1,i),idncsly(6,6,i),time,nctimestart)
      if (islputyzz.eq.1) istatus = nf90_put_var(idncsly(7,1,i),idncsly(7,6,i),time,nctimestart)      
    endif
  enddo

  ncslstart = (/1,int(mype*n2dp+1),ntdump/)
  ncslcount = (/int(n1),int(n2dp),1/)
  do i=1,size(zslice),1    
    if (slzflag.eq.1) then
      j=zslice(i)
      if (islputzu.eq.1)  istatus = nf90_put_var(idncslz(1,1,i),idncslz(1,2,i), ur(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzv.eq.1)  istatus = nf90_put_var(idncslz(2,1,i),idncslz(2,2,i), vr(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzw.eq.1)  istatus = nf90_put_var(idncslz(3,1,i),idncslz(3,2,i), wr(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzth.eq.1) istatus = nf90_put_var(idncslz(4,1,i),idncslz(4,2,i),ttr(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzzx.eq.1) istatus = nf90_put_var(idncslz(5,1,i),idncslz(5,2,i),zxr(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzzy.eq.1) istatus = nf90_put_var(idncslz(6,1,i),idncslz(6,2,i),zyr(1:n1,j,:),ncslstart,ncslcount)      
      if (islputzzz.eq.1) istatus = nf90_put_var(idncslz(7,1,i),idncslz(7,2,i),zzr(1:n1,j,:),ncslstart,ncslcount)      
      if (mype.eq.0) then 
        if (i.eq.1) print*,'Writing Z slice *.ncf at time ',time        
        if (islputzu.eq.1)  istatus = nf90_put_var(idncslz(1,1,i),idncslz(1,6,i),time,nctimestart)
        if (islputzv.eq.1)  istatus = nf90_put_var(idncslz(2,1,i),idncslz(2,6,i),time,nctimestart)
        if (islputzw.eq.1)  istatus = nf90_put_var(idncslz(3,1,i),idncslz(3,6,i),time,nctimestart)
        if (islputzth.eq.1) istatus = nf90_put_var(idncslz(4,1,i),idncslz(4,6,i),time,nctimestart)
        if (islputzzx.eq.1) istatus = nf90_put_var(idncslz(5,1,i),idncslz(5,6,i),time,nctimestart)
        if (islputzzy.eq.1) istatus = nf90_put_var(idncslz(6,1,i),idncslz(6,6,i),time,nctimestart)
        if (islputzzz.eq.1) istatus = nf90_put_var(idncslz(7,1,i),idncslz(7,6,i),time,nctimestart)
      endif
    endif    
  enddo
end subroutine ncdumpslice


subroutine ncdumprst(zx,zy,zz,tt,ntdump)
! dump restart into Zk.out.ncf
  use netcdf
  use param_netcdf
  use param
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ntdump
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  integer, dimension(1) :: nctimestart = (/1/)

  if (mype.eq.0) then
    print*,'Writing to restart files at time ',time
    istatus    = nf90_put_var(idncrst(1,1,ntdump),idncrst(1,7,ntdump),time,nctimestart)
    istatus    = nf90_put_var(idncrst(2,1,ntdump),idncrst(2,7,ntdump),time,nctimestart)
    istatus    = nf90_put_var(idncrst(3,1,ntdump),idncrst(3,7,ntdump),time,nctimestart)
    istatus    = nf90_put_var(idncrst(4,1,ntdump),idncrst(4,7,ntdump),time,nctimestart)
  endif
  ncstartk(3) = mype*iktzp+1
  ncstartk(4) = 1 ! (real part)
  istatus = nf90_put_var(idncrst(1,1,ntdump),idncrst(1,2,ntdump),real(zx),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(2,1,ntdump),idncrst(2,2,ntdump),real(zy),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(3,1,ntdump),idncrst(3,2,ntdump),real(zz),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(4,1,ntdump),idncrst(4,2,ntdump),real(tt),ncstartk,nccountk)
  ncstartk(4) = 2 ! (imag part)
  istatus = nf90_put_var(idncrst(1,1,ntdump),idncrst(1,2,ntdump),aimag(zx),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(2,1,ntdump),idncrst(2,2,ntdump),aimag(zy),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(3,1,ntdump),idncrst(3,2,ntdump),aimag(zz),ncstartk,nccountk)
  istatus = nf90_put_var(idncrst(4,1,ntdump),idncrst(4,2,ntdump),aimag(tt),ncstartk,nccountk)
end subroutine ncdumprst

! load restart files from ZXK.in.ncf,ZYK.in.ncf,ZXK.in.ncf,TTK.in.ncf
subroutine ncreadrst(zx,zy,zz,tt,wr,wi,ts)
  use param

! note: don't load param_netcdf.  This subroutine opens and reads from var.in.ncf
! The params in param_netcdf correspond to var.out.ncf

  implicit none  
  real,    intent(out) :: ts
  real :: ts1,ts2,ts3
  real,    intent(inout), dimension(iktx,ikty,iktzp) :: wr,wi
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
    
  if(mype.eq.0) print*,'Reading from restart files'
  call ncrstreadvar('RST_ZXK',ts,zx,wr,wi) 
  call ncrstreadvar('RST_ZYK',ts1,zy,wr,wi)  
  call ncrstreadvar('RST_ZZK',ts2,zz,wr,wi) 
  call ncrstreadvar('RST_TTK',ts3,tt,wr,wi)
  if(ts.ne.ts1 .or. ts.ne.ts2 .or. ts.ne.ts3) then
    print*,'Restart times ',ts,ts1,ts2,ts3
    print*,'All restart files must be from same time'
    stop
  endif
end subroutine ncreadrst

subroutine ncrstreadvar(var,ts,zx,wr,wi)    
  use param
  use netcdf
  implicit none
  include 'mpif.h'

  character(len=*),intent(in) :: var  
  integer :: idkx,idky,idkz,idri,iktx1,ikty1,iktz1,ncreststart(1),idtimesk,idnc,idvar
  integer, dimension(4) :: ncstartrk,ncstartik,nccountk
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx
  real,    intent(out) :: ts
  real,    intent(inout), dimension(iktx,ikty,iktzp) :: wr,wi
    
  ! open var.in.ncf

  istatus = nf90_open(var//'.in.ncf',IOR(NF90_NOWRITE, NF90_MPIIO),idnc,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
  if (istatus.ne.0) print*,var//'.in bad idnc'
  ! get dimension IDs
  istatus = nf90_inq_dimid(idnc,'KX',idkx)
  if (istatus.ne.0) print*,'Error getting kx ID'
  istatus = nf90_inq_dimid(idnc,'KY',idky)
  if (istatus.ne.0) print*,'Error getting ky ID'
  istatus = nf90_inq_dimid(idnc,'KZ',idkz)
  if (istatus.ne.0) print*,'Error getting kz ID'
  istatus = nf90_inq_dimid(idnc,'RI',idri)
  if (istatus.ne.0) print*,'Error getting real/imag ID'

  istatus = nf90_inquire_dimension(idnc,idkx,len=iktx1)  
  if (istatus.ne.0) print*,'Error getting idkx1'
  istatus = nf90_inquire_dimension(idnc,idky,len=ikty1)  
  if (istatus.ne.0) print*,'Error getting idky1'
  istatus = nf90_inquire_dimension(idnc,idkz,len=iktz1)  
  if (istatus.ne.0) print*,'Error getting idkz1'

  if (iktx1.ne.iktx .or. ikty1.ne.ikty .or. iktz1.ne.iktz) then
    print*,'Sorry, do not know how to change resolution.'
    stop
  endif

  istatus = nf90_inq_varid(idnc,'TIMES',idtimesk)
  if (istatus.ne.0) print*,'Error getting restart times ID '//var

  ncreststart(1) = 1
  istatus = nf90_get_var(idnc,idtimesk,ts,ncreststart)
  if (istatus.ne.0) print*,'Error reading restart time '//var

  istatus = nf90_inq_varid(idnc,var,idvar)
  if (istatus.ne.0) print*,'Error getting '//var//' ID'
  ! prep netcdf read
  ncstartrk = (/1,1,int(mype*iktzp+1),1/)
  ncstartik = (/1,1,int(mype*iktzp+1),2/)  
  nccountk = (/int(iktx),int(ikty),int(iktzp),1/)

  istatus = nf90_get_var(idnc,idvar,wr,ncstartrk,nccountk)
  if (istatus.ne.0) print*,'Error reading real '//var
  istatus = nf90_get_var(idnc,idvar,wi,ncstartik,nccountk)
  if (istatus.ne.0) print*,'Error reading imag '//var
  zx = wr + zi*wi
 
  istatus = nf90_close(idnc)
  if (istatus.ne.0) print*,'Error closing '//var//' .in.ncf'  
end subroutine ncrstreadvar



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MISC SUBROUTINES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine proj(zx,zy,zz)

! Fourier-space determination of the solenoidal part of a vector zx,y,z.

  use param

  implicit none
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zx,zy,zz

  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,k2
  complex :: c1,c2,c3

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    do iky = 1,ikty
      ky = kya(iky)
      do ikx = 1,iktx
        kx = kxa(ikx)     
        if (L(ikx,iky,ikz).eq.1) then
          k2 = max(kx*kx + ky*ky + kz*kz, 1.e-15)
          c1 =  (k2-kx*kx)*zx(ikx,iky,ikz) - kx*ky*zy(ikx,iky,ikz) - kx*kz*zz(ikx,iky,ikz)
          c2 = -ky*kx*zx(ikx,iky,ikz) + (k2-ky*ky)*zy(ikx,iky,ikz) - ky*kz*zz(ikx,iky,ikz)
          c3 = -kz*kx*zx(ikx,iky,ikz) - kz*ky*zy(ikx,iky,ikz) + (k2-kz*kz)*zz(ikx,iky,ikz)
          zx(ikx,iky,ikz) = c1 / k2
          zy(ikx,iky,ikz) = c2 / k2
          zz(ikx,iky,ikz) = c3 / k2
        endif
      enddo
    enddo
  enddo
  return
end subroutine proj


subroutine realit(zk)

! Enforces the reality condition on the plane kx=0
! by writing on modes with L(ikx,iky,ikz) = 0.

  use param
  use param_fftw
  implicit none
  complex, intent(inout) :: zk(iktx,ikty,iktzp)
  integer, parameter :: iktyh = ikty/2
  integer :: n2h,n3h
  integer :: ikx,iky,ikz,kz,inkz,ky,inky
  integer :: status(MPI_STATUS_SIZE)
  complex :: buf1(iktyh,iktzp-1),buf2(iktyh)
  integer :: nto,nfrom,nbuf1,nbuf2,nph

  nph   = npe/2
  nbuf1 = iktyh*(iktzp-1)
  nbuf2 = iktyh
  n2h = n2/2
  n3h = n3/2

! First, negative ky axis; no communication required
! Set zk(0,-ky,0) = conjg(zk(0,ky,0)) 
  if (mype.eq.0) then
     ikx = 1
     ikz = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky
        inky = n2h+1+ky
        zk(ikx,inky,ikz) = conjg( zk(ikx,iky,ikz) )           
     enddo
  endif

! Next, send kz>0,ky>=0 to kz<0,ky<=0
! Write zk(0,-ky,-kz) = conjg(zk(0,ky, kz)) 
! Have to do it in 2 parts
  if (mype.le.nph-1) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,istatus)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,istatus)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        zk(ikx,iky,inkz) = conjg( buf1(iky,ikz) )  ! do ky=0
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if ((mype.gt.0).and.(mype.le.nph-1)) then  
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,istatus)
  endif
  if (mype.gt.nph) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,istatus)
     ikx  = 1
     inkz = 1
     iky  = 1
     zk(ikx,iky,inkz) = conjg( buf2(iky) )  ! do ky=0
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif

! Finally, send kz<0,ky>=0 to kz>0,ky<=0
! Write zk(0,-ky,kz) = conjg(zk(0,ky,-kz)) 
! Have to do it in 2 parts
  if (mype.ge.nph) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,istatus)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,istatus)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if (mype.gt.nph) then
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,istatus)
  endif
  if ((mype.gt.0).and.(mype.le.nph-1)) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,istatus)
     ikx  = 1
     inkz = 1
     iky  = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif

  return
end subroutine realit

