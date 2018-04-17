module param
  use, intrinsic :: iso_c_binding
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters and global varables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -----------------------------------------------------------------
! Set Model Resolution
  integer(C_INTPTR_T), parameter :: n1=576, n2=576, n3=576
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Number of Processors (must divide n2 and iktz)
  integer, save         :: npe
! -----------------------------------------------------------------

! -----------------------------------------------------------------
  ! x,y,z planes where slices will be output
  integer(C_INTPTR_T),dimension(*) :: xslice(2) = (/ 1,int(n1)/2 /)
  integer(C_INTPTR_T),dimension(*) :: yslice(2) = (/ 1,int(n2)/2 /)
  integer(C_INTPTR_T),dimension(*) :: zslice(2) = (/ 1,int(n3)/2 /)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Misc. parameters
  real, parameter            :: twopi = 4.*asin(1.)
  real, parameter            :: sqrt2 = sqrt(2.)
  complex, parameter :: zi    = cmplx(0.,1.)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Domain Size
  real, parameter :: L1=twopi, L2=twopi, L3=twopi   
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Time
  real,    parameter   :: delt  = 0.02                   ! timestep
  real,    parameter   :: tstop = 500.                   ! length of integration
  integer, parameter :: nstop  = int(tstop/delt)         ! number of timesteps
  integer, parameter :: nout   = nstop/50.               ! do diagnostics every nout timesteps
  integer, parameter :: nslout = nstop/1                 ! output field slices every nslout timesteps 
  integer, parameter :: rsflag = 1                       ! flag for dumping real space fields
  integer, parameter :: slxflag = 0                      ! flag for dumping real space field slices in X 
  integer, parameter :: slyflag = 0                      ! flag for dumping real space field slices in Y 
  integer, parameter :: slzflag = 0                      ! flag for dumping real space field slices in Z 
  integer, parameter :: rstflag = 1                      ! flag for dumping restart files
  integer, parameter :: nbig   = nstop/1                 ! dump real space fields every nbig timesteps
  integer, parameter :: ndump  = nstop/1                 ! dump restart file every ndump timesteps
  integer, parameter :: nrsp   = rsflag*(nstop/nbig+1)   ! how many physical space outputs
  integer, parameter :: nrst   = nstop/ndump             ! how many restart outputs
  integer, parameter :: nsl    = (nstop/nslout+1)        ! how many slice outputs
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Stratification and Rotation
  real, parameter :: aj   =   0.3                         ! thermal expansivity. NOTE: N^2 = aj*bj
  real, parameter :: bj   =   0.3                         ! background theta gradient
  real, parameter :: bf2  =   aj*bj                       ! Brunt-Vaisalla frequency squared
  real, parameter :: bf   =   sqrt(bf2)                   ! Brunt-Vaisalla frequency
  real, parameter :: cor  =   1.e-15                      ! Coriolis parameter
  real, parameter :: cor2 =   cor*cor                     ! Coriolis parameter squared
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Physical and Fourier Array Dimensions
  integer(C_INTPTR_T), parameter :: n1d = n1+2, n2d=n2, n3d=n3         ! size of physical arrays (padded for in-place transforms)
  integer(C_INTPTR_T), save      :: n2dp, n2p                          ! size of local arrays with mpi, n2dp=n2d/npe, n2p=n2/npe
  integer(C_INTPTR_T), parameter :: ktx = n1/2,  kty =n2/2, ktz =n3/2  ! max integer wavenumber
  integer(C_INTPTR_T), parameter :: iktx= ktx+1, ikty=n2,   iktz=n3    ! number of wavenumbers
  integer(C_INTPTR_T), save      :: iktzp                              ! number of local wavenumbers with mpi, iktzp=iktz/npe
  integer(C_INTPTR_T), parameter :: kts = n1                      ! for spectra; should be max(ktx,ktz)
  real, parameter           :: ktrunc_x = twopi/L1*float(n1)/3.   ! dimensional truncation wavenumber (x)
  real, parameter           :: ktrunc_y = twopi/L2*float(n2)/3.   ! dimensional truncation wavenumber (y)
  real, parameter           :: ktrunc_z = twopi/L3*float(n3)/3.   ! dimensional truncation wavenumber (z)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Dissipation
  integer, parameter :: ilap = 1, ilap2 = 2*ilap          ! hyperviscosity order (ilap=1 is regular viscosity)
  real, parameter    :: visch   =   0.06e-4               ! viscosity coeff horizontal
  real, parameter    :: viscz   =   visch                 ! viscosity coeff vertical
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Wavenumbers 
  integer, dimension(:,:,:), allocatable,save :: L                    ! Mask for wavenumber truncation
  real, dimension(:), allocatable, save       :: kxa, kya, kza        ! Wavenumbers
  real, save                                  :: time
! -----------------------------------------------------------------
! MPI Stuff
  integer, save :: mype
  integer(C_INTPTR_T), save :: locz, loczstart             ! size and position of block in z   
  integer(C_INTPTR_T), save :: alloc_local                 ! local data size for malloc
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Misc Stuff
  integer, save :: istatus
! -----------------------------------------------------------------

! -----------------------------------------------------------------
  ! subroutine flags
  integer, parameter :: SLICE_MODE=1, FULLFIELD_MODE=0  
  real,    parameter :: fftnorm = float(n1*n2*n3)
! -----------------------------------------------------------------

end module param
