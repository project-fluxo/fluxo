module MOD_IDP_Vars
  implicit none
  
  public
  
  ! IDP methods:
  logical           :: IDPMathEntropy = .FALSE.
  logical           :: IDPSpecEntropy = .FALSE.
  logical           :: IDPSemiDiscEnt = .FALSE.
  logical           :: IDPPositivity  = .FALSE.
  logical           :: IDPDensityTVD  = .FALSE.
  
  ! Additional options
  integer           :: IDPMaxIter
  
  ! Internal definitions:
  logical :: IDPneedsUprev_ext = .FALSE.
  logical :: IDPneedsUbar      = .FALSE.
  logical :: IDPneedsUsafe     = .FALSE.
  logical :: IDPneedsUprev     = .FALSE.
  
  ! Variables to modify the local stencil where the bounds are computed
  integer, allocatable :: idx_p1(:)
  integer, allocatable :: idx_m1(:)
  
  ! Parameters:
  real, parameter :: alpha_maxIDP  = 1.0
  real, parameter :: NEWTON_RELTOL = 1.e-12 ! Relative tolerance to exit Newton loop
  real, parameter :: NEWTON_ABSTOL = 1.e-8  ! Absolute tolerance (with respect to the value of the entropy goal, tol = NEWTON_ABSTOL*s_goal)
  
  ! Containers
  real,allocatable    :: FFV_m_FDG(:,:,:,:,:)
  
  
  real,allocatable    :: Uprev         (:,:,:,:,:)
  real,allocatable    :: Uprev_ext     (:,:,:,:,:)
  
  real,allocatable    :: Usafe         (:,:,:,:,:)
  real,allocatable    :: p_safe          (:,:,:,:)
  real,allocatable    :: Usafe_ext     (:,:,:,:,:)
  
  real,allocatable    :: Flux_ext      (:,:,:,:,:)
  
  real,allocatable    :: EntPrev       (:,:,:,:,:)
  real,allocatable    :: EntPrev_master(:,:,:,:)
  real,allocatable    :: EntPrev_slave (:,:,:,:)
  real,allocatable    :: EntPrev_ext   (:,:,:,:,:)
  
  real,allocatable    :: Ubar_xi       (:,:,:,:,:)
  real,allocatable    :: Ubar_eta      (:,:,:,:,:)
  real,allocatable    :: Ubar_zeta     (:,:,:,:,:)
end module MOD_IDP_Vars
