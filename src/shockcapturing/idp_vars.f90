module MOD_IDP_Vars
  implicit none
  
  public
  
  ! IDP methods:
  logical           :: IDPMathEntropy = .FALSE.
  logical           :: IDPSpecEntropy = .FALSE.
  logical           :: IDPPositivity  = .FALSE.
  logical           :: IDPDensityTVD  = .FALSE.
  logical           :: IDPPressureTVD = .FALSE.
  
  ! Additional options
  logical           :: IDPForce2D
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
  real, parameter :: NEWTON_ABSTOL = 1.e-14 ! Absolute tolerance (with respect to the value of the entropy goal, tol = NEWTON_ABSTOL*s_goal)
  real            :: IDPgamma               ! User-defined: Constant for the subcell limiting of convex (nonlinear) constraints (must be IDPgamma>=2*d, where d is the number of dimensions of the problem)
  
  ! Variables for analyze bounds
#if DEBUG || IDP_CHECKBOUNDS
  integer             :: idp_bounds_num
  character(len=255)  :: idp_bounds_names(6)
  real                :: idp_bounds_delta(6)
#endif /*DEBUG || IDP_CHECKBOUNDS*/

  ! Containers
  real,allocatable    :: FFV_m_FDG(:,:,:,:,:)
  
  
  real,allocatable    :: Uprev         (:,:,:,:,:)
  real,allocatable    :: Uprev_ext     (:,:,:,:,:)
  
  real,allocatable    :: Usafe         (:,:,:,:,:)
  real,allocatable    :: p_safe          (:,:,:,:)
  real,allocatable    :: Usafe_ext     (:,:,:,:,:)
  
  real,allocatable    :: EntPrev       (:,:,:,:,:)
  real,allocatable    :: EntPrev_master(:,:,:,:)
  real,allocatable    :: EntPrev_slave (:,:,:,:)
  real,allocatable    :: EntPrev_ext   (:,:,:,:,:)
  
  real,allocatable    :: Ubar_xi       (:,:,:,:,:)
  real,allocatable    :: Ubar_eta      (:,:,:,:,:)
  real,allocatable    :: Ubar_zeta     (:,:,:,:,:)
  
  real                :: maxdt_IDP = huge(1.0)
  
  real                :: dalpha
#if LOCAL_ALPHA
  real,allocatable    :: dalpha_loc    (:,:,:)    ! Node-wise alpha correction for an element
#endif /*LOCAL_ALPHA*/

  ! Bounds
  real,allocatable    :: rho_min    (:,:,:)
  real,allocatable    :: rho_max    (:,:,:)
  real,allocatable    :: s_min      (:,:,:)
  real,allocatable    :: s_max      (:,:,:)
  real,allocatable    :: p_min      (:,:,:)
  real,allocatable    :: p_max      (:,:,:)
  
  ! Type to store parameters for newton iterations
  type IDPparam_t
    real :: bound               ! Specific bound
    real :: F_antidiff(PP_nVar) ! IDPGamma- (if needed) and mass-matrix-scaled antidiffusive flux with the right sign. Fan := ± sWGP*sJ*(F_DG - F_FV) ... (± must be consistent with the equation update)
    real :: dt                  ! current time-step size
  end type IDPparam_t
  
  abstract interface 
    ! Interface for Newton's method goal function (and its derivative)
    pure function i_sub_Goal(param,Ucurr) result(goal)
      import IDPparam_t
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: goal
    end function i_sub_Goal
    ! Interface for Newton's method initial check
    pure function i_sub_InitialCheck(param,goalFunction) result(check)
      import IDPparam_t
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: goalFunction ! Current solution
      logical                      :: check
    end function i_sub_InitialCheck
  end interface
end module MOD_IDP_Vars