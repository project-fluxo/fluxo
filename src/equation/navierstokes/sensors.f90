
!
!> Module containing some phisics-based sensors
!
module MOD_Sensors
  USE MOD_PreProc       , only: PP_Pi
  implicit none
  
  
  integer, parameter :: SENS_NUM = 3 ! 1: shock sensor, 2: thermal sensor, 3: shear sensor
  
  ! Artificial viscosity variables
  character(LEN=255), parameter :: StrArtViscNames(SENS_NUM)=(/ character(len=255) :: 'Art_bulk_Viscosity',  &
                                                                                      'Art_Conductivity',  &
                                                                                      'Art_dyn_Viscosity'  /)
  
  ! Parameters for the smooth filter
  real, parameter :: b = 100.
  real, parameter :: atanb_sPi = atan(100.) / PP_Pi
contains
  
  pure function lmax(s)
    
    implicit none
    !-arguments-----------------------------------------
    real, intent(in) :: s
    real             :: lmax
    !---------------------------------------------------
    
    lmax = s * atan(b*s)/PP_Pi + 0.5 * s - atanb_sPi + 0.5
    
  end function lmax
  
  pure function lmin(s)
    
    implicit none
    !-arguments-----------------------------------------
    real, intent(in) :: s
    real             :: lmin
    !---------------------------------------------------
    
    lmin = s - lmax(s)
    
  end function lmin
  
  pure function SmoothFilter(s,s0,smax)
    
    implicit none
    !-arguments-----------------------------------------
    real, intent(in) :: s
    real, intent(in) :: s0
    real, intent(in) :: smax
    real             :: SmoothFilter
    !---------------------------------------------------
    
    SmoothFilter = lmin(lmax(s-s0) - smax) + smax
    
  end function SmoothFilter

!===================================================================================================================================
!> Physics based sensors, see:
!> Fernandez et al. "A physics-based shock capturing method for unsteady laminar and turbulent flows and turbulent flows"
!===================================================================================================================================
  pure subroutine SensorsByFernandezEtAl(sensors,U,gradPx,gradPy,gradPz,Mh_inv,covMetrics,sigmamin_Mh,artVisc)
    USE MOD_PreProc       , only: PP_N
    use MOD_Equation_Vars , only: ConsToPrim
    use MOD_Equation_Vars , only: KappaM1, sKappaM1, Kappa, sKappaP1 !R, 
    implicit none
    !-arguments-----------------------------------------
    real                    , intent(out) :: sensors(SENS_NUM)    !< Shock sensor 
    real, dimension(PP_nVar), intent(in)  :: U
    real, dimension(PP_nVar), intent(in)  :: gradPx,gradPy,gradPz !< Gradient of the primitive variables
    real                    , intent(in)  :: Mh_inv(3,3)          !< Inverse of Möller's metric tensor
    real                    , intent(in)  :: covMetrics(3,3)      !< Covariant metric tensor
    real                    , intent(in)  :: sigmamin_Mh
    real          , intent(out) :: artVisc(SENS_NUM)    !< Artificial viscosity ! optional?? bad idea..
    !-local-variables-----------------------------------
    real :: h_beta        ! Approximate elment size along the direction of the density gradient
    real :: s_theta       ! 
    real :: s_omega       !
    real :: gradRho(3)    ! Density gradient
    real :: gradT(3)      ! Temperature gradient (Nondimensional, excluding R)
    real :: divV, divV2   ! Divergence of velocity
    real :: rotV(3)       ! Rotational of velocity
    real :: a             ! Speed of sound
    real :: a2            ! Speed of sound squared
    real :: prim(PP_nVar) ! Primitive vars
    real :: srho          ! Inverse of density
    real :: Vel2          ! Velocity squared
    real :: Vmax          ! Maximum isentropic velocity
    real :: M2            ! Mach**2
    real :: T0            ! Stagnation temperature (Nondimensional, excluding R)
    real :: L_u(3,3)      ! Shear velocity gradients: L_u = gradV - diag(gradV)
    real :: astar         ! Critical speed of sound
    real :: astar2        ! Square of the critical speed of sound
    real :: h_kappa       ! 
    real :: h_mu
    real, parameter :: eps = 1.e-8 ! TODO: Change to a measure of machine zero... Ciuca (~eps_m), Fernandez (~eps_m^2), Ducros used 1e-30 (?!!)
    real, parameter :: k_beta = 1.5
    real, parameter :: sPr_beta = 1./0.9
    real :: Cp            ! Nondimensional Cp (excluding R)
    !---------------------------------------------------
    
    Cp = Kappa*sKappaM1 ! R* TODO: Declare in equation
    
    srho = 1./U(1)
    gradRho = [gradPx(1), gradPy(1), gradPz(1)]
    
    call ConsToPrim(prim,U)
    
    Vel2 = sum(prim(2:4)**2)
    
    ! Speed of sound
    a2 = Kappa*prim(5)*srho
    astar2 = 2. * sKappaP1 * (a2 + 0.5 * KappaM1 * Vel2)
    astar  = sqrt(astar2)
    
    
    ! Mach number
    M2   = Vel2 / a2
    
    ! Divergence
    divV  = gradPx(2) + gradPy(3) + gradPz(4)
    divV2 = divV**2
    
    ! Rotational
    rotV  = [gradPy(4) - gradPz(3), gradPz(2) - gradPx(4), gradPx(3) - gradPy(2)]
    
    ! Shock sensor
    ! ************
    h_beta = 2. * norm2(gradRho) / sqrt(dot_product(gradRho,matmul(Mh_inv,gradRho)) + eps)
    
    s_theta = - h_beta * divV / (PP_N * astar)
    s_omega = divV2 / (divV2 + sum(rotV**2) + eps)
    
    sensors(1) = s_theta * s_omega
    
    ! Smooth it:
    sensors(1) = SmoothFilter(sensors(1), 0.01, 2./sqrt(Kappa**2 - 1)) ! TODO: Precompute
    
    ! Compute artificial bulk viscosity to ensure a cell Péclet of O(1)
      artVisc(1) = sensors(1) * U(1) * k_beta * h_beta * sqrt(Vel2 + astar2) / PP_N
    
    ! Thermal sensor
    ! **************
    gradT  = sRho  *[gradPx(5)-srho*prim(5)*gradPx(1), &  ! T_x = 1/(rho*R) *(p_x - p/rho*rho_x)
                     gradPy(5)-srho*prim(5)*gradPy(1), & 
                     gradPz(5)-srho*prim(5)*gradPz(1) ]  ! /R
    
    T0 = prim(5)*srho * (1. + 0.5*KappaM1*M2) ! /R
    
    sensors(2) = 2. * norm2( matmul(covMetrics,gradT) ) / (PP_N * T0)
    
    ! Smooth it:
    sensors(2) = SmoothFilter(sensors(2), 1., 2.)
    
    ! Compute artificial conductivity to ensure a cell Péclet of O(1)
      h_kappa = 2. * norm2(gradT) / sqrt(dot_product(gradT,matmul(Mh_inv,gradT)) + eps)
      artVisc(2) = artVisc(1) * Cp * sPr_beta + sensors(2) * Cp * U(1) * h_kappa * sqrt(Vel2 + astar2) / PP_N
    
    ! Shear sensor
    ! ************
    
    L_u(:,1) = [0.       , gradPx(3), gradPx(4)]
    L_u(:,2) = [gradPy(2), 0.       , gradPy(4)]
    L_u(:,3) = [gradPz(2), gradPz(3), 0.       ]
    
    L_u = matmul(L_u,covMetrics)
    
    Vmax = sqrt (Vel2 + 2. * a2 * sKappaM1)
    
    sensors(3) = 2. * sqrt(sum(L_u**2)) / (PP_N * Vmax) ! TODO: Taking Frobenius norm... Change for exact spectral norm?
    
    ! Smooth it:
    sensors(3) = SmoothFilter(sensors(3), 1., 2.)
    
    ! Compute artificial dynamic viscosity to ensure a cell Péclet of O(1)
      h_mu = 2. * sigmamin_Mh
      artVisc(3) = sensors(3) * U(1) * h_mu * sqrt(Vel2 + astar2) / PP_N
    
  end subroutine SensorsByFernandezEtAl
  
end module MOD_Sensors
