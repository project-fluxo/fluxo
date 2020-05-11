
!
!> Module containing some phisics-based sensors
!
module MOD_Sensors
  implicit none
  
  
  integer, parameter :: SENS_NUM = 3 ! 1: shock sensor, 2: thermal sensor, 3: shear sensor
  
  ! Artificial viscosity variables
  character(LEN=255), parameter :: StrArtViscNames(SENS_NUM)=(/ character(len=255) :: 'Art_bulk_Viscosity',  &
                                                                                      'Art_Conductivity',  &
                                                                                      'Art_dyn_Viscosity'  /)
  
contains

!===================================================================================================================================
!> Physics based sensors, see:
!> Fernandez et al. "A physics-based shock capturing method for unsteady laminar and turbulent flows and turbulent flows"
!===================================================================================================================================
  pure subroutine SensorsByFernandezEtAl(sensors,U,gradPx,gradPy,gradPz,Mh_inv)
    USE MOD_PreProc       , only: PP_N
    use MOD_Equation_Vars , only: ConsToPrim
    use MOD_Equation_Vars , only: R, KappaM1, sKappaM1
    implicit none
    !-arguments-----------------------------------------
    real                    , intent(out) :: sensors(SENS_NUM)    !< Shock sensor 
    real, dimension(PP_nVar), intent(in)  :: U
    real, dimension(PP_nVar), intent(in)  :: gradPx,gradPy,gradPz !< Gradient of the primitive variables
    real                    , intent(in)  :: Mh_inv(3,3)          !< Inverse of metric tensor
    !-local-variables-----------------------------------
    real :: h_beta        ! Approximate elment size along the direction of the density gradient
    real :: s_theta       ! 
    real :: s_omega       !
    real :: gradRho(3)    ! Density gradient
    real :: gradT(3)      ! Temperature gradient
    real :: divV, divV2   ! Divergence of velocity
    real :: rotV(3)       ! Rotational of velocity
    real :: a             ! Speed of sound
    real :: a2            ! Speed of sound squared
    real :: prim(PP_nVar) ! Primitive vars
    real :: srho          ! Inverse of density
    real :: Vel2          ! Velocity squared
    real :: Vmax          ! Maximum isentropic velocity
    real :: M2            ! Mach**2
    real :: T0            ! Stanation temperature
    real :: L_u(3,3)      ! Shear velocity gradients: L_u = gradV - diag(gradV)
    real, parameter :: eps = 1.e-8 ! TODO: Change to a measure of machine zero... Ciuca (~eps_m), Fernandez (~eps_m^2), Ducros used 1e-30 (?!!)
    !---------------------------------------------------
    
    srho = 1./U(1)
    gradRho = [gradPx(1), gradPy(1), gradPz(1)]
    
    call ConsToPrim(prim,U)
    
    ! Speed of sound
    a2 = prim(5)/U(1)
    a  = sqrt(a2)
    
    ! Mach number
    Vel2 = sum(prim(2:4)**2)
    M2   = Vel2 / a2
    
    ! Divergence
    divV  = gradPx(2) + gradPy(3) + gradPz(4)
    divV2 = divV**2
    
    ! Rotational
    rotV  = [gradPy(4) - gradPz(3), gradPz(2) - gradPx(4), gradPx(3) - gradPy(2)]
    
    ! Shock sensor
    ! ************
    h_beta = 2. * norm2(gradRho) / sqrt(dot_product(gradRho,matmul(Mh_inv,gradRho)) + eps)  ! TODO: This is probably NOT the right metric tensor
    
    s_theta = - h_beta * divV / (PP_N * a)          ! TODO: "Critical" speed of sound?
    s_omega = divV2 / (divV2 + sum(rotV**2) + eps)
    
    sensors(1) = s_theta * s_omega
    
    ! Thermal sensor
    ! **************
    gradT  = sRho/R*[gradPx(5)-srho*prim(5)*gradPx(1), &  ! T_x = 1/(rho*R) *(p_x - p/rho*rho_x)
                     gradPy(5)-srho*prim(5)*gradPy(1), & 
                     gradPz(5)-srho*prim(5)*gradPz(1) ]
    
    T0 = prim(5)*srho/R * (1. + 0.5*KappaM1*M2)
    
    sensors(2) = 2. * norm2( matmul(Mh_inv,gradT) ) / (PP_N * T0)
    
    ! Shear sensor
    ! ************
    
    L_u(:,1) = [0.       , gradPx(3), gradPx(4)]
    L_u(:,2) = [gradPy(2), 0.       , gradPy(4)]
    L_u(:,3) = [gradPz(2), gradPz(3), 0.       ]
    
    L_u = matmul(L_u,Mh_inv)
    
    Vmax = sqrt (Vel2 + 2. * a2 * sKappaM1)
    
    sensors(3) = 2. * sqrt(sum(L_u**2)) / (PP_N * Vmax) ! TODO: Taking Frobenius norm... Change for exact spectral norm?
    
  end subroutine SensorsByFernandezEtAl
  
end module MOD_Sensors
