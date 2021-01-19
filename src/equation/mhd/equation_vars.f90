!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================

!==================================================================================================================================
!> Contains global variables needed for the MHD equation
!==================================================================================================================================
#include "defines.h"
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcSource      !< logical to define if a source term (e.g. exactfunc) is added
INTEGER             :: IniExactFunc      !< Exact Function for initialization
INTEGER             :: IniRefState       !< RefState for initialization (case IniExactFunc=1 only)
INTEGER,PARAMETER   :: nAuxVar=8         !< number of auxiliary variables for average flux
INTEGER             :: nRefState         !< number of RefState in inifile
REAL,ALLOCATABLE    :: RefStatePrim(:,:) !< primitive reference states
REAL,ALLOCATABLE    :: RefStateCons(:,:) !< =primToCons(RefStatePrim)
! Boundary condition arrays
REAL,ALLOCATABLE    :: BCData(:,:,:,:)  !< data for steady state boundary conditions
INTEGER,ALLOCATABLE :: nBCByType(:)     !< Number of sides for each boundary
INTEGER,ALLOCATABLE :: BCSideID(:,:)    !< SideIDs for BC types
REAL                :: R                !< Gas constant
#if PARABOLIC
REAL                :: mu               !< fluid viscosity, NOT mu0 like in Navier-stokes ANYMORE!!!!
REAL                :: eta              !< Current resistivity
REAL                :: etasmu_0         !< =eta/mu0 
REAL                :: Pr               !< Prandtl number
REAL                :: KappasPr         !< =kappa/Pr
REAL                :: s23              !< (=2/3 for Navier stokes) part of stress tensor: mu*((nabla v)+(nabla v)^T-s23*div(v))
#  ifdef PP_ANISO_HEAT
REAL                :: kperp            !< perpendicular (to magnetic field) heat diffusion coefficient
REAL                :: kpar             !< parallel (to magnetic field) heat diffusion coeffcient 
#  endif /*PP_ANISO_HEAT*/             
#endif /*PARABOLIC*/                   
REAL                :: mu_0             !< magnetic permeability in vacuum
REAL                :: smu_0            !< =1/mu_0
REAL                :: s2mu_0           !< =1/(2*mu_0)
REAL                :: Kappa            !< ratio of specific heats
REAL                :: KappaM1          !< = kappa - 1
REAL                :: KappaM2          !< = kappa - 2
REAL                :: sKappaM1         !< = 1/(kappa -1)
REAL                :: KappaP1          !< = kappa + 1
REAL                :: sKappaP1         !< = 1/(kappa +1)
REAL                :: IniWavenumber(3) !< wavenumbers in 3 directions (sinus periodic with exactfunc=6)
REAL                :: IniFrequency     !< frequency for exactfunc
REAL                :: IniAmplitude     !< amplitude for exactfunc
REAL                :: IniHalfwidth     !< halfwidth for exactfunc
REAL                :: IniCenter(3)     !< center point for exactfunc
REAL                :: IniDisturbance   !< disturbance scaling for exactfunc
#ifdef PP_GLM
LOGICAL             :: GLM_init=.FALSE. !< switch set true when GLM_dtch1 is computed
REAL                :: GLM_dtch1        !< timestep for ch=1 (initialized in calctimestep)
REAL                :: GLM_ch           !< Divergence correction speed
REAL                :: GLM_scale        !< scaling of maximum divergence speed (timestep security)
REAL                :: GLM_scr          !< 1/cr. damping of divergence error, factor chi=ch^2/cp^2,cr=cp^2/ch~0.18,chi=1/cr
LOGICAL             :: divBSource       !< switch for adding source terms depending on diverngece errors
#endif /*PP_GLM*/


CHARACTER(LEN=255),DIMENSION(PP_nVar),PARAMETER :: StrVarNames(PP_nVar)=(/ CHARACTER(LEN=255) :: &
                   'Density',    &
                   'MomentumX',  &
                   'MomentumY',  &
                   'MomentumZ',  &
                   'TotalEnergy',&
                   'MagneticFieldX',  &
                   'MagneticFieldY',  &
                   'MagneticFieldZ'   &

#ifdef PP_GLM
                   ,'Psi' &
#endif 
                   /)
CHARACTER(LEN=255),DIMENSION(PP_nVar),PARAMETER :: StrVarNamesPrim(PP_nVar)=(/ CHARACTER(LEN=255) :: &
                   'Density',    &
                   'VelocityX',  &
                   'VelocityY',  &
                   'VelocityZ',  &
                   'Pressure',&
                   'MagneticFieldX',  &
                   'MagneticFieldY',  &
                   'MagneticFieldZ'   &

#ifdef PP_GLM
                   ,'Psi' &
#endif 
                   /)

integer, parameter :: nIndVar = 4
character(len=255) :: IndicatorQuantityNames(nIndVar) = (/character(len=132) :: 'Density','Pressure','DensityTimesPressure','KinPlusMagEnergy'/)

LOGICAL           :: EquationInitIsDone=.FALSE. !< Init switch  
!procedure pointers
INTEGER             :: WhichRiemannSolver       !< choice of riemann solver
#if (PP_DiscType==2)
!procedure pointers for split form DG
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
#endif /*PP_DiscType==2*/
PROCEDURE(i_sub_RiemannVolFluxAndDissipMatrices),POINTER :: RiemannVolFluxAndDissipMatrices =>Null()
PROCEDURE(i_sub_SolveRiemannProblem ),POINTER :: SolveRiemannProblem  =>Null() !< procedure pointer to riemann solver 
PROCEDURE(i_sub_VolumeFluxAverage   ),POINTER :: VolumeFluxAverage    =>Null() !< procedure pointer to 1D two-point average flux
PROCEDURE(i_sub_VolumeFluxAverageVec),POINTER :: VolumeFluxAverageVec =>Null() !< procedure pointer to 3D two-point average flux
!==================================================================================================================================
ABSTRACT INTERFACE
  SUBROUTINE i_sub_SolveRiemannProblem(ConsL,ConsR,Flux)
    REAL,DIMENSION(1:PP_nVar),INTENT(IN)  :: ConsL !<  left conservative state  
    REAL,DIMENSION(1:PP_nVar),INTENT(IN)  :: ConsR !< right conservative state
    REAL,DIMENSION(1:PP_nVar),INTENT(OUT) :: Flux  !< numerical flux
  END SUBROUTINE i_sub_SolveRiemannProblem
  
  PURE SUBROUTINE i_sub_RiemannVolFluxAndDissipMatrices(ConsL,ConsR,F,Dmatrix,Rmatrix)
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(IN)  :: ConsL    !<  left conservative state  
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(IN)  :: ConsR    !< right conservative state
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(OUT) :: F        !< Central volume flux
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(OUT) :: Dmatrix  !< Dissipation matrix
    REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: Rmatrix  !< Right-eigenvector matrix
  END SUBROUTINE i_sub_RiemannVolFluxAndDissipMatrices

  PURE SUBROUTINE i_sub_VolumeFluxAverage(UL,UR,Fstar)
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< central flux in x
  END SUBROUTINE i_sub_VolumeFluxAverage

  PURE SUBROUTINE i_sub_VolumeFluxAverageVec (UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
    REAL,DIMENSION(8),INTENT(IN)        :: UauxL          !< left auxiliary variables
    REAL,DIMENSION(8),INTENT(IN)        :: UauxR          !< right auxiliary variables
    REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
    REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed central flux
  END SUBROUTINE i_sub_VolumeFluxAverageVec
  
  pure subroutine i_indicatorFunction(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_indicatorFunction
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim 
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

!INTERFACE ConsToPrimVec
!  MODULE PROCEDURE ConsToPrimVec
!END INTERFACE

!INTERFACE PrimToConsVec
!  MODULE PROCEDURE PrimToConsVec
!END INTERFACE

INTERFACE ConsToEntropy
  MODULE PROCEDURE ConsToEntropy
END INTERFACE

!INTERFACE ConsToEntropyVec
!  MODULE PROCEDURE ConsToEntropyVec
!END INTERFACE

#if PARABOLIC
INTERFACE ConvertToGradPrim
  MODULE PROCEDURE ConvertToGradPrim
END INTERFACE

!INTERFACE ConvertToGradPrimVec
!  MODULE PROCEDURE ConvertToGradPrimVec
!END INTERFACE
#endif /*PARABOLIC*/

INTERFACE WaveSpeeds1D
  MODULE PROCEDURE WaveSpeeds1D
END INTERFACE

INTERFACE FastestWave1D
  MODULE PROCEDURE FastestWave1D
END INTERFACE

INTERFACE FastestWave3D
  MODULE PROCEDURE FastestWave3D
END INTERFACE

INTERFACE FastestWave1D_Roe
  MODULE PROCEDURE FastestWave1D_Roe
END INTERFACE

CONTAINS


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(PP_nVar) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(PP_nVar) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! rho
prim(1)=cons(1)
! velocity
prim(2:4)=cons(2:4)*sRho
! GAS PRESSURE (WITHOUT magnetic pressure)
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4))-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar)))
! B,psi
prim(6:PP_nVar)=cons(6:PP_nVar)
END SUBROUTINE ConsToPrim



!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
PURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: prim(PP_nVar) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cons(PP_nVar) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
! rho
cons(1)=prim(1)
! velocity
cons(2:4)=prim(2:4)*prim(1)
! total energy
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))+s2mu_0*(SUM(prim(6:PP_nVar)*prim(6:PP_nVar)))
! B,psi
cons(6:PP_nVar)=prim(6:PP_nVar)
END SUBROUTINE PrimToCons


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
PURE SUBROUTINE ConsToPrimVec(dim2,prim,cons)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: dim2 
REAL,INTENT(IN)     :: cons(PP_nVar,dim2) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(PP_nVar,dim2) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i
!==================================================================================================================================
DO i=1,dim2
  CALL ConsToPrim(Prim(:,i),Cons(:,i))
END DO!i
END SUBROUTINE ConsToPrimVec


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
PURE SUBROUTINE PrimToConsVec(dim2,prim,cons)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: dim2 
REAL,INTENT(IN)     :: prim(PP_nVar,dim2) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cons(PP_nVar,dim2) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i
!==================================================================================================================================
DO i=1,dim2
  CALL PrimToCons(Prim(:,i),Cons(:,i))
END DO!i
END SUBROUTINE PrimToConsVec


!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PURE FUNCTION ConsToEntropy(cons) RESULT(Entropy)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar)             :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: srho,u,v,w,v2s2,rho_sp,s
!==================================================================================================================================
srho   = 1./cons(1)
u      = cons(2)*srho
v      = cons(3)*srho
w      = cons(4)*srho
v2s2   = 0.5*(u*u+v*v+w*w)
rho_sp = cons(1)/(KappaM1*(cons(5)-cons(1)*v2s2-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar)))) ! pressure includes psi^2 if PP_nVar=9
!s      = LOG(p) - kappa*LOG(cons(1))
s      = - LOG(rho_sp*(cons(1)**kappaM1))

! Convert to entropy variables
entropy(1)         =  (kappa-s)*skappaM1 - rho_sp*v2s2  !(kappa-s)/(kappa-1)-beta*|v|^2
entropy(2)         =  rho_sp*u                  ! 2*beta*u
entropy(3)         =  rho_sp*v                  ! 2*beta*v
entropy(4)         =  rho_sp*w                  ! 2*beta*w
entropy(5)         = -rho_sp                    !-2*beta
entropy(6:PP_nVar) =  rho_sp*cons(6:PP_nVar)    ! 2*beta*B +2*beta*psi

END FUNCTION ConsToEntropy

!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
SUBROUTINE ConsToEntropyVec(dim2,Entropy,cons)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: dim2 
REAL,INTENT(IN)     :: cons(PP_nVar,dim2) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Entropy(PP_nVar,dim2) !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i
REAL                :: srho,u,v,w,v2s2,rho_sp,s
!==================================================================================================================================
DO i=1,dim2
  srho   = 1./cons(1,i)
  u      = cons(2,i)*srho
  v      = cons(3,i)*srho
  w      = cons(4,i)*srho
  v2s2   = 0.5*(u*u+v*v+w*w)
  rho_sp = cons(1,i)/(KappaM1*(cons(5,i)-cons(1,i)*v2s2-s2mu_0*SUM(cons(6:PP_nVar,i)*cons(6:PP_nVar,i)))) ! pressure includes psi^2 if PP_nVar=9
  !s      = LOG(p) - kappa*LOG(cons(1))
  s      = - LOG(rho_sp*(cons(1,i)**kappaM1))
  
  ! Convert to entropy variables
  entropy(1,i)         =  (kappa-s)*skappaM1 - rho_sp*v2s2  !(kappa-s)/(kappa-1)-beta*|v|^2
  entropy(2,i)         =  rho_sp*u                  ! 2*beta*u
  entropy(3,i)         =  rho_sp*v                  ! 2*beta*v
  entropy(4,i)         =  rho_sp*w                  ! 2*beta*w
  entropy(5,i)         = -rho_sp                    !-2*beta
  entropy(6:PP_nVar,i) =  rho_sp*cons(6:PP_nVar,i)    ! 2*beta*B +2*beta*psi
END DO!i
END SUBROUTINE ConsToEntropyVec


#if PARABOLIC
!==================================================================================================================================
!> transform gradient from conservative / primitive or entropy variables to primitive variables 
!==================================================================================================================================
FUNCTION ConvertToGradPrim(cons,grad_in) RESULT(gradP)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(PP_nVar)    !< conservative state 
REAL,INTENT(IN)     :: grad_in(PP_nVar) !< can be gradient of conservative / primivite /entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: gradP(PP_nVar) !<  gradient of primitive variables (rho,v1,v2,v3,p,B1,B2,B3,psi)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
#if (PP_Lifting_Var==1) 
REAL  :: sRho,u,v,w,gradu,gradv,gradw
#elif (PP_Lifting_Var==3) 
REAL  :: sRho,u,v,w,gradu,gradv,gradw,Ekin,p,rho_sp,p_srho
#endif /*PP_Lifting_Var*/
!==================================================================================================================================

#if (PP_Lifting_Var==1) 
  !grad_in is gradient of conservative variable
  sRho      = 1./cons(1)
      u     = cons(2)*sRho
      v     = cons(3)*sRho
      w     = cons(4)*sRho
  gradu     = sRho*(grad_in(2)-grad_in(1)*u)
  gradv     = sRho*(grad_in(3)-grad_in(1)*v)
  gradw     = sRho*(grad_in(4)-grad_in(1)*w)
  
  !density gradient
  gradP(1)  = grad_in(1)
  !velocity gradient
  gradP(2)  = gradu
  gradP(3)  = gradv
  gradP(4)  = gradw
  !pressure gradient
  gradP(5)  = KappaM1*(grad_in(5)                                        & !gradE
                       -(0.5*grad_in(1)*(u*u+v*v+w*w) + (cons(2)*gradu+cons(3)*gradv+cons(4)*gradw)) & !-grad_Ekin
                       - smu_0*SUM(cons(6:PP_nVar)*grad_in(6:PP_nVar))   ) !-grad_Emag
  !gradient of B,psi same in primitive
  gradP(6:PP_nVar)=grad_in(6:PP_nVar) 
#elif (PP_Lifting_Var==2) 
  !grad_in is gradient of primitive variable, do nothing
  gradP(:)=grad_in(:)
#elif (PP_Lifting_Var==3) 
  !grad_in is gradient of entropy variable,  entropy variables are:
  ! w(1)= (kappa-s)/(kappa-1)+(rho/p)/2*|v|^2  ,w(2:4)=(rho/p)*v, w(5)=-(rho/p), w(6:8)=(rho/p)*B(:) , w(9)=(rho/p)*psi 

  !gradient of (p/rho),  (p/rho)_x = -grad_in(5)

  sRho   = 1./cons(1)
  u      = cons(2)*sRho
  v      = cons(3)*sRho
  w      = cons(4)*sRho
  Ekin   = 0.5*(cons(2)*u+cons(3)*v+cons(4)*w)
  p      = KappaM1*(cons(5)-Ekin-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar))) ! includes psi^2 if PP_nVar=9
  rho_sp = cons(1)/p
  p_srho = p * sRho
  
  gradu  = p_sRho * (grad_in(2) +u*grad_in(5))
  gradv  = p_sRho * (grad_in(3) +v*grad_in(5))
  gradw  = p_sRho * (grad_in(4) +w*grad_in(5))
  
  !density gradient, rho_x = rho*w1_x + (rho/p)_x * (-p/(gamma-1) + 1/2*rho*|v|^2 )  + (rho/p)*(rho*v) . v_x
  gradP(1)  = cons(1)*grad_in(1) - grad_in(5)*(Ekin -p*sKappaM1) + rho_sp*(cons(2)*gradu+cons(3)*gradv+cons(4)*gradw)
  !velocity gradient
  gradP(2)  = gradu
  gradP(3)  = gradv
  gradP(4)  = gradw
  !pressure gradient, =1/(rho/p)*(rho_x-p*(rho/p)_x)
  gradP(5)  = p_srho * (gradP(1) + p *grad_in(5))
  !gradient of B,psi: 
  gradP(6:PP_nVar) = p_sRho * (grad_in(6:PP_nVar) +cons(6:PP_nVar)*grad_in(5))
#endif /*PP_Lifting_Var*/
END FUNCTION ConvertToGradPrim

!==================================================================================================================================
!> transform gradient from conservative / primitive or entropy variables to primitive variables 
!==================================================================================================================================
SUBROUTINE ConvertToGradPrimVec(dim2,cons,gradP)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: dim2 
REAL,INTENT(IN)     :: cons(PP_nVar,dim2)    !< conservative state 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: gradP(PP_nVar,dim2) !<  on intput: can be gradient of conservative / primivite /entropy variables
                                             !<  on output: gradient of primitive variables (rho,v1,v2,v3,p,B1,B2,B3,psi)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i
#if (PP_Lifting_Var==1) 
REAL  :: sRho,u,v,w,gradu,gradv,gradw,grad_in1,grad_in5
#elif (PP_Lifting_Var==3) 
REAL  :: sRho,u,v,w,gradu,gradv,gradw,Ekin,p,rho_sp,p_srho,grad_in1,grad_in5
#endif /*PP_Lifting_Var*/
!==================================================================================================================================
DO i=1,dim2
#if (PP_Lifting_Var==1) 
  !grad_in is gradient of conservative variable
  sRho      = 1./cons(1,i)
      u     = cons(2,i)*sRho
      v     = cons(3,i)*sRho
      w     = cons(4,i)*sRho
  grad_in1  = gradP(1,i)
  gradu     = sRho*(gradP(2,i)-gradP(1,i)*u)
  gradv     = sRho*(gradP(3,i)-gradP(1,i)*v)
  gradw     = sRho*(gradP(4,i)-gradP(1,i)*w)
  grad_in5  = gradP(5,i)
  
  !density gradient
  !gradP(1,i)  = gradP(1,i)
  !velocity gradient
  gradP(2,i)  = gradu
  gradP(3,i)  = gradv
  gradP(4,i)  = gradw
  !pressure gradient
  gradP(5,i)  = KappaM1*(grad_in5                                        & !gradE
                       -(0.5*grad_in1*(u*u+v*v+w*w) + (cons(2,i)*gradu+cons(3,i)*gradv+cons(4,i)*gradw)) & !-grad_Ekin
                       - smu_0*SUM(cons(6:PP_nVar,i)*gradP(6:PP_nVar,i))   ) !-grad_Emag
  !gradient of B,psi same in primitive
  !gradP(6:PP_nVar,i)=grad_P(6:PP_nVar,i) 
#elif (PP_Lifting_Var==2) 
  !grad_in is gradient of primitive variable, do nothing
  gradP(:,i)=gradP(:,i)
#elif (PP_Lifting_Var==3) 
  !grad_in is gradient of entropy variable,  entropy variables are:
  ! w(1)= (kappa-s)/(kappa-1)+(rho/p)/2*|v|^2  ,w(2:4)=(rho/p)*v, w(5)=-(rho/p), w(6:8)=(rho/p)*B(:) , w(9)=(rho/p)*psi 

  !gradient of (p/rho),  (p/rho)_x = -grad_in(5)

  sRho   = 1./cons(1,i)
  u      = cons(2,i)*sRho
  v      = cons(3,i)*sRho
  w      = cons(4,i)*sRho
  Ekin   = 0.5*(cons(2,i)*u+cons(3,i)*v+cons(4,i)*w)
  p      = KappaM1*(cons(5,i)-Ekin-s2mu_0*SUM(cons(6:PP_nVar,i)*cons(6:PP_nVar,i))) ! includes psi^2 if PP_nVar=9
  rho_sp = cons(1,i)/p
  p_srho = p * sRho
  
  grad_in1=gradP(1,i)
  gradu  = p_sRho * (gradP(2,i) +u*gradP(5,i))
  gradv  = p_sRho * (gradP(3,i) +v*gradP(5,i))
  gradw  = p_sRho * (gradP(4,i) +w*gradP(5,i))
  grad_in5=gradP(5,i)
  
  !density gradient, rho_x = rho*w1_x + (rho/p)_x * (-p/(gamma-1) + 1/2*rho*|v|^2 )  + (rho/p)*(rho*v) . v_x
  gradP(1,i)  = cons(1,i)*grad_in1 - grad_in5*(Ekin -p*sKappaM1) + rho_sp*(cons(2,i)*gradu+cons(3,i)*gradv+cons(4,i)*gradw)
  !velocity gradient
  gradP(2,i)  = gradu
  gradP(3,i)  = gradv
  gradP(4,i)  = gradw
  !pressure gradient, =1/(rho/p)*(rho_x-p*(rho/p)_x)
  gradP(5,i)  = p_srho * (gradP(1,i) + p *grad_in5)
  !gradient of B,psi: 
  gradP(6:PP_nVar,i) = p_sRho * (gradP(6:PP_nVar,i) +cons(6:PP_nVar,i)*grad_in5)
#endif /*PP_Lifting_Var*/
END DO!i
END SUBROUTINE ConvertToGradPrimVec
#endif /*PARABOLIC*/

!==================================================================================================================================
!> calculate all wave speeds 
!==================================================================================================================================
SUBROUTINE WaveSpeeds1D(Prim,ca,cs,cf)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: prim(PP_nVar) !< vector of primitive variables, rotated!!!
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: ca,cs,cf      !< alfven, sound and fast magnetosonic wave speed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho,c2,va2,ca2
REAL                :: astar
!==================================================================================================================================
sRho=1./prim(1)
c2=kappa*prim(5)*sRho
ca2=Prim(6)*Prim(6)*sRho*smu_0 !Alfen wave speed
va2=SUM(prim(7:8)*prim(7:8))*sRho*smu_0+ca2
ca=SQRT(ca2)

astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
cs=SQRT(0.5*(c2+va2-astar))
cf=SQRT(0.5*(c2+va2+astar))
END SUBROUTINE WaveSpeeds1D


!==================================================================================================================================
!> calculate fastest wave speed 
!==================================================================================================================================
SUBROUTINE FastestWave1D(Prim,cf)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: prim(PP_nVar) !< vector of primitive variables, rotated!!!
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cf            !< fast magnetosonic wave speed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho,c2,va2,ca2
REAL                :: astar
!==================================================================================================================================
sRho=1./prim(1)
c2=kappa*prim(5)*sRho
ca2=Prim(6)*Prim(6)*sRho*smu_0 !Alfen wave speed
va2=SUM(prim(7:8)*prim(7:8))*sRho*smu_0+ca2

astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
cf=SQRT(0.5*(c2+va2+astar))
END SUBROUTINE FastestWave1D


!==================================================================================================================================
!> calculate fastest wave speed in 3D 
!==================================================================================================================================
SUBROUTINE FastestWave3D(Prim,cf)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: prim(PP_nVar) !< vector of primitive variables, rotated!!!
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cf            !< fast magnetosonic wave speed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
! sRho=1./prim(1)
! c2=kappa*prim(5)*sRho
! va2=SUM(prim(6:8)*prim(6:8))*sRho*smu_0
! astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)

! cf max in direction of ca2=0 => astar=c2+va2 => cf2=c2+va2
cf=SQRT((kappa*prim(5)+SUM(prim(6:8)*prim(6:8))*smu_0)/prim(1))
END SUBROUTINE FastestWave3D


!==================================================================================================================================
!> calculate fastest roe wave speed 
!> use Roe mean wavespeeds from  Roe meanvalues 
!>   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)
!==================================================================================================================================
SUBROUTINE FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: ConsL(PP_nVar) !< vector of left conserved variables, rotated!!!
REAL,INTENT(IN)     :: ConsR(PP_nVar) !< vector of right conserved variables, rotated!!!
REAL,INTENT(IN)     :: primL(PP_nVar) !< vector of left primitive variables, rotated!!!
REAL,INTENT(IN)     :: primR(PP_nVar) !< vector of right primitive variables, rotated!!!
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: RoeVelx        !< roe velocity 
REAL,INTENT(OUT)    :: cf_Roe         !< fastest wave
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL             :: SqrtRho_L,SqrtRho_R,sSqrtRoe_LR
REAL             :: Roe_L,Roe_R
REAL             :: ptot_L,ptot_R
REAL             :: sRho_Roe,RoeVel(3),RoeB(3)
REAL             :: H_L,H_R,RoeH
REAL             :: XX, va2_Roe,c2_Roe,ca2_Roe,astar_Roe
!==================================================================================================================================
SqrtRho_L   = SQRT(PrimL(1))
SqrtRho_R   = SQRT(PrimR(1))
sSqrtRoe_LR = 1./(SqrtRho_L+SqrtRho_R)
Roe_L       = SqrtRho_L*sSqrtRoe_LR
Roe_R       = SqrtRho_R*sSqrtRoe_LR
ptot_L      = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!
ptot_R      = PrimR(5)+s2mu_0*SUM(PrimR(6:8)**2) !Total presssure!

sRho_Roe    = 1./(SqrtRho_L*SqrtRho_R)
RoeVel(1:3) = Roe_L*PrimL(2:4)+Roe_R*PrimR(2:4)
RoeB(1:3)   = Roe_L*PrimR(6:8)+Roe_R*PrimL(6:8)

H_L  = (ConsL(5)+ptot_L)/PrimL(1)
H_R  = (ConsR(5)+ptot_R)/PrimR(1)
RoeH = Roe_L*H_L+Roe_R*H_R

XX   = 0.5*SUM((PrimL(6:8)-PrimR(6:8))**2)*(sSqrtRoe_LR**2)

va2_Roe   = SUM(RoeB(:)**2)*smu_0*sRho_Roe
c2_Roe    = (2.-Kappa)*XX+ (Kappa-1)*(RoeH-0.5*SUM(RoeVel(:)**2)-va2_Roe)
ca2_Roe   = RoeB(1)**2*smu_0*sRho_Roe
astar_Roe = SQRT((c2_Roe+va2_Roe)**2-4.*c2_Roe*ca2_Roe)
cf_Roe    = SQRT(0.5*(c2_Roe+va2_Roe+astar_Roe))
RoeVelx   = RoeVel(1)
END SUBROUTINE FastestWave1D_Roe

!===================================================================================================================================
!> Subroutine to initialize a procedure pointer to a specific indicator function
!===================================================================================================================================
subroutine SetIndicatorFunction(ind_id,IndicatorFunc)
  use MOD_Globals
  implicit none
  !-arguments---------------------------------------
  character(len=255)        , intent(in) :: ind_id
  procedure(i_indicatorFunction),pointer :: IndicatorFunc
  !-------------------------------------------------
  
  select case (trim(ind_id))
    case('Density')              ; IndicatorFunc => Get_Density
    case('Pressure')             ; IndicatorFunc => Get_Pressure
    case('DensityTimesPressure') ; IndicatorFunc => Get_DensityTimesPressure
    case('KinPlusMagEnergy')     ; IndicatorFunc => Get_KinPlusMagEnergy
    case default
      CALL abort(__STAMP__,'Indicator quantity "'//trim(ind_id)//'" is not defined for this equation!',999,999.)
      RETURN
  end select
  
end subroutine SetIndicatorFunction

!===================================================================================================================================
!> Returns the density
!===================================================================================================================================
pure subroutine Get_Density(U,rho)
  implicit none
  !-arguments---------------------------------------
  real,intent(in)  :: U(PP_nVar)
  real,intent(out) :: rho
  !-------------------------------------------------
  
  rho = U(1)
end subroutine Get_Density
!===================================================================================================================================
!> Returns the pressure
!===================================================================================================================================
pure subroutine Get_Pressure(U,p)
  implicit none
  !-arguments---------------------------------------
  real,intent(in)  :: U(PP_nVar)
  real,intent(out) :: p
  !-------------------------------------------------
  
  p = KappaM1*(U(5)-0.5*(SUM(U(2:4)*U(2:4))/U(1))-s2mu_0*SUM(U(6:PP_nVar)*U(6:PP_nVar)))
  
end subroutine Get_Pressure
!============================================================================================================================
!> Get Density Times Pressure
!============================================================================================================================
pure subroutine Get_DensityTimesPressure(U,rhop)
  implicit none
  !-arguments---------------------------------------
  real, intent(in)  :: U(PP_nVar)
  real, intent(out) :: rhop
  !-------------------------------------------------
  real :: p
  !-------------------------------------------------
  
  call Get_Pressure(U,p)
  rhop = U(1) * p
end subroutine Get_DensityTimesPressure
!============================================================================================================================
!> Get kinetic energy
!============================================================================================================================
  pure subroutine Get_KinPlusMagEnergy(U,kinen)
    implicit none
    !-arguments---------------------------------------
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: kinen
    !-------------------------------------------------
    
    kinen = 0.5*SUM(U(2:4)*U(2:4))/U(1) + s2mu_0*SUM(U(6:8)*U(6:8))
  end subroutine Get_KinPlusMagEnergy

!============================================================================================================================
!> Get dp/du
!============================================================================================================================
  pure subroutine Get_dpdU(U,dpdu)
    implicit none
    !-arguments---------------------------------------
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: dpdU(PP_nVar)
    !-local-variables---------------------------------
    real :: vel(3)
    !-------------------------------------------------
    
    vel = U(2:4)/U(1)
    
    dpdu(1)         = 0.5*sum(vel**2)
    dpdu(2:4)       = -vel
    dpdu(5)         = 1.0
    dpdu(6:PP_nVar) = -smu_0*U(6:PP_nVar)
    
    dpdu = dpdu*KappaM1
    
  end subroutine Get_dpdU
  
END MODULE MOD_Equation_Vars
