!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
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
#include "defines.h"
!==================================================================================================================================
!> Contains parameters for the Navier stokes equation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
#if PP_N == N
USE MOD_PreProc,ONLY:PP_N
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcSource        !< logical to define if a source term (e.g. exactfunc) is added
INTEGER             :: IniExactFunc        !< Exact Function for initialization
INTEGER             :: IniRefState         !< RefState for initialization (case IniExactFunc=1 only)
REAL,ALLOCATABLE    :: RefStatePrim(:,:)   !< primitive reference states
REAL,ALLOCATABLE    :: RefStateCons(:,:)   !< =primToCons(RefStatePrim)
INTEGER,PARAMETER   :: nAuxVar=6           !< number of auxiliary variables for average flux
! Boundary condition arrays
REAL,ALLOCATABLE    :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE :: nBCByType(:)        !< Number of sides for each boundary
INTEGER,ALLOCATABLE :: BCSideID(:,:)       !< SideIDs for BC types
#if PARABOLIC
REAL                :: mu0                 !< viscosity
REAL                :: Pr                  !< Prandtl number
REAL                :: KappasPr            !< Kappa/Pr
#if PP_VISC==1
REAL                :: sTref,ExpoSuth      !< Parameters used in muSuth sutherland visosity law
REAL                :: Ts,cSuth            !< Parameters used in muSuth sutherland visosity law
#endif
#if PP_VISC==2
REAL                :: ExpoPow             !< Exponent used in  power law viscosity
#endif
#endif /*PARABOLIC*/
REAL                :: Kappa               !< ratio of specific heats
REAL                :: sKappa              !< = 1/kappa
REAL                :: KappaM1             !< = kappa - 1
REAL                :: sKappaM1            !< = 1/(kappa -1)
REAL                :: KappaP1             !< = kappa + 1
REAL                :: sKappaP1            !< = 1/(kappa +1)
REAL                :: R                   !< Gas constant
REAL                :: s23
REAL                :: MachShock           !< Shoch Mach speed for ExactFunction = 6 (shock)
REAL                :: PreShockDens        !< Pre-shock density for ExactFunction = 6 (shock)
REAL                :: AdvVel(3)           !< Advection Velocity for the test cases
REAL                :: IniCenter(3)        !< for Iniexactfunc, center point
REAL                :: IniAxis(3)          !< for Iniexactfunc, center axis
REAL                :: IniWaveNumber(3)    !< for Iniexactfunc, wave numbers in xyz
REAL                :: IniFrequency        !< for Iniexactfunc, Frequeny       
REAL                :: IniAmplitude        !< for Iniexactfunc, Amplitude
REAL                :: IniHalfwidth        !< for Iniexactfunc, Halfwidth

CHARACTER(LEN=255),DIMENSION(5),PARAMETER :: StrVarNames(5)=(/ CHARACTER(LEN=255) :: 'Density',    &
                                                                                     'MomentumX',  &
                                                                                     'MomentumY',  &
                                                                                     'MomentumZ',  &
                                                                                     'EnergyStagnationDensity'/)

CHARACTER(LEN=255),DIMENSION(PP_nVar),PARAMETER :: StrVarNamesPrim(PP_nVar)=(/ CHARACTER(LEN=255) :: &
                   'Density',    &
                   'VelocityX',  &
                   'VelocityY',  &
                   'VelocityZ',  &
                   'Pressure'/)
                   
integer, parameter :: nIndVar = 4
character(len=255) :: IndicatorQuantityNames(nIndVar) = (/character(len=132) :: 'Density','Pressure','DensityTimesPressure','KinEnergy'/)
                   
LOGICAL           :: EquationInitIsDone=.FALSE. !< Init switch  
INTEGER             :: WhichRiemannSolver       !< choose riemann solver
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
PROCEDURE(i_sub_RiemannVolFluxAndDissipMatrices),POINTER :: RiemannVolFluxAndDissipMatrices =>Null()
PROCEDURE(i_sub_SolveRiemannProblem ),POINTER :: SolveRiemannProblem  =>Null() !< procedure pointer to riemann solver 
PROCEDURE(i_sub_VolumeFluxAverage   ),POINTER :: VolumeFluxAverage    =>Null() !< procedure pointer to 1D two-point average flux
PROCEDURE(i_sub_VolumeFluxAverageVec),POINTER :: VolumeFluxAverageVec =>Null() !< procedure pointer to 3D two-point average flux
!==================================================================================================================================
ABSTRACT INTERFACE
  SUBROUTINE i_sub_SolveRiemannProblem(F,U_LL,U_RR)
#if PP_N == N
    IMPORT PP_N
#endif
    REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
    REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
    REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N),INTENT(INOUT) :: F     !< numerical flux
  END SUBROUTINE i_sub_SolveRiemannProblem
  
  SUBROUTINE i_sub_RiemannVolFluxAndDissipMatrices(ConsL,ConsR,F,Dmatrix,Rmatrix)
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(IN)  :: ConsL !<  left conservative state  
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(IN)  :: ConsR !< right conservative state
    REAL,DIMENSION(1:PP_nVar)      ,INTENT(OUT) :: F     !< Volume flux
    REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: Dmatrix,Rmatrix  !< numerical flux
  END SUBROUTINE i_sub_RiemannVolFluxAndDissipMatrices
  
  PURE SUBROUTINE i_sub_VolumeFluxAverage(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< central flux in x
    REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
  END SUBROUTINE i_sub_VolumeFluxAverage

  PURE SUBROUTINE i_sub_VolumeFluxAverageVec (UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
    REAL,DIMENSION(6),INTENT(IN)        :: UauxL          !< left auxiliary variables
    REAL,DIMENSION(6),INTENT(IN)        :: UauxR          !< right auxiliary variables
    REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
    REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed central flux
  END SUBROUTINE i_sub_VolumeFluxAverageVec
  
  pure subroutine i_indicatorFunction(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_indicatorFunction
END INTERFACE

INTERFACE ConsToPrim_aux
  MODULE PROCEDURE ConsToPrim_aux
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

INTERFACE SoundSpeed2
  MODULE PROCEDURE SoundSpeed2
END INTERFACE

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

#if PP_VISC==1
INTERFACE muSuth
  MODULE PROCEDURE muSuth
END INTERFACE
#endif /*PP_VISC==1*/

#endif /*PARABOLIC*/

CONTAINS

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables, also
!> provides the sound speed and the energy and total energy 
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_aux(prim,cons)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(8) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! rho
prim(1)=cons(1)
! vel/rho
prim(2:4)=cons(2:4)*sRho
! pressure
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! Additional information
prim(6)=SQRT(Kappa*prim(5)*sRho) ! soundspeed
prim(7)=cons(5)*sRho ! e
prim(8)=prim(7)+prim(5)*sRho ! e+p/rho
END SUBROUTINE ConsToPrim_aux



!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho    !< 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! rho
prim(1)=cons(1)
! vel/rho
prim(2:4)=cons(2:4)*sRho
!pressure
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
END SUBROUTINE ConsToPrim



!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
PURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
! rho
cons(1)=prim(1)
! vel/rho
cons(2:4)=prim(2:4)*prim(1)
! inner energy
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))
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
!> calculate soundspeed^2 ,c^2 = kappa*p/rho 
!==================================================================================================================================
FUNCTION SoundSpeed2(Cons) RESULT(cs2)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(PP_nVar) !< vector of conservative variables, already rotated!
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: cs2      !< soundspeed^2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
! p*rho=(kappa-1)*(E - 1/2*rho*|v|^2)*rho = (kappa-1)*(E*rho - 1/2(rho^2*|v|^2)
cs2=kappa*kappaM1*(cons(5)*cons(1) -  0.5*SUM(cons(2:4)*cons(2:4)))/(cons(1)*cons(1)) !=(kappa*p*rho)/(rho^2)
END FUNCTION SoundSpeed2


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
rho_sp = cons(1)/(KappaM1*(cons(5)-cons(1)*v2s2))
!s      = LOG(p) - kappa*LOG(cons(1))
s      = - LOG(rho_sp*(cons(1)**kappaM1))

! Convert to entropy variables
entropy(1)   =  (kappa-s)*skappaM1 - rho_sp*v2s2  !(kappa-s)/(kappa-1)-beta*|v|^2
entropy(2)   =  rho_sp*u  ! 2*beta*v
entropy(3)   =  rho_sp*v  ! 2*beta*v
entropy(4)   =  rho_sp*w  ! 2*beta*v
entropy(5)   = -rho_sp    !-2*beta
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
REAL,INTENT(IN)     :: cons(PP_nVar,dim2) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Entropy(PP_nVar,dim2) !< vector of conservative variables
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
  rho_sp = cons(1,i)/(KappaM1*(cons(5,i)-cons(1,i)*v2s2))
  !s      = LOG(p) - kappa*LOG(cons(1,i))
  s      = - LOG(rho_sp*(cons(1,i)**kappaM1))
  
  ! Convert to entropy variables
  entropy(1,i)   =  (kappa-s)*skappaM1 - rho_sp*v2s2  !(kappa-s)/(kappa-1)-beta*|v|^2
  entropy(2,i)   =  rho_sp*u  ! 2*beta*v
  entropy(3,i)   =  rho_sp*v  ! 2*beta*v
  entropy(4,i)   =  rho_sp*w  ! 2*beta*v
  entropy(5,i)   = -rho_sp    !-2*beta
END DO!i
END SUBROUTINE ConsToEntropyVec


#if PARABOLIC
!==================================================================================================================================
!> transform gradient from conservative / primitive or entropy variables to primitive variables 
!==================================================================================================================================
PURE FUNCTION ConvertToGradPrim(cons,grad_in) RESULT(gradP)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: cons(PP_nVar)    !< conservative state 
REAL,INTENT(IN)     :: grad_in(PP_nVar) !< can be gradient of conservative / primivite /entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: gradP(PP_nVar) !<  gradient of primitive variables (rho,v1,v2,v3,p)
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
                       -(0.5*grad_in(1)*(u*u+v*v+w*w)+ (cons(2)*gradu+cons(3)*gradv+cons(4)*gradw)) ) !-grad_Ekin
#elif (PP_Lifting_Var==2) 
  !grad_in is gradient of primitive variable, do nothing
  gradP(:)=grad_in(:)
#elif (PP_Lifting_Var==3) 
  !grad_in is gradient of entropy variable,  entropy variables are:
  ! w(1)= (kappa-s)/(kappa-1)+(rho/p)/2*|v|^2  ,w(2:4)=(rho/p)*v, w(5)=-(rho/p)

  !gradient of (p/rho),  (p/rho)_x = -grad_in(5)

  sRho   = 1./cons(1)
  u      = cons(2)*sRho
  v      = cons(3)*sRho
  w      = cons(4)*sRho
  Ekin   = 0.5*(cons(2)*u+cons(3)*v+cons(4)*w)
  p      = KappaM1*(cons(5)-Ekin) 
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
                                             !<  on output: gradient of primitive variables (rho,v1,v2,v3,p)
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
                       -(0.5*grad_in1*(u*u+v*v+w*w)+ (cons(2,i)*gradu+cons(3,i)*gradv+cons(4,i)*gradw)) ) !-grad_Ekin
#elif (PP_Lifting_Var==2) 
  !grad_in is gradient of primitive variable, do nothing
  gradP(:,i)=gradP(:,i)
#elif (PP_Lifting_Var==3) 
  !grad_in is gradient of entropy variable,  entropy variables are:
  ! w(1)= (kappa-s)/(kappa-1)+(rho/p)/2*|v|^2  ,w(2:4)=(rho/p)*v, w(5)=-(rho/p)

  !gradient of (p/rho),  (p/rho)_x = -grad_in(5)

  sRho   = 1./cons(1,i)
  u      = cons(2,i)*sRho
  v      = cons(3,i)*sRho
  w      = cons(4,i)*sRho
  Ekin   = 0.5*(cons(2,i)*u+cons(3,i)*v+cons(4,i)*w)
  p      = KappaM1*(cons(5,i)-Ekin)
  rho_sp = cons(1,i)/p
  p_srho = p * sRho
  
  grad_in1=gradP(1,i)
  gradu   = p_sRho * (gradP(2,i) +u*gradP(5,i))
  gradv   = p_sRho * (gradP(3,i) +v*gradP(5,i))
  gradw   = p_sRho * (gradP(4,i) +w*gradP(5,i))
  grad_in5=gradP(5,i)
  
  !density gradient, rho_x = rho*w1_x + (rho/p)_x * (-p/(gamma-1) + 1/2*rho*|v|^2 )  + (rho/p)*(rho*v) . v_x
  gradP(1,i)  = cons(1,i)*grad_in1 - grad_in5*(Ekin -p*sKappaM1) + rho_sp*(cons(2,i)*gradu+cons(3,i)*gradv+cons(4,i)*gradw)
  !velocity gradient
  gradP(2,i)  = gradu
  gradP(3,i)  = gradv
  gradP(4,i)  = gradw
  !pressure gradient, =1/(rho/p)*(rho_x-p*(rho/p)_x)
  gradP(5,i)  = p_srho * (gradP(1,i) + p *grad_in5)
#endif /*PP_Lifting_Var*/
END DO!i
END SUBROUTINE ConvertToGradPrimVec
#endif /*PARABOLIC*/


#if PARABOLIC
#if PP_VISC==1
!==================================================================================================================================
!> Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature 
!>
!> Initialization of mu0, Ts, Tref, ExpoSuth and cSuth takes place in SUBROUTINE IniEquation 
!>
!> Temperatures above the Sutherlands Temperature Ts are computed according to (1)
!> 1) T >= Ts:    mu = mu0 * (T/Tref)^(expo) *  (Tref+TS)/(T+TS)
!> Example values would be Ts=110.4K, Tref 280K and mu0=mu(Tref)=1.735E-5Kg/ms
!> and expo = 3/2
!>
!> below Ts a linear dependence is assumed, (2)
!> 2) T < Ts:    mu = mu0*T/Tref*c
!>
!> with c = (Ts/Tref)^exp*(1+(Ts/Tref))/(2(Ts/Tref)²) for steady transition from (1) to (2) at T = Ts.
!>
!> This is only valid for Temperatures in the range 0 < T < 555 K
!> For further informration check out the HALOWiki and Babucke's Diss. and Code. 'NS3D'
!> ATTENTION!!!!! The global variable sTref=1./Tref and Ts=Ts/Tref !!!!!
!==================================================================================================================================
FUNCTION muSuth(T)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: T !< Temperature
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: muSuth !< sutherland viscosity
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: TnoDim
!==================================================================================================================================
TnoDim=T*sTref ! sTref=1./Tref !
IF(TnoDim .GE. TS)THEN  ! Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  muSuth=mu0*TnoDim**ExpoSuth*(1+Ts)/(TnoDim+Ts)  ! Ts=Ts/Tref !
ELSE
  muSuth=mu0*TnoDim*cSuth
END IF
END FUNCTION muSuth
#endif
#endif /*PARABOLIC*/

!==================================================================================================================================
!> calculate the sound speed
!==================================================================================================================================
SUBROUTINE FastestWave3D(Prim,c)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: prim(PP_nVar) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: c             !< sound speed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
c=SQRT(kappa*prim(5)/prim(1))
END SUBROUTINE FastestWave3D

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
    case('KinEnergy')            ; IndicatorFunc => Get_KinEnergy
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
  
  p = KappaM1*(U(5)-0.5*(SUM(U(2:4)*U(2:4))/U(1)))
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
  pure subroutine Get_KinEnergy(U,kinen)
    implicit none
    !-arguments---------------------------------------
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: kinen
    !-------------------------------------------------
    
    kinen = SUM(U(2:4)*U(2:4))/U(1)
  end subroutine Get_KinEnergy

END MODULE MOD_Equation_Vars
