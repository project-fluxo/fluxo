!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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

LOGICAL           :: EquationInitIsDone=.FALSE. !< Init switch  
!procedure pointers
INTEGER             :: WhichRiemannSolver       !< choice of riemann solver
PROCEDURE(),POINTER :: SolveRiemannProblem      !< pointer to riemann solver routine (depends on WhichRiemannSolver)
#if (PP_DiscType==2)
!procedure pointers for split form DG
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
PROCEDURE(),POINTER :: VolumeFluxAverageVec     !< procedure pointer to two-point average flux
#endif /*PP_DiscType==2*/
!==================================================================================================================================

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
FUNCTION ConsToEntropy(cons) RESULT(Entropy)
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
REAL                                :: p,v(3),v2,rho_sp,s
!==================================================================================================================================
v(:)   = cons(2:4)/cons(1)
v2     = SUM(v*v)
p      = KappaM1*(cons(5)-0.5*cons(1)*v2-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar))) ! includes psi^2 if PP_nVar=9
s      = LOG(p) - kappa*LOG(cons(1))
rho_sp = cons(1)/p

! Convert to entropy variables
entropy(1)         =  (kappa-s)*skappaM1 - 0.5*rho_sp*v2  !(kappa-s)/(kappa-1)-beta*|v|^2
entropy(2:4)       =  rho_sp*v(:)                  ! 2*beta*v
entropy(5)         = -rho_sp                       !-2*beta
entropy(6:PP_nVar) =  rho_sp*cons(6:PP_nVar)       ! 2*beta*B +2*beta*psi

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
!==================================================================================================================================
DO i=1,dim2
  Entropy(:,i)=ConsToEntropy(Cons(:,i))
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
REAL  :: sRho,v(3),gradv(3)
#elif (PP_Lifting_Var==3) 
REAL  :: sRho,v(3),gradv(3),Ekin,p,rho_sp,p_srho
#endif /*PP_Lifting_Var*/
!==================================================================================================================================

#if (PP_Lifting_Var==1) 
  !grad_in is gradient of conservative variable
  sRho      = 1./cons(1)
      v(:)  = cons(2:4)*sRho
  gradv(:)  = sRho*(grad_in(2:4)-grad_in(1)*v(:))
  
  !density gradient
  gradP(1)  = grad_in(1)
  !velocity gradient
  gradP(2:4)= gradv(:)
  !pressure gradient
  gradP(5)  = KappaM1*(grad_in(5)                                        & !gradE
                       -(0.5*grad_in(1)*SUM(v*v) + SUM(cons(2:4)*gradv)) & !-grad_Ekin
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
  v(:)   = cons(2:4)*sRho
  Ekin   = 0.5*SUM(cons(2:4)*v)
  p      = KappaM1*(cons(5)-Ekin-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar))) ! includes psi^2 if PP_nVar=9
  rho_sp = cons(1)/p
  p_srho = p * sRho
  
  gradv  = p_sRho * (grad_in(2:4) +v(:)*grad_in(5))
  
  !density gradient, rho_x = rho*w1_x + (rho/p)_x * (-p/(gamma-1) + 1/2*rho*|v|^2 )  + (rho/p)*(rho*v) . v_x
  gradP(1)  = cons(1)*grad_in(1) - grad_in(5)*(Ekin -p*sKappaM1) + rho_sp*SUM(cons(2:4)*gradv(:))
  !velocity gradient
  gradP(2:4)= gradv(:)
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
!==================================================================================================================================
DO i=1,dim2
  gradP(:,i)=ConvertToGradPrim(cons(:,i),gradP(:,i))
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


END MODULE MOD_Equation_Vars
