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

!==================================================================================================================================
!> Contains parameters for the Navier stokes equation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
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
LOGICAL           :: EquationInitIsDone=.FALSE. !< Init switch  
INTEGER             :: WhichRiemannSolver       !< choose riemann solver
PROCEDURE(),POINTER :: SolveRiemannProblem      !< procedure pointer to riemann solver 
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
PROCEDURE(),POINTER :: VolumeFluxAverage        !< procedure pointer to 1D two-point average flux
PROCEDURE(),POINTER :: VolumeFluxAverageVec     !< procedure pointer to 3D two-point average flux
!==================================================================================================================================

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
REAL                                :: v(3),v2,rho_sp,p,s
!==================================================================================================================================
v(:)   = cons(2:4)/cons(1)
v2     = SUM(v*v)
p      = KappaM1*(cons(5)-0.5*cons(1)*v2)
s      = LOG(p) - kappa*LOG(cons(1))
rho_sp = cons(1)/p   !=2*beta

! Convert to entropy variables
entropy(1)   =  (kappa-s)*skappaM1 - 0.5*rho_sp*v2  !(kappa-s)/(kappa-1)-beta*|v|^2
entropy(2:4) =  rho_sp*v(:)  ! 2*beta*v
entropy(5)   = -rho_sp       !-2*beta
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
                       -(0.5*grad_in(1)*SUM(v*v) + SUM(cons(2:4)*gradv)) ) !-grad_Ekin
#elif (PP_Lifting_Var==2) 
  !grad_in is gradient of primitive variable, do nothing
  gradP(:)=grad_in(:)
#elif (PP_Lifting_Var==3) 
  !grad_in is gradient of entropy variable,  entropy variables are:
  ! w(1)= (kappa-s)/(kappa-1)+(rho/p)/2*|v|^2  ,w(2:4)=(rho/p)*v, w(5)=-(rho/p)

  !gradient of (p/rho),  (p/rho)_x = -grad_in(5)

  sRho   = 1./cons(1)
  v(:)   = cons(2:4)*sRho
  Ekin   = 0.5*SUM(cons(2:4)*v)
  p      = KappaM1*(cons(5)-Ekin) ! includes psi^2 if PP_nVar=9
  rho_sp = cons(1)/p
  p_srho = p * sRho
  
  gradv  = p_sRho * (grad_in(2:4) +v(:)*grad_in(5))
  
  !density gradient, rho_x = rho*w1_x + (rho/p)_x * (-p/(gamma-1) + 1/2*rho*|v|^2 )  + (rho/p)*(rho*v) . v_x
  gradP(1)  = cons(1)*grad_in(1) - grad_in(5)*(Ekin -p*sKappaM1) + rho_sp*SUM(cons(2:4)*gradv(:))
  !velocity gradient
  gradP(2:4)= gradv(:)
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

END MODULE MOD_Equation_Vars
