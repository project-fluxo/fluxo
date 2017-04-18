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
INTEGER             :: IniExactFunc        !< Exact Function for initialization
INTEGER             :: IniRefState         !< RefState for initialization (case IniExactFunc=1 only)
REAL,ALLOCATABLE    :: RefStatePrim(:,:)
REAL,ALLOCATABLE    :: RefStateCons(:,:)
! Boundary condition arrays
REAL,ALLOCATABLE    :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE :: nBCByType(:)        !< Number of sides for each boundary
INTEGER,ALLOCATABLE :: BCSideID(:,:)       !< SideIDs for BC types
#ifdef PARABOLIC
REAL                :: mu0                 !< viscosity
REAL                :: Pr                  !< Prandtl number
REAL                :: KappasPr            !< Kappa/Pr
#if PP_VISC==1
REAL                :: Ts,cSuth            !< Parameters used in muSuth
#endif
#if (PP_VISC==1) || (PP_VISC==2)
REAL                :: Tref,ExpoSuth       !< Parameters used in muSuth and power law
#endif
#endif /*PARABOLIC*/
REAL                :: Kappa               !< ratio of specific heats
REAL                :: sKappa              !< = 1/kappa
REAL                :: KappaM1             !< = kappa - 1
REAL                :: sKappaM1            !< = 1/(kappa -1)
REAL                :: KappaP1             !< = kappa + 1
REAL                :: sKappaP1            !< = 1/(kappa +1)
REAL                :: R                   !< Gas constant
REAL                :: s43,s23
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
#if (PP_DiscType==2)
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
INTEGER,PARAMETER   :: nAuxVar=6 
PROCEDURE(),POINTER :: VolumeFluxAverage        !< procedure pointer to 1D two-point average flux
PROCEDURE(),POINTER :: VolumeFluxAverageVec     !< procedure pointer to 3D two-point average flux
#endif /*PP_DiscType==2*/
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

INTERFACE ConsToEntropy
  MODULE PROCEDURE ConsToEntropy
END INTERFACE

#ifdef PARABOLIC
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
REAL,INTENT(INOUT)  :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(8) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho    ! 1/Rho
!==================================================================================================================================
cons(1)=MAX(0.0000001,abs(cons(1)))
sRho=1./cons(1)
! conversion
prim(1)=cons(1)
! rho
prim(2:4)=cons(2:4)*sRho
! vel/rho
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! pressure
! pressure must not be negative
prim(5)=MAX(0.0000001,abs(prim(5)))
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
REAL,INTENT(INOUT)  :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: sRho    !< 1/Rho
!==================================================================================================================================
cons(1)=MAX(0.0000001,abs(cons(1)))
sRho=1./cons(1)
! conversion
prim(1)=cons(1)
! rho
prim(2:4)=cons(2:4)*sRho
! vel/rho
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! pressure
! pressure must not be negative
prim(5)=MAX(0.0000001,abs(prim(5)))
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
REAL,INTENT(INOUT)  :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
prim(1)=MAX(0.0000001,abs(prim(1)))
! conversion
cons(1)=prim(1)
! rho
cons(2:4)=prim(2:4)*prim(1)
! vel/rho
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))
! inner energy
END SUBROUTINE PrimToCons


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables a la Ismail and Roe
!==================================================================================================================================
SUBROUTINE ConsToEntropy(entropy,cons)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(5),INTENT(IN)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(5),INTENT(OUT) :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: sRho ! 1/rho
REAL                                :: s,p,u,v,w,rho,rho_p
!==================================================================================================================================
! Pull the conservative variables and extra constants

sRho  = 1./cons(1)
rho   = cons(1)
u     = cons(2)*sRho
v     = cons(3)*sRho
w     = cons(4)*sRho
p     = KappaM1*(cons(5) - 0.5*rho*(u*u + v*v + w*w))
s     = LOG(p) - kappa*LOG(rho)
rho_p = rho/p

! Convert to entropy variables
entropy(1) =  (kappa-s)/kappaM1 - 0.5*rho_p*(u*u + v*v * w*w)
entropy(2) =  rho_p*u
entropy(3) =  rho_p*v
entropy(4) =  rho_p*w
entropy(5) = -rho_p
END SUBROUTINE ConsToEntropy


#ifdef PARABOLIC
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
!> ATTENTION!!!!! The global variable Tref=1./Tref and Ts=Ts/Tref !!!!!
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
TnoDim=T*Tref ! Tref=1./Tref !
IF(TnoDim .GE. TS)THEN  ! Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  muSuth=mu0*TnoDim**ExpoSuth*(1+Ts)/(TnoDim+Ts)  ! Ts=Ts/Tref !
ELSE
  muSuth=mu0*TnoDim*cSuth
END IF
END FUNCTION muSuth
#endif
#endif /*PARABOLIC*/

END MODULE MOD_Equation_Vars
