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
#include "defines.h"

!==================================================================================================================================
!> contains riemann solvers for MHD
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE RiemannSolverByHLL
  MODULE PROCEDURE RiemannSolverByHLL
END INTERFACE

INTERFACE RiemannSolverByHLLC
  MODULE PROCEDURE RiemannSolverByHLLC
END INTERFACE

INTERFACE RiemannSolverByHLLD
  MODULE PROCEDURE RiemannSolverByHLLD
END INTERFACE

INTERFACE RiemannSolverByRoe
  MODULE PROCEDURE RiemannSolverByRoe
END INTERFACE

INTERFACE RiemannSolverByRusanov
  MODULE PROCEDURE RiemannSolverByRusanov
END INTERFACE


PUBLIC :: Riemann
PUBLIC :: RiemannSolverByHLL
PUBLIC :: RiemannSolverByHLLC
PUBLIC :: RiemannSolverByHLLD
PUBLIC :: RiemannSolverByRoe
PUBLIC :: RiemannSolverByRusanov
!==================================================================================================================================


CONTAINS
!==================================================================================================================================
!> main routine calling different riemann solvers
!==================================================================================================================================
SUBROUTINE Riemann(F,UL,UR,                                                                       &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R,                         &
#endif
                   nv,t1,t2)
USE MOD_PreProc
USE MOD_Equation_vars,ONLY: SolveRiemannProblem
#if PARABOLIC
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: UL(      PP_nVar,0:PP_N,0:PP_N) !<  left state on face
REAL,INTENT(IN) :: UR(      PP_nVar,0:PP_N,0:PP_N) !< right state on face
#if PARABOLIC                                                 
REAL,INTENT(IN) :: gradUx_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in x 
REAL,INTENT(IN) :: gradUx_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in x 
REAL,INTENT(IN) :: gradUy_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in y 
REAL,INTENT(IN) :: gradUy_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in y 
REAL,INTENT(IN) :: gradUz_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in z 
REAL,INTENT(IN) :: gradUz_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in z 
#endif /*PARABOLIC*/
REAL,INTENT(IN) :: nv(            3,0:PP_N,0:PP_N) !< normal vector of face
REAL,INTENT(IN) :: t1(            3,0:PP_N,0:PP_N) !< 1st tangential vector of face
REAL,INTENT(IN) :: t2(            3,0:PP_N,0:PP_N) !< 2nd tangential vector of face
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT):: F(       PP_nVar,0:PP_N,0:PP_N) !< numerical flux on face
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: ConsL(1:PP_nVar)
REAL             :: ConsR(1:PP_nVar)
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N) :: k_L,g_L,j_L,k_R,g_R,j_R
INTEGER          :: iVar
#endif /*PARABOLIC*/
INTEGER          :: i,j
!==================================================================================================================================
DO j = 0,PP_N
  DO i = 0,PP_N
    !rotate fields
    ConsL(:) =     UL(:  ,i,j)
    ConsL(2) = SUM(UL(2:4,i,j)*nv(:,i,j))
    ConsL(3) = SUM(UL(2:4,i,j)*t1(:,i,j))
    ConsL(4) = SUM(UL(2:4,i,j)*t2(:,i,j))
    ConsL(6) = SUM(UL(6:8,i,j)*nv(:,i,j))
    ConsL(7) = SUM(UL(6:8,i,j)*t1(:,i,j))
    ConsL(8) = SUM(UL(6:8,i,j)*t2(:,i,j))

    ConsR(:) =     UR(:,i,j)
    ConsR(2) = SUM(UR(2:4,i,j)*nv(:,i,j))
    ConsR(3) = SUM(UR(2:4,i,j)*t1(:,i,j))
    ConsR(4) = SUM(UR(2:4,i,j)*t2(:,i,j))
    ConsR(6) = SUM(UR(6:8,i,j)*nv(:,i,j))
    ConsR(7) = SUM(UR(6:8,i,j)*t1(:,i,j))
    ConsR(8) = SUM(UR(6:8,i,j)*t2(:,i,j))


    CALL SolveRiemannProblem(ConsL(:),ConsR(:),F(:,i,j))

    ! Back Rotate the normal flux into Cartesian direction
    F(2:4,i,j) =   nv(:,i,j)*F(2,i,j) &
                 + t1(:,i,j)*F(3,i,j) &
                 + t2(:,i,j)*F(4,i,j)
    F(6:8,i,j) =   nv(:,i,j)*F(6,i,j) &
                 + t1(:,i,j)*F(7,i,j) &
                 + t2(:,i,j)*F(8,i,j)

  END DO !i
END DO !j
#if PARABOLIC
!! Don#t forget the diffusion contribution, my young padawan
!! Compute MHD Diffusion flux
CALL EvalDiffFlux3D(k_L,g_L,j_L,UL,gradUx_L,gradUy_L,gradUz_L)
CALL EvalDiffFlux3D(k_R,g_R,j_R,UR,gradUx_R,gradUy_R,gradUz_R)
!
! !BR1/BR2 uses arithmetic mean of the fluxes
DO iVar=2,PP_nVar
  F(iVar,:,:)=F(iVar,:,:)+0.5*( nv(1,:,:)*(k_L(iVar,:,:)+k_R(iVar,:,:))  & 
                               +nv(2,:,:)*(g_L(iVar,:,:)+g_R(iVar,:,:))  &
                               +nv(3,:,:)*(j_L(iVar,:,:)+j_R(iVar,:,:))  )
END DO
#endif /* PARABOLIC */

END SUBROUTINE Riemann



!==================================================================================================================================
!> Lax-friedrichs / rusanov flux
!==================================================================================================================================
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,Flux)
USE MOD_PreProc
USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
USE MOD_Equation_vars, ONLY: FastestWave1D
USE MOD_Equation_vars, ONLY: ConsToPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: ConsL(1:PP_nVar) !<  left conservative state  
REAL,INTENT(IN)  :: ConsR(1:PP_nVar) !< right conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar) !<numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
!----------------------------------------------------------------------------------------------------------------------------------
REAL             :: FluxL(1:PP_nVar), FluxR(1:PP_nVar)
REAL             :: PrimL(1:PP_nVar), PrimR(1:PP_nVar) !< left/right primitive state
REAL             :: LambdaMax
REAL             :: cf_L 
REAL             :: cf_R 
!==================================================================================================================================

CALL EvalAdvectionFlux1D(ConsL(1:PP_nVar),FluxL(1:PP_nVar))
CALL EvalAdvectionFlux1D(ConsR(1:PP_nVar),FluxR(1:PP_nVar))
CALL ConsToPrim(PrimL(:),ConsL(:))
CALL ConsToPrim(PrimR(:),ConsR(:))
CALL FastestWave1D(PrimL(1:PP_nVar),cf_L)
CALL FastestWave1D(PrimR(1:PP_nVar),cf_R)

LambdaMax = MAX(ABS(PrimL(2))+cf_L,ABS(PrimR(2))+cf_R)

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

END SUBROUTINE RiemannSolverByRusanov


!==================================================================================================================================
!> HLL solver following the paper of Shentai Li: "An HLLC Riemann Solver for Magnetohydrodynamics"
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLL(ConsL,ConsR,Flux)
USE MOD_PreProc
USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
USE MOD_Equation_vars, ONLY: FastestWave1D
USE MOD_Equation_vars, ONLY: FastestWave1D_Roe
USE MOD_Equation_vars, ONLY: ConsToPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: ConsL(1:PP_nVar) !<  left conservative state  
REAL,INTENT(IN)  :: ConsR(1:PP_nVar) !< right conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar) !<numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
!----------------------------------------------------------------------------------------------------------------------------------
REAL             :: FluxL(1:PP_nVar), FluxR(1:PP_nVar)
REAL             :: PrimL(1:PP_nVar), PrimR(1:PP_nVar) !< left/right primitive state
REAL             :: SL, SR
REAL             :: cf_L,cf_R,cf_Roe,RoeVelx 
!==================================================================================================================================

CALL ConsToPrim(PrimL(:),ConsL(:))
CALL ConsToPrim(PrimR(:),ConsR(:))

CALL FastestWave1D(PrimL,cf_L)
CALL FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)
SL = MIN(PrimL(2)-cf_L, RoeVelx-cf_Roe   )


IF (SL .GT. 0.0)  THEN
  CALL EvalAdvectionFlux1D(ConsL,Flux)
  RETURN
END IF

CALL FastestWave1D(PrimR,cf_R)
SR = MAX(RoeVelx+cf_Roe, PrimR(2)+cf_R   )

IF (SR .LT. 0.0) THEN
  CALL EvalAdvectionFlux1D(ConsR,Flux)
  RETURN
END IF

! NOW SL<0 and SR>0
CALL EvalAdvectionFlux1D(ConsL,FluxL)
CALL EvalAdvectionFlux1D(ConsR,FluxR)
Flux = (SR*FluxL - SL*FluxR + SL*SR*(ConsR-ConsL))/(SR-SL)

END SUBROUTINE RiemannSolverByHLL



!==================================================================================================================================
!> HLLC solver following the paper of Shentai Li: "An HLLC Riemann Solver for Magnetohydrodynamics"
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLLC(ConsL,ConsR,Flux)
USE MOD_PreProc
USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
USE MOD_Equation_vars, ONLY: FastestWave1D
USE MOD_Equation_vars, ONLY: FastestWave1D_Roe
USE MOD_Equation_vars, ONLY: smu_0,s2mu_0
USE MOD_Equation_vars, ONLY: ConsToPrim
#ifdef PP_GLM
USE MOD_Equation_vars, ONLY: GLM_ch
#endif /*PP_GLM*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: ConsL(1:PP_nVar) !<  left conservative state  
REAL,INTENT(IN)  :: ConsR(1:PP_nVar) !< right conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar) !<numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL             :: FluxL(PP_nVar), FluxR(PP_nVar)
REAL             :: PrimL(1:PP_nVar), PrimR(1:PP_nVar) !< left/right primitive state
REAL             :: U_HLL(PP_nVar)
REAL             :: U_star_L(PP_nVar), U_star_R(PP_nVar)
REAL             :: SL, SR, SM, s_SL_SM, s_SR_SM
REAL             :: ptot_L, p_star
REAL             :: rho_L, rho_R
REAL             :: srho_HLL
REAL             :: Bx_star,By_star,Bz_star,E_star
REAL             :: vx_L, vx_R
REAL             :: Bx_L, Bx_R
REAL             :: cf_L,cf_R,cf_Roe,RoeVelx 
!==================================================================================================================================
CALL ConsToPrim(PrimL(:),ConsL(:))
CALL ConsToPrim(PrimR(:),ConsR(:))


!--------------------------------------------------------------------------------
!use Roe mean wavespeeds from  Roe meanvalues 
!   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)
CALL FastestWave1D(PrimL,cf_L)
CALL FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)



SL = MIN(PrimL(2)-cf_L, RoeVelx-cf_Roe   )


IF (SL .GT. 0.0)  THEN
  CALL EvalAdvectionFlux1D(ConsL,Flux)
  RETURN
END IF

CALL FastestWave1D(PrimR,cf_R)
SR = MAX(RoeVelx+cf_Roe, PrimR(2)+cf_R   )

IF (SR .LT. 0.0) THEN
  CALL EvalAdvectionFlux1D(ConsR,Flux)
  RETURN
END IF

! NOW SL<0 and SR>0
CALL EvalAdvectionFlux1D(ConsL,FluxL)
CALL EvalAdvectionFlux1D(ConsR,FluxR)
CALL EvalHLLState(ConsL,ConsR,SL,SR,FluxL,FluxR,U_HLL)

rho_L = PrimL(1)
vx_L  = PrimL(2)
Bx_L  = PrimL(6) 

rho_R = PrimR(1)
vx_R  = PrimR(2)
Bx_R  = PrimR(6)

ptot_L   = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!

!SM = ((rho_L*vx_L*(vx_L-SL)+ptot_L-smu_0*Bx_L*Bx_L)-(rho_R*vx_R*(vx_R-SR)+ptot_R-smu_0*Bx_R*Bx_R)) &
!     /(rho_L*(vx_L-SL)-rho_R*(vx_R-SR))
sRho_HLL = 1./U_HLL(1)
SM       = U_HLL(2)*sRho_HLL

Bx_star  = U_HLL(6)
By_star  = U_HLL(7)
Bz_star  = U_HLL(8)
p_star   = ptot_L + rho_L*(vx_L-SL)*(vx_L-SM)-s2mu_0*(Bx_L*Bx_L-Bx_star*Bx_star)
E_star   = p_star*SM + smu_0*(-Bx_star*SUM(U_HLL(6:8)*U_HLL(2:4))*sRho_HLL  &
#ifdef PP_GLM
                              +U_HLL(9)*(GLM_ch*Bx_star-0.5*U_HLL(9)*SM  )&
#endif /* PP_GLM */
                             )


!-------------------------------------------------


IF ((SL .LE. 0.0) .AND. (0.0 .LT. SM)) THEN
  s_SL_SM=1./(SL-SM)
  U_star_L(1) = rho_L*(SL-vx_L)*s_SL_SM
  U_star_L(2) = U_star_L(1)*SM
  U_star_L(3) = (ConsL(3)*SL-FluxL(3) -smu_0*Bx_star*By_star)*s_SL_SM
  U_star_L(4) = (ConsL(4)*SL-FluxL(4) -smu_0*Bx_star*Bz_star)*s_SL_SM
  U_star_L(5) = (E_star+ ConsL(5)*SL-FluxL(5) )*s_SL_SM
  U_star_L(6) = Bx_star 
  U_star_L(7) = By_star
  U_star_L(8) = Bz_star
#ifdef PP_GLM
  U_star_L(9) = ConsL(9)
#endif /* PP_GLM */
  Flux = FluxL+SL*(U_star_L-ConsL)
ELSE
  s_SR_SM=1./(SR-SM)
  U_star_R(1)  = rho_R*(SR-vx_R)*s_SR_SM
  U_star_R(2) = U_star_R(1)*SM
  U_star_R(3) = (ConsR(3)*SR-FluxR(3) -smu_0*Bx_star*By_star)*s_SR_SM
  U_star_R(4) = (ConsR(4)*SR-FluxR(4) -smu_0*Bx_star*Bz_star)*s_SR_SM
  U_star_R(5) = (E_star+ ConsR(5)*SR-FluxR(5) )*s_SR_SM
  U_star_R(6) = Bx_star 
  U_star_R(7) = By_star
  U_star_R(8) = Bz_star
#ifdef PP_GLM
  U_star_R(9) = ConsR(9)!0.0
#endif /* PP_GLM */
  Flux = FluxR+SR*(U_star_R-ConsR)
END IF


END SUBROUTINE RiemannSolverByHLLC



!==================================================================================================================================
!> HLLD solver following "A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics"
!> Takahiro Miyoshi  Kanya Kusano, JCP 208 (2005). Follow implementation of Elise Estibals
!> Input state already rotated to normal system, and 
!> ONLY WORKS FOR mu_0=1!!!
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLLD(ConsL_in,ConsR_in,Flux)
USE MOD_PreProc
USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
USE MOD_Equation_vars, ONLY: FastestWave1D
USE MOD_Equation_vars, ONLY: PrimToCons
USE MOD_Equation_vars, ONLY: ConsToPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: ConsL_in(1:PP_nVar) !<  left conservative state ,not used here 
REAL,INTENT(IN)  :: ConsR_in(1:PP_nVar) !< right conservative state ,not used here
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar) !<numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL                 :: ConsL(1:PP_nVar), ConsR(1:PP_nVar)
REAL                 :: PrimL(1:PP_nVar), PrimR(1:PP_nVar)
REAL                 :: FluxL(PP_nVar), FluxR(PP_nVar)
REAL                 :: SL, SR, SM
REAL                 :: SLS, SRS
REAL                 :: FastestWave 
REAL                 :: sSL_SM,sSR_SM
REAL                 :: cf_L,cf_R
REAL                 :: temp,stemp
REAL                 :: signBn,Bn,ptotS 
REAL                 :: rhoLS,sqrtRhoLS,rhoLSS,uB,uBS
REAL                 :: rhoRS,sqrtRhoRS,rhoRSS
REAL                 :: pmagL,ptotL,eLS,eLSS !,rhoL,pL,unL,eL
REAL                 :: pmagR,ptotR,eRS,eRSS !,rhoR,pR,unR,eR
REAL,DIMENSION(1:3)  :: uLS,utL,utLS  !,uL
REAL,DIMENSION(1:3)  :: uRS,utR,utRS  !,uR
REAL,DIMENSION(1:3)  :: BLS,BtL,BtLS  !,BL
REAL,DIMENSION(1:3)  :: BRS,BtR,BtRS  !,BR
REAL,DIMENSION(1:3)  :: uSS,BSS,utSS,BtSS
REAL,DIMENSION(1:PP_nVar)  :: U_LS, U_LSS
REAL,DIMENSION(1:PP_nVar)  :: U_RS, U_RSS
REAL, PARAMETER      :: hlld_eps = 1.0e-8
REAL, PARAMETER      :: hlld_small_eps = 1.0e-12
!==================================================================================================================================
CALL ConsToPrim(PrimL(:),ConsL_in(:))
CALL ConsToPrim(PrimR(:),ConsR_in(:))
! make B_n continuous 
Bn=0.5*(PrimL(6)+PrimR(6))
PrimR(6)=Bn

CALL PrimToCons(PrimL(:),ConsL(:))
CALL PrimToCons(PrimR(:),ConsR(:))
!--------------------------------------------------------------------------------
CALL FastestWave1D(PrimL,cf_L)
CALL FastestWave1D(PrimR,cf_R)

FastestWave=MAX(cf_L,cf_R)

SL = MIN(PrimL(2),PrimR(2))-FastestWave


IF (SL .GT. 0.0)  THEN
  CALL EvalAdvectionFlux1D(ConsL,Flux)
  RETURN
END IF

SR = MAX(PrimL(2),PrimR(2))+FastestWave

IF (SR .LT. 0.0) THEN
  CALL EvalAdvectionFlux1D(ConsR,Flux)
  RETURN
END IF

! NOW SL<0 and SR>0
CALL EvalAdvectionFlux1D(ConsL,FluxL)
CALL EvalAdvectionFlux1D(ConsR,FluxR)

!-------------------------------------------------

!Left and Right variables
ASSOCIATE( rhoL => PrimL(1)          , rhoR => PrimR(1)          ,  &  
           unL  => PrimL(2)          , unR  => PrimR(2)          ,  & 
           pL   => PrimL(5)          , pR   => PrimR(5)          ,  & 
           eL   => ConsL(5)          , eR   => ConsR(5)          ,  & 
           uL   => PrimL(2:4)        , uR   => PrimR(2:4)        ,  & 
           BL   => PrimL(6:8)        , BR   => PrimR(6:8)           )

utL(1:3)  = (/0.,uL(2),uL(3)/) ; utR(1:3)  = (/0.,uR(2),uR(3)/) 
BtL(1:3)  = (/0.,BL(2),BL(3)/) ; BtR(1:3)  = (/0.,BR(2),BR(3)/) 

!magnetic and total pressures
pmagL = 0.5*(BL(1)**2 + BL(2)**2 + BL(3)**2) ; ptotL = pL + pmagL
pmagR = 0.5*(BR(1)**2 + BR(2)**2 + BR(3)**2) ; ptotR = pR + pmagR

!Compute middle slope
SM  = (rhoR*unR*(SR-unR) - rhoL * unL*(SL-unL) - ptotR + ptotL) / &
     & (rhoR*(SR-unR) - rhoL*(SL-unL))

sSL_SM=1./(SL-SM)
sSR_SM=1./(SR-SM)
rhoLS = rhoL * (SL-unL) * sSL_SM
rhoRS = rhoR * (SR-unR) * sSR_SM 

sqrtRhoLS=SQRT(rhoLS)
sqrtRhoRS=SQRT(rhoRS)

SLS = SM - ABS(Bn) / SqrtRhoLS
SRS = SM + ABS(Bn) / SqrtRhoRS

!Compute intermediate states
ptotS = ptotL + rhoL * (SL - unL) * (SM - unL)

!Left * state
temp = rhoL * (SL-unL) * (SL - SM) - Bn*Bn
IF (ABS(temp)<hlld_eps) THEN
  utLS = utL
  BtLS = BtL
ELSE
  stemp=1./temp
  utLS = utL - (Bn*(SM-unL)*stemp)*BtL
  BtLS = ((rhoL*(SL-unL)*(SL-unL)-Bn*Bn)*stemp)* BtL
END IF
BLS = (/Bn,0.,0./) + BtLS
uLS = (/SM,0.,0./) + utLS

uBS = uLS(1)*BLS(1) + uLS(2)*BLS(2) + uLS(3)*BLS(3)
uB  = uL(1) *BL(1)  + uL(2) *BL(2)  + uL(3) *BL(3)

eLS = ( eL*(SL-unL) +  ptotS*SM - ptotL*unL + Bn*(uB-uBS) ) *sSL_SM

U_LS(1)   = rhoLS       
U_LS(2:4) = rhoLS * uLS 
U_LS(5)   = eLS         
U_LS(6:8) = BLS         
#ifdef PP_GLM
U_LS(9) = ConsL(9)
#endif /* PP_GLM */

!Right * state
temp = rhoR * (SR-unR) * (SR - SM) - Bn*Bn
IF (ABS(temp)<hlld_eps) THEN
   utRS = utR
   BtRS = BtR
ELSE
   stemp=1./temp
   utRS = utR - (Bn*(SM-unR)*stemp)*BtR
   BtRS = ((rhoR*(SR-unR)*(SR-unR)-Bn*Bn)*stemp) * BtR
END IF
BRS = (/Bn,0.,0./) + BtRS
uRS = (/SM,0.,0./) + utRS

uBS = uRS(1)*BRS(1) + uRS(2)*BRS(2) + uRS(3)*BRS(3)
uB  = uR(1) *BR(1)  + uR(2) *BR(2)  + uR(3) *BR(3)

eRS = (eR*(SR-unR) - ptotR*unR + ptotS*SM  + Bn*(uB-uBS) ) *sSR_SM

U_RS(1)   = rhoRS
U_RS(2:4) = rhoRS * uRS
U_RS(5)   = eRS
U_RS(6:8) = BRS
#ifdef PP_GLM
U_RS(9) = ConsR(9)
#endif /* PP_GLM */

IF (ABS(Bn)<hlld_small_eps) THEN
   !If Bn=0 four states reduce to two
   U_LSS = U_LS
   U_RSS = U_RS
ELSE
   rhoLSS = rhoLS
   rhoRSS = rhoRS

   temp = SqrtRhoLS + SqrtRhoRS

   IF (Bn<0.0) THEN
      signBn = -1.0d0
   ELSE
      signBn =  1.0d0
   END IF
   stemp=1./temp 
   utSS = (SqrtRhoLS*utLS + SqrtRhoRS*utRS + signBn * (BtRS - BtLS)) * stemp
   BtSS = (SqrtRhoLS*BtRS + SqrtRhoRS*BtLS + signBn * (utRS - utLS) * sqrtRhoLS*sqrtRhoRS) *stemp
   
   uSS  = (/SM,0.,0./) + utSS
   BSS  = (/Bn,0.,0./) + BtSS
   uBS  = utSS(1)*BtSS(1) + utSS(2)*BtSS(2) + utSS(3)*BtSS(3)
   
   uB   = utLS(1)*BtLS(1) + utLS(2)*BtLS(2) + utLS(3)*BLS(3)
   eLSS = eLS - signBn * sqrtRhoLS * (uB - uBS)
  
   uB   = uRS(1)*BtRS(1) + utRS(2)*BtRS(2) + utRS(3)*BtRS(3)
   eRSS = eRS + signBn * sqrtRhoRS * (uB - uBS)

   U_LSS(1)   = rhoLSS          ; U_RSS(1)   = rhoRSS
   U_LSS(2:4) = rhoLSS*uSS      ; U_RSS(2:4) = rhoRSS*uSS
   U_LSS(5)   = eLSS            ; U_RSS(5)   = eRSS
   U_LSS(6:8) = BSS             ; U_RSS(6:8) = BSS
#ifdef PP_GLM
   U_LSS(9) = ConsL(9)          ; U_RSS(9)   = consR(9)
#endif /* PP_GLM */
END IF !|Bn|<0

END ASSOCIATE !rhoL,rhoR,...

!Fluxes
IF (0.0<=SLS) THEN
  !U=ULS
  Flux(:) = FluxL(:) + SL*(U_LS(:) - ConsL(:))
ELSE IF (0.0<=SM) THEN
  !U=ULSS
  Flux(:) = FluxL(:) + SLS*U_LSS(:) - (SLS-SL)*U_LS(:) - SL*ConsL(:)
ELSE IF (0.0<=SRS) THEN
  !U=URSS
  Flux(:) = FluxR(:) + SRS*U_RSS(:) - (SRS-SR)*U_RS(:) - SR*ConsR(:)
ELSE
  !U=URS
  Flux(:) = FluxR(:) + SR*(U_RS(:) - ConsR(:))
END IF

END SUBROUTINE RiemannSolverByHLLD



SUBROUTINE EvalHLLState(ConsL,ConsR,SL,SR,FluxL,FluxR,U_HLL)
!==================================================================================================================================
! Calculates the HLL state for use with the MHD HLLC Riemann solver
!==================================================================================================================================
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)     :: ConsL(PP_nVar)
REAL,INTENT(IN)     :: ConsR(PP_nVar)
REAL,INTENT(IN)     :: SL
REAL,INTENT(IN)     :: SR
REAL,INTENT(IN)     :: FluxL(PP_nVar)
REAL,INTENT(IN)     :: FluxR(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(OUT)    :: U_HLL(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

U_HLL=(SR*ConsR-SL*ConsL-FluxR+FluxL)/(SR-SL)

END SUBROUTINE EvalHLLState



!==================================================================================================================================
!> Roe solver following the paper of Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997
!==================================================================================================================================
SUBROUTINE RiemannSolverByRoe(ConsL,ConsR,Flux)
USE MOD_PreProc
USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
USE MOD_Equation_vars, ONLY: WaveSpeeds1D
USE MOD_Equation_vars, ONLY: Kappa,KappaM1,KappaM2,sKappaM1,smu_0,s2mu_0
USE MOD_Equation_vars, ONLY: ConsToPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: ConsL(1:PP_nVar) !<  left conservative state  
REAL,INTENT(IN)  :: ConsR(1:PP_nVar) !< right conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar) !<numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: i
REAL             :: FluxL(PP_nVar), FluxR(PP_nVar)
REAL             :: PrimL(PP_nVar), PrimR(PP_nVar)
REAL             :: ptot_L, ptot_R
REAL             :: SqrtRho_L,SqrtRho_R,sSqrtRoe_LR
REAL             :: Roe_L,Roe_R
REAL             :: Rho_Roe,sRho_Roe,SqrtRho_Roe,sSqrtRho_Roe,RoeVel(3),RoeB(3),RoeVel2,RoeB2
REAL             :: H_L,H_R,RoeH
REAL             :: XX, va2_Roe,c2_Roe,sc2_Roe,ca2_Roe,cs2_Roe,cf2_Roe,astar_Roe
REAL             :: c_Roe,ca_Roe,cs_Roe,cf_Roe
REAL             :: lambda(1:7), eta(1:7), EVr(1:8,1:7)
REAL             :: delta_rho,delta_v(3),delta_B(3),delta_p
REAL             :: beta_y,beta_z,beta_v,beta_dv,beta_dB,BtMag,sBx 
REAL             :: alpha_f,alpha_s,scf2_cs2,as_cs_s,as_c_ssqRho,af_cf_s,af_c_ssqRho
REAL             :: efix_delta,sefix_delta
!==================================================================================================================================
efix_delta=1.0E-06
sefix_delta=1.0E+06

CALL ConsToPrim(PrimL(:),ConsL(:))
CALL ConsToPrim(PrimR(:),ConsR(:))
!--------------------------------------------------------------------------------
!use Roe mean wavespeeds from  Roe meanvalues 
!   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)
SqrtRho_L   = SQRT(PrimL(1))
SqrtRho_R   = SQRT(PrimR(1))
sSqrtRoe_LR = 1./(SqrtRho_L+SqrtRho_R)
Roe_L       = SqrtRho_L*sSqrtRoe_LR
Roe_R       = SqrtRho_R*sSqrtRoe_LR
ptot_L      = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!
ptot_R      = PrimR(5)+s2mu_0*SUM(PrimR(6:8)**2) !Total presssure!

Rho_Roe     = (SqrtRho_L*SqrtRho_R)
sRho_Roe    = 1./Rho_Roe
SqrtRho_Roe = SQRT(Rho_Roe)
sSqrtRho_Roe= 1./SqrtRho_Roe
RoeVel(1:3) = Roe_L*PrimL(2:4)+Roe_R*PrimR(2:4)
RoeVel2     = SUM(RoeVel(:)**2)
RoeB(1:3)   = Roe_L*PrimR(6:8)+Roe_R*PrimL(6:8)
RoeB2       = SUM(RoeB(:)**2)

H_L  = (ConsL(5)+ptot_L)/PrimL(1)
H_R  = (ConsR(5)+ptot_R)/PrimR(1)
RoeH = Roe_L*H_L+Roe_R*H_R

XX   = 0.5*SUM((PrimL(6:8)-PrimR(6:8))**2)*(sSqrtRoe_LR**2)

va2_Roe   = RoeB2*smu_0*sRho_Roe
c2_Roe    = (2.-Kappa)*XX+ KappaM1*(RoeH-0.5*RoeVel2-va2_Roe)
sc2_Roe   = 1./c2_Roe
c_Roe     = SQRT(c2_Roe) 
ca2_Roe   = RoeB(1)**2*smu_0*sRho_Roe
ca_Roe    = SQRT(ca2_Roe) 
astar_Roe = SQRT((c2_Roe+va2_Roe)**2-4.*c2_Roe*ca2_Roe) 
cf2_Roe   = 0.5*(c2_Roe+va2_Roe+astar_Roe)
cf_Roe    = SQRT(cf2_Roe)
cs2_Roe   = 0.5*(c2_Roe+va2_Roe-astar_Roe)
cs_Roe    = SQRT(cs2_Roe)
!--------------------------------------------------------------------------------

CALL EvalAdvectionFlux1D(ConsL(1:PP_nVar),FluxL(1:PP_nVar))
CALL EvalAdvectionFlux1D(ConsR(1:PP_nVar),FluxR(1:PP_nVar))

! the following implementation is based on the astrophysics code Pluto (roe.c), which implements the Cargo & Gallice Roe solver
sBx  = SIGN(1.,RoeB(1))
Btmag=SQRT(RoeB(2)*RoeB(2)+RoeB(3)*RoeB(3))
IF(Btmag.GT. 1.0E-9)THEN
  beta_y=RoeB(2)/BtMag
  beta_z=RoeB(3)/BtMag
ELSE
  beta_y=SQRT(0.5)
  beta_z=beta_y
END IF

IF(cf2_Roe .EQ. cs2_Roe) THEN
  alpha_f  = 1.
  alpha_s  = 0.
ELSEIF(c2_Roe .LE. cs2_Roe)THEN
  alpha_f  = 0.
  alpha_s  = 1.
ELSEIF(cf2_Roe .LE. c2_Roe)THEN
  alpha_f  = 1.
  alpha_s  = 0.
ELSE
  scf2_cs2 = 1./(cf2_Roe - cs2_Roe)
  alpha_f  = (c2_Roe  - cs2_Roe)*scf2_cs2
  alpha_s  = (cf2_Roe -  c2_Roe)*scf2_cs2
  alpha_f  = SQRT(MAX(0., alpha_f))
  alpha_s  = SQRT(MAX(0., alpha_s))
END IF

! calculate differences in primitive variables
Delta_rho    = PrimR(1) - PrimL(1)
Delta_v(1:3) = PrimR(2:4) - PrimL(2:4)
Delta_B(1:3) = PrimR(6:8) - PrimL(6:8)
Delta_p   = KappaM1*( (0.5*RoeVel2-XX)*delta_rho -SUM(RoeVel*(ConsR(2:4)-ConsL(2:4)))+(ConsR(5)-ConsL(5))-smu_0*SUM(RoeB*delta_B)) 

! mean eigenvalues
lambda(1)    = RoeVel(1)-cf_Roe
lambda(2)    = RoeVel(1)-ca_Roe
lambda(3)    = RoeVel(1)-cs_Roe
lambda(4)    = RoeVel(1)
lambda(5)    = RoeVel(1)+cs_Roe
lambda(6)    = RoeVel(1)+ca_Roe
lambda(7)    = RoeVel(1)+cf_Roe

! mean eigenvectors
! u-c_f
beta_v  = beta_y* RoeVel(2) + beta_z* RoeVel(3)
beta_dv = beta_y*delta_v(2) + beta_z*delta_v(3)
beta_dB = beta_y*delta_B(2) + beta_z*delta_B(3)

as_cs_s    = alpha_s*cs_Roe*sBx
as_c_ssqRho= alpha_s*c_Roe*sSqrtRho_Roe

EVr(1,1) = alpha_f
EVr(2,1) = alpha_f*lambda(1)
EVr(3,1) = alpha_f*RoeVel(2)+as_cs_s*beta_y
EVr(4,1) = alpha_f*RoeVel(3)+as_cs_s*beta_z
EVr(5,1) = alpha_f*(RoeH-va2_Roe-RoeVel(1)*cf_Roe)+as_cs_s*beta_v+as_c_ssqRho*Btmag
EVr(6,1) = 0.
EVr(7,1) = as_c_ssqRho*beta_y 
EVr(8,1) = as_c_ssqRho*beta_z


eta(1) = 0.5*sc2_Roe*(  alpha_f*(XX*delta_rho+delta_p)  & 
                      + Rho_Roe*(as_cs_s*beta_dv-alpha_f*cf_Roe*delta_v(1) + as_c_ssqRho*beta_dB))

! u+c_f

EVr(1,7) = alpha_f
EVr(2,7) = alpha_f*lambda(7)
EVr(3,7) = alpha_f*RoeVel(2)-as_cs_s*beta_y
EVr(4,7) = alpha_f*RoeVel(3)-as_cs_s*beta_z
EVr(5,7) = alpha_f*(RoeH-va2_Roe+RoeVel(1)*cf_Roe)-as_cs_s*beta_v+as_c_ssqRho*Btmag
EVr(6,7) = 0.
EVr(7,7) = EVr(7,1)
EVr(8,7) = EVr(8,1)


eta(7) = 0.5*sc2_Roe*(  alpha_f*(XX*delta_rho+delta_p)  &
                      - Rho_Roe*(as_cs_s*beta_dv-alpha_f*cf_Roe*delta_v(1) - as_c_ssqRho*beta_dB))

! u -c_s
af_cf_s=alpha_f*cf_Roe*sBx
af_c_ssqRho=alpha_f*c_Roe*sSqrtRho_Roe

EVr(1,3) = alpha_s
EVr(2,3) = alpha_s*lambda(3)
EVr(3,3) = alpha_s*RoeVel(2)-af_cf_s*beta_y
EVr(4,3) = alpha_s*RoeVel(3)-af_cf_s*beta_z
EVr(5,3) = alpha_s*(RoeH-va2_Roe-RoeVel(1)*cs_Roe)-af_cf_s*beta_v-af_c_ssqRho*Btmag
EVr(6,3) = 0.
EVr(7,3) = -af_c_ssqRho*beta_y 
EVr(8,3) = -af_c_ssqRho*beta_z

eta(3) = 0.5*sc2_Roe*(  alpha_s*(XX*delta_rho+delta_p) &
                      - Rho_Roe*(af_cf_s*beta_dv+alpha_s*cs_Roe*delta_v(1) + af_c_ssqRho*beta_dB))

! u +c_s

EVr(1,5) = alpha_s
EVr(2,5) = alpha_s*lambda(5)
EVr(3,5) = alpha_s*RoeVel(2)+af_cf_s*beta_y
EVr(4,5) = alpha_s*RoeVel(3)+af_cf_s*beta_z
EVr(5,5) = alpha_s*(RoeH-va2_Roe+RoeVel(1)*cs_Roe)+af_cf_s*beta_v-af_c_ssqRho*Btmag
EVr(6,5) = 0.
EVr(7,5) = EVr(7,3)
EVr(8,5) = EVr(8,3) 

eta(5) = 0.5*sc2_Roe*(  alpha_s*(XX*delta_rho+delta_p) &
                      + Rho_Roe*(af_cf_s*beta_dv+alpha_s*cs_Roe*delta_v(1) - af_c_ssqRho*beta_dB))

! u-c_a
beta_v  = beta_z* RoeVel(2) - beta_y* RoeVel(3) !overwrite!
beta_dv = beta_z*delta_v(2) - beta_y*delta_v(3) !overwrite!
beta_dB = beta_z*delta_B(2) - beta_y*delta_B(3) !overwrite!

EVr(1,2)=0.
EVr(2,2)=0.
EVr(3,2)=-Rho_Roe*beta_z
EVr(4,2)= Rho_Roe*beta_y
EVr(5,2)=-Rho_Roe*beta_v
EVr(6,2)=0. 
EVr(7,2)=-sBx*SqrtRho_roe*beta_z 
EVr(8,2)= sBx*SqrtRho_roe*beta_y 

eta(2) = 0.5*( -beta_dv  - sBx*sSqrtRho_Roe*beta_dB)

! u+c_a

EVr(1,6)=0.
EVr(2,6)=0.
EVr(3,6)=-EVr(3,2)
EVr(4,6)=-EVr(4,2)
EVr(5,6)=-EVr(5,2)
EVr(6,6)=0. 
EVr(7,6)=EVr(7,2)
EVr(8,6)=EVr(8,2)

eta(6) = 0.5*( beta_dv  - sBx*sSqrtRho_Roe*beta_dB)


!u 

Evr(1,4)   = 1.
Evr(2,4)   = RoeVel(1)
Evr(3,4)   = RoeVel(2)
Evr(4,4)   = RoeVel(3)
Evr(5,4)   = 0.5*RoeVel2+ KappaM2*sKappaM1*XX
Evr(6:8,4) = 0.

eta(4) =  sc2_Roe*((c2_Roe - XX)*delta_rho - delta_p)


!entropy fix

IF(ABS(lambda(1)) .LT. 0.5*efix_delta) lambda(1)= lambda(1)*lambda(1)*sefix_delta+0.25*efix_delta
IF(ABS(lambda(7)) .LT. 0.5*efix_delta) lambda(7)= lambda(7)*lambda(7)*sefix_delta+0.25*efix_delta
IF(ABS(lambda(3)) .LT. 0.5*efix_delta) lambda(3)= lambda(3)*lambda(3)*sefix_delta+0.25*efix_delta
IF(ABS(lambda(6)) .LT. 0.5*efix_delta) lambda(6)= lambda(6)*lambda(6)*sefix_delta+0.25*efix_delta

!CHECK Roe Matrix condition FR-FL = A*(UR-UL) = R*lambda*eta
!Flux(:)=FluxR-FluxL
!DO i=1,7
!  Flux(1:8)=Flux(1:8) - eta(i)*lambda(i)*EVr(1:8,i)
!END DO
!DO i=1,8
!  IF(ABS(Flux(i))>1.0E-08)THEN
!    WRITE(*,*)'Roe matrix not satisfied, var=',i,'Flux(i)=',Flux(i)
!  END IF 
!END DO
!IF(MAXVAL(ABS(Flux))>1.0E-08) STOP


! assemble Roe flux
Flux(:)=FluxL(:)+FluxR(:)
DO i=1,7
  Flux(1:8)=Flux(1:8) - eta(i)*ABS(lambda(i))*EVr(1:8,i)
END DO
Flux(:)=0.5*Flux(:)

END SUBROUTINE RiemannSolverByRoe


END MODULE MOD_Riemann
