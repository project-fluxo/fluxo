!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2021 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2021 Andrés Rueda
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

!INTERFACE AddNonConsFlux
!  MODULE PROCEDURE AddNonConsFlux
!END INTERFACE

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

INTERFACE EntropyStableByLLF
  MODULE PROCEDURE EntropyStableByLLF
END INTERFACE

PUBLIC :: Riemann, AdvRiemann
#if NONCONS
PUBLIC :: AddNonConsFlux
#endif /*NONCONS*/
PUBLIC :: RiemannSolverByHLL
PUBLIC :: RiemannSolverByHLLC
PUBLIC :: RiemannSolverByHLLD
PUBLIC :: RiemannSolverByRoe
PUBLIC :: RiemannSolverByRusanov
PUBLIC :: EntropyStableByLLF
#ifdef PP_GLM
PUBLIC :: EntropyStableDerigsFlux
PUBLIC :: EntropyStableDerigsFlux_VolFluxAndDissipMatrices
PUBLIC :: EntropyStableFloGorFlux
PUBLIC :: EntropyStableFloGorFlux_VolFluxAndDissipMatrices
#endif
PUBLIC :: RotateState
PUBLIC :: RotateFluxBack
!==================================================================================================================================


CONTAINS
!==================================================================================================================================
!> main routine calling different riemann solvers
!==================================================================================================================================
SUBROUTINE Riemann(F,UL,UR,                                                                       &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R,                         &
#endif /*PARABOLIC*/
                   nv,t1,t2)
USE MOD_PreProc
USE MOD_Equation_vars,ONLY: SolveRiemannProblem
#if PARABOLIC
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
#endif /*PARABOLIC*/
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
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N) :: k_L,g_L,j_L,k_R,g_R,j_R
INTEGER                                 :: iVar
#endif /*PARABOLIC*/
!==================================================================================================================================

call AdvRiemann(F,UL,UR,nv,t1,t2)

#if PARABOLIC
!! Don#t forget the diffusion contribution, my young padawan
!! Compute MHD Diffusion flux

CALL EvalDiffFlux3D(k_L,g_L,j_L,UL,gradUx_L,gradUy_L,gradUz_L)

CALL EvalDiffFlux3D(k_R,g_R,j_R,UR,gradUx_R,gradUy_R,gradUz_R)

!
! !BR1/BR2 uses arithmetic mean of the fluxes
! ATTENTION: This is done from iVar=2 because the first component of the viscous fluxes is  F(1,:,:)=0
!            ... Change to iVar=1 if that is no longer the case
DO iVar=2,PP_nVar
  F(iVar,:,:)=F(iVar,:,:)+0.5*( nv(1,:,:)*(k_L(iVar,:,:)+k_R(iVar,:,:))  &
                               +nv(2,:,:)*(g_L(iVar,:,:)+g_R(iVar,:,:))  &
                               +nv(3,:,:)*(j_L(iVar,:,:)+j_R(iVar,:,:))  )
END DO

#endif /* PARABOLIC */

END SUBROUTINE Riemann

!==================================================================================================================================
!> Advective Riemann solver
!==================================================================================================================================
pure SUBROUTINE AdvRiemann(F,UL,UR,nv,t1,t2)
USE MOD_PreProc
USE MOD_Equation_vars,ONLY: SolveRiemannProblem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: UL(      PP_nVar,0:PP_N,0:PP_N) !<  left state on face
REAL,INTENT(IN) :: UR(      PP_nVar,0:PP_N,0:PP_N) !< right state on face
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
INTEGER          :: i,j
!==================================================================================================================================

DO j = 0,PP_N
  DO i = 0,PP_N
    !rotate fields
    ConsL = RotateState(UL(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))
    ConsR = RotateState(UR(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))

    CALL SolveRiemannProblem(ConsL(:),ConsR(:),F(:,i,j))

    ! Back Rotate the normal flux into Cartesian direction
    call RotateFluxBack(F(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))

  END DO !i
END DO !j
END SUBROUTINE AdvRiemann
!==================================================================================================================================
!> Rotate the state to the normal frame of reference
!==================================================================================================================================
pure function RotateState(U,nv,t1,t2) result(rotU)
  implicit none
  real, intent(in) :: U(PP_nVar)
  real, intent(in) :: nv(3)
  real, intent(in) :: t1(3)
  real, intent(in) :: t2(3)
  real             :: rotU(PP_nVar)

  rotU(:) =     U(:  )
  rotU(2) = SUM(U(2:4)*nv(:))
  rotU(3) = SUM(U(2:4)*t1(:))
  rotU(4) = SUM(U(2:4)*t2(:))
  rotU(6) = SUM(U(6:8)*nv(:))
  rotU(7) = SUM(U(6:8)*t1(:))
  rotU(8) = SUM(U(6:8)*t2(:))
end function RotateState
!==================================================================================================================================
!> Rotate the flux from the normal frame of reference back to the physical framework
!==================================================================================================================================
pure subroutine RotateFluxBack(F,nv,t1,t2)
  implicit none
  real, intent(inout) :: F(PP_nVar)
  real, intent(in)    :: nv(3)
  real, intent(in)    :: t1(3)
  real, intent(in)    :: t2(3)

  F(2:4) =   nv(:)*F(2) &
           + t1(:)*F(3) &
           + t2(:)*F(4)
  F(6:8) =   nv(:)*F(6) &
           + t1(:)*F(7) &
           + t2(:)*F(8)
end subroutine RotateFluxBack
#if NONCONS
!==================================================================================================================================
!> Nonconservative flux on a side: 
!> ({B}·nvec)*phi_L^MHD + {psi}*phi_L^GLM·nvec
!> 
!> phi^MHD is the Powell, Brackbill or Janhunen nonconservative term:
!> * Powell:    phi^MHD = (0,B_1,B_2,B_3,v·B,v_1,v_2,v_3,0)
!> * Brackbill: phi^MHD = (0,B_1,B_2,B_3,0,0,0,0)
!> * Janhunen:  phi^MHD = (0,0,0,0,0,v_1,v_2,v_3,0)
!> phi^GLM is the GLM nonconservative term
!> * phi^GLM = (0,0,0,0,psi*v,0,0,0,v)
!==================================================================================================================================
pure SUBROUTINE AddNonConsFlux(FL,UL,UR,nv,t1,t2)
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:nTotal_Face
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: UL(      PP_nVar,nTotal_Face) !<  left state on face
REAL,INTENT(IN) :: UR(      PP_nVar,nTotal_Face) !< right state on face
REAL,INTENT(IN) :: nv(            3,nTotal_Face) !< normal vector of face
REAL,INTENT(IN) :: t1(            3,nTotal_Face) !< 1st tangential vector of face
REAL,INTENT(IN) :: t2(            3,nTotal_Face) !< 2nd tangential vector of face
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT):: FL(       PP_nVar,nTotal_Face) !< ADDING TO nonconservative flux on UL side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: v_L(3)
!==================================================================================================================================
DO i=1,nTotal_Face
  v_L=UL(2:4,i)/UL(1,i)
#if NONCONS==1 /*Powell*/
  FL(2:8,i)=FL(2:8,i) +(0.5*SUM((UL(6:8,i)+UR(6:8,i))*nv(:,i)))*(/UL(6:8,i),SUM(UL(6:8,i)*v_L(1:3)),v_L(1:3)/)
#elif NONCONS==2 /*Brackbill*/
  FL(2:4,i)=FL(2:4,i) +(0.5*SUM((UL(6:8,i)+UR(6:8,i))*nv(:,i)))*UL(6:8,i)
#elif NONCONS==3 /*Janhunen*/
  FL(6:8,i)=FL(6:8,i) +(0.5*SUM((UL(6:8,i)+UR(6:8,i))*nv(:,i)))*v_L(1:3)
#endif /*NONCONSTYPE*/


#if defined (PP_GLM) && defined (PP_NC_GLM)
  !nonconservative term to restore Galilean invariance for GLM term
  FL((/5,9/),i)=FL((/5,9/),i) +(0.5*SUM(v_L(:)*nv(:,i)))*(/UL(9,i)*(UR(9,i)+UL(9,i)),UR(9,i)+UL(9,i)/)
#endif /*PP_GLM and PP_NC_GLM*/

END DO !i=1,nTotal_Face

END SUBROUTINE AddNonConsFlux
#endif /*NONCONS*/

!==================================================================================================================================
!> Lax-friedrichs / rusanov flux
!==================================================================================================================================
pure SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,Flux)
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
REAL,INTENT(OUT) :: Flux(1:PP_nVar) !<numerical flux
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
pure SUBROUTINE RiemannSolverByHLL(ConsL,ConsR,Flux)
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
REAL,INTENT(OUT) :: Flux(1:PP_nVar) !<numerical flux
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
pure SUBROUTINE RiemannSolverByHLLC(ConsL,ConsR,Flux)
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
REAL,INTENT(OUT) :: Flux(1:PP_nVar) !<numerical flux
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
pure SUBROUTINE RiemannSolverByHLLD(ConsL_in,ConsR_in,Flux)
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
REAL,INTENT(OUT) :: Flux(1:PP_nVar) !<numerical flux
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
PrimL(6)=Bn

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



pure SUBROUTINE EvalHLLState(ConsL,ConsR,SL,SR,FluxL,FluxR,U_HLL)
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
pure SUBROUTINE RiemannSolverByRoe(ConsL,ConsR,Flux)
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
REAL,INTENT(OUT) :: Flux(1:PP_nVar) !<numerical flux
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


pure SUBROUTINE EntropyStableByLLF(UL,UR,Fstar)
!==================================================================================================================================
! entropy conservation for MHD, kinetric Energy conservation only in the Euler case
! following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system"
! mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t)
!
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Flux_Average, ONLY:LN_MEAN
USE MOD_Equation_Vars,ONLY:kappa,kappaM1,skappaM1,smu_0,s2mu_0,consToEntropy
USE MOD_Equation_Vars,ONLY:VolumeFluxAverage !pointer to EC routine
!##ifdef PP_GLM
!#USE MOD_Equation_Vars,ONLY:GLM_ch
!##endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: sbetaLN,beta_R,beta_L
REAL            :: srho_L,srho_R,rhoLN,srhoLN
REAL            :: B2_L,B2_R,B2Avg
REAL            :: u2_L,u2_R,u2Avg,uAvg2
REAL            :: p_L,p_R,pAvg,pLN !,pTilde
REAL            :: a2Avg,va2Avg,ca2Avg,cfAvg,LambdaMax_s2
REAL            :: u_L(3),u_R(3)
REAL            :: BAvg(3),uAvg(3)
!REAL            :: u1_B2Avg,uB_Avg
REAL            :: Hmatrix(5,5),tau,Eavg
REAL            :: V_jump(PP_nVar)
#ifdef PP_GLM
REAL            :: psiAvg
#endif
!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
              E_L =>UL(5)-s2mu_0*UL(9)**2, E_R =>UR(5)-s2mu_0*UR(9)**2, &
            psi_L =>UL(9)   ,  psi_R =>UR(9), &
#else
              E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
              B_L => UL(6:8),    B_R => UR(6:8)  )
! Get the inverse density, velocity, and pressure on left and right
srho_L=1./rho_L
srho_R=1./rho_R
u_L = rhoU_L(:)*srho_L
u_R = rhoU_R(:)*srho_R

u2_L = SUM(u_L(:)*u_L(:))
u2_R = SUM(u_R(:)*u_R(:))
B2_L = SUM(B_L(:)*B_L(:))
B2_R = SUM(B_R(:)*B_R(:))

!beta=rho/(2*p)
p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*B2_L))
p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*B2_R))
beta_L = 0.5*rho_L/p_L
beta_R = 0.5*rho_R/P_R

! Get the averages for the numerical flux

rhoLN      = LN_MEAN( rho_L, rho_R)
sbetaLN    = 1./LN_MEAN(beta_L,beta_R)
uAvg       = 0.5 * ( u_L +  u_R)
u2Avg      = 0.5 * (u2_L + u2_R)
BAvg       = 0.5 * ( B_L +  B_R)
B2Avg      = 0.5 * (B2_L + B2_R)

pAvg       = 0.5*(rho_L+rho_R)/(beta_L+beta_R) !rhoMEAN/(2*betaMEAN)
#ifdef PP_GLM
psiAvg     = 0.5*(psi_L+psi_R)
#endif

! Entropy conserving and kinetic energy conserving flux, pointer to EC routine
CALL VolumeFluxAverage(UL,UR,Fstar)

!u1_B2Avg   = 0.5 * (u_L(1)*B2_L       + u_R(1)*B2_R)
!uB_Avg     = 0.5 * (SUM(u_L(:)*B_L(:))+ SUM(u_R(:)*B_R(:)))
!pTilde     = pAvg+ s2mu_0*B2Avg !+1/(2mu_0)({{|B|^2}}...)
!Fstar(1) = rhoLN*uAvg(1)
!Fstar(2) = Fstar(1)*uAvg(1) - smu_0*Bavg(1)*BAvg(1) + pTilde
!Fstar(3) = Fstar(1)*uAvg(2) - smu_0*Bavg(1)*BAvg(2)
!Fstar(4) = Fstar(1)*uAvg(3) - smu_0*Bavg(1)*BAvg(3)
!Fstar(7) = uAvg(1)*Bavg(2) - BAvg(1)*uAvg(2)
!Fstar(8) = uAvg(1)*Bavg(3) - BAvg(1)*uAvg(3)
!#ifdef PP_GLM
!Fstar(6) = GLM_ch*psiAvg
!Fstar(9) = GLM_ch*BAvg(1)
!#endif
!
!Fstar(5) = Fstar(1)*0.5*(skappaM1*sbetaLN - u2Avg)  &
!           + SUM(uAvg(:)*Fstar(2:4)) &
!           +smu_0*( SUM(BAvg(:)*Fstar(6:8)) &
!                   -0.5*u1_B2Avg +BAvg(1)*uB_Avg &
!#ifdef PP_GLM
!                   +Fstar(9)*psiAvg-GLM_ch*0.5*(B_L(1)*psi_L+B_R(1)*psi_R)     &
!#endif
!                   )
! MHD wavespeed
pLN = 0.5*rhoLN*sbetaLN

srhoLN = 1./rhoLN
a2Avg  = kappa*pAvg*srhoLN
va2Avg = B2Avg*sRhoLN
ca2Avg = BAvg(1)*BAvg(1)*srhoLN

cfAvg  = SQRT(0.5 * ((a2Avg+va2Avg) + SQRT((a2Avg+va2Avg)**2 - 4.0*a2Avg*ca2Avg)))

LambdaMax_s2 =0.5*(ABS(uAvg(1))+cfAvg)

! H matrix stabilization

uAvg2 = SUM(uAvg(:)*uAvg(:))
Eavg = pLN*sKappaM1 + 0.5*rhoLN*(2.0*uAvg2-u2Avg)

! Euler components as derived by Dominik, April 18-20, 2016
! MATRIX IS SYMMETRIC!!
! Hmatrix = 0.0
Hmatrix(1:5,1)=(/rhoLN       ,rhoLN*uAvg(1)            ,rhoLN*uAvg(2)            ,rhoLN*uAvg(3)            , Eavg               /)
!                -------------     ||                       ||                         ||
Hmatrix(1:5,2)=(/Hmatrix(2,1),Hmatrix(2,1)*uAvg(1)+pAvg,Hmatrix(3,1)*uAvg(1)     ,Hmatrix(4,1)*uAvg(1)     ,(Eavg+pAvg)*uAvg(1) /)
!                            ---------------------------    ||                         ||
Hmatrix(1:5,3)=(/Hmatrix(3,1),    Hmatrix(3,2)         ,Hmatrix(3,1)*uAvg(2)+pAvg,Hmatrix(4,1)*uAvg(2)     ,(Eavg+pAvg)*uAvg(2) /)
!                                                      ---------------------------     ||
Hmatrix(1:5,4)=(/Hmatrix(4,1),    Hmatrix(4,2)         ,    Hmatrix(4,3)         ,Hmatrix(4,1)*uAvg(3)+pAvg,(Eavg+pAvg)*uAvg(3) /)
!                                                                                ---------------------------
Hmatrix(1:4,5)=(/Hmatrix(5,1),    Hmatrix(5,2)         ,    Hmatrix(5,3)         ,     Hmatrix(5,4)       /)   !Hmatrix(5,5)

tau = 1./(beta_L+beta_R) !0.5/betaA
!Hmatrix(5,5) = 0.25*rhoA/(betaL*betaR)*sKappaM1 + Eavg*Eavg/rhoLN + (u1A*u1A+v1A*v1A+w1A*w1A)*pAvg + tau*(B11A*B11A + B21A*B21A + B31A*B31A)
!Special averaging such that RSRT-H = 0 in the Euler case
Hmatrix(5,5) = (pLN*pLN*sKappaM1 + Eavg*Eavg)*srhoLN + uAvg2*pAvg + tau*( SUM(BAvg(:)*Bavg(:)) &
#ifdef PP_GLM
                                                                         + psiAvg*psiAvg &
#endif
                                                                        )

!! Magnetic field components
!Hmatrix(6:8,5) = tau * Bavg(:)
!Hmatrix(5,6:8) = Hmatrix(5,6:8)
!
!Hmatrix(6,6) = tau
!Hmatrix(7,7) = tau
!Hmatrix(8,8) = tau
!
!#ifdef PP_GLM
!! Psi components
!Hmatrix(9,5) = tau*psiAvg
!Hmatrix(5,9) = Hmatrix(9,5)
!Hmatrix(9,9) = tau
!#endif /*PP_GLM*/

!V_jump=consToEntropy(UR)-consToEntropy(UL) !better use already computed values!!

!V(1)=(kappa-s)/(kappa-1)-beta*|v|^2 , s = LOG(p) - kappa*LOG(cons(1)), VR-VL, using Log(R)-Log(L)=LOG(R/L)
V_jump(1)         = (kappa*(LOG(rho_R*srho_L))-LOG(p_R/p_L))*skappaM1 &
                       - (beta_R*u2_R         -beta_L*u2_L  )
V_jump(2:4)       =  2.0*(beta_R*u_R(:)       -beta_L*u_L(:))        ! 2*beta*v
V_jump(5)         = -2.0*(beta_R              -beta_L       )        !-2*beta
V_jump(6:PP_nVar) =  2.0*(beta_R*UR(6:PP_nVar)-beta_L*UL(6:PP_nVar)) ! 2*beta*B

!Fstar(:) = Fstar -LambdaMax_s2*MATMUL(Hmatrix,(V_jump))

!Hmatrix has only non-zero entries (1:5,1:5),(5,6:9),(6:9,5),diag(6:9,6:9)
Fstar(1:4) = Fstar(1:4) - LambdaMax_s2       *MATMUL(Hmatrix(1:4,1:5),V_jump(1:5))
Fstar(  5) = Fstar(  5) - LambdaMax_s2       *(  SUM(Hmatrix(1:5,  5)*V_jump(1:5))  &
                                                +tau*  (SUM(Bavg(1:3)*V_jump(6:8)) &
#ifdef PP_GLM
                                                        +      psiAvg*V_jump(  9)  &
#endif
                                                       ))
Fstar(6:8) = Fstar(6:8) - (LambdaMax_s2*tau)*(Bavg(1:3)*V_jump(5)   + V_jump(6:8))
#ifdef PP_GLM
Fstar(  9) = Fstar(  9) - (LambdaMax_s2*tau)*(   psiAvg*V_jump(5)   + V_jump(  9))
#endif

END ASSOCIATE
END SUBROUTINE EntropyStableByLLF



#ifdef PP_GLM
!==================================================================================================================================
!> The following three routines are for the entropy stable 9wave solver.
!>    See: Derigs et al. (2018). "Ideal GLM-MHD: About the entropy consistent nine-wave magnetic field divergence diminishing ideal
!>         magnetohydrodynamics equations. Journal of Computational Physics, 364, 420-467."
!> ATTENTION: 1) Central flux is entropy conservaing for MHD, kinetic Energy conservation only in the Euler case
!>            2) mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t)
!>            3) Dissipation for each characteristic wave seperately
!==================================================================================================================================
pure SUBROUTINE EntropyStableDerigsFlux(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Flux_Average, ONLY:LN_MEAN
USE MOD_Equation_Vars,ONLY:kappa,kappaM1,skappaM1,smu_0,s2mu_0,consToEntropy
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: beta_R,beta_L
REAL            :: srho_L,srho_R
REAL            :: u2_L,u2_R
REAL            :: p_L,p_R
REAL            :: u_L(3),u_R(3)
REAL            :: V_jump(PP_nVar)
REAL            :: Dmatrix(PP_nVar),Rmatrix(PP_nVar,PP_nVar),RT(PP_nVar,PP_nVar)
!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
!
           rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
              E_L =>UL(5)-s2mu_0*UL(9)**2, E_R =>UR(5)-s2mu_0*UR(9)**2, &
#else
              E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
              B_L => UL(6:8),    B_R => UR(6:8)  )

! Get auxiliar variables
srho_L=1./rho_L
srho_R=1./rho_R
u_L = rhoU_L(:)*srho_L
u_R = rhoU_R(:)*srho_R
u2_L = SUM(u_L(:)*u_L(:))
u2_R = SUM(u_R(:)*u_R(:))
p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*SUM(B_L(:)*B_L(:))))
p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*SUM(B_R(:)*B_R(:))))
beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
beta_R = 0.5*rho_R/P_R

! Jump in entropy vars
V_jump(1)         = (kappa*(LOG(rho_R*srho_L))-LOG(p_R/p_L))*skappaM1 &
                       - (beta_R*u2_R         -beta_L*u2_L  )
V_jump(2:4)       =  2.0*(beta_R*u_R(:)       -beta_L*u_L(:))        ! 2*beta*v
V_jump(5)         = -2.0*(beta_R              -beta_L       )        !-2*beta
V_jump(6:PP_nVar) =  2.0*(beta_R*UR(6:PP_nVar)-beta_L*UL(6:PP_nVar)) ! 2*beta*B

call EntropyStableDerigsFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)
RT = TRANSPOSE(Rmatrix)

! Compute entropy-stable fluxes
Fstar = Fstar - 0.5*MATMUL(Rmatrix,Dmatrix*MATMUL(RT,V_jump))

END ASSOCIATE
END SUBROUTINE EntropyStableDerigsFlux
!==================================================================================================================================
!> Computation of the entropy conserving flux and the dissipation matrices for the entropy stable 9wave solver
!==================================================================================================================================
pure subroutine EntropyStableDerigsFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)
  USE MOD_Equation_Vars,ONLY:kappaM1, smu_0, s2mu_0, sKappaM1, Kappa
  USE MOD_Flux_Average, ONLY:LN_MEAN
#ifdef PP_GLM
  USE MOD_Equation_Vars,ONLY:GLM_ch
#endif
  implicit none
  !-arguments-----------------------------------------------
  real, intent(in)  :: UL(PP_nVar)      !< left state
  real, intent(in)  :: UR(PP_nVar)      !< right state
  real, intent(out) :: Fstar  (PP_nVar)
  real, intent(out) :: Dmatrix(PP_nVar)
  real, intent(out) :: Rmatrix(PP_nVar,PP_nVar)
  !-local-variables-----------------------------------------
  ! Mean states
  real    :: sbetaLN,betaAvg,rhoLN,pAvg,psiAvg,BAvg(3),uAvg(3)
  real    :: p_L,p_R,u_L(3), u_R(3), pTilde
  ! Additional
  REAL :: B2_L,B2_R,B2Avg
  REAL :: u2_L,u2_R,u2Avg
  REAL :: beta_R,beta_L
  REAL :: srho_L,srho_R
  REAL :: u1_B2Avg, uB_Avg
  real :: Tmatrix(PP_nVar)
  REAL :: aA,aLN,abeta,bb1A,bb2A,bb3A,bbA,aabbA,ca,xx,xxx,cf,cs,bperpA,beta2A,beta3A,alphaf,alphas,sgnb1
  REAL :: psiSplus,psiSminus,psiFplus,psiFminus, srhoLN, rhoAvg, pLN
  real :: LambdaMax, phi
  !---------------------------------------------------------

  ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
             rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
                E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
              psi_L =>UL(9)   ,  psi_R =>UR(9), &
#else
                E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
                B_L => UL(6:8),    B_R => UR(6:8)  )
  ! Get auxiliar variables
  srho_L=1./rho_L
  srho_R=1./rho_R
  u_L = rhoU_L(:)*srho_L
  u_R = rhoU_R(:)*srho_R
  u2_L = SUM(u_L(:)*u_L(:))
  u2_R = SUM(u_R(:)*u_R(:))
  B2_L = SUM(B_L(:)*B_L(:))
  B2_R = SUM(B_R(:)*B_R(:))
  p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*B2_L))
  p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*B2_R))
  beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
  beta_R = 0.5*rho_R/P_R

  ! Get the averages for the numerical flux
  rhoLN      = LN_MEAN( rho_L, rho_R)
  sbetaLN    = 1./LN_MEAN(beta_L,beta_R)
  uAvg       = 0.5 * ( u_L +  u_R)
  u2Avg      = 0.5 * (u2_L + u2_R)
  BAvg       = 0.5 * ( B_L +  B_R)
  B2Avg      = 0.5 * (B2_L + B2_R)
  u1_B2Avg   = 0.5 * (u_L(1)*B2_L       + u_R(1)*B2_R)
  uB_Avg     = 0.5 * (SUM(u_L(:)*B_L(:))+ SUM(u_R(:)*B_R(:)))
  betaAvg    = 0.5 * (beta_L + beta_R)
  pAvg       = 0.5*(rho_L+rho_R)/(beta_L+beta_R) !rhoMEAN/(2*betaMEAN)
#ifdef PP_GLM
  psiAvg     = 0.5*(psi_L+psi_R)
#endif
  pTilde     = pAvg+ s2mu_0*B2Avg !+1/(2mu_0)({{|B|^2}}...)

  ! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rhoLN*uAvg(1)
Fstar(2) = Fstar(1)*uAvg(1) - smu_0*Bavg(1)*BAvg(1) + pTilde
Fstar(3) = Fstar(1)*uAvg(2) - smu_0*Bavg(1)*BAvg(2)
Fstar(4) = Fstar(1)*uAvg(3) - smu_0*Bavg(1)*BAvg(3)
Fstar(7) = uAvg(1)*Bavg(2) - BAvg(1)*uAvg(2)
Fstar(8) = uAvg(1)*Bavg(3) - BAvg(1)*uAvg(3)
#ifdef PP_GLM
Fstar(6) = GLM_ch*psiAvg
Fstar(9) = GLM_ch*BAvg(1)
#endif

Fstar(5) = Fstar(1)*0.5*(skappaM1*sbetaLN - u2Avg)  &
           + SUM(uAvg(:)*Fstar(2:4)) &
           +smu_0*( SUM(BAvg(:)*Fstar(6:8)) &
                   -0.5*u1_B2Avg +BAvg(1)*uB_Avg &
#ifdef PP_GLM
                   +Fstar(9)*psiAvg-GLM_ch*0.5*(B_L(1)*psi_L+B_R(1)*psi_R)     &
#endif
                   )
  ! MHD waves selective dissipation

    ! Compute additional averages
    rhoAvg = 0.5*(rho_L + rho_R)
    pLN = 0.5*rhoLN*sbetaLN
    u2Avg = u_L(1)*u_R(1) + u_L(2)*u_R(2) + u_L(3)*u_R(3)  ! Redefinition of u2Avg -> \bar{\norm{u}^2}
    srhoLN = 1./rhoLN
    aA = SQRT(kappa*pAvg/rhoLN)
    aLN = SQRT(kappa*pLN/rhoLN)
    abeta = SQRT(0.5*kappa/betaAvg)
    bb1A = BAvg(1)/SQRT(rhoLN)
    bb2A = BAvg(2)/SQRT(rhoLN)
    bb3A = BAvg(3)/SQRT(rhoLN)
    bbA = SQRT(bb1A*bb1A + bb2A*bb2A + bb3A*bb3A)
    aabbA = aA*aA + bbA*bbA

    ! Alfven speed
    ! Derigs et al. (2018), (4.63)
    ca = ABS(bb1A)

    ! Control round-off errors
    xx = aabbA*aabbA - 4.0*aA*aA*bb1A*bb1A
    if(xx .lt. 0.0) xx = 0.0
    xxx = aabbA + SQRT(xx)

    ! Fast magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cf = SQRT(0.5 * xxx)

    ! Control round-off errors
    xxx = aabbA - SQRT(xx)
    if(xxx .lt. 0.0) xxx = 0.0

    ! Slow magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cs = SQRT(0.5 * xxx)

    bperpA = SQRT(bb2A*bb2A + bb3A*bb3A)
    ! In case of very small bperpA, the values of betaA_{2,3}
    ! are indeterminable so we make them pairwise orthogonal
    ! Derigs et al. (2018), (4.63)
    if(bperpA .gt. 1.0E-14) then
      beta2A = bb2A/bperpA
      beta3A = bb3A/bperpA
    else
      bperpA = 0.0
      beta2A = 1.0/sqrt(2.0)
      beta3A = 1.0/sqrt(2.0)
    endif

    ! Avoid negative round-off errors when computing alphaf (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (aA*aA-cs*cs)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphaf = SQRT(xx)
    else
      alphaf = 0.0
    end if

    ! Avoid negative round-off errors when computing alphas (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (cf*cf-aA*aA)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphas = SQRT(xx)
    else
      alphas = 0.0
    end if

    sgnb1 = sign(1.,bb1A)

    ! Derigs et al. (2018), (4.63)
    psiSplus =  0.5*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 + alphas*cs*rhoLN*uAvg(1) + &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiSminus = 0.5*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 - alphas*cs*rhoLN*uAvg(1) - &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    psiFplus =  0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 + alphaf*cf*rhoLN*uAvg(1) - &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiFminus = 0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 - alphaf*cf*rhoLN*uAvg(1) + &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    ! + fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,1) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) + cf), &
                      rhoLN*(alphaf*uAvg(2) - alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) - alphas*cs*beta3A*sgnb1), &
                      psiFplus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,2) = (/ 0.0, &
                      0.0,  &
                      rhoLN*SQRT(rhoAvg)*beta3A, &
                      -rhoLN*SQRT(rhoAvg)*beta2A, &
                      -rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! + slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,3) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) + cs), &
                      rhoLN*(alphas*uAvg(2) + alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) + alphaf*cf*beta3A*sgnb1), &
                      psiSplus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,4) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) + psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      1.0 /)

    ! Entropy wave
    ! Derigs et al. (2018), (4.66)
    Rmatrix(:,5) = (/ 1.0, &
                      uAvg(1), &
                      uAvg(2), &
                      uAvg(3), &
                      0.5*u2Avg, &
                      0.0, &
                      0.0, &
                      0.0, &
                      0.0 /)

    ! - GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,6) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) - psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      -1.0 /)

    ! - slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,7) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) - cs), &
                      rhoLN*(alphas*uAvg(2) - alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) - alphaf*cf*beta3A*sgnb1), &
                      psiSminus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! - Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,8) = (/ 0.0, &
                      0.0, &
                      -rhoLN*SQRT(rhoAvg)*beta3A, &
                      rhoLN*SQRT(rhoAvg)*beta2A, &
                      rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! - fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,9) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) - cf), &
                      rhoLN*(alphaf*uAvg(2) + alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) + alphas*cs*beta3A*sgnb1), &
                      psiFminus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! A blend of the 9waves solver and LLF
      phi = sqrt(ABS(1.0-(p_R*rho_R/(p_L*rho_L)))/(1.0+(p_R*rho_R/(p_L*rho_L))))
      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))

      Dmatrix(1) = (1.-phi)*ABS( uAvg(1) + cf ) + phi*LambdaMax
      Dmatrix(2) = (1.-phi)*ABS( uAvg(1) + ca ) + phi*LambdaMax
      Dmatrix(3) = (1.-phi)*ABS( uAvg(1) + cs ) + phi*LambdaMax
      Dmatrix(4) = (1.-phi)*ABS( uAvg(1) + GLM_ch ) + phi*LambdaMax
      Dmatrix(5) = (1.-phi)*ABS( uAvg(1)      ) + phi*LambdaMax
      Dmatrix(6) = (1.-phi)*ABS( uAvg(1) - GLM_ch ) + phi*LambdaMax
      Dmatrix(7) = (1.-phi)*ABS( uAvg(1) - cs ) + phi*LambdaMax
      Dmatrix(8) = (1.-phi)*ABS( uAvg(1) - ca ) + phi*LambdaMax
      Dmatrix(9) = (1.-phi)*ABS( uAvg(1) - cf ) + phi*LambdaMax

    ! Pure 9waves solver (no LLF)
!#
!#      Dmatrix(1) = ABS( uAvg(1) + cf ) ! + fast magnetoacoustic wave
!#      Dmatrix(2) = ABS( uAvg(1) + ca ) ! + Alfven wave
!#      Dmatrix(3) = ABS( uAvg(1) + cs ) ! + slow magnetoacoustic wave
!#      Dmatrix(4) = ABS( uAvg(1) + GLM_ch ) ! + GLM wave
!#      Dmatrix(5) = ABS( uAvg(1)      ) ! / Entropy wave
!#      Dmatrix(6) = ABS( uAvg(1) - GLM_ch ) ! - GLM wave
!#      Dmatrix(7) = ABS( uAvg(1) - cs ) ! - slow magnetoacoustic wave
!#      Dmatrix(8) = ABS( uAvg(1) - ca ) ! - Alfven wave
!#      Dmatrix(9) = ABS( uAvg(1) - cf ) ! - fast magnetoacoustic wave

    ! LLF solver
!#      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))
!#
!#      Dmatrix(1) = LambdaMax
!#      Dmatrix(2) = LambdaMax
!#      Dmatrix(3) = LambdaMax
!#      Dmatrix(4) = LambdaMax
!#      Dmatrix(5) = LambdaMax
!#      Dmatrix(6) = LambdaMax
!#      Dmatrix(7) = LambdaMax
!#      Dmatrix(8) = LambdaMax
!#      Dmatrix(9) = LambdaMax

    ! Diagonal scaling matrix as described in Winters et al., eq. (4.15)
    Tmatrix = 0.0
    Tmatrix(1) = 0.5/kappa/rhoLN ! + f
    Tmatrix(2) = 0.25/betaAvg/rhoLN/rhoLN ! + a
    Tmatrix(3) = Tmatrix(1) ! + s
    Tmatrix(4) = 0.25 / betaAvg ! + GLM
    Tmatrix(5) = rhoLN*kappaM1/kappa ! E
    Tmatrix(6) = Tmatrix(4) ! - GLM
    Tmatrix(7) = Tmatrix(1) ! - s
    Tmatrix(8) = Tmatrix(2) ! - a
    Tmatrix(9) = Tmatrix(1) ! - f

    ! Scale D matrix
    Dmatrix = Dmatrix*Tmatrix
end ASSOCIATE
end subroutine EntropyStableDerigsFlux_VolFluxAndDissipMatrices
#endif


#ifdef PP_GLM
!==================================================================================================================================
!> The following three routines are for an entropy stable 9wave solver.
!> ATTENTION: 1) Central flux is the FloGor flux that is entropy conservaing for MHD, kinetic Energy conservation only in the Euler case
!>            2) mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t)
!>            3) Dissipation for each characteristic wave seperately
!==================================================================================================================================
PURE SUBROUTINE EntropyStableFloGorFlux(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1,skappaM1,smu_0,s2mu_0
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: beta_R,beta_L
REAL            :: srho_L,srho_R
REAL            :: u2_L,u2_R
REAL            :: p_L,p_R
REAL            :: u_L(3),u_R(3)
REAL            :: V_jump(PP_nVar)
REAL            :: Dmatrix(PP_nVar),Rmatrix(PP_nVar,PP_nVar),RT(PP_nVar,PP_nVar)
!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
!
           rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
              E_L =>UL(5)-s2mu_0*UL(9)**2, E_R =>UR(5)-s2mu_0*UR(9)**2, &
#else
              E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
              B_L => UL(6:8),    B_R => UR(6:8)  )

! Get auxiliar variables
srho_L=1./rho_L
srho_R=1./rho_R
u_L = rhoU_L(:)*srho_L
u_R = rhoU_R(:)*srho_R
u2_L = SUM(u_L(:)*u_L(:))
u2_R = SUM(u_R(:)*u_R(:))
p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*SUM(B_L(:)*B_L(:))))
p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*SUM(B_R(:)*B_R(:))))
beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
beta_R = 0.5*rho_R/P_R

! Jump in entropy vars
V_jump(1)         = (kappa*(LOG(rho_R*srho_L))-LOG(p_R/p_L))*skappaM1 &
                       - (beta_R*u2_R         -beta_L*u2_L  )
V_jump(2:4)       =  2.0*(beta_R*u_R(:)       -beta_L*u_L(:))        ! 2*beta*v
V_jump(5)         = -2.0*(beta_R              -beta_L       )        !-2*beta
V_jump(6:PP_nVar) =  2.0*(beta_R*UR(6:PP_nVar)-beta_L*UL(6:PP_nVar)) ! 2*beta*B

call EntropyStableFloGorFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)
RT = TRANSPOSE(Rmatrix)

! Compute entropy-stable fluxes
Fstar = Fstar - 0.5*MATMUL(Rmatrix,Dmatrix*MATMUL(RT,V_jump))

END ASSOCIATE
END SUBROUTINE EntropyStableFloGorFlux
!==================================================================================================================================
!> Computation of the entropy conserving flux and the dissipation matrices for the entropy stable FloGor solver
!==================================================================================================================================
PURE SUBROUTINE EntropyStableFloGorFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)
  USE MOD_Equation_Vars,ONLY:kappaM1, smu_0, s2mu_0, sKappaM1, Kappa
  USE MOD_Flux_Average, ONLY:LN_MEAN
#ifdef PP_GLM
  USE MOD_Equation_Vars,ONLY:GLM_ch
#endif
  implicit none
  !-arguments-----------------------------------------------
  real, intent(in)  :: UL(PP_nVar)      !< left state
  real, intent(in)  :: UR(PP_nVar)      !< right state
  real, intent(out) :: Fstar  (PP_nVar)
  real, intent(out) :: Dmatrix(PP_nVar)
  real, intent(out) :: Rmatrix(PP_nVar,PP_nVar)
  !-local-variables-----------------------------------------
  ! Mean states
  real    :: sbetaLN,betaAvg,rhoLN,pAvg,psiAvg,BAvg(3),uAvg(3)
  real    :: p_L,p_R,u_L(3), u_R(3), pTilde
  ! Additional
  REAL :: B2_L,B2_R,B2Avg,B2_ZIP
  REAL :: u2_L,u2_R,u2Avg,u2_ZIP
  REAL :: beta_R,beta_L
  REAL :: srho_L,srho_R
  REAL :: u1_B2Avg, uB_Avg
  real :: Tmatrix(PP_nVar)
  REAL :: aA,aLN,abeta,bb1A,bb2A,bb3A,bbA,aabbA,ca,xx,xxx,cf,cs,bperpA,beta2A,beta3A,alphaf,alphas,sgnb1
  REAL :: psiSplus,psiSminus,psiFplus,psiFminus, srhoLN, rhoAvg, pLN
  REAL :: in_e_L,in_e_R
  real :: LambdaMax, phi
  !---------------------------------------------------------

  ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
             rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
                E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
              psi_L =>UL(9)   ,  psi_R =>UR(9), &
#else
                E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
                B_L => UL(6:8),    B_R => UR(6:8)  )

! Get the inverse density, velocity, and pressure on left and right
srho_L = 1./rho_L
srho_R = 1./rho_R
u_L = rhoU_L(:)*srho_L
u_R = rhoU_R(:)*srho_R

u2_L = SUM(u_L(:)*u_L(:))
u2_R = SUM(u_R(:)*u_R(:))
B2_L = SUM(B_L(:)*B_L(:))
B2_R = SUM(B_R(:)*B_R(:))

p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*B2_L))
p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*B2_R))

! Specific inner energy
in_e_L = p_L*srho_L*skappaM1
in_e_R = p_R*srho_R*skappaM1

! Get the averages for the numerical flux
rhoLN  = LN_MEAN( rho_L, rho_R)
uAvg   = 0.5 * ( u_L +  u_R)
u2_ZIP = (u_L(1)*u_R(1)+u_L(2)*u_R(2)+u_L(3)*u_R(3))
BAvg   = 0.5 * ( B_L +  B_R)
B2_ZIP = (B_L(1)*B_R(1)+B_L(2)*B_R(2)+B_L(3)*B_R(3))
pAvg   = 0.5*(p_L+p_R)

#define ZIP(a,b,c,d) 0.5*(a*d+b*c)
! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rhoLN*uAvg(1)

Fstar(2) = Fstar(1)*uAvg(1)+pAvg+s2mu_0*B2_ZIP- smu_0*ZIP(B_L(1),B_R(1),B_L(1),B_R(1))
Fstar(3) = Fstar(1)*uAvg(2)                    - smu_0*ZIP(B_L(1),B_R(1),B_L(2),B_R(2))
Fstar(4) = Fstar(1)*uAvg(3)                    - smu_0*ZIP(B_L(1),B_R(1),B_L(3),B_R(3))

Fstar(5) = Fstar(1)*(0.5*u2_ZIP+in_e_L*in_e_R/LN_MEAN(in_e_L,in_e_R))+ZIP(p_L,p_R,u_L(1),u_R(1))  &
           + smu_0*(ZIP(u_L(1)*B_L(2),u_R(1)*B_R(2),B_L(2),B_R(2))-ZIP(u_L(2)*B_L(1),u_R(2)*B_R(1),B_L(2),B_R(2)) &
           +        ZIP(u_L(1)*B_L(3),u_R(1)*B_R(3),B_L(3),B_R(3))-ZIP(u_L(3)*B_L(1),u_R(3)*B_R(1),B_L(3),B_R(3)) &
#ifdef PP_GLM
           + GLM_ch*ZIP(B_L(1),B_R(1),psi_L,psi_R)  &
#endif /*PP_GLM*/
              )

#ifdef PP_GLM
Fstar(6) = GLM_ch*0.5*(psi_L+psi_R)
Fstar(9) = GLM_ch*BAvg(1)
#else
Fstar(6) = 0.
#endif /*PP_GLM*/
Fstar(7) = 0.5* ((u_L(1)*B_L(2)-u_L(2)*B_L(1)) + (u_R(1)*B_R(2)-u_R(2)*B_R(1)))
Fstar(8) = 0.5* ((u_L(1)*B_L(3)-u_L(3)*B_L(1)) + (u_R(1)*B_R(3)-u_R(3)*B_R(1)))

#undef ZIP

  !!
  ! MHD waves selective dissipation opertor components
  !!
    ! Compute additional averages
    beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
    beta_R = 0.5*rho_R/P_R
    sbetaLN = 1./LN_MEAN(beta_L,beta_R)
    betaAvg = 0.5 * (beta_L + beta_R)
    rhoAvg = 0.5*(rho_L + rho_R)
    pAvg = 0.5*(rho_L+rho_R)/(beta_L+beta_R) ! Redefinition of pAvg -> rhoAvg/(2*betaAvg)
    pLN = 0.5*rhoLN*sbetaLN
    u2Avg = u2_ZIP  ! Redefinition of u2Avg -> \bar{\norm{u}^2}
    srhoLN = 1./rhoLN
    psiAvg = 0.5*(psi_L+psi_R)
    aA = SQRT(kappa*pAvg/rhoLN)
    aLN = SQRT(kappa*pLN/rhoLN)
    abeta = SQRT(0.5*kappa/betaAvg)
    bb1A = BAvg(1)/SQRT(rhoLN)
    bb2A = BAvg(2)/SQRT(rhoLN)
    bb3A = BAvg(3)/SQRT(rhoLN)
    bbA = SQRT(bb1A*bb1A + bb2A*bb2A + bb3A*bb3A)
    aabbA = aA*aA + bbA*bbA

    ! Alfven speed
    ! Derigs et al. (2018), (4.63)
    ca = ABS(bb1A)

    ! Control round-off errors
    xx = aabbA*aabbA - 4.0*aA*aA*bb1A*bb1A
    if(xx .lt. 0.0) xx = 0.0
    xxx = aabbA + SQRT(xx)

    ! Fast magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cf = SQRT(0.5 * xxx)

    ! Control round-off errors
    xxx = aabbA - SQRT(xx)
    if(xxx .lt. 0.0) xxx = 0.0

    ! Slow magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cs = SQRT(0.5 * xxx)

    bperpA = SQRT(bb2A*bb2A + bb3A*bb3A)
    ! In case of very small bperpA, the values of betaA_{2,3}
    ! are indeterminable so we make them pairwise orthogonal
    ! Derigs et al. (2018), (4.63)
    if(bperpA .gt. 1.0E-14) then
      beta2A = bb2A/bperpA
      beta3A = bb3A/bperpA
    else
      bperpA = 0.0
      beta2A = 1.0/sqrt(2.0)
      beta3A = 1.0/sqrt(2.0)
    endif

    ! Avoid negative round-off errors when computing alphaf (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (aA*aA-cs*cs)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphaf = SQRT(xx)
    else
      alphaf = 0.0
    end if

    ! Avoid negative round-off errors when computing alphas (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (cf*cf-aA*aA)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphas = SQRT(xx)
    else
      alphas = 0.0
    end if

    sgnb1 = sign(1.,bb1A)

    ! Derigs et al. (2018), (4.63)
    psiSplus =  0.5*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 + alphas*cs*rhoLN*uAvg(1) + &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiSminus = 0.5*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 - alphas*cs*rhoLN*uAvg(1) - &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    psiFplus =  0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 + alphaf*cf*rhoLN*uAvg(1) - &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiFminus = 0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 - alphaf*cf*rhoLN*uAvg(1) + &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    ! + fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,1) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) + cf), &
                      rhoLN*(alphaf*uAvg(2) - alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) - alphas*cs*beta3A*sgnb1), &
                      psiFplus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,2) = (/ 0.0, &
                      0.0,  &
                      rhoLN*SQRT(rhoAvg)*beta3A, &
                      -rhoLN*SQRT(rhoAvg)*beta2A, &
                      -rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! + slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,3) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) + cs), &
                      rhoLN*(alphas*uAvg(2) + alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) + alphaf*cf*beta3A*sgnb1), &
                      psiSplus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,4) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) + psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      1.0 /)

    ! Entropy wave
    ! Derigs et al. (2018), (4.66)
    Rmatrix(:,5) = (/ 1.0, &
                      uAvg(1), &
                      uAvg(2), &
                      uAvg(3), &
                      0.5*u2Avg, &
                      0.0, &
                      0.0, &
                      0.0, &
                      0.0 /)

    ! - GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,6) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) - psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      -1.0 /)

    ! - slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,7) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) - cs), &
                      rhoLN*(alphas*uAvg(2) - alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) - alphaf*cf*beta3A*sgnb1), &
                      psiSminus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! - Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,8) = (/ 0.0, &
                      0.0, &
                      -rhoLN*SQRT(rhoAvg)*beta3A, &
                      rhoLN*SQRT(rhoAvg)*beta2A, &
                      rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! - fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,9) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) - cf), &
                      rhoLN*(alphaf*uAvg(2) + alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) + alphas*cs*beta3A*sgnb1), &
                      psiFminus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! A blend of the 9waves solver and LLF
      phi = sqrt(ABS(1.0-(p_R*rho_R/(p_L*rho_L)))/(1.0+(p_R*rho_R/(p_L*rho_L))))
      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))

      Dmatrix(1) = (1.-phi)*ABS( uAvg(1) + cf ) + phi*LambdaMax
      Dmatrix(2) = (1.-phi)*ABS( uAvg(1) + ca ) + phi*LambdaMax
      Dmatrix(3) = (1.-phi)*ABS( uAvg(1) + cs ) + phi*LambdaMax
      Dmatrix(4) = (1.-phi)*ABS( uAvg(1) + GLM_ch ) + phi*LambdaMax
      Dmatrix(5) = (1.-phi)*ABS( uAvg(1)      ) + phi*LambdaMax
      Dmatrix(6) = (1.-phi)*ABS( uAvg(1) - GLM_ch ) + phi*LambdaMax
      Dmatrix(7) = (1.-phi)*ABS( uAvg(1) - cs ) + phi*LambdaMax
      Dmatrix(8) = (1.-phi)*ABS( uAvg(1) - ca ) + phi*LambdaMax
      Dmatrix(9) = (1.-phi)*ABS( uAvg(1) - cf ) + phi*LambdaMax

    ! Pure 9waves solver (no LLF)
!#
!#      Dmatrix(1) = ABS( uAvg(1) + cf ) ! + fast magnetoacoustic wave
!#      Dmatrix(2) = ABS( uAvg(1) + ca ) ! + Alfven wave
!#      Dmatrix(3) = ABS( uAvg(1) + cs ) ! + slow magnetoacoustic wave
!#      Dmatrix(4) = ABS( uAvg(1) + GLM_ch ) ! + GLM wave
!#      Dmatrix(5) = ABS( uAvg(1)      ) ! / Entropy wave
!#      Dmatrix(6) = ABS( uAvg(1) - GLM_ch ) ! - GLM wave
!#      Dmatrix(7) = ABS( uAvg(1) - cs ) ! - slow magnetoacoustic wave
!#      Dmatrix(8) = ABS( uAvg(1) - ca ) ! - Alfven wave
!#      Dmatrix(9) = ABS( uAvg(1) - cf ) ! - fast magnetoacoustic wave

    ! LLF solver
!#      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))
!#
!#      Dmatrix(1) = LambdaMax
!#      Dmatrix(2) = LambdaMax
!#      Dmatrix(3) = LambdaMax
!#      Dmatrix(4) = LambdaMax
!#      Dmatrix(5) = LambdaMax
!#      Dmatrix(6) = LambdaMax
!#      Dmatrix(7) = LambdaMax
!#      Dmatrix(8) = LambdaMax
!#      Dmatrix(9) = LambdaMax

    ! Diagonal scaling matrix as described in Winters et al., eq. (4.15)
    Tmatrix = 0.0
    Tmatrix(1) = 0.5/kappa/rhoLN ! + f
    Tmatrix(2) = 0.25/betaAvg/rhoLN/rhoLN ! + a
    Tmatrix(3) = Tmatrix(1) ! + s
    Tmatrix(4) = 0.25 / betaAvg ! + GLM
    Tmatrix(5) = rhoLN*kappaM1/kappa ! E
    Tmatrix(6) = Tmatrix(4) ! - GLM
    Tmatrix(7) = Tmatrix(1) ! - s
    Tmatrix(8) = Tmatrix(2) ! - a
    Tmatrix(9) = Tmatrix(1) ! - f

    ! Scale D matrix
    Dmatrix = Dmatrix*Tmatrix
END ASSOCIATE
END SUBROUTINE EntropyStableFloGorFlux_VolFluxAndDissipMatrices
#endif

END MODULE MOD_Riemann
