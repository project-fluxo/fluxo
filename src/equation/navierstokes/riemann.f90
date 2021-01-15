!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!> Contains riemann solvers for Euler equations, include call to diffusion fluxes
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE RiemannSolverCentral
  MODULE PROCEDURE RiemannSolverCentral
END INTERFACE

INTERFACE RiemannSolverByRusanov
  MODULE PROCEDURE RiemannSolverByRusanov
END INTERFACE

INTERFACE RiemannSolverByHLL
  MODULE PROCEDURE RiemannSolverByHLL
END INTERFACE

INTERFACE RiemannSolverByHLLC
  MODULE PROCEDURE RiemannSolverByHLLC
END INTERFACE

INTERFACE RiemannSolverByRoe
  MODULE PROCEDURE RiemannSolverByRoe
END INTERFACE

INTERFACE RiemannSolver_EntropyStable
  MODULE PROCEDURE RiemannSolver_EntropyStable
END INTERFACE

INTERFACE RiemannSolver_VolumeFluxAverage
  MODULE PROCEDURE RiemannSolver_VolumeFluxAverage
END INTERFACE

INTERFACE RiemannSolver_EC_LLF
  MODULE PROCEDURE RiemannSolver_EC_LLF
END INTERFACE

INTERFACE RiemannSolver_VolumeFluxAverage_LLF
  MODULE PROCEDURE RiemannSolver_VolumeFluxAverage_LLF
END INTERFACE

INTERFACE RiemannSolver_ECKEP_LLF
  MODULE PROCEDURE RiemannSolver_ECKEP_LLF
END INTERFACE

PUBLIC:: Riemann, AdvRiemann
PUBLIC:: RiemannSolverCentral
PUBLIC:: RiemannSolverByRusanov
PUBLIC:: RiemannSolverByHLL
PUBLIC:: RiemannSolverByHLLC
PUBLIC:: RiemannSolverByRoe
PUBLIC:: RiemannSolver_EntropyStable
PUBLIC:: RiemannSolver_EntropyStable_VolFluxAndDissipMatrices
PUBLIC:: RiemannSolver_VolumeFluxAverage
PUBLIC:: RiemannSolver_EC_LLF
PUBLIC:: RiemannSolver_VolumeFluxAverage_LLF
PUBLIC:: RiemannSolver_ECKEP_LLF
PUBLIC:: RotateState
PUBLIC:: RotateFluxBack
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> main routine calling different riemann solvers
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann(F,U_L,U_R, &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif /*PARABOLIC*/
                   nv,t1,t2)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars   ,ONLY:SolveRiemannProblem
USE MOD_Flux_Average
!USE MOD_Analyze_Vars  ,ONLY: wGPSurf ! ECMORTAR
!USE MOD_Equation_Vars  ,ONLY: ConsToPrim, ConsToEntropy, kappa ! ECMORTAR
!USE MOD_Equation_Vars ,ONLY:kappaM1 ! ECMORTAR
#if PARABOLIC
USE MOD_Flux            ,ONLY:EvalDiffFlux3D    ! and the NSE diffusion fluxes in all directions to approximate the numerical flux
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: U_L(     PP_nVar,0:PP_N,0:PP_N) !<  left state on face, not rotated
REAL,INTENT(IN) :: U_R(     PP_nVar,0:PP_N,0:PP_N) !< right state on face, not rotated
#if PARABOLIC                                                 
REAL,INTENT(IN) :: gradUx_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in x, not rotated 
REAL,INTENT(IN) :: gradUx_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in x, not rotated 
REAL,INTENT(IN) :: gradUy_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in y, not rotated 
REAL,INTENT(IN) :: gradUy_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in y, not rotated 
REAL,INTENT(IN) :: gradUz_L(PP_nVar,0:PP_N,0:PP_N) !<  left state gradient in z, not rotated 
REAL,INTENT(IN) :: gradUz_R(PP_nVar,0:PP_N,0:PP_N) !< right state gradient in z, not rotated 
#endif /*PARABOLIC*/
REAL,INTENT(IN) :: nv(            3,0:PP_N,0:PP_N) !< normal vector of face
REAL,INTENT(IN) :: t1(            3,0:PP_N,0:PP_N) !< 1st tangential vector of face
REAL,INTENT(IN) :: t2(            3,0:PP_N,0:PP_N) !< 2nd tangential vector of face
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT):: F(       PP_nVar,0:PP_N,0:PP_N) !< numerical flux on face
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                       :: i,j,iVar
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N)         :: U_LL,U_RR
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N)         :: g_L,g_R,k_L,k_R,j_L,j_R
#endif
!REAL                                          :: Sum1, sum2, sum3, sum3_1, sum3_2, suma, sumb ! ECMORTAR
!REAL                                          :: math_entropy_R, math_entropy_L, prim(PP_nVar), sRho, pres,v1,v2,v3 ! ECMORTAR
!REAL                                          :: Flux_F_R(PP_nVar), Flux_F_L(PP_nVar) ! ECMORTAR
!REAL                                          :: entropy_vars_R(PP_nVar), entropy_vars_L(PP_nVar) ! ECMORTAR
!REAL                                          :: entropy_flux_R, entropy_flux_L ! ECMORTAR
!==================================================================================================================================

call AdvRiemann(F,U_L,U_R,nv,t1,t2)

#if PARABOLIC
!! Don#t forget the diffusion contribution, my young padawan
!! Compute NSE Diffusion flux in cartesian coordinates
CALL EvalDiffFlux3D(k_L,g_L,j_L,U_L,gradUx_L,gradUy_L,gradUz_L)

CALL EvalDiffFlux3D(k_R,g_R,j_R,U_R,gradUx_R,gradUy_R,gradUz_R)
!
! !BR1/BR2 uses arithmetic mean of the fluxes
DO iVar=2,PP_nVar
  F(iVar,:,:)=F(iVar,:,:)+0.5*( nv(1,:,:)*(k_L(iVar,:,:)+k_R(iVar,:,:)) &
                               +nv(2,:,:)*(g_L(iVar,:,:)+g_R(iVar,:,:)) &
                               +nv(3,:,:)*(j_L(iVar,:,:)+j_R(iVar,:,:)))
END DO
#endif /* PARABOLIC */
END SUBROUTINE Riemann

!==================================================================================================================================
!> Advective Riemann solver
!==================================================================================================================================
SUBROUTINE AdvRiemann(F,U_L,U_R,nv,t1,t2)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars   ,ONLY:SolveRiemannProblem
USE MOD_Flux_Average
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: U_L(     PP_nVar,0:PP_N,0:PP_N) !<  left state on face, not rotated
REAL,INTENT(IN) :: U_R(     PP_nVar,0:PP_N,0:PP_N) !< right state on face, not rotated
REAL,INTENT(IN) :: nv(            3,0:PP_N,0:PP_N) !< normal vector of face
REAL,INTENT(IN) :: t1(            3,0:PP_N,0:PP_N) !< 1st tangential vector of face
REAL,INTENT(IN) :: t2(            3,0:PP_N,0:PP_N) !< 2nd tangential vector of face
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT):: F(       PP_nVar,0:PP_N,0:PP_N) !< numerical flux on face
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                       :: i,j,iVar
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N)         :: U_LL,U_RR
!==================================================================================================================================
! Momentum has to be rotated using the normal system individual for each
! Gauss point i,j
DO j=0,PP_N
  DO i=0,PP_N
    U_LL(:,i,j) = RotateState(U_L(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))
    U_RR(:,i,j) = RotateState(U_R(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))
  END DO ! i 
END DO ! j

CALL SolveRiemannProblem(F,U_LL,U_RR)

! Back Rotate the normal flux into Cartesian direction
DO j=0,PP_N
  DO i=0,PP_N
    call RotateFluxBack(F(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j))
  END DO ! i
END DO ! j

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
  
  rotU(1) =     U(1  )
  rotU(2) = SUM(U(2:4)*nv(:))
  rotU(3) = SUM(U(2:4)*t1(:))
  rotU(4) = SUM(U(2:4)*t2(:))
  rotU(5) =     U(5  )
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
  
end subroutine RotateFluxBack
!==================================================================================================================================
!> Central / Average Euler flux
!==================================================================================================================================
SUBROUTINE RiemannSolverCentral(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Flux         ,ONLY:EvalEulerFlux1D   ! we use the Euler fluxes in normal direction to approximate the numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
!==================================================================================================================================
! Euler Fluxes
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)

F=0.5*(F_L+F_R)

END SUBROUTINE RiemannSolverCentral


!==================================================================================================================================
!> Rusanov / lax-Friedrichs Riemann solver
!==================================================================================================================================
SUBROUTINE RiemannSolverByRusanov(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:SoundSpeed2
USE MOD_Flux         ,ONLY:EvalEulerFlux1D   ! we use the Euler fluxes in normal direction to approximate the numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER :: iVar,i,j
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
REAL    :: LambdaMax(0:PP_N,0:PP_N)
!==================================================================================================================================
DO j=0,PP_N;  DO i=0,PP_N
  LambdaMax(i,j)=MAX(ABS(U_LL(2,i,j)/U_LL(1,i,j)),ABS(U_RR(2,i,j)/U_RR(1,i,j))) &
                 +SQRT(MAX(SoundSpeed2(U_LL(:,i,j)),SoundSpeed2(U_RR(:,i,j))) )
END DO; END DO

! Euler Fluxes
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)
! Rusanov /Lax-Friedrichs flux
DO iVar=1,PP_nVar
  ! compute flux
  ! f=0.5*(f(U_l)+f(U_r))-0.5*lambdamax*(U_r-U_l)
  F(iVar,:,:)=0.5*((F_L(iVar,:,:)+F_R(iVar,:,:))-LambdaMax(:,:)*(U_RR(iVar,:,:)-U_LL(iVar,:,:)))
END DO
END SUBROUTINE RiemannSolverByRusanov



!==================================================================================================================================
!> HLL Riemann solver
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLL(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
USE MOD_Flux         ,ONLY:EvalEulerFlux1D   ! we use the Euler fluxes in normal direction to approximate the numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER :: i,j
REAL    :: ssl,ssr,SStar,sMu_L,sMu_R,U_Star(1:5)
REAL    :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R, H_L, H_R, HTilde, vTilde, aTilde
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
!==================================================================================================================================
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)
DO j=0,PP_N
  DO i=0,PP_N
    sRho_L=1./U_LL(1,i,j)
    sRho_R=1./U_RR(1,i,j)
    Vel_L =U_LL(2:4,i,j)*sRho_L
    Vel_R =U_RR(2:4,i,j)*sRho_R
    p_L   =KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R   =KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
!   
!   Davis (not recommended according to Toro)
!   -----------------------------------------
!~     Ssl = Vel_L(1) - SQRT(kappa*p_L*sRho_L)
!~     Ssr = Vel_R(1) + SQRT(kappa*p_R*sRho_R)
!
!   Davis and Einfeldt
!   ------------------
    !! enthalpy
    H_L = (U_LL(5,i,j) + p_L)/U_LL(1,i,j)
    H_R = (U_RR(5,i,j) + p_R)/U_RR(1,i,j)

    !! Roe averages
    HTilde = (SQRT(U_LL(1,i,j))*H_L      + SQRT(U_RR(1,i,j))*H_R     )/(SQRT(U_LL(1,i,j)) + SQRT(U_RR(1,i,j)))
    vTilde = (SQRT(U_LL(1,i,j))*Vel_L(1) + SQRT(U_RR(1,i,j))*Vel_R(1))/(SQRT(U_LL(1,i,j)) + SQRT(U_RR(1,i,j)))

    aTilde = SQRT((kappa-1.)*(HTilde - 0.5*vTilde**2))

    Ssl = vTilde - aTilde
    Ssr = vTilde + aTilde
    
    ! positive supersonic speed
    IF(Ssl .GE. 0. .and. Ssr .GT. 0.)THEN
      F(:,i,j)=F_L(:,i,j)
    ! negative supersonic speed
    ELSEIF(Ssr .LE. 0. .and. Ssl .LT. 0.)THEN
      F(:,i,j)=F_R(:,i,j)
    ! subsonic case
    ELSE
      F(:,i,j) = (ssR*F_L(:,i,j) - ssL*F_R(:,i,j) + ssL*ssR*(U_RR(:,i,j)-U_LL(:,i,j)))/(ssR-ssL)
    END IF ! subsonic case
  END DO ! i 
END DO ! j
END SUBROUTINE RiemannSolverByHLL


!==================================================================================================================================
!> HLLC Riemann solver
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLLC(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
USE MOD_Flux         ,ONLY:EvalEulerFlux1D   ! we use the Euler fluxes in normal direction to approximate the numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER :: i,j
REAL    :: ssl,ssr,SStar,sMu_L,sMu_R,U_Star(1:5)
REAL    :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
!==================================================================================================================================
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)
DO j=0,PP_N
  DO i=0,PP_N
    sRho_L=1./U_LL(1,i,j)
    sRho_R=1./U_RR(1,i,j)
    Vel_L =U_LL(2:4,i,j)*sRho_L
    Vel_R =U_RR(2:4,i,j)*sRho_R
    p_L   =KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R   =KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
    Ssl = Vel_L(1) - SQRT(kappa*p_L*sRho_L)
    Ssr = Vel_R(1) + SQRT(kappa*p_R*sRho_R)
    ! positive supersonic speed
    IF(Ssl .GE. 0.)THEN
      F(:,i,j)=F_L(:,i,j)
    ! negative supersonic speed
    ELSEIF(Ssr .LE. 0.)THEN
      F(:,i,j)=F_R(:,i,j)
    ! subsonic case
    ELSE
      sMu_L = Ssl - Vel_L(1)
      sMu_R = Ssr - Vel_R(1)
      SStar = (p_R - p_L +                    &
               U_LL(2,i,j)*sMu_L - U_RR(2,i,j)*sMu_R) / &
              (U_LL(1,i,j)*sMu_L - U_RR(1,i,j)*sMu_R)
      IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
        U_Star(:) = U_LL(1,i,j) * sMu_L/(Ssl-SStar)*                  &
          (/ 1., SStar, U_LL(3:4,i,j)*sRho_L,                         &
             U_LL(PP_nVar,i,j)*sRho_L + (SStar-Vel_L(1))*             &
            (SStar + p_L*sRho_L/sMu_L)                              /)
        F(:,i,j)=F_L(:,i,j)+Ssl*(U_Star(:)-U_LL(:,i,j))
      ELSE
        U_Star(:) = U_RR(1,i,j) * sMu_R/(Ssr-SStar)*                  &
          (/ 1., SStar, U_RR(3:4,i,j)*sRho_R,                         &
             U_RR(PP_nVar,i,j)*sRho_R + (SStar-Vel_R(1))*             &
            (SStar + p_R*sRho_R/sMu_R)                              /)
        F(:,i,j)=F_R(:,i,j)+Ssr*(U_Star(:)-U_RR(:,i,j))
      END IF
    END IF ! subsonic case
  END DO ! i 
END DO ! j
END SUBROUTINE RiemannSolverByHLLC


!==================================================================================================================================
!> HLLC Riemann solver with Carbuncle fix
!==================================================================================================================================
SUBROUTINE RiemannSolverByHLLCnoCarbuncle(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
USE MOD_Flux         ,ONLY:EvalEulerFlux1D   ! we use the Euler fluxes in normal direction to approximate the numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER :: i,j
REAL    :: ssl,ssr,SStar,sMu_L,sMu_R,U_Star(1:5)
REAL    :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
real,DIMENSION(1:5) :: F_HLLE, F_HLLC
real :: smooth, a_L, a_R, H_L, H_R, HTilde, vTilde, aTilde
!==================================================================================================================================
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)
DO j=0,PP_N
  DO i=0,PP_N
    sRho_L=1./U_LL(1,i,j)
    sRho_R=1./U_RR(1,i,j)
    Vel_L =U_LL(2:4,i,j)*sRho_L
    Vel_R =U_RR(2:4,i,j)*sRho_R
    p_L   =KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R   =KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
    
!   ***************************    
!   Direct wave speed estimates
!   ***************************
    
!   Davis estimates:
!   ----------------   ...not recommended for practical computations [Toro,2009]
!~     a_L   =SQRT(kappa*p_L*sRho_L)
!~     a_R   =SQRT(kappa*p_R*sRho_R)
!~     !(1)
!~     Ssl = Vel_L(1) - a_L
!~     Ssr = Vel_R(1) + a_R
!~     !(2)
!~     Ssl = min( Vel_L(1) - a_L, Vel_R(1) - a_R )
!~     Ssr = max( Vel_L(1) + a_L, Vel_R(1) + a_R )
    
!   Davis (1988) and Einfeldt (1988):
!   ---------------------------------
    !! enthalpy
    H_L = (U_LL(5,i,j) + p_L)/U_LL(1,i,j)
    H_R = (U_RR(5,i,j) + p_R)/U_RR(1,i,j)

    !! Roe averages
    HTilde = (SQRT(U_LL(1,i,j))*H_L      + SQRT(U_RR(1,i,j))*H_R     )/(SQRT(U_LL(1,i,j)) + SQRT(U_RR(1,i,j)))
    vTilde = (SQRT(U_LL(1,i,j))*Vel_L(1) + SQRT(U_RR(1,i,j))*Vel_R(1))/(SQRT(U_LL(1,i,j)) + SQRT(U_RR(1,i,j)))

    aTilde = SQRT((kappa-1.)*(HTilde - 0.5*vTilde**2))

    Ssl = vTilde - aTilde
    Ssr = vTilde + aTilde
    
!   ***************
!   Riemann solver!
!   ***************
    
    ! positive supersonic speed
    IF(Ssl .GE. 0. .and. Ssr .GT. 0.)THEN
      F(:,i,j)=F_L(:,i,j)
    ! negative supersonic speed
    ELSEIF(Ssr .LE. 0. .and. Ssl .LT. 0.)THEN
      F(:,i,j)=F_R(:,i,j)
    ! subsonic case
    ELSE
      sMu_L = Ssl - Vel_L(1)
      sMu_R = Ssr - Vel_R(1)
      SStar = (p_R - p_L +                    &
               U_LL(2,i,j)*sMu_L - U_RR(2,i,j)*sMu_R) / &
              (U_LL(1,i,j)*sMu_L - U_RR(1,i,j)*sMu_R)
      IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
        U_Star(:) = U_LL(1,i,j) * sMu_L/(Ssl-SStar)*                  &
          (/ 1., SStar, U_LL(3:4,i,j)*sRho_L,                         &
             U_LL(PP_nVar,i,j)*sRho_L + (SStar-Vel_L(1))*             &
            (SStar + p_L*sRho_L/sMu_L)                              /)
        
        ! Original HLLC:
        !F(:,i,j)=F_L(:,i,j)+Ssl*(U_Star(:)-U_LL(:,i,j))
        
        ! For Carbuncle:
        F_HLLE = (Ssr*F_L(:,i,j) - Ssl*F_R(:,i,j) + Ssl*Ssr*(U_RR(:,i,j) - U_LL(:,i,j))) / (Ssr - Ssl)
        F_HLLC = F_L(:,i,j) + Ssl*(U_Star(:) - U_LL(:,i,j))

        smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))

        F(:,i,j) = smooth*F_HLLE + (1.-smooth)*F_HLLC

        
      ELSE
        U_Star(:) = U_RR(1,i,j) * sMu_R/(Ssr-SStar)*                  &
          (/ 1., SStar, U_RR(3:4,i,j)*sRho_R,                         &
             U_RR(PP_nVar,i,j)*sRho_R + (SStar-Vel_R(1))*             &
            (SStar + p_R*sRho_R/sMu_R)                              /)
        
        ! Original HLLC:
        !F(:,i,j)=F_R(:,i,j)+Ssr*(U_Star(:)-U_RR(:,i,j))
        
        ! For Carbuncle:
        F_HLLE = (Ssr*F_L(:,i,j) - Ssl*F_R(:,i,j) + Ssl*Ssr*(U_RR(:,i,j) - U_LL(:,i,j))) / (Ssr - Ssl)
        F_HLLC = F_R(:,i,j) + Ssr*(U_Star(:) - U_RR(:,i,j))

        smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))

        F(:,i,j) = smooth*F_HLLE + (1.-smooth)*F_HLLC
      END IF
    END IF ! subsonic case
  END DO ! i 
END DO ! j
END SUBROUTINE RiemannSolverByHLLCnoCarbuncle

!==================================================================================================================================
!> Roe riemann solver, depending on "whichvolumeflux" global parameter, consistent dissipation is added
!==================================================================================================================================
SUBROUTINE RiemannSolverByRoe(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: kappaM1,WhichVolumeFlux
USE MOD_Flux         , ONLY: EvalEulerFlux1D 
USE MOD_Flux_Average , ONLY: KennedyAndGruberFlux1,DucrosFlux, KennedyAndGruberFlux2
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: i,j
REAL                 :: sRho_L,sRho_R,Vel_L(3),Vel_R(3)
REAL                 :: F_c(5),uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
REAL                 :: H_L,H_R,SqrtRho_L,SqrtRho_R,sSqrtRho,RoeVel(1:3),RoeH,RoeC
REAL                 :: a1,a2,a3,a4,a5
REAL                 :: alpha1,alpha2,alpha3,alpha4,alpha5
REAL                 :: Delta_U(1:6)
REAL,DIMENSION(1:5)  :: r1,r2,r3,r4,r5
REAL,DIMENSION(1:5,0:PP_N,0:PP_N) :: F_L,F_R
!==================================================================================================================================
CALL EvalEulerFlux1D(U_LL,F_L)
CALL EvalEulerFlux1D(U_RR,F_R)
DO j=0,PP_N
  DO i=0,PP_N
    sRho_L    = 1./U_LL(1,i,j)
    sRho_R    = 1./U_RR(1,i,j)
    Vel_L     = U_LL(2:4,i,j)*sRho_L
    Vel_R     = U_RR(2:4,i,j)*sRho_R
    H_L       = (U_LL(PP_nVar,i,j)+KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3))))*sRho_L
    H_R       = (U_RR(PP_nVar,i,j)+KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3))))*sRho_R
    SqrtRho_L = SQRT(U_LL(1,i,j))
    SqrtRho_R = SQRT(U_RR(1,i,j))
    sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
    ! Roe mean values
    RoeVel(:) = (SqrtRho_R*Vel_R(:) + SqrtRho_L*Vel_L(:)) * sSqrtRho
    RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
    Roec      = SQRT(kappaM1*(RoeH-0.5*SUM(RoeVel(:)*RoeVel(:))))
    ! mean eigenvalues
    a1    = RoeVel(1)-Roec
    a2    = RoeVel(1)
    a3    = RoeVel(1)
    a4    = RoeVel(1)
    a5    = RoeVel(1)+Roec
    ! mean eigenvectors
    r1(1) = 1.
    r1(2) = a1
    r1(3) = RoeVel(2)
    r1(4) = RoeVel(3)
    r1(5) = RoeH-RoeVel(1)*Roec

    r2(1) = 1.
    r2(2) = RoeVel(1)
    r2(3) = RoeVel(2)
    r2(4) = RoeVel(3)
    r2(5) = 0.5*SUM(RoeVel(:)*RoeVel(:))

    r3(1) = 0.
    r3(2) = 0.
    r3(3) = 1.
    r3(4) = 0.
    r3(5) = RoeVel(2)

    r4(1) = 0.
    r4(2) = 0.
    r4(3) = 0.
    r4(4) = 1.
    r4(5) = RoeVel(3)

    r5(1) = 1.
    r5(2) = a5
    r5(3) = RoeVel(2)
    r5(4) = RoeVel(3)
    r5(5) = RoeH+RoeVel(1)*Roec
    ! calculate differences
    Delta_U(1:5) = U_RR(1:5,i,j) - U_LL(1:5,i,j)
    Delta_U(6)   = Delta_U(5)-(Delta_U(3)-RoeVel(2)*Delta_U(1))*RoeVel(2) - (Delta_U(4)-RoeVel(3)*Delta_U(1))*RoeVel(3)
    ! calculate factors
    Alpha3 = Delta_U(3) - RoeVel(2)*Delta_U(1)
    Alpha4 = Delta_U(4) - RoeVel(3)*Delta_U(1)
    Alpha2 = kappaM1/(Roec*Roec) * (Delta_U(1)*(RoeH-RoeVel(1)*RoeVel(1)) - Delta_U(6) + RoeVel(1)*Delta_U(2))
    Alpha1 = 0.5/Roec * (Delta_U(1)*(RoeVel(1)+Roec) - Delta_U(2) - Roec*Alpha2)
    Alpha5 = Delta_U(1) - Alpha1 - Alpha2
    ! assemble Roe flux
    ! Get the baseline flux needed for the standard Roe solver
    SELECT CASE(WhichVolumeFlux)
    CASE DEFAULT
      F_c=0.5*(F_L(:,i,j)+F_R(:,i,j))
    CASE(2)
      CALL KennedyAndGruberFlux1(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                                 uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
      a1 = MAX(ABS(RoeVel(1)-Roec),ABS(RoeVel(1)+Roec)) ! ensures consistent KE dissipation
      a5 = MAX(ABS(RoeVel(1)-Roec),ABS(RoeVel(1)+Roec))
    CASE(3)
      CALL DucrosFlux(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                      uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
    CASE(4)
      CALL KennedyAndGruberFlux2(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                                 uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
      a1 = MAX(ABS(RoeVel(1)-Roec),ABS(RoeVel(1)+Roec)) ! ensures consistent KE dissipation
      a5 = MAX(ABS(RoeVel(1)-Roec),ABS(RoeVel(1)+Roec))
    END SELECT   
    F(:,i,j) = F_c - 0.5*(Alpha1*ABS(a1)*r1(:) + &
                                      Alpha2*ABS(a2)*r2(:) + &
                                      Alpha3*ABS(a3)*r3(:) + &
                                      Alpha4*ABS(a4)*r4(:) + &
                                      Alpha5*ABS(a5)*r5(:))
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolverByRoe

!==================================================================================================================================
!> Entropy stable Riemann solver with full matrix dissipation (uses TwoPoint Entropy Conserving flux)
!==================================================================================================================================
SUBROUTINE RiemannSolver_EntropyStable(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY: kappa,KappaM1
!USE MOD_Flux_Average ,ONLY: TwoPointEntropyConservingFlux
USE MOD_Equation_Vars,ONLY:VolumeFluxAverage
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: j,i,iVar
REAL                 :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R
REAL                 :: F_c(5),uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
REAL                 :: sR,sL,rho_pR,rho_pL 
REAL                 :: Rhat(5,5)
REAL,DIMENSION(1:5)  :: Dhat,vR,vL,vjump,diss
!==================================================================================================================================
DO j=0,PP_N
  DO i=0,PP_N
! Compute entropy conserving flux and the dissipation matrices
    call RiemannSolver_EntropyStable_VolFluxAndDissipMatrices(U_LL(:,i,j),U_RR(:,i,j),F_c,Dhat,Rhat)
    
! Get the entropy variables locally
    sRho_L = 1./U_LL(1,i,j)
    sRho_R = 1./U_RR(1,i,j)
    Vel_L  = U_LL(2:4,i,j)*sRho_L
    Vel_R  = U_RR(2:4,i,j)*sRho_R
    p_L    = KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R    = KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
    sR     =  LOG(p_R) - kappa*LOG(U_RR(1,i,j))
    sL     =  LOG(p_L) - kappa*LOG(U_LL(1,i,j))
    rho_pR =  U_RR(1,i,j)/p_R
    rho_pL =  U_LL(1,i,j)/p_L
    vL(1)  =  (kappa-sL)/(kappaM1) - 0.5*rho_pL*(SUM(Vel_L(:)*Vel_L(:)))
    vR(1)  =  (kappa-sR)/(kappaM1) - 0.5*rho_pR*(SUM(Vel_R(:)*Vel_R(:)))
    vL(2:4)=  rho_pL*Vel_L(1:3)
    vR(2:4)=  rho_pR*Vel_R(1:3)
    vR(5)  = -rho_pR
    vL(5)  = -rho_pL
! Compute the dissipation term RHat*DHat*RHat^T*[v]
    vJump = vR - vL
    diss  = RHat(1,:)*vJump(1) + RHat(2,:)*vJump(2) + RHat(3,:)*vJump(3) + RHat(4,:)*vJump(4) + RHat(5,:)*vJump(5)
    DO iVar = 1,5
      diss(iVar) = DHat(iVar)*diss(iVar)
    END DO
    diss = RHat(:,1)*diss(1) + RHat(:,2)*diss(2) + RHat(:,3)*diss(3) + RHat(:,4)*diss(4) + RHat(:,5)*diss(5)
! Compute entropy stable numerical flux
    F(:,i,j) = F_c - 0.5*diss
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolver_EntropyStable
!==================================================================================================================================
!> Volume flux and dissipation matrices evaluation for entropy stable Riemann solver with full matrix dissipation 
!> (uses TwoPoint Entropy Conserving flux)
!==================================================================================================================================
PURE SUBROUTINE RiemannSolver_EntropyStable_VolFluxAndDissipMatrices(U_LL,U_RR,F,Dhat,Rhat)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY: kappa,KappaM1
!USE MOD_Flux_Average ,ONLY: TwoPointEntropyConservingFlux
USE MOD_Equation_Vars,ONLY:VolumeFluxAverage
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5)    ,INTENT(OUT)  :: F         !< numerical flux
REAL,DIMENSION(1:5)    ,INTENT(OUT)  :: Dhat      !< Dissipation matrix
REAL,DIMENSION(1:5,1:5),INTENT(OUT)  :: Rhat      !< Right-eigenvector matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
!==================================================================================================================================

CALL VolumeFluxAverage(F,U_LL,U_RR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)

! Matrix of right eigenvectors
RHat(1,:) = (/ 1.               , 1.                                      , 0.   ,  0.  , 1.               /)
RHat(2,:) = (/ uHat - aHat      , uHat                                    , 0.   ,  0.  , uHat + aHat      /)
RHat(3,:) = (/ vHat             , vHat                                    , 1.   ,  0.  , vHat             /)
RHat(4,:) = (/ wHat             , wHat                                    , 0.   ,  1.  , wHat             /)
RHat(5,:) = (/ HHat - uHat*aHat , 0.5*(uHat*uHat + vHat*vHat + wHat*wHat) , vHat , wHat , HHat + uHat*aHat /)
! Diagonal scaling matrix where DHat = ABS(\Lambda)S
DHat = 0.0
DHat(1) = 0.5*ABS(uHat - aHat)*rhoHat/kappa
DHat(2) = ABS(uHat)*(kappaM1/kappa)*rhoHat!*rhoHat*rhoHat
DHat(3) = ABS(uHat)*p1Hat!*rhoHat
DHat(4) = DHat(3)
DHat(5) = 0.5*ABS(uHat + aHat)*rhoHat/kappa
    
!----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE RiemannSolver_EntropyStable_VolFluxAndDissipMatrices

!==================================================================================================================================
!> Riemann solver, using simply the two point average chosen by the global parameter whichVolumeFlux
!==================================================================================================================================
SUBROUTINE RiemannSolver_VolumeFluxAverage(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:VolumeFluxAverage
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: i,j
REAL                 :: uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat !dummy
!==================================================================================================================================
DO j=0,PP_N
  DO i=0,PP_N
! Compute entropy conserving flux
    CALL VolumeFluxAverage(F(:,i,j),U_LL(:,i,j),U_RR(:,i,j),&
                                       uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolver_VolumeFluxAverage


!==================================================================================================================================
!> Entropy conservative Ismali and Roe with LLF diss
!==================================================================================================================================
SUBROUTINE RiemannSolver_EC_LLF(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1,SoundSpeed2
USE MOD_Flux_Average,ONLY:TwoPointEntropyConservingFlux
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: j,i
REAL                 :: F_c(5),uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
REAL                 :: sL,sR,rhoa,v1a,v2a,v3a,pa,aa,rhoea,Ha
REAL                 :: vJump(5),H(5,5),diss(5)
REAL                 :: lambdaMaxMax
REAL                 :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R
!==================================================================================================================================
LambdaMaxMax=0.
DO j=0,PP_N;  DO i=0,PP_N
  LambdaMaxMax= MAX(LambdaMaxMax, &
                    MAX(ABS(U_LL(2,i,j)/U_LL(1,i,j)),ABS(U_RR(2,i,j)/U_RR(1,i,j))) &
                    +SQRT(MAX(SoundSpeed2(U_LL(:,i,j)),SoundSpeed2(U_RR(:,i,j))))) 
END DO; END DO

DO j=0,PP_N
  DO i=0,PP_N
    CALL TwoPointEntropyConservingFlux(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                                       uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
    sRho_L=1./U_LL(1,i,j)
    sRho_R=1./U_RR(1,i,j)
    Vel_L =U_LL(2:4,i,j)*sRho_L
    Vel_R =U_RR(2:4,i,j)*sRho_R
    p_L   =KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R   =KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
! Compute necessary values for the entropy Jacobian H and the entropy variable
! jump
! Left
    sL    = LOG(p_L) - kappa*LOG(U_LL(1,i,j))

    vJump(1) = -((kappa-sL)/kappaM1 - 0.5*U_LL(1,i,j)*(SUM(Vel_L(:)**2))/p_L)
    vJump(2:4) = -U_LL(1,i,j)*vel_L(1:3)/p_L
    vJump(5) =  U_LL(1,i,j)/p_L

! Right
    sR    = LOG(p_R) - kappa*LOG(U_RR(1,i,j))

    vJump(1)  = vJump(1)   + ((kappa-sR)/kappaM1 - 0.5*U_RR(1,i,j)*(SUM(Vel_R(:)**2))/p_R)
    vJump(2:4)= vJump(2:4) + U_RR(1,i,j)*vel_R(1:3)/p_R
    vJump(5)  = vJump(5)   - U_RR(1,i,j)/p_R

! H, using local average mean values
    rhoa  = rhoHat
    v1a   = uHat
    v2a   = vHat
    v3a   = wHat
    pa    = p1Hat
    aa    = aHat
    rhoea = pa/kappaM1 + 0.5*rhoa*(v1a*v1a+v2a*v2a+v3a*v3a)
    Ha    = HHat
! create the entropy Jacobian matrix H
    H      = 0. ! initialize to zero
    H(1,1) = rhoa
    H(2,2) = rhoa*v1a*v1a+pa
    H(3,3) = rhoa*v2a*v2a+pa
    H(4,4) = rhoa*v3a*v3a+pa
    H(5,5) = rhoa*Ha*Ha - aa*aa*pa/kappaM1

    H(1,2) = rhoa*v1a
    H(2,1) = rhoa*v1a
    H(1,3) = rhoa*v2a
    H(3,1) = rhoa*v2a
    H(1,4) = rhoa*v3a
    H(4,1) = rhoa*v3a
    H(1,5) = rhoea
    H(5,1) = rhoea

    H(2,3) = rhoa*v1a*v2a
    H(3,2) = rhoa*v1a*v2a
    H(2,4) = rhoa*v1a*v3a
    H(4,2) = rhoa*v1a*v3a
    H(2,5) = rhoa*v1a*Ha
    H(5,2) = rhoa*v1a*Ha

    H(3,4) = rhoa*v2a*v3a
    H(4,3) = rhoa*v2a*v3a
    H(3,5) = rhoa*v2a*Ha
    H(5,3) = rhoa*v2a*Ha

    H(4,5) = rhoa*v3a*Ha
    H(5,4) = rhoa*v3a*Ha

    diss = H(:,1)*vJump(1) + H(:,2)*vJump(2) + H(:,3)*vJump(3) + H(:,4)*vJump(4) + H(:,5)*vJump(5)
    F(:,i,j) = F_c - 0.5*LambdaMaxMax*diss
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolver_EC_LLF


!==================================================================================================================================
!> Kennedy Gruber /Decros with LLF diss, EC+KEP - pressure with LLF diss
!==================================================================================================================================
SUBROUTINE RiemannSolver_VolumeFluxAverage_LLF(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:VolumeFluxAverage,SoundSpeed2
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: j,i
REAL                 :: F_c(5),uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
REAL                 :: lambdamaxmax
!==================================================================================================================================
LambdaMaxMax=0.
DO j=0,PP_N;  DO i=0,PP_N
  LambdaMaxMax= MAX(LambdaMaxMax, &
                    MAX(ABS(U_LL(2,i,j)/U_LL(1,i,j)),ABS(U_RR(2,i,j)/U_RR(1,i,j))) &
                    +SQRT(MAX(SoundSpeed2(U_LL(:,i,j)),SoundSpeed2(U_RR(:,i,j))))) 
END DO; END DO
DO j=0,PP_N
  DO i=0,PP_N
! Compute entropy conserving flux
    CALL VolumeFluxAverage(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                               uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
    F(:,i,j) = F_c - 0.5*LambdaMaxMax*(U_RR(:,i,j) - U_LL(:,i,j))
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolver_VolumeFluxAverage_LLF


!==================================================================================================================================
!> Kennedy Gruber /Decros with LLF diss
!==================================================================================================================================
SUBROUTINE RiemannSolver_ECKEP_LLF(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,SoundSpeed2
USE MOD_Flux_Average,ONLY:EntropyAndEnergyConservingFlux,LN_MEAN
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_LL  !< rotated conservative state left
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(IN)    :: U_RR  !< rotated conservative state right
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:5,0:PP_N,0:PP_N),INTENT(INOUT) :: F    !< numerical flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                                                               !
INTEGER              :: j,i
REAL                 :: F_c(5),uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat
REAL                 :: beta_L,beta_R,rhoa,betaLN,lambdamaxmax 
REAL                 :: sRho_L,sRho_R,Vel_L(3),Vel_R(3),p_L,p_R
!==================================================================================================================================
LambdaMaxMax=0.
DO j=0,PP_N;  DO i=0,PP_N
  LambdaMaxMax= MAX(LambdaMaxMax, &
                    MAX(ABS(U_LL(2,i,j)/U_LL(1,i,j)),ABS(U_RR(2,i,j)/U_RR(1,i,j))) &
                    +SQRT(MAX(SoundSpeed2(U_LL(:,i,j)),SoundSpeed2(U_RR(:,i,j))))) 
END DO; END DO

DO j=0,PP_N
  DO i=0,PP_N
    CALL EntropyAndEnergyConservingFlux(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                                        uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
    sRho_L=1./U_LL(1,i,j)
    sRho_R=1./U_RR(1,i,j)
    Vel_L =U_LL(2:4,i,j)*sRho_L
    Vel_R =U_RR(2:4,i,j)*sRho_R
    p_L   =KappaM1*(U_LL(5,i,j)-0.5*SUM(U_LL(2:4,i,j)*Vel_L(1:3)))
    p_R   =KappaM1*(U_RR(5,i,j)-0.5*SUM(U_RR(2:4,i,j)*Vel_R(1:3)))
    beta_L = 0.5*U_LL(1,i,j)/p_L
    beta_R = 0.5*U_RR(1,i,j)/p_R
! averages for the energy dissipation, already have rhoHat,uHat,vHat,wHat
    rhoa   = 0.5*(U_RR(1,i,j)+U_LL(1,i,j))
    betaLN = LN_MEAN(beta_L,beta_R)
    F(1:4,i,j) = F_c(1:4) - 0.5*LambdaMaxMax*(U_RR(1:4,i,j) - U_LL(1:4,i,j))
    F(5  ,i,j) = F_c(5)   - 0.5*LambdaMaxMax*(                                                                        &
                                    ( 0.5/(kappaM1*betaLN) +0.5*(SUM(vel_L(:)*vel_R(:))))*(U_RR(1,i,j) - U_LL(1,i,j)) &
                                  + rhoa*( uHat*(vel_R(1)-vel_L(1))                                                   &
                                          +vHat*(vel_R(2)-vel_L(2))                                                   &
                                          +wHat*(vel_R(3)-vel_L(3))                                                   &
                                          + 0.5/kappaM1*(1./beta_R - 1./beta_L)))
  END DO ! i
END DO ! j
END SUBROUTINE RiemannSolver_ECKEP_LLF


END MODULE MOD_Riemann
