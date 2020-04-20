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

INTERFACE RiemannSolverByRusanov
  MODULE PROCEDURE RiemannSolverByRusanov
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

INTERFACE RiemannSolver_ESM
  MODULE PROCEDURE RiemannSolver_ESM
END INTERFACE

INTERFACE RiemannSolver_VolumeFluxAverage_LLF
  MODULE PROCEDURE RiemannSolver_VolumeFluxAverage_LLF
END INTERFACE

INTERFACE RiemannSolver_ECKEP_LLF
  MODULE PROCEDURE RiemannSolver_ECKEP_LLF
END INTERFACE

PUBLIC:: Riemann
PUBLIC:: RiemannSolverByRusanov
PUBLIC:: RiemannSolverByHLLC
PUBLIC:: RiemannSolverByRoe
PUBLIC:: RiemannSolver_EntropyStable
PUBLIC:: RiemannSolver_VolumeFluxAverage
PUBLIC:: RiemannSolver_EC_LLF
PUBLIC:: RiemannSolver_ESM
PUBLIC:: RiemannSolver_VolumeFluxAverage_LLF
PUBLIC:: RiemannSolver_ECKEP_LLF
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
#endif
                   nv,t1,t2)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars   ,ONLY:SolveRiemannProblem
USE MOD_Flux_Average
#if PARABOLIC
USE MOD_Flux            ,ONLY:EvalDiffFlux3D    ! and the NSE diffusion fluxes in all directions to approximate the numerical flux
#endif
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
!==================================================================================================================================
! Momentum has to be rotatet using the normal system individual for each
! Gauss point i,j
DO j=0,PP_N
  DO i=0,PP_N
    !LEFT
    U_LL(1,i,j)=U_L(1,i,j)
    ! rotate momentum
    U_LL(2,i,j)=SUM(U_L(2:4,i,j)*nv(:,i,j))
    U_LL(3,i,j)=SUM(U_L(2:4,i,j)*t1(:,i,j))
    U_LL(4,i,j)=SUM(U_L(2:4,i,j)*t2(:,i,j))
    U_LL(5,i,j)=U_L(5,i,j)
    !right
    U_RR(1,i,j)=U_R(1,i,j)
    ! rotate momentum
    U_RR(2,i,j)=SUM(U_R(2:4,i,j)*nv(:,i,j))
    U_RR(3,i,j)=SUM(U_R(2:4,i,j)*t1(:,i,j))
    U_RR(4,i,j)=SUM(U_R(2:4,i,j)*t2(:,i,j))
    U_RR(5,i,j)=U_R(5,i,j)
  END DO ! i 
END DO ! j


CALL SolveRiemannProblem(F,U_LL,U_RR)


! Back Rotate the normal flux into Cartesian direction
DO j=0,PP_N
  DO i=0,PP_N
    F(2:4,i,j)= nv(:,i,j)*F(2,i,j) &
               +t1(:,i,j)*F(3,i,j) &
               +t2(:,i,j)*F(4,i,j)
  END DO ! i
END DO ! j
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
!> Entropy stable Riemann solver for Mortar Sides, uses TwoPoint Entropy Conserving flux. Ratio 2-to-1 is assumed
!==================================================================================================================================

SUBROUTINE RiemannSolver_ESM(Uface_master,Uface_slave,Flux_master,Flux_slave,doMPISides, weak)
  USE MOD_Preproc
  USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2
  USE MOD_Mortar_Vars, ONLY: M_1_0,M_2_0
  USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo, FS2M,nElems,nMortarSides
  USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
  USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
  USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,lastSlaveSide
  USE MOD_Mesh_Vars,   ONLY: NormVec,TangVec1,TangVec2,FS2M,nSides, SurfElem
  USE MOD_Flux_Average ,ONLY: TwoPointEntropyConservingFlux
  USE MOD_Interpolation_Vars  ,ONLY: wGP
  USE MOD_Analyze_Vars  ,ONLY: wGPSurf
  USE MOD_Equation_Vars  ,ONLY: ConsToPrim, ConsToEntropy, kappa
  USE MOD_Equation_Vars ,ONLY:kappaM1
  
  IMPLICIT NONE
  !----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT/OUTPUT VARIABLES
  REAL,INTENT(INOUT) :: Uface_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
  REAL,INTENT(INOUT) :: Uface_slave( 1:PP_nVar,0:PP_N,0:PP_N,FirstSlaveSide:LastSlaveSide) !< (INOUT) can be U or Grad_Ux/y/z_master
  REAL,INTENT(INOUT)   :: Flux_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< on input: has flux from small mortar sides 
                                                                      !< on output: flux on big mortar sides filled
  REAL,INTENT(INOUT   )   :: Flux_slave(1:PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide) !<has flux from small mortar sides,
                                                                      !< set -F_slave in call if surfint is weak (dg.f90)
                                                                      !< set +F_slave in call if surfint is strong (lifting)
  LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
  LOGICAL,INTENT(IN) :: weak                                          !< flag whether strong or weak form is used

  !----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  INTEGER      :: p,q,l, m, s,t, i,j, iSide
INTEGER      :: iMortar,nMortars,iVar,iVar1
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID(4),locSide,flip(4)
REAL         :: Flux_l( PP_nVar,0:PP_N,0:PP_N,0:1,0:1) ! For small mortar sides
REAL         :: Flux_R( PP_nVar,0:PP_N,0:PP_N) !Big mortar Side
REAL         :: U_LL( PP_nVar,0:PP_N,0:PP_N,0:1, 0:1) !For small mortar sides
REAL         :: U_RR( PP_nVar,0:PP_N,0:PP_N) !Big mortar Side
REAL         :: U_L( PP_nVar,0:PP_N,0:PP_N,0:1, 0:1) !For small mortar sides
REAL         :: U_R( PP_nVar,0:PP_N,0:PP_N) !Big mortar Side
REAL         :: SUM2, SUM3, sRho, pres,v1,v2,v3, SumA, SumB, uint_r( PP_nVar), uint_l( PP_nVar)
REAL         :: ent_loc, FluxF(PP_nVar), prim(PP_nVar)

REAL         :: F_c( PP_nVar) !Returned Value from TwoPointFlux
! REAL         :: U_tmp( PP_nVar,0:PP_N,0:PP_N,1:4)
! REAL         :: U_tmp2(PP_nVar,0:PP_N,0:PP_N,1:2)
REAL,POINTER :: M1(:,:),M2(:,:)
REAL         ::PR2L(0:1,0:PP_N,0:PP_N),PL2R(0:1,0:PP_N,0:PP_N)
REAL                  :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
REAL         :: nv(            3,0:PP_N,0:PP_N) !< normal vector of face
REAL         :: t1(            3,0:PP_N,0:PP_N) !< 1st tangential vector of face
REAL         :: t2(            3,0:PP_N,0:PP_N) !< 2nd tangential vector of face
REAL         :: SUM1, pr2l_upper(PP_N+1,PP_N+1), pr2l_lower(PP_N+1,PP_N+1) , pl2r_upper(PP_N+1,PP_N+1) , pl2r_lower(PP_N+1,PP_N+1)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
  lastMortarSideID =  lastMortarMPISide
ELSE
  firstMortarSideID = firstMortarInnerSide
  lastMortarSideID =  lastMortarInnerSide
END IF !doMPISides

! First compute F^{L_st}
PR2L(0,:,:)=M_0_1(:,:); PR2L(1,:,:)=M_0_2(:,:);

PL2R(0,:,:)=M_1_0(:,:); Pl2R(1,:,:)=M_2_0(:,:);


!
DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  locSide=MortarType(2,MortarSideID)
  
! For debugging, store MortarSideID
  iVar1=MortarSideID

  DO iMortar=1,nMortars
    SideID(iMortar)= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip(iMortar)  = MortarInfo(MI_FLIP,iMortar,locSide)
  ENDDO   

  
  DO q=0,PP_N; DO p=0,PP_N
    U_L(:,p,q,0,0) = Uface_slave(:,FS2M(1,p,q,flip(1)),FS2M(2,p,q,flip(1)),SideID(0+1))
    U_L(:,p,q,1,0) = Uface_slave(:,FS2M(1,p,q,flip(2)),FS2M(2,p,q,flip(2)),SideID(0+2))
    U_L(:,p,q,0,1) = Uface_slave(:,FS2M(1,p,q,flip(3)),FS2M(2,p,q,flip(3)),SideID(0+3))
    U_L(:,p,q,1,1) = Uface_slave(:,FS2M(1,p,q,flip(4)),FS2M(2,p,q,flip(4)),SideID(0+4))

    ! U_L(:,p,q,0,0) = Uface_slave(:,p,q,SideID(0+1))
    ! U_L(:,p,q,1,0) = Uface_slave(:,p,q,SideID(0+2))
    ! U_L(:,p,q,0,1) = Uface_slave(:,p,q,SideID(0+3))
    ! U_L(:,p,q,1,1) = Uface_slave(:,p,q,SideID(0+4))
  END DO; END DO ! q, p
  ! DO iMortar=1,nMortars
    ! SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    ! flip   = MortarInfo(MI_FLIP,iMortar,iSide)
  
  DO S = 0, 1; DO T = 0,1;
    ! PRINT *, "Riemann SideID(S + T*2 + 1)", SideID(S + T*2 + 1), S + T*2 + 1
    
    nv = NormVec(:,:,:,SideID(S + T*2 + 1));
    t1 = TangVec1(:,:,:,SideID(S + T*2 + 1));
    t2 = TangVec2(:,:,:,SideID(S + T*2 + 1));
    
    
   DO j=0,PP_N
    DO i=0,PP_N
      !LEFT
      U_LL(1,i,j,S,T)=U_L(1,i,j,S,T)
      ! rotate momentum
      U_LL(2,i,j,S,T)=SUM(U_L(2:4,i,j,S,T)*nv(:,i,j))
      U_LL(3,i,j,S,T)=SUM(U_L(2:4,i,j,S,T)*t1(:,i,j))
      U_LL(4,i,j,S,T)=SUM(U_L(2:4,i,j,S,T)*t2(:,i,j))
      U_LL(5,i,j,S,T)=U_L(5,i,j,S,T)
    ENDDO
  ENDDO
ENDDO; ENDDO;


U_R(:,:,:) = Uface_master(:,:,:,MortarSideID)
nv =  NormVec(:,:,:,  MortarSideID);
t1 = TangVec1(:,:,:, MortarSideID);
t2 = TangVec2(:,:,:, MortarSideID);
  DO j=0,PP_N
    DO i=0,PP_N
  
      !right
      U_RR(1,i,j)=U_R(1,i,j)
      ! rotate momentum
      U_RR(2,i,j)=SUM(U_R(2:4,i,j)*nv(:,i,j))
      U_RR(3,i,j)=SUM(U_R(2:4,i,j)*t1(:,i,j))
      U_RR(4,i,j)=SUM(U_R(2:4,i,j)*t2(:,i,j))
      U_RR(5,i,j)=U_R(5,i,j)
    END DO ! i 
  END DO ! j
  ! PRINT *, "Riemann 1"
  Flux_l = 0;
    DO S = 0, 1; DO T = 0,1;
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     nv = NormVec(:,:,:,SideID(S + T*2 + 1));
     t1 = TangVec1(:,:,:,SideID(S + T*2 + 1));
     t2 = TangVec2(:,:,:,SideID(S + T*2 + 1));
      print *, "s = ", s, "     t = ", t
      print *, "nv = ", nv(:, 1, 1)
      print *, "t1 = ", t1(:, 1, 1)
      print *, "t2 = ", t2(:, 1, 1)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  end do; end do

  DO i = 0,PP_N; DO j = 0,PP_N;
    DO S = 0, 1; DO T = 0,1;
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     nv = NormVec(:,:,:,SideID(S + T*2 + 1));
     t1 = TangVec1(:,:,:,SideID(S + T*2 + 1));
     t2 = TangVec2(:,:,:,SideID(S + T*2 + 1));
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      DO l = 0,PP_N; DO m = 0,PP_N;
        CALL TwoPointEntropyConservingFlux(F_c,U_LL(:,i,j,s,t),U_RR(:,l,m),&
        uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        F_c(2:4)= nv(:,i,j)*F_c(2) &
        +t1(:,i,j)*F_c(3) &
        +t2(:,i,j)*F_c(4)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        ! U_LL is U_slave (S,T)
        ! U_RR is U_mortar
       Flux_l(:,i,j,S,T) = Flux_l(:,i,j,S,T) + PR2L(s,L,I) * PR2L(t,M,J) * F_c(:); !1(a)
        
      END DO; END DO !l,m
    END DO; END DO !S,T
  END DO; END DO !i ,j 
  
  
  ! Back Rotate the normal flux into Cartesian direction
  !DO S = 0, 1; DO T = 0,1;
  !  nv = NormVec(:,:,:,SideID(S + T*2 + 1));
  !  t1 = TangVec1(:,:,:,SideID(S + T*2 + 1));
  !  t2 = TangVec2(:,:,:,SideID(S + T*2 + 1));
  !  !print *, "S,T=",S,T
  !  !print *, "nv=",nv
  !  !print *, "t1=",t1
  !  !print *, "t2=",t2
  !  DO j=0,PP_N
  !    DO i=0,PP_N
  !      Flux_l(2:4,i,j,S,T)= nv(:,i,j)*Flux_l(2,i,j,S,T) &
  !      +t1(:,i,j)*Flux_l(3,i,j,S,T) &
  !      +t2(:,i,j)*Flux_l(4,i,j,S,T)
  !      !Flux_l(2,i,j,S,T)= -Flux_l(2,i,j,S,T)
  !      !Flux_l(3,i,j,S,T)=  Flux_l(3,i,j,S,T)
  !      !Flux_l(4,i,j,S,T)= -Flux_l(4,i,j,S,T)
  !    END DO ! i
  !  END DO ! j
  !END DO; END DO !S,T
  
  DO S = 0, 1; DO T = 0,1;
    DO q=0,PP_N
      DO p=0,PP_N
        
        Flux_master(:,p,q,SideID((S + T*2 + 1))) = Flux_l(:,p,q,S,T)
        Flux_slave(:,p,q,SideID((S + T*2 + 1))) = Flux_l(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),S,T)
        ! Flux_slave(:,p,q,SideID((S + T*2 + 1))) = Flux_l(:,p,q,S,T)
        ! Flux_slave(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),SideID((S + T*2 + 1))) = Flux_l(:,p,q,S,T)

        ! Flux_slave(:,:,:,SideID(0+2)) = Flux_l(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),1,0) 
        ! Flux_slave(:,:,:,SideID(0+3)) = Flux_l(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),0,1) 
        ! Flux_slave(:,:,:,SideID(0+4)) = Flux_l(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),1,1) 
        
      ENDDO; ENDDO;
      
    END DO; END DO !S,T
    ! Flux_master(:,:,:,MortarSideID) = Flux_slave(:,:,:,SideID(1))
  ! PRINT *, " Flux_slave(:,1,1,SideID(1))",  Flux_slave(:,1,1,SideID(1))!, Flux_slave(:,1,1,SideID(1))
! PRINT *, " Flux_slave(:,1,1,SideID(1))",  Flux_slave(:,1,1,SideID(1))
DO i = 1,4
  DO q=0,PP_N
    DO p=0,PP_N
      IF(weak)THEN
        Flux_slave(:,p,q,SideID(i))= Flux_slave(:,p,q,SideID(i))*SurfElem(p,q,SideID(i))
      else
        Flux_slave(:,p,q,SideID(i))= -Flux_slave(:,p,q,SideID(i))*SurfElem(p,q,SideID(i))
      ENDIF
    END DO
  END DO
  
END DO
! PRINT *, "NACHDEM Flux_slave(:,1,1,SideID(1))",  Flux_slave(:,1,1,SideID(1))!, Flux_slave(:,1,1,SideID(1))

! PRINT *, "Riemann 2"
! NOw BIG Master Side
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  nv =  NormVec(:,:,:,  MortarSideID);
  t1 = TangVec1(:,:,:, MortarSideID);
  t2 = TangVec2(:,:,:, MortarSideID);
  print *, "nv = ", nv(:, 1, 1)
  print *, "t1 = ", t1(:, 1, 1)
  print *, "t2 = ", t2(:, 1, 1)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Flux_R(:,:,:) = 0;

DO i = 0,PP_N; DO j = 0,PP_N;
  ! s = 0; t = 0;
  DO S = 0, 1; DO T = 0,1;
    
    DO l = 0,PP_N; DO m = 0,PP_N;
      CALL TwoPointEntropyConservingFlux(F_c,U_LL(:,l,m,s,t),U_RR(:,i,j),&
      uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        F_c(2:4)= nv(:,i,j)*F_c(2) &
        +t1(:,i,j)*F_c(3) &
        +t2(:,i,j)*F_c(4)
  !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! U_LL is U_slave (S,T)
      ! U_RR is U_mortar
      Flux_R(:,i,j) =  Flux_R(:,i,j) + PL2R(s,L,I) * PL2R(t,M,J) * F_c(:);
     
    END DO; END DO ! l, m 
  END DO; END DO !S,T,
END DO; END DO !i,j


!nv =  NormVec(:,:,:,  MortarSideID);
!t1 = TangVec1(:,:,:, MortarSideID);
!t2 = TangVec2(:,:,:, MortarSideID);
!!print *, "nv=",nv
!!print *, "t1=",t1
!!print *, "t2=",t2
!DO j=0,PP_N
!  DO i=0,PP_N
!    Flux_R(2:4,i,j)= nv(:,i,j)*Flux_R(2,i,j) &
!    +t1(:,i,j)*Flux_R(3,i,j) &
!    +t2(:,i,j)*Flux_R(4,i,j)
!    !Flux_R(2,i,j)= -Flux_R(2,i,j)
!    !Flux_R(3,i,j)=  Flux_R(3,i,j)
!    !Flux_R(4,i,j)= -Flux_R(4,i,j)
!  END DO ! i
!END DO ! j
!cALL EXIT()

Flux_master(:,:,:,MortarSideID) = 0.25* Flux_R(:,:,:)

! Flux_slave(:,p,q,SideID((S + T*2 + 1))) = Flux_l(:,FS2M(1,p,q,flip(S + T*2 + 1)),FS2M(2,p,q,flip(S + T*2 + 1)),S,T)

! PRINT *, "Flux_master(:,1,1,MortarSideID)", Flux_master(:,1,1,MortarSideID)
DO q=0,PP_N
  DO p=0,PP_N
    ! Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID)*SurfElem(p,q,MortarSideID)
    Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID)*SurfElem(p,q,MortarSideID)
  END DO
END DO
    ! PRINT *, "Riemann 3"
  ! PRINT *, "Flux_master(:,1,1,MortarSideID)", Flux_master(:,1,1,MortarSideID)

! PRINT *, "NACHDEM Flux_master(:,1,1,MortarSideID)", Flux_master(:,1,1,MortarSideID)
!  CALL EXIT()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Remove after tests
  print *, "nElems=",nElems
  print *, "nMortarSides=",nMortarSides

  Flux_R = Flux_R/4.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PRIMARY CONSERVATION - SECOND TERM
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  uint_l = 0.
  S=0; T=0;
  DO i = 0,PP_N; DO j = 0,PP_N;
    DO S = 0, 1; DO T = 0,1;
        DO iVar = 1, PP_nVar
          uint_l(iVar) = uint_l(iVar) + wGPSurf(i,j) * Flux_L(iVar, i,j,S,T)
        ENDDO
      END DO; END DO !S,T,
  END DO; END DO !i,j
  uint_l = uint_l /4.
  PRINT *, "IU_t second part = ", uint_l

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PRIMARY CONSERVATION - FIRST TERM
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  uint_r = 0.
  DO i = 0,PP_N; DO j = 0,PP_N;
        DO iVar = 1, PP_nVar
          uint_r(iVar) = uint_r(iVar) + wGPSurf(i,j) * Flux_R(iVar, i,j)
        ENDDO
    SUM1 = SUM1 + wGPSurf(i,j)*(SUM2-sum3)
  END DO; END DO !i,j
  PRINT *, "IU_t  first part = ", uint_r
  PRINT *, "Diff = ", uint_r - uint_l
  PRINT *
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUM1 = 0
  SUM2 = 0.
  S=0; T=0;
  DO i = 0,PP_N; DO j = 0,PP_N;
    ! s = 0; t = 0;
    SUM3=0.
    SUM2=0.
    DO S = 0, 1; DO T = 0,1;
        CALL ConsToPrim(prim,U_LL(:,i,j,S,T))
        ent_loc  = -prim(1)*(LOG(prim(5))-kappa*LOG(prim(1)))/kappaM1
        F_c = ConsToEntropy(U_LL(:,i,j,S,T)) ! '\NU^T
        ASSOCIATE(rho   =>U_LL(1,i,j,S,T), &
          rhov1 =>U_LL(2,i,j,S,T), &
          rhov2 =>U_LL(3,i,j,S,T), &
          rhov3 =>U_LL(4,i,j,S,T), &
          rhoE  =>U_LL(5,i,j,S,T)   )
          
          srho = 1./rho
          v1   = rhov1*srho 
          v2   = rhov2*srho 
          v3   = rhov3*srho 
          pres    = kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))

       
        FluxF(1) = rhov1
        FluxF(2) = rhov1*v1+pres
        FluxF(3) = rhov1*v2    
        FluxF(4) = rhov1*v3    
        FluxF(5) = (rhoE+pres)*v1
        END ASSOCIATE !v_x/y/z...
       
        DO iVar = 1, PP_nVar
          SUM2 = SUM2 + F_C(iVar) * Flux_L(iVar, i,j,S,T)
        ENDDO
      
        DO iVar = 1, PP_nVar
          SUM3 = SUM3 + F_C(iVar) * FluxF(iVar)
        ENDDO
        SUM3 = SUM3 - ent_loc * U_LL(2,i,j,S,T)/U_LL(1,i,j,S,T)
      
      END DO; END DO !S,T,
      SUM1 = SUM1 + wGPSurf(i,j)*(SUM2-sum3)
  END DO; END DO !i,j
  SumB = Sum1 /4.
  PRINT *, "iVar1 = ", iVar1, "= IS_t second part = ", SumB

  
  SUM1 = 0.
  SUM2 = 0.
  DO i = 0,PP_N; DO j = 0,PP_N;
        CALL ConsToPrim(prim,U_RR(:,i,j))
        ent_loc  = -prim(1)*(LOG(prim(5))-kappa*LOG(prim(1)))/kappaM1
        F_c = ConsToEntropy(U_RR(:,i,j)) ! '\NU^T
        ASSOCIATE(rho   =>U_RR(1,i,j), &
          rhov1 =>U_RR(2,i,j), &
          rhov2 =>U_RR(3,i,j), &
          rhov3 =>U_RR(4,i,j), &
          rhoE  =>U_RR(5,i,j)   )
          
          srho = 1./rho
          v1   = rhov1*srho 
          v2   = rhov2*srho 
          v3   = rhov3*srho 
          pres    = kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))

       
        FluxF(1) = rhov1
        FluxF(2) = rhov1*v1+pres
        FluxF(3) = rhov1*v2    
        FluxF(4) = rhov1*v3    
        FluxF(5) = (rhoE+pres)*v1
        END ASSOCIATE !v_x/y/z...

        SUM2=0.
        DO iVar = 1, PP_nVar
          SUM2 = SUM2 + F_C(iVar) * Flux_R(iVar, i,j)
        ENDDO

        SUM3=0.
        DO iVar = 1, PP_nVar
          SUM3 = SUM3 + F_C(iVar) * FluxF(iVar)
        ENDDO
        SUM3 = SUM3 - ent_loc * U_RR(2,i,j)/U_RR(1,i,j)
      
    SUM1 = SUM1 + wGPSurf(i,j)*(SUM2-sum3)
  END DO; END DO !i,j
  SumA = Sum1
  PRINT *, "iVar1 = ", iVar1, "= IS_t  first part = ", SumA
  PRINT *, "iVar1 = ", iVar1, "= Diff = ", SumA - SumB

  PRINT *, "TEST FINISH"
  CALL EXIT()
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if(iVar1 .eq. 32) then
if(iVar1 .eq. 24) then
  call exit()
end if
ENDDO !DO MortarSideID=firstMortarSideID,lastMortarSideID
  ! PRINT *, "Riamann ENd"

END SUBROUTINE RiemannSolver_ESM

!==================================================================================================================================
!> Entropy stable Riemann solver, uses TwoPoint Entropy Conserving flux
!==================================================================================================================================
SUBROUTINE RiemannSolver_EntropyStable(F,U_LL,U_RR)
!MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY: kappa,KappaM1
USE MOD_Flux_Average ,ONLY: TwoPointEntropyConservingFlux
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
! Compute entropy conserving flux
    CALL TwoPointEntropyConservingFlux(F_c,U_LL(:,i,j),U_RR(:,i,j),&
                                       uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
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
! Matrix of right eigenvectors
    RHat(1,:) = (/ 1.               , 1.                                      , 0.   ,  0.  , 1.               /)
    RHat(2,:) = (/ uHat - aHat      , uHat                                    , 0.   ,  0.  , uHat + aHat      /)
    RHat(3,:) = (/ vHat             , vHat                                    , 1.   ,  0.  , vHat             /)
    RHat(4,:) = (/ wHat             , wHat                                    , 0.   ,  1.  , wHat             /)
    RHat(5,:) = (/ HHat - uHat*aHat , 0.5*(uHat*uHat + vHat*vHat + wHat*wHat) , vHat , wHat , HHat + uHat*aHat /)
! Diagonal scaling matrix where DHat = ABS(\Lambda)S
    DHat(1) = 0.5*ABS(uHat - aHat)*rhoHat/kappa
    DHat(2) = ABS(uHat)*(kappaM1/kappa)*rhoHat!*rhoHat*rhoHat
    DHat(3) = ABS(uHat)*p1Hat!*rhoHat
    DHat(4) = DHat(3)
    DHat(5) = 0.5*ABS(uHat + aHat)*rhoHat/kappa
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
