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
!> Module for the shock capturing routines
!==================================================================================================================================

MODULE MOD_ShockCapturing
! MODULES
IMPLICIT NONE
PRIVATE
! ----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersShockCapturing
   MODULE PROCEDURE DefineParametersShockCapturing
END INTERFACE

INTERFACE InitShockCapturing
   MODULE PROCEDURE InitShockCapturing
END INTERFACE

#if SHOCKCAPTURE
INTERFACE CalcArtificialViscosity
   MODULE PROCEDURE CalcArtificialViscosity
END INTERFACE
#endif /*SHOCKCAPTURE*/

INTERFACE FinalizeShockCapturing
   MODULE PROCEDURE FinalizeShockCapturing
END INTERFACE

PUBLIC :: DefineParametersShockCapturing
PUBLIC :: InitShockCapturing
#if SHOCKCAPTURE
PUBLIC :: CalcArtificialViscosity
#endif /*SHOCKCAPTURE*/
PUBLIC :: FinalizeShockCapturing
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersShockCapturing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("ShockCapturing")
END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars         ,ONLY: nElems
USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================
IF (ShockCapturingInitIsDone.OR.(.NOT.InterpolationInitIsDone)) THEN
  SWRITE(*,*) "InitShockCapturing not ready to be called or already called."
  RETURN
END IF
IF (PP_N.LT.2) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF
! shock caturing parameters
ALLOCATE(nu(nElems))
nu     = 0.
nu_max = 0.

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING...'
CALL InitBasisTrans(PP_N,xGP)
#if (PP_Indicator_Var==1)
SWRITE(UNIT_StdOut,'(A)') '    USING DENSITY AS SHOCK INDICATOR!'
#elif (PP_Indicator_Var==2)
SWRITE(UNIT_StdOut,'(A)') '    USING PRESSURE AS SHOCK INDICATOR!'
#elif (PP_Indicator_Var==3)
SWRITE(UNIT_StdOut,'(A)') '    USING PRESSURE TIMES DENSITY AS SHOCK INDICATOR!'
#endif
ShockCapturingInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SHOCKCAPTURING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitShockCapturing


SUBROUTINE InitBasisTrans(N_in,xGP)
!===================================================================================================================================
!> Initialize Vandermodematrix for basis transformation
!===================================================================================================================================
! MODULES
USE MOD_ShockCapturing_Vars,ONLY:sVdm_Leg
USE MOD_Basis, ONLY :BuildLegendreVdm
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP
REAL,DIMENSION(0:N_in,0:N_in)              :: Vdm_Leg
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!  NODAL <--> MODAL
! Compute the 1D Vandermondematrix, needed to tranform the nodal basis into a modal (Legendre) basis
ALLOCATE(sVdm_Leg(0:N_in,0:N_in))
CALL BuildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
END SUBROUTINE InitBasisTrans

#if SHOCKCAPTURE
SUBROUTINE CalcArtificialViscosity(U)
!===================================================================================================================================
!> Use framework of Persson and Peraire to measure shocks with DOF energy indicator and calculate artificial viscosity, if necessary
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ShockCapturing_Vars, ONLY: sVdm_Leg,nu,nu_max
USE MOD_Equation_Vars      , ONLY: KappaM1,kappa
#ifdef mhd
USE MOD_Equation_Vars      , ONLY: s2mu_0
#endif /*mhd*/
USE MOD_TimeDisc_Vars      , ONLY: dt
USE MOD_Mesh_Vars          , ONLY: nElems,sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars      , ONLY: ConsToPrim
USE MOD_Equation_Vars      , ONLY: FastestWave3D
USE MOD_ChangeBasis        , ONLY: ChangeBasis3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems),INTENT(IN) :: U
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
REAL,DIMENSION(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
REAL                                     :: LU,LUM1,LUM2,LU_N,LU_NM1,eta_dof,eta_min,eta_max,eps0,p
REAL                                     :: v(3),Prim(1:PP_nVar),cf,Max_Lambda(6),lambda_max,h,lambda_max2
INTEGER                                  :: l,ind,i,j,k

nu    =0.
nu_max=0.
DO l=1,nElems

  ! Select a shock indicator: density = 1, pressure = 2
  SELECT CASE(PP_Indicator_Var)
    CASE(1)
    Uind(1,:,:,:) = U(1,:,:,:,l)
    CASE(2)
    Uind(1,:,:,:) = KappaM1*(U(5,:,:,:,l)-0.5*(SUM(U(2:4,:,:,:,l)*U(2:4,:,:,:,l))/U(1,:,:,:,l)))
#ifdef mhd
    Uind(1,:,:,:) = Uind(1,:,:,:)-KappaM1*s2mu_0*SUM(U(6:8,:,:,:,l)*U(6:8,:,:,:,l))
#ifdef PP_GLM
    Uind(1,:,:,:) = Uind(1,:,:,:)-0.5*KappaM1*U(9,:,:,:,l)*U(9,:,:,:,l)
#endif /*PP_GLM*/
#endif /*mhd*/
    CASE(3)
    Uind(1,:,:,:) = KappaM1*(U(5,:,:,:,l)-0.5*(SUM(U(2:4,:,:,:,l)*U(2:4,:,:,:,l))/U(1,:,:,:,l)))
#ifdef mhd
    Uind(1,:,:,:) = Uind(1,:,:,:)-KappaM1*s2mu_0*SUM(U(6:8,:,:,:,l)*U(6:8,:,:,:,l))
#ifdef PP_GLM
    Uind(1,:,:,:) = Uind(1,:,:,:)-0.5*KappaM1*U(9,:,:,:,l)*U(9,:,:,:,l)
#endif /*PP_GLM*/
#endif /*mhd*/
    Uind(1,:,:,:) = Uind(1,:,:,:)*U(1,:,:,:,l)
  END SELECT
  
  ! Transform Uind into modal Legendre interpolant Umod
  CALL ChangeBasis3D(1,PP_N,PP_N,sVdm_Leg,Uind,Umod)

  ! Compute (truncated) error norms
  LU     = SUM(Umod(1,:,:,:)**2)
  LUM1   = SUM(Umod(1,0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
  LUM2   = SUM(Umod(1,0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
  LU_N   = LU-LUM1
  LU_NM1 = LUM1-LUM2

  ! DOF energy indicator
  eta_dof = LOG10(MAX(LU_N/LU,LU_NM1/LUM1))

  ! Artificial Viscosity
  eta_min = -5.0
  eta_max = -3.0
  eps0 = 0.05
  IF (eta_dof.GE.eta_max) THEN
    nu(l) = eps0
  ELSE IF (eta_dof.LE.eta_min) THEN
    nu(l) = 0.
  ELSE
    nu(l) = 0.5*eps0*(1.0+SIN(PP_Pi*(eta_dof-0.5*(eta_max+eta_min))/(eta_max-eta_min)))
  END IF
   
  ! Save max artificial viscosity for DFL timestepping
  nu_max = MAX(nu_max,nu(l))

  ! Get (transformed) max eigenvalue:
  Max_Lambda = 0.0
  DO i=1,PP_N
    DO j=1,PP_N
      DO k=1,PP_N
        CALL ConsToPrim(Prim,U(:,i,j,k,l))
        CALL FastestWave3D(Prim,cf)
        v(:)=Prim(2:4) 
        Max_Lambda(1)=MAX(Max_Lambda(1),sJ(i,j,k,l)*(ABS(SUM(Metrics_fTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_fTilde(:,i,j,k,l)*Metrics_fTilde(:,i,j,k,l)))))
        Max_Lambda(2)=MAX(Max_Lambda(2),sJ(i,j,k,l)*(ABS(SUM(Metrics_gTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_gTilde(:,i,j,k,l)*Metrics_gTilde(:,i,j,k,l)))))
        Max_Lambda(3)=MAX(Max_Lambda(3),sJ(i,j,k,l)*(ABS(SUM(Metrics_hTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_hTilde(:,i,j,k,l)*Metrics_hTilde(:,i,j,k,l)))))
        Max_Lambda(4)=MAX(Max_Lambda(4),ABS(v(1)) + cf)
        Max_Lambda(5)=MAX(Max_Lambda(5),ABS(v(2)) + cf)
        Max_Lambda(6)=MAX(Max_Lambda(6),ABS(v(3)) + cf)
      END DO
    END DO
  END DO

  v(1) = MAXVAL(Metrics_fTilde(:,:,:,:,l))
  v(2) = MAXVAL(Metrics_gTilde(:,:,:,:,l))
  v(3) = MAXVAL(Metrics_hTilde(:,:,:,:,l))
  eps0 = 1.0/SQRT(MAXVAL(v))

!  h=2.0*eps0/MINVAL(sJ(:,:,:,l))
!  lambda_max=MAXVAL(Max_Lambda(1:3))*h
  h=(8.0/MINVAL(sJ(:,:,:,l)))**(1.0/3.0)
  lambda_max2=MAXVAL(Max_Lambda(4:6))*h

  ! Scaling of artificial viscosity
  nu(l) = nu(l)*lambda_max2/(REAL(PP_N))

END DO ! l

END SUBROUTINE CalcArtificialViscosity
#endif /*SHOCKCAPTURE*/

SUBROUTINE FinalizeShockCapturing()
!============================================================================================================================
!> Deallocate all global shock capturing variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ShockCapturing_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================
IF (.NOT.ShockCapturingInitIsDone) THEN
  WRITE(UNIT_stdOut,*) "InitShockCapturing was not called before."
  RETURN
END IF
ShockCapturingInitIsDone = .FALSE.
nu_max = 0.
SDEALLOCATE(sVdm_Leg)
SDEALLOCATE(nu)
END SUBROUTINE FinalizeShockCapturing

END MODULE MOD_ShockCapturing
