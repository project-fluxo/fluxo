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

INTERFACE CalcArtificialViscosity
   MODULE PROCEDURE CalcArtificialViscosity
END INTERFACE

INTERFACE FinalizeShockCapturing
   MODULE PROCEDURE FinalizeShockCapturing
END INTERFACE

PUBLIC::DefineParametersShockCapturing,InitShockCapturing,CalcArtificialViscosity,FinalizeShockCapturing
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
CALL prms%CreateLogicalOption("CaptureShocks", "Shock Capturing on/off")
END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
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
  CALL abort(__STAMP__,'InitShockCapturing not ready to be called or already called.',999,999.)
  RETURN
END IF
CaptureShocks=GETLOGICAL('CaptureShocks','.FALSE.') 
IF (PP_N.LT.2.AND.CaptureShocks) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING...'
CALL InitBasisTrans(PP_N,xGP)
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
USE MOD_Basis, ONLY :buildLegendreVdm
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
CALL buildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
END SUBROUTINE InitBasisTrans


SUBROUTINE CalcArtificialViscosity(U)
!===================================================================================================================================
!> Use framework of Persson and Peraire to measure shocks with DOF energy indicator and calculate artificial viscosity, if necessary
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ShockCapturing_Vars, ONLY: CaptureShocks,sVdm_Leg
USE MOD_DG_Vars            , ONLY: nu,nu_max!,lambda_max
USE MOD_Equation_Vars      , ONLY: KappaM1,s2mu_0,kappa
USE MOD_TimeDisc_Vars      , ONLY: dt
USE MOD_Mesh_Vars	   , ONLY: PP_nElems=>nElems,sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars	   , ONLY: ConsToPrim
USE MOD_Equation_Vars	   , ONLY: FastestWave3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems),INTENT(IN)  :: U
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
REAL                                 :: LU,LUM1,LUM2,LU_N,LU_NM1,eta_dof,eta_min,eta_max,eps0!,c2,ca2,va2,astar,cf,lambda
REAL                                 :: v(3),Prim(1:PP_nVar),cf,Max_Lambda(6),lambda_max,h,lambda_max2
INTEGER                              :: l,ind,i,j,k

! Check if shock capturing is turned on
IF(.NOT.CaptureShocks) RETURN

nu=0.
nu_max=0.
DO l=1,PP_nElems

  ! Choose indicator (density: ind=1, pressure: ind=2)
  ind = 2
  SELECT CASE(ind)
    CASE(1)
    Uind = U(1,:,:,:,l)
    CASE(2)
    Uind = KappaM1*(U(5,:,:,:,l)-0.5*(SUM(U(2:4,:,:,:,l)*U(2:4,:,:,:,l))/U(1,:,:,:,l))-s2mu_0*SUM(U(6:9,:,:,:,l)*U(6:9,:,:,:,l)))
  END SELECT
  
  ! Transform Uind into modal Legendre interpolant Umod
  CALL ChangeBasis3D(PP_N,PP_N,sVdm_Leg,Uind,Umod)

  ! Compute (truncated) error norms
  LU     = SUM(Umod(:,:,:)**2)
  LUM1   = SUM(Umod(0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
  LUM2   = SUM(Umod(0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
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

  h=2.0*eps0/MINVAL(sJ(:,:,:,l))
  lambda_max=MAXVAL(Max_Lambda(1:3))*h
  h=(8.0/MINVAL(sJ(:,:,:,l)))**(1.0/3.0)
  lambda_max2=MAXVAL(Max_Lambda(4:6))*h

  ! Scaling of artificial viscosity
  nu(l) = nu(l)*lambda_max2/(REAL(PP_N))

END DO ! l

END SUBROUTINE CalcArtificialViscosity


SUBROUTINE ChangeBasis3D(N_In,N_Out,Vdm,X3D_In,X3D_Out)
!===================================================================================================================================
!> interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
!> to another 3D tensor product node positions (number of nodes N_out+1)
!> defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!>  xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: N_In,N_Out
REAL,INTENT(IN)     :: X3D_In(0:N_In,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(0:N_Out,0:N_Out,0:N_Out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(0:N_Out,0:N_In,0:N_In)  !< first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(0:N_Out,0:N_Out,0:N_In) !< second intermediate results from 1D interpolations
!===================================================================================================================================
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_Out
        X3D_Buf1(iN_Out,jN_In,kN_In)=X3D_Buf1(iN_Out,jN_In,kN_In)+Vdm(iN_Out,iN_In)*X3D_In(iN_In,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Buf2(iN_Out,jN_Out,kN_In)=X3D_Buf2(iN_Out,jN_Out,kN_In)+Vdm(jN_Out,jN_In)*X3D_Buf1(iN_Out,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  DO kN_Out=0,N_Out
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Out(iN_Out,jN_Out,kN_Out)=X3D_Out(iN_Out,jN_Out,kN_Out)+Vdm(kN_Out,kN_In)*X3D_Buf2(iN_Out,jN_Out,kN_In)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D


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
END SUBROUTINE FinalizeShockCapturing

END MODULE MOD_ShockCapturing
