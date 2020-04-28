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
!> Contains analyze routines specific to the navierstokes advection equation
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC:: InitAnalyzeEquation
PUBLIC:: AnalyzeEquation
PUBLIC:: FinalizeAnalyzeEquation
PUBLIC:: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS
!==================================================================================================================================
!> Define parameters for analyze navierstokes 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcBodyForces'  , "Set true to compute integral pressure and viscous forces on wall boundaries"  &
                                                , '.FALSE.')
CALL prms%CreateLogicalOption('CalcBulkVelocity',"Set true to compute mean velocity,momentum and density integrated over domain" &
                                                , '.FALSE.')
CALL prms%CreateLogicalOption('CalcWallVelocity', "Set true to compute min/max/mean velocity at each wall boundary separately"   &
                                                , '.FALSE.')
CALL prms%CreateLogicalOption('CalcEntropy', "Set true to compute the integrated entropy"&
           , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeEquation

!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_Mesh_Vars,          ONLY: BoundaryName,nBCs,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i
!==================================================================================================================================
! Get the various analysis/output variables 
doCalcBodyForces    =GETLOGICAL('CalcBodyForces','.FALSE.')
doCalcBulkVelocity  =GETLOGICAL('CalcBulkVelocity','.FALSE.')
doCalcWallVelocity  =GETLOGICAL('CalcWallVelocity','.FALSE.')
doCalcEntropy       = GETLOGICAL('CalcEntropy'   ,'.FALSE.')

! Initialize eval routines
IF(doCalcWallVelocity) ALLOCATE(meanV(nBCs),maxV(nBCs),minV(nBCs))

IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doCalcBodyForces)THEN
    A2F_iVar=A2F_iVar+3
    A2F_VarNames(A2F_iVar-2:A2F_iVar)=(/'"BodyForce_X"','"BodyForce_Y"','"BodyForce_Z"'/)
    A2F_iVar=A2F_iVar+3
    A2F_VarNames(A2F_iVar-2:A2F_iVar)=(/'"BodyForceP_X"','"BodyForceP_Y"','"BodyForceP_Z"'/)
    A2F_iVar=A2F_iVar+3
    A2F_VarNames(A2F_iVar-2:A2F_iVar)=(/'"BodyForceV_X"','"BodyForceV_Y"','"BodyForceV_Z"'/)
  END IF  !doCalcDivergence
  IF(doCalcBulkVelocity)THEN
    A2F_iVar=A2F_iVar+3
    A2F_VarNames(A2F_iVar-2:A2F_iVar)=(/'"BulkVel_X"','"BulkVel_Y"','"BulkVel_Z"'/)
    A2F_iVar=A2F_iVar+3
    A2F_VarNames(A2F_iVar-2:A2F_iVar)=(/'"BulkMom_X"','"BulkMom_Y"','"BulkMom_Z"'/)
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"BulkDensity"'
  END IF ! doCalcBulk
  IF(doCalcWallVelocity)THEN
    DO i=1,nBCs
      IF((BoundaryType(i,BC_TYPE).NE.4).AND.((BoundaryType(i,BC_TYPE).NE.9))) CYCLE
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"WallVel_'//TRIM(BoundaryName(i))//'_mean"'
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"WallVel_'//TRIM(BoundaryName(i))//'_min"'
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"WallVel_'//TRIM(BoundaryName(i))//'_max"'
    END DO !i=1,nBCs
  END IF !doCalcEnergy
  IF(doCalcEntropy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"Entropy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"dSdU_Ut"'
  END IF !doCalcEntropy
END IF !MPIroot & doAnalyzeToFile
END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> execute the analyze steps 
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_Restart_Vars,       ONLY: RestartTime
USE MOD_Mesh_Vars,          ONLY: BoundaryName,nBCs,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=40)               :: formatStr
REAL                            :: Fv(3),Fp(3),BodyForce(3) ! Viscous force, pressure force and surface area, resulting body force
REAL                            :: bulkVel(3),bulkDensity,bulkMom(3)
REAL                            :: tmp(1),dSdU_Ut 
INTEGER                         :: i
!==================================================================================================================================
! Attention: during the initialization phase no face data / gradients available!
IF(doCalcBodyforces)THEN
  CALL CalcBodyforces(Time,Bodyforce,Fp,Fv)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',3,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' BodyForce  : ',BodyForce
    WRITE(UNIT_StdOut,formatStr)' BodyForceP : ',Fp
    WRITE(UNIT_StdOut,formatStr)' BodyForceV : ',Fv
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=BodyForce
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=Fp
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=Fv
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF  !(doCalcBodyforces)

! Calculate bulk velocity 
IF(doCalcBulkVelocity)THEN
  CALL CalcBulkVelocity(bulkVel,bulkMom,bulkDensity)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',3,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' Bulk Vel   : ',bulkVel
    WRITE(UNIT_StdOut,formatStr)' Bulk Mom   : ',bulkMom
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',1,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' Bulk Dens  : ',bulkDensity
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=bulkVel
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=bulkMom
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=bulkDensity
    END IF !doAnalyzeToFile
  END IF
END IF

! Calculate wall velocities / only terminal output 
IF(doCalcWallVelocity)THEN
  CALL CalcWallVelocity(maxV,minV,meanV)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',3,'ES16.7)'
    DO i=1,nBCs
      IF((BoundaryType(i,BC_TYPE).NE.4).AND.((BoundaryType(i,BC_TYPE).NE.9))) CYCLE
      WRITE(UNIT_StdOut,*)'Wall velocities for ',TRIM(BoundaryName(i)),' (mean/min/max) : '
      WRITE(UNIT_StdOut,formatStr)'              ',meanV(i), minV(i), maxV(i)
      IF(doAnalyzeToFile)THEN
        A2F_iVar=A2F_iVar+1
        A2F_Data(A2F_iVar)=meanV(i)
        A2F_iVar=A2F_iVar+1
        A2F_Data(A2F_iVar)=minV(i)
        A2F_iVar=A2F_iVar+1
        A2F_Data(A2F_iVar)=maxV(i)
      END IF !doAnalyzeToFile
    END DO !i=1,nBCs
  END IF
END IF  !(doCalcWallVelocity)

IF(doCalcEntropy)THEN
  tmp(1)=Entropy
  CALL CalcEntropy(Entropy,dSdU_Ut)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A)')'(A21,ES21.12)'
    WRITE(UNIT_StdOut,formatStr)  ' Entropy      : ',Entropy
    IF((time-RestartTime).GE.Analyze_dt)THEN
      WRITE(UNIT_StdOut,formatStr)' dEntropy/dt  : ',(Entropy-tmp(1))/Analyze_dt
    END IF !time>Analyze_dt
    WRITE(UNIT_StdOut,formatStr)  ' dSdU*Ut         : ',dSdU_Ut
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Entropy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=dSdU_Ut
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEntropy
!
END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Calculates bulk velocities over whole domain
!==================================================================================================================================
SUBROUTINE CalcBulkVelocity(BulkVel,BulkMom,BulkDensity)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                :: BulkVel(3),BulkDensity,BulkMom(3)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: IntegrationWeight,Moment(3)
INTEGER                         :: iElem,i,j,k
#if MPI
REAL                            :: box(7)
#endif
!==================================================================================================================================
! Needed for the computation of the forcing term for the channel flow
BulkMom=0.
BulkVel=0.
BulkDensity  =0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    Moment          =U(2:4,i,j,k,iElem)*IntegrationWeight
    BulkMom          =BulkMom+Moment
    BulkVel          =BulkVel+Moment/U(1,i,j,k,iElem)
    BulkDensity      =BulkDensity+U(1,i,j,k,iElem)*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
Box(1:3)=BulkMom; Box(4:6)=BulkVel; Box(7)=BulkDensity
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,7,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  BulkMom=Box(1:3); BulkVel=Box(4:6); BulkDensity=Box(7)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,7,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/

BulkVel    =BulkVel/Vol
BulkMom    =BulkMom/Vol
BulkDensity=BulkDensity/Vol

END SUBROUTINE CalcBulkVelocity 


!==================================================================================================================================
!>
!==================================================================================================================================

SUBROUTINE CalcWallVelocity(maxV,minV,meanV)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: U_Master
USE MOD_Mesh_Vars,         ONLY: SurfElem
USE MOD_Mesh_Vars,         ONLY: nBCSides,BC,BoundaryType,nBCs
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: maxV(nBCs),minV(nBCs),meanV(nBCs)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA,Vel(3),locV
INTEGER                        :: SideID,i,j,iBC
!==================================================================================================================================
minV =  1.e14
maxV = -1.e14
meanV= 0.
DO SideID=1,nBCSides
  iBC=BC(SideID)
  IF((BoundaryType(iBC,BC_TYPE).EQ.4).OR.(BoundaryType(iBC,BC_TYPE).EQ.9))THEN
    DO j=0,PP_N; DO i=0,PP_N
      Vel=U_Master(2:4,i,j,SideID)/U_Master(1,i,j,SideID)
      ! Calculate velocity magnitude
      locV=SQRT(Vel(1)*Vel(1)+Vel(2)*Vel(2)+Vel(3)*Vel(3))
      maxV(iBC)=MAX(maxV(iBC),locV)
      minV(iBC)=MIN(minV(iBC),locV)
      dA=wGPSurf(i,j)*SurfElem(i,j,SideID)
      meanV(iBC)=meanV(iBC)+locV*dA
    END DO; END DO
  END IF
END DO

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,maxV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,minV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,meanV,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(maxV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(minV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(meanV       ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
DO iBC=1,nBCs
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
  MeanV(iBC)=MeanV(iBC)/Surf(iBC)
END DO

END SUBROUTINE CalcWallVelocity


!==================================================================================================================================
!> Integrate the forces on the wall boundary condition
!==================================================================================================================================
SUBROUTINE CalcBodyForces(Time,BodyForce,Fp,Fv)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,         ONLY: U_Master
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY: gradPx_Master,gradPy_Master,gradPz_Master
#endif /*PARABOLIC*/
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nBCSides,BC,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)                :: Time 
REAL,INTENT(OUT)               :: Fv(3),Fp(3),BodyForce(3) ! Viscous force, pressure force and surface area, resulting body force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Fp_loc(3),Fv_loc(3)
INTEGER                        :: BCType
INTEGER                        :: SideID
#if MPI
REAL                           :: Box(6)
#endif /*MPI*/
!==================================================================================================================================
! Calculate body forces  ! Attention: during the initialization phase no face data / gradients available!
Fp=0.
Fv=0.
BodyForce=0.
DO SideID=1,nBCSides
  BCType=Boundarytype(BC(SideID),BC_TYPE)
  ! Calculate pressure force (Euler wall / Navier-Stokes wall)
  IF((BCType .EQ. 9) .OR. (BCType .EQ. 4))THEN
    CALL CalcPressureForce(Fp_loc,U_Master(:,:,:,SideID),SurfElem(:,:,SideID),NormVec(:,:,:,SideID))
    Fp=Fp+Fp_loc
  END IF
#if PARABOLIC
  ! Calculate viscous force (Navier-Stokes wall)
  IF(BCType .EQ. 4)THEN
    CALL CalcViscousForce(Fv_loc,                     &
                          U_Master(:,:,:,SideID),      &
                          gradPx_Master(:,:,:,SideID), &
                          gradPy_Master(:,:,:,SideID), &
                          gradPz_Master(:,:,:,SideID), &
                          SurfElem(:,:,SideID),       &
                          NormVec(:,:,:,SideID))
    Fv=Fv+Fv_loc
  END IF
#endif /*PARABOLIC*/
END DO

#if MPI
Box(1:3)=Fv; Box(4:6)=Fp
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Box,6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  Fv=Box(1:3); Fp=Box(4:6)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/


END SUBROUTINE CalcBodyForces



!==================================================================================================================================
!> Integrate pressure force on wall boundary condition
!==================================================================================================================================
SUBROUTINE CalcPressureForce(Fp,U_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,     ONLY:KappaM1
USE MOD_Analyze_Vars,      ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: U_Face(PP_nVar,0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)              :: Fp(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: p,dA
INTEGER                        :: i, j
!==================================================================================================================================
Fp=0.
DO j=0,PP_N; DO i=0,PP_N
  ! Calculate pressure
  p =KappaM1* (U_Face(5,i,j)-0.5*SUM(U_Face(2:4,i,j)*U_Face(2:4,i,j))/U_Face(1,i,j))
  dA=wGPSurf(i,j)*SurfElem(i,j)
  Fp=Fp+p*NormVec(:,i,j)*dA
END DO; END DO
END SUBROUTINE CalcPressureForce


#if PARABOLIC
!==================================================================================================================================
!> Integrate the viscous forces on the wall boundary condition
!==================================================================================================================================
SUBROUTINE CalcViscousForce(Fv,U_Face,gradPx_Face,gradPy_Face,gradPz_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:s23
USE MOD_Analyze_Vars  ,ONLY:wGPSurf
#if PP_VISC == 0
USE MOD_Equation_Vars ,ONLY:mu0
#endif
#if PP_VISC == 1
USE MOD_Equation_Vars ,ONLY:KappaM1,R,muSuth
#endif
#if PP_VISC == 2
USE MOD_Equation_Vars ,ONLY:KappaM1,R,mu0,ExpoPow
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: U_Face(PP_nVar,0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: gradPx_Face(PP_nVar,0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: gradPy_Face(PP_nVar,0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: gradPz_Face(PP_nVar,0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_N)
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)              :: Fv(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: srho                      ! 1/rho
REAL                           :: tau(3,3)                  ! Viscous stress tensor
REAL                           :: muS
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                           :: p,T
#endif
REAL                           :: GradV(3,3),vel(3),DivV
INTEGER                        :: i, j
!==================================================================================================================================
Fv       =0.

DO j=0,PP_N; DO i=0,PP_N
  srho=1./U_Face(1,i,j)
! Constant mu
#if PP_VISC == 0
  muS=mu0
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
  ! Calculate pressure
  p=KappaM1*(U_Face(5,i,j)-0.5*srho*SUM(U_Face(2:4,i,j)*U_Face(2:4,i,j)))
  ! Calculate temperature
  T=srho*p/R
#endif
! mu-Sutherland
#if PP_VISC == 1
  muS=muSuth(T)
#endif
! mu power-law
#if PP_VISC == 2
  muS=mu0*T**ExpoPow  ! mu0=mu0/T0^n
#endif
  ! Compute velocity derivatives via product rule (a*b)'=a'*b+a*b'
  ! GradV(:,1)=x-gradient of V
  vel = U_Face(2:4,i,j)*srho
  GradV(1:3,1)=gradPx_Face(2:4,i,j)
  ! GradV(:,2)=y-gradient of V
  GradV(1:3,2)=gradPy_Face(2:4,i,j)
  ! GradV(:,3)=z-gradient of V
  GradV(1:3,3)=gradPz_Face(2:4,i,j)
  ! Velocity divergence
  DivV=GradV(1,1)+GradV(2,2)+GradV(3,3)
  ! Calculate shear stress tensor
  tau=muS*(GradV + TRANSPOSE(GradV))
  tau(1,1)=tau(1,1)-s23*muS*DivV
  tau(2,2)=tau(2,2)-s23*muS*DivV
  tau(3,3)=tau(3,3)-s23*muS*DivV
  ! Calculate viscous force vector
  Fv=Fv+MATMUL(tau,NormVec(:,i,j))*wGPSurf(i,j)*SurfElem(i,j)
END DO; END DO

Fv=-Fv  ! Change direction to get the force acting on the wall
END SUBROUTINE CalcViscousForce
#endif /*PARABOLIC*/

!==================================================================================================================================
!> Calculates  Entropy over whole domain Entropy=-rho*s/(kappa-1), s=ln(p rho^(-kappa))=ln(p)-kappa*ln(rho) 
!> and integral of q^T*Ut, which is the semi-discrete entropy update
!==================================================================================================================================
SUBROUTINE CalcEntropy(Entropy,dSdU_Ut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U,Ut
USE MOD_Equation_Vars,      ONLY: kappa,sKappaM1,ConsToPrim,ConsToEntropy
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)            :: Entropy !< Entropy
REAL,INTENT(OUT)            :: dSdU_Ut !< dS/dU * U_t
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: iElem,i,j,k
REAL                        :: ent_loc,prim(PP_nVar),dSdU(PP_nVar)
#if MPI
REAL                        :: box(2)
#endif 
!==================================================================================================================================
Entropy=0.
dSdU_Ut=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim(prim,U(:,i,j,k,iElem))
    ent_loc  = -prim(1)*(LOG(prim(5))-kappa*LOG(prim(1)))
    Entropy  = Entropy+ent_loc*wGPVol(i,j,k)/sJ(i,j,k,iElem)
    dSdU(:)  = ConsToEntropy(U(:,i,j,k,iElem))
    ent_loc  = SUM(dSdU(:)*Ut(:,i,j,k,iElem))
    dSdU_Ut  = dSdU_Ut+ent_loc*wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
box(1) = Entropy
box(2) = dSdU_Ut
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  Entropy = box(1)
  dSdU_Ut = box(2)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
Entropy=Entropy*sKappaM1


END SUBROUTINE CalcEntropy


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(minV)
SDEALLOCATE(maxV)
SDEALLOCATE(meanV)
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
