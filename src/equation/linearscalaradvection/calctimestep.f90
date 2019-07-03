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
!> Contains a routine to calculate the maximum explicit time step.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Calculate the time step for the current update of U for the Linear Scalar Advection Equation du/dt + a du/dx = 0
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems,sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP
USE MOD_Equation_Vars,ONLY:AdvVel
USE MOD_TimeDisc_Vars,ONLY:CFLScale,ViscousTimeStep,dtElem
USE MOD_PreProc
#if PARABOLIC
USE MOD_Equation_Vars,ONLY:DiffC
USE MOD_TimeDisc_Vars,ONLY:DFLScale
#endif /*PARABOLIC*/
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep  !< Smallest permitted time step
INTEGER,INTENT(OUT)          :: errType       !< Error code
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: maxLambda1,maxLambda2,maxLambda3
REAL                         :: TimeStepConv, TimeStepVisc,TimeStepViscElem
#if PARABOLIC
REAL                         :: maxLambda_v1,maxLambda_v2,maxLambda_v3
#endif /*PARABOLIC*/
#if MPI
REAL                         :: buf(3)
#endif /*MPI*/
!==================================================================================================================================
errType=0
TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
DO iElem=1,nElems
  maxLambda1=1.0E-12
  maxLambda2=1.0E-12
  maxLambda3=1.0E-12
#if PARABOLIC
  maxLambda_v1=1.0E-12
  maxLambda_v2=1.0E-12
  maxLambda_v3=1.0E-12
#endif /*PARABOLIC*/
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        maxLambda1=MAX(maxLambda1,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_fTilde(:,i,j,k,iElem)*AdvVel(:)))))
        maxLambda2=MAX(maxLambda2,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_gTilde(:,i,j,k,iElem)*AdvVel(:)))))
        maxLambda3=MAX(maxLambda3,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_hTilde(:,i,j,k,iElem)*AdvVel(:)))))
#if PARABOLIC
        maxLambda_v1=MAX(maxLambda_v1,DiffC*(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        maxLambda_v2=MAX(maxLambda_v2,DiffC*(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        maxLambda_v3=MAX(maxLambda_v3,DiffC*(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
#endif /* PARABOLIC*/
      END DO ! i
    END DO ! j
  END DO ! k          
  dtElem(iElem)=CFLScale*2./(maxLambda1+maxLambda2+maxLambda3)
  TimeStepConv=MIN(TimeStepConv,dtElem(iElem))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Convective timestep NaN on proc',myRank,' for element: ',iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_conv=',TimeStepConv,' dt_visc=',TimeStepVisc
    errType=2
  END IF
#if PARABOLIC
  TimeStepViscElem= DFLScale*4./(MaxLambda_v1+MaxLambda_v2+MaxLambda_v3)
  TimeStepVisc= MIN(TimeStepVisc, TimeStepViscElem)
  dtElem(iElem)=MIN(dtElem(iElem),TimeStepViscElem)
  IF(IEEE_IS_NAN(TimeStepVisc))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Viscous timestep NaN on proc ',myRank,' for element: ', iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_visc=',TimeStepVisc,' dt_conv=',TimeStepConv
    errType=3
  END IF
#endif /* PARABOLIC*/
END DO ! iElem=1,nElems
#if MPI
buf(1)=TimeStepConv
buf(2)=TimeStepVisc
buf(3)=-REAL(errType)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,buf,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
TimeStepConv=buf(1)
TimeStepVisc=buf(2)
errType=NINT(-buf(3))
#endif /*MPI*/
ViscousTimeStep=(TimeStepVisc .LT. TimeStepConv)
CalcTimeStep=MIN(TimeStepConv,TimeStepVisc)
END FUNCTION CalcTimeStep

END MODULE MOD_CalcTimeStep
