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
!> Computes the time step to be used in the explicit Runge-Kutta schemes
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
!> Calculate the time step for the current update of U for the Maxwell equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_Equation_Vars,ONLY:c,c_corr
USE MOD_TimeDisc_Vars,ONLY:CFLScale
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: CalcTimeStep !< computed time step
INTEGER,INTENT(OUT)            :: errType      !< integer flag for the type of error that occured
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k,iElem
REAL                           :: Max_Lambda1,Max_Lambda2,Max_Lambda3
REAL                           :: TimeStepConv
REAL                           :: buf(2)
!==================================================================================================================================
TimeStepConv=HUGE(1.)
errType=0
DO iElem=1,nElems
  Max_Lambda1=0.
  Max_Lambda2=0.
  Max_Lambda3=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ! Convective Eigenvalues
        Max_Lambda1=MAX(Max_Lambda1,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)))))
        Max_Lambda2=MAX(Max_Lambda2,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)))))
        Max_Lambda3=MAX(Max_Lambda3,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)))))
      END DO ! i
    END DO ! j
  END DO ! k
  TimeStepConv=MIN(TimeStepConv,CFLScale*2./(Max_Lambda1+Max_Lambda2+Max_Lambda3))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ERRWRITE(*,'(A,I0,A,I0,A,I0,A,I0)')'Convective timestep NaN on proc ',myRank,' at global position (iElem): ',iElem
    ERRWRITE(*,*)'dt_conv=',TimeStepConv
    errType=1
  END IF
END DO ! iElem

IF(IEEE_IS_NAN(TimeStepConv))THEN
  errType=1
END IF
#if MPI
buf(1)=TimeStepConv
buf(2)=-errType ! reduce with timestep, minus due to MPI_MIN
CALL MPI_ALLREDUCE(MPI_IN_PLACE,buf,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
TimeStepConv=buf(1)
errType=INT(-buf(2))
#endif /*MPI*/

CalcTimeStep=TimeStepConv

END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
