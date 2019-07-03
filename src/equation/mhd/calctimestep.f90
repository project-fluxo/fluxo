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
!> Contains a routine to calculate the maximum explicit time step and estimate the wave speed ch for GLM from  current timestep
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
!> Calculate the time step for the current update of U for the MHD-Equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,nElems
USE MOD_Equation_Vars,ONLY: ConsToPrim
USE MOD_Equation_Vars,ONLY: FastestWave3D
#ifdef PP_GLM
USE MOD_Equation_Vars,  ONLY: GLM_init,GLM_dtch1,GLM_ch,GLM_scale
#endif /*PP_GLM*/
#if PARABOLIC
USE MOD_Equation_Vars,ONLY:mu,KappasPr,etasmu_0
USE MOD_TimeDisc_Vars,ONLY:DFLScale
#ifdef PP_ANISO_HEAT
USE MOD_Equation_Vars,ONLY:KappaM1,kperp,kpar
#endif /*PP_ANISO_HEAT*/
#endif /*PARABOLIC*/
USE MOD_TimeDisc_Vars,ONLY:CFLScale,ViscousTimeStep,dtElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep  !< Smallest permitted time step
INTEGER,INTENT(OUT)          :: errType       !< Error code
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: sRho,v(3),Prim(1:PP_nVar),cf
REAL                         :: TimeStepConv, TimeStepVisc, TimeStep(2)
REAL                         :: Max_Lambda(3)
#if PARABOLIC
REAL                         :: Max_Lambda_v
REAL                         :: Lambda_v(3)
REAL                         :: muKappasPr_max,diffC_max
#endif /*PARABOLIC*/
#if MPI
REAL                         :: buf(3)
#endif /*MPI*/
!==================================================================================================================================
errType=0
#if PARABOLIC
muKappasPr_max=mu*MAX(4./3.,KappasPr)
diffC_max=etasmu_0
Max_Lambda_v=0.  ! Viscous
#endif /*PARABOLIC*/
#ifdef PP_GLM
IF(.NOT.GLM_init) CALL InitTimeStep_GLM() !called only once!
#endif /*PP_GLM*/
TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
DO iElem=1,nElems
  Max_Lambda=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ! Convective Eigenvalues
        IF(IEEE_IS_NAN(U(1,i,j,k,iElem)))THEN
          ERRWRITE(*,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=1
        END IF
        sRho=1./U(1,i,j,k,iElem)
        CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
        IF(prim(5).LE.0.)THEN
          ERRWRITE(*,'(A,3ES16.7)')'Pressure Negative, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=2
        END IF
        CALL FastestWave3D(Prim,cf)
        v(:)=Prim(2:4) 
        Max_Lambda(1)=MAX(Max_Lambda(1),sJ(i,j,k,iElem)*(ABS(SUM(Metrics_fTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)))))
        Max_Lambda(2)=MAX(Max_Lambda(2),sJ(i,j,k,iElem)*(ABS(SUM(Metrics_gTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)))))
        Max_Lambda(3)=MAX(Max_Lambda(3),sJ(i,j,k,iElem)*(ABS(SUM(Metrics_hTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)))))
#if PARABOLIC
        ! Viscous Eigenvalues, isotropic part
        Lambda_v(1)=(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Lambda_v(2)=(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Lambda_v(3)=(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Max_Lambda_v=MAX(Max_Lambda_v,MAX(diffC_max,sRho*muKappasPR_max)*SUM(lambda_v(:)))
#ifdef PP_ANISO_HEAT
        !isotropic part of anisotropic heat flux I*kperp
        Max_Lambda_v=MAX(Max_Lambda_v,sRho*KappaM1*kperp*SUM(lambda_v(:)))
        !anisotropic part: (Bn Bn^T)(kpar-kperp)
        !compute contravariant components of BnCon and ^2
        ASSOCIATE(bvec=> U(6:8,i,j,k,iElem))
        !bvec(:)=bvec(:)/SQRT(SUM(bvec(:)*bvec(:))) !normalization, ^2 added below
        Lambda_v(1)=(SUM(bvec(:)*Metrics_fTilde(:,i,j,k,iElem))*sJ(i,j,k,iElem))**2
        Lambda_v(2)=(SUM(bvec(:)*Metrics_gTilde(:,i,j,k,iElem))*sJ(i,j,k,iElem))**2
        Lambda_v(3)=(SUM(bvec(:)*Metrics_hTilde(:,i,j,k,iElem))*sJ(i,j,k,iElem))**2
        Max_Lambda_v=MAX(Max_Lambda_v,((kpar-kperp)*KappaM1*srho/SUM(bvec(:)*bvec(:)))*SUM(lambda_v(:)))
        END ASSOCIATE !bvec
#endif /*PP_ANISO_HEAT*/
#endif /* PARABOLIC*/
      END DO ! i
    END DO ! j
  END DO ! k
  dtElem(iElem)=CFLScale*2./SUM(Max_Lambda)
  TimeStepConv=MIN(TimeStepConv,dtElem(iElem))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Convective timestep NaN on proc',myRank,' for element: ',iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_conv=',TimeStepConv,' dt_visc=',TimeStepVisc
    errType=3
  END IF
#if PARABOLIC
  IF(Max_Lambda_v.GT.0.)THEN
    dtElem(iElem)=MIN(dtElem(iElem),DFLScale*4./Max_Lambda_v)
    TimeStepVisc= MIN(TimeStepVisc, DFLScale*4./Max_Lambda_v)
  END IF
  IF(IEEE_IS_NAN(TimeStepVisc))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Viscous timestep NaN on proc ',myRank,' for element: ', iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_visc=',TimeStepVisc,' dt_conv=',TimeStepConv
    errType=4
  END IF
  Max_Lambda_v=0.  ! Viscous
#endif /* PARABOLIC*/
  IF(errType.NE.0)EXIT
END DO ! iElem=1,nElems
TimeStep(1)=TimeStepConv
TimeStep(2)=TimeStepVisc
#if MPI
buf(1)=TimeStep(1)
buf(2)=TimeStep(2)
buf(3)=-REAL(errType)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,buf,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
TimeStep(1)=buf(1)
TimeStep(2)=buf(2)
errType=NINT(-buf(3))
#endif /*MPI*/
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep)

#ifdef PP_GLM
GLM_ch=GLM_scale*(GLM_dtch1/CalcTimeStep) ! dt~1/ch -> dt/dtch1 = 1/ch -> ch=dtch1/dt
#endif /*PP_GLM*/

END FUNCTION CALCTIMESTEP

#ifdef PP_GLM
!==================================================================================================================================
!> Calculate the time step for the GLM wave system, assuming ch=1. 
!> This is computed only once, then ch is computed frm the current timestep:
!> dt~1/ch -> dt/dtch1 = 1/ch -> ch=dtch1/dt
!==================================================================================================================================
SUBROUTINE InitTimeStep_GLM()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_TimeDisc_Vars,ONLY:CFLScale
USE MOD_Equation_Vars,  ONLY: GLM_init,GLM_dtch1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: Max_Lambda(3)
!==================================================================================================================================
GLM_dtch1=HUGE(1.)
DO iElem=1,nElems
  Max_Lambda=0.
  DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
        Max_Lambda(1)=MAX(Max_Lambda(1),sJ(i,j,k,iElem) &
                      *SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem))))
        Max_Lambda(2)=MAX(Max_Lambda(2),sJ(i,j,k,iElem) &
                      *SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem))))
        Max_Lambda(3)=MAX(Max_Lambda(3),sJ(i,j,k,iElem) &
                          *SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem))))
  END DO; END DO; END DO ! i,j,k
  GLM_dtch1=MIN(GLM_dtch1,CFLScale*2./SUM(Max_Lambda))
END DO !iElem
#if MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,GLM_dtch1,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(A)')       '  GLM correction speed from timestep, (GLM_ch=GLM_scale*GLM_dtch1/dt)'
SWRITE(UNIT_StdOut,'(A,ES16.7)')'    GLM_dt for ch=1 : ', GLM_dtch1

GLM_init=.TRUE.

END SUBROUTINE InitTimeStep_GLM
#endif /*PP_GLM*/

END MODULE MOD_CalcTimeStep
