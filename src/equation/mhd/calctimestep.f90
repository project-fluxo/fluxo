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


PUBLIC :: CALCTIMESTEP, InitTimeStep
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Calculate the time step for the current update of U for the MHD-Equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,nElems
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale
#endif /*PARABOLIC*/
USE MOD_TimeDisc_Vars,ONLY:CFLScale,CFLScale_usr,ViscousTimeStep,dtElem, FVTimeStep
#if FV_TIMESTEP
USE MOD_NFVSE_Vars   ,ONLY:sWGP
#if NFVSE_CORR
USE MOD_IDP_Vars           ,ONLY: maxdt_IDP
#endif /*NFVSE_CORR*/
#endif /*FV_TIMESTEP*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep  !< Smallest permitted time step
INTEGER,INTENT(OUT)          :: errType       !< Error code
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: sRho,v(3),Prim(1:PP_nVar),cf
REAL                         :: TimeStepConv, TimeStepVisc, TimeStepFV, TimeStep(3)
REAL                         :: Max_Lambda(3), Lambda(3)
#if FV_TIMESTEP
REAL                         :: TimeStepFVElem
#endif /*FV_TIMESTEP*/
#if PARABOLIC
REAL                         :: Max_Lambda_v
REAL                         :: Lambda_v(3)
REAL                         :: muKappasPr_max,diffC_max
#endif /*PARABOLIC*/
#if MPI
REAL                         :: buf(4)
#endif /*MPI*/
!==================================================================================================================================
errType=0
#if PARABOLIC
#if PP_NumComponents==1
! For single-component flows, KappasPr is constant in the domain
muKappasPr_max=mu*MAX(4./3.,KappasPr)
#endif /*PP_NumComponents==1*/
diffC_max=etasmu_0
Max_Lambda_v=0.  ! Viscous
#endif /*PARABOLIC*/

TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
TimeStepFV  =HUGE(1.)
DO iElem=1,nElems
  Max_Lambda=0.
#if FV_TIMESTEP
  TimeStepFVElem =-HUGE(1.) ! Initialize inverse of time-step size
#endif /*FV_TIMESTEP*/
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
#if PP_NumComponents>1 && PARABOLIC
        ! For multi-component flows, KappasPr is NOT constant in the domain
        muKappasPr_max=mu*MAX(4./3.,totalKappa(U(:,i,j,k,iElem))/Pr)
#endif /*PP_NumComponents>1*/
        ! Convective Eigenvalues
        IF(IEEE_IS_NAN(U(IRHO1,i,j,k,iElem)))THEN
          ERRWRITE(*,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=1
        END IF
        sRho=1./U(IRHO1,i,j,k,iElem)
        CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
        IF(prim(IP).LE.0.)THEN
          ERRWRITE(*,'(A,3ES16.7)')'Pressure Negative, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=2
        END IF
        CALL FastestWave3D(Prim,cf)
        v(:)=Prim(IU:IW) 
        Lambda(1) = sJ(i,j,k,iElem)*(ABS(SUM(Metrics_fTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem))))
        Max_Lambda(1)=MAX(Max_Lambda(1),Lambda(1))
        Lambda(2) = sJ(i,j,k,iElem)*(ABS(SUM(Metrics_gTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem))))
        Max_Lambda(2)=MAX(Max_Lambda(2),Lambda(2))
        Lambda(3) = sJ(i,j,k,iElem)*(ABS(SUM(Metrics_hTilde(:,i,j,k,iElem)*v)) + &
                        cf*SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem))))
        Max_Lambda(3)=MAX(Max_Lambda(3),Lambda(3))
#if FV_TIMESTEP
        ! first compute the inverse of the time-step
        TimeStepFVElem = max (TimeStepFVElem, Lambda(1) * sWGP(i), Lambda(2) * sWGP(j), Lambda(3) * sWGP(k)) ! Subcell-local approximation of low-order CFL condition
#endif /*FV_TIMESTEP*/
#if PARABOLIC
        ! Viscous Eigenvalues, isotropic part
        Lambda_v(1)=(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Lambda_v(2)=(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Lambda_v(3)=(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2))
        Max_Lambda_v=MAX(Max_Lambda_v,MAX(diffC_max,sRho*muKappasPr_max)*SUM(lambda_v(:)))
#ifdef PP_ANISO_HEAT
        !isotropic part of anisotropic heat flux I*kperp
        Max_Lambda_v=MAX(Max_Lambda_v,sRho*KappaM1*kperp*SUM(lambda_v(:)))
        !anisotropic part: (Bn Bn^T)(kpar-kperp)
        !compute contravariant components of BnCon and ^2
        ASSOCIATE(bvec=> U(IB1:IB3,i,j,k,iElem))
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
#if FV_TIMESTEP
  TimeStepFVElem=CFLScale_usr*0.5/TimeStepFVElem
  TimeStepFV=MIN(TimeStepFV,TimeStepFVElem)
  dtElem(iElem)=MIN(dtElem(iElem),TimeStepFVElem)
  IF(IEEE_IS_NAN(TimeStepFV))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'FV timestep NaN on proc ',myRank,' for element: ', iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_FV=',TimeStepFV,' dt_conv=',TimeStepConv,'dt_visc=',TimeStepVisc
    errType=5
  END IF
#endif /*FV_TIMESTEP*/
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
TimeStep(3)=TimeStepFV
#if MPI
buf(1:3)=TimeStep(1:3)
buf(4)=-REAL(errType)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,buf,4,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
TimeStep(1:3)=buf(1:3)
errType=NINT(-buf(4))
#endif /*MPI*/
! Correct FVTimeStep with the one computed from IDP LLF method
#if NFVSE_CORR && FV_TIMESTEP
TimeStepFV = min(TimeStepFV,CFLScale_usr*maxdt_IDP)
#endif /*NFVSE_CORR && FV_TIMESTEP*/

ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1)) .and. (TimeStep(2) .LT. TimeStep(3))
FVTimeStep=(TimeStep(3) .LT. TimeStep(1)) .and. (TimeStep(3) .LT. TimeStep(2))

CalcTimeStep=MINVAL(TimeStep) !(1:2)

#ifdef PP_GLM
GLM_ch=GLM_scale*(GLM_dtch1/CalcTimeStep) ! dt~1/ch -> dt/dtch1 = 1/ch -> ch=dtch1/dt
#endif /*PP_GLM*/

END FUNCTION CALCTIMESTEP


!==================================================================================================================================
!> Calculate the time step for the GLM wave system, assuming ch=1. 
!> This is computed only once, then ch is computed frm the current timestep:
!> dt~1/ch -> dt/dtch1 = 1/ch -> ch=dtch1/dt
!==================================================================================================================================
SUBROUTINE InitTimeStep(dt)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_TimeDisc_Vars,ONLY:CFLScale, CFLScale_usr
#ifdef PP_GLM
USE MOD_Equation_Vars,  ONLY: GLM_init,GLM_dtch1,GLM_scale,GLM_ch
#endif /*def PP_GLM*/
#if FV_TIMESTEP
USE MOD_NFVSE_Vars   ,ONLY:sWGP
#endif /*FV_TIMESTEP*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, OPTIONAL :: dt
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: Max_Lambda(3), Lambda(3)
#if FV_TIMESTEP
REAL                         :: TimeStepFVElem
#endif /*FV_TIMESTEP*/
!==================================================================================================================================

#ifdef PP_GLM
GLM_dtch1=HUGE(1.)
DO iElem=1,nElems
  Max_Lambda=0.
#if FV_TIMESTEP
  TimeStepFVElem=-HUGE(1.) ! Initialize inverse of time-step size
#endif /*FV_TIMESTEP*/
  DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
        Lambda(1)=sJ(i,j,k,iElem)*SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)))
        Max_Lambda(1)=MAX(Max_Lambda(1),Lambda(1))
        Lambda(2)=sJ(i,j,k,iElem)*SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)))
        Max_Lambda(2)=MAX(Max_Lambda(2),Lambda(2))
        Lambda(3)=sJ(i,j,k,iElem)*SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)))
        Max_Lambda(3)=MAX(Max_Lambda(3),Lambda(3))
#if FV_TIMESTEP
        TimeStepFVElem = max (TimeStepFVElem, Lambda(1) * sWGP(i), Lambda(2) * sWGP(j),Lambda(3) * sWGP(k)) ! Subcell-local approximation of low-order CFL condition
#endif /*FV_TIMESTEP*/
  END DO; END DO; END DO ! i,j,k
  GLM_dtch1=MIN(GLM_dtch1,CFLScale*2./SUM(Max_Lambda))
#if FV_TIMESTEP
  GLM_dtch1=MIN(GLM_dtch1,CFLScale_usr*0.5/TimeStepFVElem)
#endif /*FV_TIMESTEP*/
END DO !iElem
#if MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,GLM_dtch1,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(A,ES16.7)')'  GLM correction speed from timestep, (GLM_ch=GLM_scale*GLM_dtch1/dt). GLM_dt for ch=1 : ', GLM_dtch1

GLM_init=.TRUE.

if (present(dt)) then
  GLM_ch=GLM_scale*(GLM_dtch1/dt) ! dt~1/ch -> dt/dtch1 = 1/ch -> ch=dtch1/dt
end if
#endif /*PP_GLM*/

END SUBROUTINE InitTimeStep


END MODULE MOD_CalcTimeStep
