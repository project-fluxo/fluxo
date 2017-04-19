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
!> Calculate the time step for the current update of U for the Navierstokes-Equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Mesh_Vars     ,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP
USE MOD_Equation_Vars ,ONLY:kappa,kappaM1
USE MOD_TimeDisc_Vars ,ONLY:CFLScale,ViscousTimeStep,dtElem
#if PARABOLIC
USE MOD_Equation_Vars ,ONLY:KappasPr
USE MOD_TimeDisc_Vars ,ONLY:DFLScale
#if PP_VISC == 0
USE MOD_Equation_Vars ,ONLY:mu0
#endif
#if PP_VISC == 1
USE MOD_Equation_Vars ,ONLY:muSuth,R
#endif
#if PP_VISC == 2
USE MOD_Equation_Vars ,ONLY:mu0,ExpoSuth,R
#endif
#endif /*PARABOLIC*/
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep  !< Minimum time step
INTEGER,INTENT(OUT)          :: errType       !< Error code
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: sRho,v(3),p,c
REAL                         :: TimeStepConv, TimeStepVisc,TimeStepViscElem
REAL                         :: maxLambda1,maxLambda2,maxLambda3
#if PARABOLIC
REAL                         :: maxLambda_v1,maxLambda_v2,maxLambda_v3
REAL                         :: muX,KappasPr_max
#endif /*PARABOLIC*/
#if MPI
REAL                         :: buf(3)
#endif /*MPI*/
!==================================================================================================================================
errType=0
errMsg=""
#if PARABOLIC
KappasPr_max=MAX(4./3.,KappasPr)
#endif /*PARABOLIC*/

TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
DO iElem=1,PP_nElems
  maxLambda1=1.0E-12
  maxLambda2=1.0E-12
  maxLambda3=1.0E-12
#if PARABOLIC
  MaxLambda_v1=1.0E-12  ! Viscous
  MaxLambda_v2=1.0E-12  ! Viscous
  MaxLambda_v3=1.0E-12  ! Viscous
#endif /*PARABOLIC*/
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ! Convective Eigenvalues
        IF(IEEE_IS_NAN(U(1,i,j,k,iElem)))THEN
          WRITE(errMsg,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
          ERRWRITE(*,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=1
        END IF
        sRho=1./U(1,i,j,k,iElem)
        v=U(2:4,i,j,k,iElem)*sRho
        p=kappaM1*(U(5,i,j,k,iElem)-0.5*U(1,i,j,k,iElem)*SUM(v*v))
        IF(p.LE.0.)THEN
          WRITE(ErrMsg,'(A,3ES16.7)')'Pressure Negative, Position= ',Elem_xGP(:,i,j,k,iElem)
          ERRWRITE(*,'(A,3ES16.7)')'Pressure Negative, Position= ',Elem_xGP(:,i,j,k,iElem)
          errType=2
        END IF
        c=SQRT(kappa*p*sRho)
        MaxLambda1=MAX(MaxLambda1,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_fTilde(:,i,j,k,iElem)*v)) + &
                        c*SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)))))
        MaxLambda2=MAX(MaxLambda2,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_gTilde(:,i,j,k,iElem)*v)) + &
                        c*SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)))))
        MaxLambda3=MAX(MaxLambda3,sJ(i,j,k,iElem)*(ABS(SUM(Metrics_hTilde(:,i,j,k,iElem)*v)) + &
                        c*SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)))))
#ifdef PARABOLIC
        ! Viscous Eigenvalues
#if   PP_VISC == 0
        muX=sRho*KappasPr_max*mu0   ! Constant mu
#elif PP_VISC == 1
        muX=sRho*KappasPr_max*muSuth(srho*p/R)         ! compute viscosity with Sutherlands law
#elif PP_VISC == 2
        muX=sRho*KappasPr_max*mu0*(srho*p/R)**ExpoSuth ! compute vsicosity using the power-law
#endif /*PP_VISC*/
        MaxLambda_v1=MAX(MaxLambda_v1,muX*(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        MaxLambda_v2=MAX(MaxLambda_v2,muX*(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        MaxLambda_v3=MAX(MaxLambda_v3,muX*(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
#endif /* PARABOLIC*/
      END DO ! i
    END DO ! j
  END DO ! k
  dtElem(iElem)=CFLScale*2./(maxLambda1+maxLambda2+maxLambda3)
  TimeStepConv=MIN(TimeStepConv,dtElem(iElem))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ErrMsg='Convective timestep NaN '
    ERRWRITE(*,'(A,I0,A,I0)')'Convective timestep NaN on proc',myRank,' for element: ',iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_conv=',TimeStepConv,' dt_visc=',TimeStepVisc
    errType=3
  END IF
#if PARABOLIC
  TimeStepViscElem= DFLScale*4./(maxLambda_v1+maxLambda_v2+maxLambda_v3)
  TimeStepVisc= MIN(TimeStepVisc, TimeStepViscElem)
  dtElem(iElem)=MIN(dtElem(iElem),TimeStepViscElem)
  IF(IEEE_IS_NAN(TimeStepVisc))THEN
    ErrMsg='Viscous timestep NaN '
    ERRWRITE(*,'(A,I0,A,I0)')'Viscous timestep NaN on proc ',myRank,' for element: ', iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_visc=',TimeStepVisc,' dt_conv=',TimeStepConv
    errType=4
  END IF
  Max_Lambda_v=0.  ! Viscous
#endif /* PARABOLIC*/
  IF(errType.NE.0)EXIT
END DO ! iElem=1,PP_nElems

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
END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
