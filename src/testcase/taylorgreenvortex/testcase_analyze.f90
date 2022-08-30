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
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitAnalyzeTestCase
  MODULE PROCEDURE InitAnalyzeTestCase
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE AnalyzeTestCase
END INTERFACE

PUBLIC:: DefineParametersAnalyzeTestcase
PUBLIC:: InitAnalyzeTestCase
PUBLIC:: AnalyzeTestCase

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeTestcase()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("AnalyzeTestcase")
CALL prms%CreateLogicalOption('doTCanalyze', "switch off/on TC_analyze" , '.FALSE.')
CALL prms%CreateLogicalOption('doPlotGrads', "switch off/on grad(dens) and vorticity plots" , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeTestcase


!==================================================================================================================================
!> Initialize Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE InitAnalyzeTestcase()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,        ONLY: nElems
USE MOD_Analyze_Vars,     ONLY: doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_Testcase_Vars,    ONLY: doTCanalyze,last_Ekin_comp,last_Emag_comp,doPlotGrads,grad2plot
USE MOD_ReadInTools,      ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
!prepare AnalyzeToFile

doTCanalyze   = GETLOGICAL('doTCanalyze','.TRUE.')
doPlotGrads   = GETLOGICAL('doPlotGrads','.FALSE.')

last_Ekin_comp = 0.
last_Emag_comp = 0.

IF(MPIroot.AND.doAnalyzeToFile) THEN
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Dissipation Rate Incompressible"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Dissipation Rate Compressible"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Ekin incomp"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Ekin comp"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Enstrophy comp"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"DR_u"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"DR_S"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"DR_Sd"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"DR_p"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Maximum Vorticity"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Mean Temperature"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"uprime"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Entropy comp"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"-dEkin/dt"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Emag comp"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"-dEmag/dt"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Maximum Current"'
  A2F_iVar=A2F_iVar+1
  A2F_VarNames(A2F_iVar)='"Magnetic Enstrophy"'
END IF !MPIroot & doAnalyzeToFile

! Variables to plot
if (doPlotGrads) then
  ALLOCATE( grad2plot    (7,0:PP_N ,0:PP_N ,0:PP_N ,nElems) )
end if

END SUBROUTINE InitAnalyzeTestcase


!===================================================================================================================================
!> Analyze for insulating or conducting Taylor-Green Vortex configurations from "Ideal evolution of MHD turbulence when imposing 
!> Taylor-Green symmetries" by Brachet et al. (2012)
!===================================================================================================================================
SUBROUTINE AnalyzeTestCase(Time)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,        ONLY: U
#if PARABOLIC
USE MOD_Lifting_Vars,   ONLY: GradPx,GradPy,GradPz
#endif
USE MOD_Analyze_Vars,   ONLY: wGPVol,Vol
USE MOD_Analyze_Vars,   ONLY: A2F_iVar,A2F_Data
USE MOD_Analyze_Vars,   ONLY: doAnalyzeToFile,Analyze_dt
USE MOD_Restart_Vars,   ONLY: RestartTime
USE MOD_Testcase_Vars,  ONLY: doTCanalyze,last_Ekin_comp,last_Emag_comp, doPlotGrads, grad2plot
USE MOD_Mesh_Vars,      ONLY: sJ,nElems
USE MOD_Equation_Vars
USE MOD_TimeDisc_Vars,  ONLY: doWriteData
USE MOD_Mesh_Vars,      ONLY: Elem_xGP
USE MOD_Output,         ONLY: VisualizeAny
USE MOD_Restart_Vars   ,ONLY: DoRestart,RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,i,j,k,p,q
REAL                            :: Vel(1:3),GradVel(1:3,1:3),Ekin,uprime,BField(1:3),GradB(1:3,1:3)
REAL                            :: Vorticity(1:3),max_Vorticity,mean_temperature,temperature,curlB(1:3),max_Current
REAL                            :: S(1:3,1:3)                    !< Strain rate tensor S (symmetric)
REAL                            :: Sd(1:3,1:3)                   !< Deviatoric part of the strain rate tensor S 
REAL                            :: divU                          !< Divergence of velocity vector
REAL                            :: kE                            !< Integrand: 0.5*v*v
REAL                            :: kE_comp                       !< Integrand: 0.5*rho*v*v
REAL                            :: magE_comp                     !< Integrand: 0.5*B*B
REAL                            :: ens                           !< Integrand: 0.5*rho*omega*omega
REAL                            :: ent                           !< Integrand: -rho*s/kappaM1, s = log(p) - kappa*log(rho)
REAL                            :: eps3                          !< Integrand: p*(div u)
REAL                            :: mag_ens                       !< Integrand: 0.5*j*j, j = curl(B)
REAL                            :: u_tens, s_tens, sd_tens       !< matrix : matrix product, integrands of Gradvel, S, Sd
REAL                            :: Intfactor                     !< Integrationweights with Jacobian
REAL                            :: rho,srho                      !< rho,1/rho
REAL                            :: Ekin_comp,Enstrophy_comp,Entropy_comp,Emag_comp,mag_comp,Mag_Enstrophy
REAL                            :: DR_u,DR_S,DR_Sd,DR_p          !< Contributions to dissipation rate
REAL                            :: Pressure,rho0,negdEkindt,negdEmagdt
#ifdef MPI
REAL                            :: buf(11)
#endif
#if PP_NumComponents>1
REAL :: Kappa
#endif /*PP_NumComponents>1*/
CHARACTER(LEN=255)            :: FileTypeStr
CHARACTER(LEN=255)            :: VarNamesIn(7)
LOGICAL :: doPlot
!===================================================================================================================================
IF(.NOT.doTCanalyze)RETURN
Ekin=0.
Ekin_comp=0.
Emag_comp=0.
Enstrophy_comp=0.
Entropy_comp=0.
Mag_Enstrophy=0.

DR_u=0.;DR_S=0.;DR_Sd=0.;DR_p=0.
max_Vorticity=-1.
max_Current=-1.
mean_Temperature=0.

DO iElem=1,nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        rho=sum(U(IRHO1:PP_NumComponents,i,j,k,iElem))
        srho=1./rho
        ! compute primitive gradients (of u,v,w) at each GP
        Vel(1:3)    =U(IRHOU:IRHOW,i,j,k,iElem)*srho
        BField(1:3) =U(IB1:IB3,i,j,k,iElem)
#if PARABOLIC
        GradVel(:,1)=GradPx(IU:IW,i,j,k,iElem)
        GradVel(:,2)=GradPy(IU:IW,i,j,k,iElem)
        GradVel(:,3)=GradPz(IU:IW,i,j,k,iElem)
        GradB(:,1)=GradPx(IB1:IB3,i,j,k,iElem)
        GradB(:,2)=GradPy(IB1:IB3,i,j,k,iElem)
        GradB(:,3)=GradPz(IB1:IB3,i,j,k,iElem)
        if (doPlotGrads) then
          grad2plot(1,i,j,k,iElem) = sum(GradPx(IRHO1:PP_NumComponents,i,j,k,iElem))
          grad2plot(2,i,j,k,iElem) = sum(GradPy(IRHO1:PP_NumComponents,i,j,k,iElem))
          grad2plot(3,i,j,k,iElem) = sum(GradPz(IRHO1:PP_NumComponents,i,j,k,iElem))
        end if
#else
        GradVel=0.
        GradB=0.
        if (doPlotGrads) grad2plot=0.
#endif
        ! Pressure
        call Get_Pressure(U(:,i,j,k,iElem),Pressure)
        
        ! compute divergence of velocity
        divU=GradVel(1,1)+GradVel(2,2)+GradVel(3,3)
        ! compute tensor of velocity gradients
        S=0.5*(Gradvel+TRANSPOSE(GradVel))
        ! deviatoric part of strain tensor
        Sd=S
        DO p=1,3
          Sd(p,p)=Sd(p,p)-1./3.*divU
        END DO
        ! compute kinetic energy integrand (incomp)
        kE=0.5*SUM(Vel(1:3)*Vel(1:3))
        ! compute kinetic energy integrand (compr)
        kE_comp=rho*kE
        ! compute the magnetic energy integrand (incomp or comp)
        mag_comp=s2mu_0*SUM(BField(1:3)*BField(1:3))
        ! compute vorticity and max(vorticity)
        Vorticity(1)=GradVel(3,2) - GradVel(2,3)
        Vorticity(2)=GradVel(1,3) - GradVel(3,1)
        Vorticity(3)=GradVel(2,1) - GradVel(1,2)
        max_Vorticity=MAX(max_Vorticity,SQRT(SUM(Vorticity(:)*Vorticity(:))))
        ! Store variable to plot
        if (doPlotGrads) then
          grad2plot(4:6,i,j,k,iElem) = Vorticity
          
#if PP_NumComponents>1
          Kappa = totalKappa(U(:,i,j,k,iElem))
#endif /*PP_NumComponents>1*/
          grad2plot(7,i,j,k,iElem) = sqrt(sum(Vel**2))/sqrt(Kappa*Pressure*srho)
        end if
        
        ! compute enstrophy integrand
        ens=0.5*rho*SUM(Vorticity(1:3)*Vorticity(1:3))
        ! compute curl(B) (the current) and max(current)
        curlB(1)=GradB(3,2) - GradB(2,3)
        curlB(2)=GradB(1,3) - GradB(3,1)
        curlB(3)=GradB(2,1) - GradB(1,2)
        max_Current=MAX(max_Current,SQRT(SUM(curlB(:)*curlB(:))))
        ! compute magnetic enstrophy integrand
        mag_ens=s2mu_0*SUM(curlB(1:3)*curlB(1:3))
        ! compute the entropy integrand
        ent=Get_MathEntropy(U(:,i,j,k,iElem))
        ! compute integrand for epsilon3, pressure contribution to dissipation (compressiblity effect)
        eps3=Pressure*divU
        ! Matrix : Matrix product for velocity gradient tensor, S:S and Sd:Sd
        u_tens=0.;s_tens=0.;sd_tens=0.
        DO p=1,3
          DO q=1,3
            u_tens=u_tens+GradVel(p,q)*GradVel(p,q)
            s_tens=s_tens+S(p,q)*S(p,q)
            sd_tens=sd_tens+Sd(p,q)*Sd(p,q)
          END DO
        END DO
        Intfactor=wGPVol(i,j,k)/sJ(i,j,k,iElem)
        ! compute integrals:
          ! Kinetic Energy incompressible
        Ekin=Ekin+kE*Intfactor
          ! Kinetic Energy compressible
        Ekin_comp=Ekin_comp+kE_comp*IntFactor
          ! Magnetic Energy compressible
        Emag_comp=Emag_comp+mag_comp*IntFactor
          ! Enstrophy compressible
        Enstrophy_comp=Enstrophy_comp+ens*IntFactor
          ! Magnetic enstrophy (norm of the current)
        Mag_Enstrophy=Mag_Enstrophy+mag_ens*IntFactor
          ! Entropy compressible
        Entropy_comp=Entropy_comp+ent*IntFactor
          ! dissipation rate epsilon incomp from velocity gradient tensor (Diss Fauconnier)
        DR_u=DR_u+u_tens*IntFactor
          ! dissipation rate epsilon incomp from strain rate tensor (incomp) (Sagaut)
        DR_S=DR_S+S_tens*IntFactor
          ! dissipation rate epsilon 1 from deviatoric part of strain rate tensor Sd (compressible)
        DR_SD=DR_SD+sd_tens*Intfactor
          ! dissipation rate epsilon 3 from pressure times div u (compressible)
        DR_p=DR_p+eps3*Intfactor
          ! compute mean temperature
        Temperature=Get_Temperature(U(:,i,j,k,iElem))
        mean_temperature=mean_temperature+Temperature*Intfactor       
      END DO
    END DO
  END DO
END DO

#ifdef MPI
! Communicate quantities to be maximized
buf( 1) = max_Vorticity
buf( 2) = max_Current
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,buf,2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(buf  ,0         ,2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
max_Vorticity    = buf( 1)
max_Current      = buf( 2)

! Communicate quantities to be summed
buf( 1) = DR_u
buf( 2) = DR_S
buf( 3) = DR_Sd
buf( 4) = DR_p
buf( 5) = Ekin
buf( 6) = Ekin_comp
buf( 7) = Enstrophy_comp
buf( 8) = Entropy_comp
buf( 9) = mean_temperature
buf(10) = Emag_comp
buf(11) = Mag_Enstrophy
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,buf  ,11,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(buf  ,0           ,11,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
DR_u             = buf( 1)
DR_S             = buf( 2)
DR_Sd            = buf( 3)
DR_p             = buf( 4)
Ekin             = buf( 5)
Ekin_comp        = buf( 6)
Enstrophy_comp   = buf( 7)
Entropy_comp     = buf( 8)
mean_temperature = buf( 9)
Emag_comp        = buf(10)
Mag_Enstrophy    = buf(11)
#endif /*MPI*/
IF(MPIroot)THEN
  ! some turbulent quantities
  uprime=SQRT(Ekin/Vol*2./3.)
  ! for TGV = nu = mu, since rho = 1 = const
  !lambda=SQRT(15*mu0*uprime**2/(DR_u*mu0/(Vol)))
  !nu=(mu0**3./((DR_u*mu0/Vol)))**0.25
  !tnu=SQRT(mu0/(DR_u*mu0/Vol))
  !Rlambda=uprime*lambda/mu0
  
  ! now do the normalization of integrals
  ! warning, rho0=1 for the TGV runs, not in general case!
  rho0=1.
  Ekin=Ekin/Vol
  Ekin_comp=Ekin_comp/(rho0*Vol)
  Emag_comp=Emag_comp/(rho0*Vol)

  Enstrophy_comp=Enstrophy_comp/(rho0*Vol)
  Mag_Enstrophy=Mag_Enstrophy/(rho0*Vol)
  Entropy_comp=Entropy_comp/(rho0*Vol)
  DR_u=DR_u*mu_0/Vol
  DR_S=DR_S*2.*mu_0/(rho0*Vol)
  DR_SD=DR_SD*2.*mu_0/(rho0*Vol)
  DR_p=-DR_p/(rho0*Vol)
  
  mean_temperature=mean_temperature/Vol
  
  WRITE(UNIT_StdOut,'(A,E21.11)')  ' TGV Analyze    Ekin_comp   : ',Ekin_comp
  WRITE(UNIT_StdOut,'(A,E21.11)')  ' TGV Analyze    Emag_comp   : ',Emag_comp
  negdEkindt = 0.
  negdEmagdt = 0.
  IF((time-RestartTime).GE.Analyze_dt)THEN
    negdEkindt = -(Ekin_comp - last_Ekin_comp)/Analyze_dt
    last_Ekin_comp=Ekin_comp
    WRITE(UNIT_StdOut,'(A,E21.11)')' TGV Analyze -dEkin_comp/dt : ',negdEkindt
    negdEmagdt = -(Emag_comp - last_Emag_comp)/Analyze_dt
    last_Emag_comp=Emag_comp
    WRITE(UNIT_StdOut,'(A,E21.11)')' TGV Analyze -dEmag_comp/dt : ',negdEmagdt
  END IF !time>Analyze_dt

  IF(doAnalyzeToFile) THEN
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=DR_S             !"Dissipation Rate Incompressible"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=DR_Sd+DR_p       !"Dissipation Rate Compressible"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=Ekin             !"Ekin incomp"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=Ekin_comp        !"Ekin comp"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=Enstrophy_comp   !"Enstrophy comp"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=DR_u             !"DR_u"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=DR_S             !"DR_S"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=DR_Sd            !"DR_Sd"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=DR_p             !"DR_p"'
    A2F_iVar=A2F_iVar+1                     
    A2F_data(A2F_iVar)=max_Vorticity    !"Maximum Vorticity"'
    A2F_iVar=A2F_iVar+1           
    A2F_data(A2F_iVar)=mean_temperature !"Mean Temperature"'
    A2F_iVar=A2F_iVar+1           
    A2F_data(A2F_iVar)=uprime           !"uprime"'
    A2F_iVar=A2F_iVar+1           
    A2F_data(A2F_iVar)=Entropy_comp     !"Entropy comp"'
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=negdEkindt       !'"-dEkin/dt"'
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=Emag_comp        !'"Emag comp"'
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=negdEmagdt       !'"-dEmag/dt"'
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=max_Current      !'"Maximum Current"'
    A2F_iVar=A2F_iVar+1
    A2F_data(A2F_iVar)=Mag_Enstrophy    !'"Magnetic Enstrophy"'
  END IF !doAnalyzeToFile
END IF !MPIroot

! Plot variables
! --------------
IF (.NOT. doPlotGrads) RETURN
! check if we need to plot
doPlot = .FALSE.
IF(.NOT.DoRestart)THEN
  doPlot = (Time==0.0)
ELSE
  doPlot = (Time==RestartTime)
END IF
IF(.NOT.doPlot) doPlot = doWriteData

! Plot
FileTypeStr='Gradients'
VarNamesIn=(/'DensGradX','DensGradY','DensGradZ','VorticityX','VorticityY','VorticityZ','Mach'/)    
IF(doPlot) CALL VisualizeAny(Time,7,PP_N,.TRUE.,Elem_xGP,grad2plot,FileTypeStr,VarNamesIn)

END SUBROUTINE AnalyzeTestCase

!TODO: FinalizeAnalyzeTestCase

END MODULE MOD_Testcase_Analyze
