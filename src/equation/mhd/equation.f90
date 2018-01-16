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
!> Module for the MHD equations
!==================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE

INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc 
END INTERFACE

INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC:: InitEquation
PUBLIC:: FillIni
PUBLIC:: ExactFunc
PUBLIC:: CalcSource
PUBLIC:: FinalizeEquation
PUBLIC:: DefineParametersEquation 
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(     "IniExactFunc"  , " Specifies exactfunc to be used for initialization ")
CALL prms%CreateIntOption(     "IniRefState"   , " Specifies exactfunc to be used for initialization ")
CALL prms%CreateRealArrayOption("IniWaveNumber", " For exactfunction: wavenumber of solution.")
CALL prms%CreateRealArrayOption("IniCenter"    , " For exactfunction: center coordinates.","0.,0.,0.")
CALL prms%CreateRealOption(   "IniAmplitude"   , " For exactfunction: amplitude","0.1")
CALL prms%CreateRealOption(   "IniFrequency"   , " For exactfunction: Frequency","1.0")
CALL prms%CreateRealOption(   "IniHalfwidth"   , " For exactfunction: Halfwidth","0.1")
CALL prms%CreateRealOption(   "IniDisturbance" , " For exactfunction: Strength of initial disturbance","0.")
CALL prms%CreateRealOption(   "kappa"          , " Ratio of specific heats","1.6666666666666667")
CALL prms%CreateRealOption(     'R'            , " Gas constant.","287.058")
CALL prms%CreateRealOption(   "mu_0"           , " Magnetic permeability in vacuum","1.0")
#if PARABOLIC
CALL prms%CreateRealOption(   "eta"            , " Magnetic resistivity","0.")
CALL prms%CreateRealOption(   "mu"             , " Fluid viscosity","0.")
CALL prms%CreateRealOption(   "s23"            , " Stress tensor scaling (normally 2/3)","0.6666666666666667")
CALL prms%CreateRealOption(   "Pr"             , " Prandtl number","0.72")
#ifdef PP_ANISO_HEAT
CALL prms%CreateRealOption(   "kpar"  , " If anisotropic heat terms enabled: diffusion parallel to magnetic field","0.")
CALL prms%CreateRealOption(   "kperp" , " If anisotropic heat terms enabled: diffusion perpendicular to magnetic field","0.")
#endif /*PP_ANISO_HEAT*/
#endif /*PARABOLIC*/

#ifdef PP_GLM
CALL prms%CreateRealOption(   "GLM_scale", "MHD with GLM option: save ch from timestep <1","0.5")
CALL prms%CreateRealOption(   "GLM_scr", "MHD with GLM option: damping term of GLM variable 1/cr=5.555 (cr=0.18),"//&
                                         "set zero for no damping.","0.0")
CALL prms%CreateLogicalOption("DivBsource" , "Set true to add divB-error dependent source terms.",&
                                                 ".FALSE.")
#endif /*PP_GLM*/
CALL prms%CreateRealArrayOption(   "RefState", "primitive constant reference state, used for exactfunction/initialization" &
                                ,multiple=.TRUE.)
CALL prms%CreateIntOption(     "Riemann",  " Specifies Riemann solver:"//&
                                           "1: Lax-Friedrichs, "//&
                                           "2: HLLC, "//&
                                           "3: Roe, "//&
                                           "4: HLL, "//&
                                           "5: HLLD (only with mu_0=1), " )

#if (PP_DiscType==2)
CALL prms%CreateIntOption(     "VolumeFlux",  " Specifies the two-point flux to be used in the flux of the split-form "//&
                                              "DG volume integral "//&
                                              "0: Standard DG Flux"//&
                                              "1: standard DG Flux with metric dealiasing" &
                            ,"0")
#endif /*PP_DiscType==2*/
END SUBROUTINE DefineParametersEquation


!==================================================================================================================================
!> Get parameters needed by MHD equation modules and initialize equations
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools       ,ONLY:COUNTOPTION,GETINT,GETREAL,GETREALARRAY,GETINTARRAY,GETLOGICAL
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY:MeshInitIsDone,nBCSides,BC,BoundaryType
USE MOD_Equation_Vars 
USE MOD_Riemann
USE MOD_Flux_Average
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: MaxBCState,locType,locState,nRefState
!==================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.EquationInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "InitEquation not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MHD...'
doCalcSource=.TRUE.

IniRefState=-1
! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
IniRefState   = GETINT('IniRefState')
IniWavenumber = GETREALARRAY('IniWaveNumber',3,'1.,1.,1.')
IniCenter     = GETREALARRAY('IniCenter',3,'0.,0.,0.')
IniAmplitude  = GETREAL('IniAmplitude','0.1')
IniFrequency  = GETREAL('IniFrequency','1.0')
IniHalfwidth  = GETREAL('IniHalfwidth','0.1')
IniDisturbance= GETREAL('IniDisturbance','0.')



! Gas constants
Kappa    =GETREAL('kappa','1.4')
KappaM1  =Kappa-1.
KappaM2  =Kappa-2.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)
!permeability
mu_0    =GETREAL('mu_0','1.')
smu_0  = 1./(mu_0)
s2mu_0  =0.5*smu_0
#if PARABOLIC
!resistivity
eta    =GETREAL('eta','0.') 
etasmu_0 = eta*smu_0 
! Viscosity
mu        =GETREAL('mu','0.')
s23      = GETREAL('s23','0.6666666666666667')
Pr       =GETREAL('Pr','0.72')
R        =GETREAL('R','287.058')
KappasPr =Kappa/Pr
!heat diffusion coefficients 
#ifdef PP_ANISO_HEAT
  SWRITE(UNIT_StdOut,'(A)') "ANISOTROPIC HEAT COEFFICIENTS ENABLED."
  kperp = GETREAL('kperp','0.')
  kpar  = GETREAL('kpar','0.')
#endif /*PP_ANISO_HEAT*/
#endif /*PARABOLIC*/

#ifdef PP_GLM
GLM_scale = GETREAL('GLM_scale','0.5')
GLM_scr    = GETREAL('GLM_scr','0.0') !damping 1/cr, cr=0.18, 1/cr = 5.555
DivBSource = GETLOGICAL('DivBSource','.FALSE.')
!compute timestep for ch=1, then compute ch from timestep:
!dt ~ 1/ch -> dt/dtch1=1/ch -> ch=dtch1/dt
GLM_dtch1=0. !must be initialized to correct value in first call of calctimestep
#endif /*PP_GLM*/


! Read Boundary information / RefStates / perform sanity check
nRefState=COUNTOPTION('RefState')

! determine max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  IF((locType.NE.22).AND.(locType.NE.220)) MaxBCState = MAX(MaxBCState,locState) !BCType=22: special BC with different exactfuncs
END DO
#if MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCState,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
#endif /*MPI*/

! Sanity check for BCs
IF(IniRefState.GT.nRefState)&
  CALL abort(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))
IF(MaxBCState.GT.nRefState)&
  CALL abort(__STAMP__,&
    'ERROR: Boundary RefState not defined! (MaxBCState,nRefState):',MaxBCState,REAL(nRefState))

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefStatePrim(nRefState,PP_nVar))
  ALLOCATE(RefStateCons(nRefState,PP_nVar))
  DO i=1,nRefState
    RefStatePrim(i,:)  = GETREALARRAY('RefState',8)
#ifdef PP_GLM
    RefStatePrim(i,9)  =0.
#endif 
    CALL PrimToCons(RefStatePrim(i,:),RefStateCons(i,:))
    IF(RefStateCons(i,5).LT.0.)THEN
      CALL abort(__STAMP__, &
          "Refstate has negative energy",i,RefStateCons(i,5))
    END IF !neg. Energy
  END DO
END IF

WhichRiemannSolver = GETINT('Riemann','1')
SELECT CASE(WhichRiemannSolver)
CASE(1)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Lax-Friedrichs'
  SolveRiemannProblem => RiemannSolverByRusanov
CASE(2)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLC'
  SolveRiemannProblem => RiemannSolverByHLLC 
CASE(3)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Roe'
  SolveRiemannProblem => RiemannSolverByRoe
CASE(4)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLL'
  SolveRiemannProblem => RiemannSolverByHLL  
CASE(5)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLD'
  IF(ABS(mu_0-1.).GT.1.0E-12) &
   CALL abort(__STAMP__,&
  'HLLD solver only for mu_0=1 implemented!')
  SolveRiemannProblem => RiemannSolverByHLLD  
CASE DEFAULT
  CALL ABORT(__STAMP__,&
       "Riemann solver not implemented")
END SELECT

#if (PP_DiscType==2)
WhichVolumeFlux = GETINT('VolumeFlux','0')
SELECT CASE(WhichVolumeFlux)
CASE(0)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Standard DG'
  VolumeFluxAverageVec => StandardDGFluxVec
CASE(1)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Standard DG with Metrics Dealiasing'
  VolumeFluxAverageVec => StandardDGFluxDealiasedMetricVec
CASE DEFAULT
  CALL ABORT(__STAMP__,&
         "volume flux not implemented")
END SELECT
#endif /*PP_DiscType==2*/

IF(MPIroot) CALL CheckFluxes()

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MHD DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation


!==================================================================================================================================
!> Fill initial solution with IniExactFunc
!==================================================================================================================================
SUBROUTINE FillIni(IniExactFunc_in,U_in)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:Elem_xGP,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: IniExactFunc_in  !< handle to specify exactfunction
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_in(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< initialized DG solution 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================
! Determine Size of the Loops, i.e. the number of grid cells in the
! corresponding directions
DO iElem=1,nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        CALL ExactFunc(IniExactFunc_in,0.,Elem_xGP(1:3,i,j,k,iElem),U_in(1:PP_nVar,i,j,k,iElem))
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iElem=1,nElems
END SUBROUTINE FillIni



!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu) 
! MODULES
USE MOD_Globals,ONLY:Abort,CROSS
USE MOD_Preproc
USE MOD_Equation_Vars,ONLY:Kappa,sKappaM1,RefStateCons,RefStatePrim,IniRefState
USE MOD_Equation_Vars,ONLY:smu_0,mu_0
USE MOD_Equation_Vars,ONLY:IniCenter,IniFrequency,IniAmplitude,IniHalfwidth,IniWaveNumber
USE MOD_Equation_Vars,ONLY:IniDisturbance
USE MOD_Equation_Vars,ONLY:PrimToCons
USE MOD_TestCase_ExactFunc,ONLY: TestcaseExactFunc
USE MOD_TimeDisc_vars,ONLY:dt,CurrentStage,FullBoundaryOrder,RKc,RKb,t
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
REAL,INTENT(IN)                 :: x(3)             !< evaluation position
REAL,INTENT(IN)                 :: tIn              !< evaluation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar)      ! state in conservative variables
INTEGER                         :: i,j
REAL                            :: Omega,a
REAL                            :: Prim(1:PP_nVar) 
REAL                            :: r, e, nx,ny,sqr,va,phi_alv
REAL                            :: r2(1:16),Bphi,dp
REAL                            :: q0,q1,Lz
REAL                            :: B_R,r0,B_tor,PsiN,psi_a
REAL                            :: b0(3),xc(3)
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

! Determine the value, the first and the second time derivative
SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL TestcaseExactFunc(ExactFunction,tEval,x,Resu,Resu_t,Resu_tt)
CASE(0)        
  CALL TestcaseExactFunc(ExactFunction,tEval,x,Resu,Resu_t,Resu_tt)
CASE(1) ! constant
  Resu = RefStateCons(IniRefState,:) ! prim=(/1.,0.3,0.,0.,0.71428571/)
  !Resu(1)  = Prim(1)
  !Resu(2:4)= Resu(1)*Prim(2:4)
  !Resu(5)  = Prim(5)*sKappaM1 + 0.5*Resu(1)*SUM(Prim(2:4)*Prim(2:4))
CASE(2) ! non-divergence-free magnetic field,diss. Altmann
  Resu(1)=1.0
  Resu(2:4)=0.
  Resu(5)=6.0
  Resu(6)=IniAmplitude*EXP(-(SUM(((x(:)-IniCenter(:))/IniHalfwidth)**2)))
  Resu(7:PP_nVar)=0.
CASE(3) ! alfven wave 
  Omega=PP_Pi*IniFrequency
  ! r: lenght-variable = lenght of computational domain
  r=2.
  ! e: epsilon = 0.2
  e=0.2
  nx  = 1./SQRT(r**2+1.)
  ny  = r/SQRT(r**2+1.)
  sqr = 1. !2*SQRT(PP_Pi)
  Va  = omega/(ny*sqr)
  phi_alv = omega/ny*(nx*(x(1)-0.5*r) + ny*(x(2)-0.5*r)) - Va*tEval
  Resu=0.
  Resu(1) = 1.
  Resu(2) = -e*ny*COS(phi_alv)
  Resu(3) =  e*nx*COS(phi_alv)
  Resu(4) =  e*SIN(phi_alv)
  Resu(5) =  sKappaM1+0.5*SUM(Resu(2:4)*Resu(2:4)) !p=1, rho=1
  Resu(6) = nx -Resu(2)*sqr
  Resu(7) = ny -Resu(3)*sqr
  Resu(8) =    -Resu(4)*sqr
  ! g'(t)
  Resu_t     =0.
  Resu_t(1)  =0.
  Resu_t(2)  =-Va*e*ny*SIN(phi_alv) 
  Resu_t(3)  = Va*e*nx*SIN(phi_alv) 
  Resu_t(4)  =-Va*e*COS(phi_alv)
  Resu_t(5)  = SUM(Resu(2:4)*Resu_t(2:4))
  Resu_t(6)  = -Resu_t(2)*sqr
  Resu_t(7)  = -Resu_t(3)*sqr
  Resu_t(8)  = -Resu_t(4)*sqr
  Resu_tt    =0.
  Resu_tt(1) =0.
  Resu_tt(2) =-Va*Va*resu(2)
  Resu_tt(3) =-Va*Va*resu(3)
  Resu_tt(4) =-Va*Va*resu(4)
  Resu_tt(5) = SUM(Resu(2:4)*Resu_tt(2:4)+Resu_t(2:4)*Resu_t(2:4))
  Resu_tt(6) = -Resu_tt(2)*sqr
  Resu_tt(7) = -Resu_tt(3)*sqr
  Resu_tt(8) = -Resu_tt(4)*sqr

CASE(31,32,33) ! linear shear alfven wave , linearized MHD,|B|>=1 , p,rho from inirefstate 
         !IniWavenumber=(k_x,k_yk_z): k_parallel=k_x*e_x+k_y*e_y, k_perp=k_z*e_z
         !IniAmplitude should be small compare to density and pressure (1e-08)
         !31: backward moving, 32: forward moving, 33: standing (31+32!)
  b0(3)=0.
  IF(IniWavenumber(2).LT.0.01) THEN
    b0(1:2)=(/1.,0./)
  ELSEIF(IniWavenumber(1).LT.0.01) THEN
    b0(1:2)=(/0.,1./)
  ELSE
    b0(1:2)=(/MIN(1.,IniWavenumber(1)/IniWavenumber(2)),MIN(1.,IniWavenumber(2)/IniWavenumber(1)) /)
  END IF
  ASSOCIATE(rho_0=>RefStatePrim(IniRefState,1),p_0=>RefStatePrim(IniRefState,5))
  q0=SQRT(SUM(b0(:)**2))        !=|B_0|
  a=SQRT(mu_0*rho_0)  !=|B_0|/va = sqrt(mu_0*rho_0)
  IF(Exactfunction.EQ.32) a=-a ! case(32) -va!
  va=q0/a      !(+va)=|B_0|/sqrt(mu_0*rho_0)
  IF(Exactfunction.EQ.33)THEN
    r0=0. !switch for standing wave
  ELSE
    r0=1.
  END IF
  xc(1:3)=x(1:3)-r0*b0(1:3)/a*tEval ! b_0/a = B_0/|B_0|*va
  e=IniAmplitude*SIN(2.0*PP_Pi*SUM(xc(:)*IniWavenumber(:)))
  Prim=0.
  Prim(1)  =rho_0
  Prim(2:3)=(/-b0(2),b0(1)/)*(e/q0) !vperp
  Prim(5)  =p_0
  Prim(6:7)=b0(1:2)-r0*Prim(2:3)*a !-B0/(+va)=sqrt(mu0*rho_0) 
  CALL PrimToCons(Prim,Resu)
  !second time derivative
  Resu_tt     =0.
  e=e*(-va*2.0*PP_Pi*SUM(b0(:)*IniWavenumber(:)))**2
  Resu_tt(2:3)=(rho_0*e/q0)*(/-b0(2),b0(1)/)
  Resu_tt(6:7)=-Resu_tt(2:3)*r0*a
  !first time derivative
  Resu_t     =0.
  e=IniAmplitude*COS(2.0*PP_Pi*SUM(xc(:)*IniWavenumber(:)))*(-va*2.0*PP_Pi*SUM(b0(:)*IniWavenumber(:)))
  Resu_t(2:3)=(rho_0*e/q0)*(/-b0(2),b0(1)/)
  Resu_t(6:7)=-Resu_t(2:3)*r0*a

  END ASSOCIATE !rho_0,p_0
CASE(4) ! navierstokes exact function
  Omega=PP_Pi*IniFrequency
  a=RefStatePrim(IniRefState,2)*2.*PP_Pi

  ! g(t)
  Resu(1:4)=2.+ IniAmplitude*sin(Omega*SUM(x) - a*tEval)
  Resu(5)=Resu(1)*Resu(1)
  Resu(6:PP_nVar)=0.
  ! g'(t)
  Resu_t(1:4)=(-a)*IniAmplitude*cos(Omega*SUM(x) - a*tEval)
  Resu_t(5)=2.*Resu(1)*Resu_t(1)
  Resu_t(6:PP_nVar)=0.
  ! g''(t)
  Resu_tt(1:4)=-a*a*IniAmplitude*sin(Omega*SUM(x) - a*tEval)
  Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  Resu_tt(6:PP_nVar)=0.
CASE(5) ! mhd exact function (KAPPA=2., mu_0=1)
  IF(.NOT.((kappa.EQ.2.0).AND.(smu_0.EQ.1.0)))THEN
    CALL abort(__STAMP__,&
               'Exactfuntion 5 works only with kappa=2 and mu_0=1 !')
  END IF
  Omega=PP_Pi*IniFrequency
  ! g(t)
  Resu(1:3)         = 2. + IniAmplitude*SIN(Omega*(SUM(x) - tEval))
  Resu(4)           = 0.
  Resu(5)           = 2*Resu(1)*Resu(1)+resu(1)
  Resu(6)           = Resu(1)
  Resu(7)           =-Resu(1)
  Resu(8:PP_nVar)   = 0.
  ! g'(t)
  Resu_t(1:3)       = -IniAmplitude*omega*COS(Omega*(SUM(x) - tEval))
  Resu_t(4)         = 0.
  Resu_t(5)         = 4.*Resu(1)*Resu_t(1) +Resu_t(1)
  Resu_t(6)         = Resu_t(1)
  Resu_t(7)         =-Resu_t(1)
  Resu_t(8:PP_nVar) = 0.
  ! g''(t)
  Resu_tt(1:3)      =-IniAmplitude*omega*omega*SIN(Omega*(SUM(x) - tEval))
  Resu_tt(4)        = 0.
  Resu_tt(5)        = 4.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))+Resu_tt(1)
  Resu_tt(6)        = Resu_tt(1)
  Resu_tt(7)        =-Resu_tt(1)
  Resu_tt(8:PP_nVar)= 0.
CASE(6) ! case 5 rotated
  IF(.NOT.((kappa.EQ.2.0).AND.(smu_0.EQ.1.0)))THEN
    CALL abort(__STAMP__,&
               'Exactfuntion 5 works only with kappa=2 and mu_0=1 !')
  END IF
  Omega=PP_Pi*IniFrequency
  ! g(t)
  Resu              = 0.
  Resu(1:2)         = 2. + IniAmplitude*SIN(Omega*(SUM(x) - tEval))
  !Resu(3)           = 0.
  Resu(4)           = Resu(1)
  Resu(5)           = 2*Resu(1)*Resu(1)+resu(1)
  Resu(6)           = Resu(1)
  !Resu(7)           = 0.
  Resu(8)           =-Resu(1)
  ! g'(t)
  Resu_t            = 0.
  Resu_t(1:2)       = -IniAmplitude*omega*COS(Omega*(SUM(x) - tEval))
  !Resu_t(3)         = 0.
  Resu_t(4)         = Resu_t(1)
  Resu_t(5)         = 4.*Resu(1)*Resu_t(1) +Resu_t(1)
  Resu_t(6)         = Resu_t(1)
  !Resu_t(7)         = 0.
  Resu_t(8)         =-Resu_t(1)
  ! g''(t)
  Resu_tt           = 0.
  Resu_tt(1:2)      =-IniAmplitude*omega*omega*SIN(Omega*(SUM(x) - tEval))
  !Resu_tt(3)        = 0.
  Resu_tt(4)        = Resu_tt(1)
  Resu_tt(5)        = 4.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))+Resu_tt(1)
  Resu_tt(6)        = Resu_tt(1)
  !Resu_tt(7)        = 0.
  Resu_tt(8)        =-Resu_tt(1)
CASE(10) ! mhd exact equilibrium, from potential A=(0,0,A3), A3=IniAmplitude*PRODUCT(sin(omega*x(:)))
         !domain should be a cube [0,1]^2, boundary conditions either periodic of perfectly conducting wall
  Prim(:)= RefStatePrim(IniRefState,:)
  Omega=2*PP_Pi*IniFrequency
  a=SQRT(IniAmplitude*Prim(5))/omega !IniAmplitude is related to the change of pressure (IniAmplitude=0.1: 10% change) 
  Prim(6)= a*omega*SIN(Omega*x(1))*COS(Omega*x(2))
  Prim(7)=-a*omega*COS(Omega*x(1))*SIN(Omega*x(2))
  Prim(5)=Prim(5)*(1+ IniAmplitude*(SIN(Omega*x(1))*SIN(Omega*x(2)))**2) !a^2omega^2=p0*IniAmplitude

  CALL PrimToCons(Prim,Resu)
CASE(11) ! mhd exact equilibrium, Psi=a(x^2+y^2), B_x= dPsi/dy=2ay, B_y=-dPsi/dx=-2ax p=int(-Laplace(Psi),Psi)=-4a Psi
         ! domain |x|,|y|<1, Dirichlet BC needed
  Prim(:)= RefStatePrim(1,:)
  a      = SQRT(0.25*(0.5*IniAmplitude)*Prim(5)) !IniAmplitude is related to the change of pressure 
                                           ! at x,y=1 (x^2+y^2=2) (IniAmplitude=0.1: 10% change) 
  Prim(5)= Prim(5)-4*a*a*SUM(x(1:2)**2)
  Prim(6)= 2*a*x(2)
  Prim(7)=-2*a*x(1)
  Prim(:)= Prim(:)+RefStatePrim(2,:)*IniDisturbance*PRODUCT(SIN(2*PP_Pi*x(1:2)))

  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 11 !',999,prim(5))
  CALL PrimToCons(Prim,Resu)
CASE(12) ! mhd exact equilibrium, Psi=a*exp(-(x^2+y^2)/H^2), B_x= dPsi/dy=-2x/H^2*Psi, B_y=-dPsi/dx=2y/H^2*Psi
         ! p=int(-Laplace(Psi),Psi)=1/H^2(2*ln(Psi/a)+1)Psi^2
  Prim(:)= RefStatePrim(IniRefState,:)
  a      = SQRT(IniAmplitude*prim(5))*IniHalfwidth !IniAmplitude is related to the change of pressure (at Psi=Psi_max=a)
  psi_a  = a*EXP(-SUM(((x(1:2)-IniCenter(1:2))/IniHalfwidth)**2)) 
  Prim(5)= Prim(5)+(2.*LOG(Psi_a/a)+1.)*(Psi_a/IniHalfwidth)**2
  Prim(6)=-2*(x(2)-IniCenter(2))*Psi_a/IniHalfwidth**2
  Prim(7)= 2*(x(1)-IniCenter(1))*Psi_a/IniHalfwidth**2

  CALL PrimToCons(Prim,Resu)
CASE(13) ! 3D mhd exact equilibrium with constant pressure, 
         ! Potential A_1= (y^3*z-y*z^3), A_2=x^3*z-x*z^3, A_3= x^3*y-x*y^3
         ! B = [3*x*z^2-3*x*y^2, -3*y*z^2+2*y^3-3*x^2*y, 3*x^2*z-3*y^2*z]
         ! J=curl(B) = (0,0,0) -> dp=0, constant pressure
  Prim(:)= RefStatePrim(IniRefState,:)
  Prim(2:4)=0.
  Prim(6)=3.*x(1)*(x(3)**2-x(2)**2)
  Prim(7)=   x(2)*(2.*x(2)**2-3.*(x(3)**2+x(1)**2))
  Prim(8)=3.*x(3)*(x(1)**2-x(2)**2)

  CALL PrimToCons(Prim,Resu)
CASE(14) ! 3D mhd exact equilibrium with constant pressure & constant magnetic field
         ! Potential A = 1/2 [B_2*z-B_3*y ,B_3*x-B_1*z,, B_1*y-B_2*x] 
         ! B = [B_1,B_2,B_3] =const
         ! J=curl(B) = (0,0,0) -> dp=0, constant pressure
  Prim(:)= RefStatePrim(IniRefState,:)
  Prim(2:4)=0.
  Prim(6)=0.15
  Prim(7)=0.3
  Prim(8)=1.0

  CALL PrimToCons(Prim,Resu)

CASE(60) !TEST for MagneticEnergyModes computation
  Prim=0.
  Prim(1)=1.
  Prim(5)=1.
  Prim(6)=1.+0.1*COS(x(1)/3.) +0.2*COS(2./3.*x(1)) -0.3*COS(x(1)) &
            -0.1*SIN(x(1)/3.) -0.2*SIN(2./3.*x(1)) +0.3*SIN(x(1)) 
  CALL PrimToCons(Prim,Resu)

CASE(70) !Tearing mode instability, of paper Landi et al. , domain [0,6*pi]x[-pi/2,pi/2]x[0:2Pi]
        ! "Three-dimensional simulations of compressible tearing instability"
        ! rho_0=1, p0=0.5*beta (choose with refstate)
        ! Re_eta=5000, mu=0.,kappa=5/3 1/delta=0.1(=IniHalfwidth)  IniAmplitude=1.0E-04
  Prim=0.
  Prim(6)=TANH(x(2)/IniHalfwidth) !tanh(y*delta) delta=10.
  Prim(8)=SQRT(1-Prim(6)*Prim(6))
  
  Prim(1)=RefStatePrim(IniRefState,1) 
  Prim(3)=IniAmplitude*Prim(6)*Prim(8)*SIN(x(1)/3.*IniWaveNumber(1)+ x(3)*IniWaveNumber(3))
  Prim(5)=RefStatePrim(IniRefState,5)
  !Prim(6:8)=sSqrt4pi*Prim(6:8) ! scaling with sqrt(4pi)!?!
  CALL PrimToCons(Prim,Resu)
CASE(71) !Tearing mode instability, of paper Landi et al. , domain [0,6*pi]x[-pi/2,pi/2]x[0:2Pi]
        ! "Three-dimensional simulations of compressible tearing instability"
        ! rho_0=1, p0=0.5*beta (choose with refstate)
        ! Re_eta=5000, mu=0.,kappa=5/3 1/delta=0.1(=IniHalfwidth)  IniAmplitude=1.0E-04
  Prim=0.
  Prim(6)=TANH((x(2)-0.01)/IniHalfwidth) !tanh(y*delta) delta=10. !NOT FULLY CENTERED
  Prim(8)=SQRT(1-Prim(6)*Prim(6))
  
  Prim(1)=RefStatePrim(IniRefState,1) 
  DO j=0,NINT(IniWaveNumber(3))
    DO i=0,NINT(IniWaveNumber(1))
      a=REAL(1+0.8*i+0.9*j)/(1+0.8*IniWaveNumber(1)+0.9*IniWaveNumber(3))
      Prim(3)=Prim(3)+SIN(x(1)/3.*i+ x(3)*j+2*PP_Pi*a)
    END DO
  END DO
  Prim(3)=IniDisturbance*Prim(6)*Prim(8)*Prim(3)
  Prim(5)=RefStatePrim(IniRefState,5)
  !Prim(6:8)=sSqrt4pi*Prim(6:8) ! scaling with sqrt(4pi)!?!
  CALL PrimToCons(Prim,Resu)

CASE(75) !2D tearing mode instability, domain [0,1]x[0,4]
        ! from L.Chacon "A non-staggered, conservative, ∇ · B = 0, finite-volume scheme for 3D implicit extended MHD..",2004
        ! assuming rho_0=p_0=1
        ! eta=1.0E-02, mu=1.0E-03,kappa=5/3
  Prim=0.
  Prim(1)=1.+1.0E-03*COS(2*PP_Pi*x(1))*SIN(0.5*PP_Pi*x(2))
  Prim(5)=1.
  Prim(7)=1.-2./(1+EXP(2*5*(x(1)-0.5))) !tanh((x(1)-0.5)/lambda) lambda=0.2
  Prim(8)=SQRT(1-Prim(7)*Prim(7))
  CALL PrimToCons(Prim,Resu)

CASE(80) ! 2D island coalesence domain [-1,1]^2
        ! from L.Chacon "An optimal, parallel, fully implicit Newton-Krylov solver for 3D viscoresitive MHD",2008
        ! assuming rho_0=p_0=1
        ! eta=1.0E-02, mu=1.0E-03,kappa=5/3, lambda=2/(4*pi), eps=0.2
        ! Psi = -lambda*LN(COSH(x/lambda)+EPS*COS(y/lambda))
        ! deltaPsi = 1.0E-03*SIN(pi/2*x)*COS(pi*y)
        ! Bx = -Psi_y = -Eps*SIN(y/lambda) / (COSH(x/lambda)+EPS*COS(y/lambda)) 
        ! By =  Psi_x = -SINH(x/lambda) / (COSH(x/lambda)+EPS*COS(y/lambda)) 
        ! deltaBx = -deltaPsi_y = 1.0E-03*pi*SIN(pi/2*x)*SIN(pi*y) 
        ! deltaBy =  deltaPsi_x = 1.0E-03*pi/2*COS(pi/2*x)*COS(pi*y) 

        ! p  = p0 + 0.5*(1-eps^2)* (COSH(x/lambda)+EPS*COS(y/lambda))^-2
  nx=2*PP_Pi*x(1)
  ny=2*PP_Pi*x(2)
  r = 1./(COSH(nx)+0.2*COS(ny))
  Prim=0.
  Prim(1)=1.
  Prim(5)=1. + 0.5*(1-0.2*0.2)*r*r
  Prim(6)=-r*(0.2*SIN(ny)) + 1.0E-03*PP_Pi*SIN(0.25*nx)*SIN(0.5*ny)
  Prim(7)=-r*(SINH(nx))    + 1.0E-03*0.5*PP_Pi*COS(0.25*nx)*COS(0.5*ny)
  Prim(8)=0.
  CALL PrimToCons(Prim,Resu)

CASE(90) !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=j1*(1-r^4)-j2*(1-r^2)^8, j1=0.2, j2=0.266
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  Prim(:)=0.
  Prim(1)=1. 
  Prim(5)=0.025/9. ! =pmax, pmin=p(r=1)=pmax-0.0025 and pmax/pmin ~=10
  Prim(8)=SQRT(2*mu_0*Prim(5)*1.0E4) !beta=p/(0.5*smu0 B^2) =1.0E-04
  Prim(5)=Prim(5)+0.01  ! BETA CHANGED!!! ~ 1E-02  
  r2(1)=SUM(x(1:2)**2)
  r=SQRT(r2(1))
  DO i=2,16
    r2(i)=r2(i-1)*r2(1)
  END DO  
  Bphi=-mu_0*(133*r2(8)-1197*r2(7)+4788*r2(6)-11172*r2(5)+16758*r2(4)-16758*r2(3)+11472*r2(2)-4788*r2(1)+297)/9.0E+3 !*r
  dp  =-r2(1)*(      607086480.*r2(16)   -10965499545.*r2(15)   +93572262784.*r2(14)  -501279979200.*r2(13) &
                +1889439921600.*r2(12) -5321922445840.*r2(11)+11614288882560.*r2(10)-20096857548768.*r2( 9) &
               +27971578034400.*r2( 8)-31525903338930.*r2( 7)+28714629715200.*r2( 6)-20859299741280.*r2( 5) &
               +11747978409024.*r2( 4) -4854306857400.*r2( 3) +1285266977280.*r2( 2)  -138278780640.*r2( 1) &
               +5718295440.)/5.25096E+12
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  Prim(5)=Prim(5)+smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
              
  CALL PrimToCons(Prim,Resu)
CASE(91) !like 90, BUT WITH REAL BETA!!
         !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=j1*(1-r^4)-j2*(1-r^2)^8, j1=0.2, j2=0.266
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  Prim(:)=0.
  Prim(1)=1. 
  Prim(5)=0.025/9. ! =pmax, pmin=p(r=1)=pmax-0.0025 and pmax/pmin ~=10
  Prim(8)=SQRT(2*mu_0*Prim(5)*1.0E4) !beta=p/(0.5*smu0 B^2) =1.0E-04

  r2(1)=SUM(x(1:2)**2)
  r=SQRT(r2(1))
  DO i=2,16
    r2(i)=r2(i-1)*r2(1)
  END DO  
  Bphi=-mu_0*(133*r2(8)-1197*r2(7)+4788*r2(6)-11172*r2(5)+16758*r2(4)-16758*r2(3)+11472*r2(2)-4788*r2(1)+297)/9.0E+3 !*r
  dp  =-r2(1)*(      607086480.*r2(16)   -10965499545.*r2(15)   +93572262784.*r2(14)  -501279979200.*r2(13) &
                +1889439921600.*r2(12) -5321922445840.*r2(11)+11614288882560.*r2(10)-20096857548768.*r2( 9) &
               +27971578034400.*r2( 8)-31525903338930.*r2( 7)+28714629715200.*r2( 6)-20859299741280.*r2( 5) &
               +11747978409024.*r2( 4) -4854306857400.*r2( 3) +1285266977280.*r2( 2)  -138278780640.*r2( 1) &
               +5718295440.)/5.25096E+12
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  Prim(5)=Prim(5)+smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
  !disturbance of the velocity in xy plane  (scaled with (1-r^2)^2 and IniDisturbance) 
  Prim(2)=IniDisturbance*(1.-r2(1))**2
  Prim(3)=x(2)*Prim(2)
              
  CALL PrimToCons(Prim,Resu)
CASE(92) !cylindrical equilibrium for ideal MHD for internal kink (Jorek paper Huysmans),q profile given, Bz=1
         ! cylindrical domain r[0,1], z[0,100]
         ! q(r)=0.73*(1-r^2)+1.6*r^2 => Bphi(r)=2*pi*r*Bz/(L0*q(r))
         ! pressure difference from gradp=J x B:
         ! dp(r) = smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
         ! dp(r) = smu_0*2*pi^2*r^2/Lz^2*((q1-q0)*r^2+2*q0)/(q0*(((q1-q0)*r^2+2*q0)*(q1-q0)*r^2 + q0^2))
         ! p=p0-dp, p0=dp(1)/0.98, dp(1) = smu_0*2*pi^2/Lz^2*(q0+q1)/(q0*q1^2), so that pmax/pmin = 50
         ! small perturbation u_z=1.0E-12*r^2*(1-r^2)*cos(2*pi*z/Lz)
         !   in the first mode of the z velocity z \in[0,100] !! 
  r2(1)=SUM(x(1:2)**2)
  q0=0.73  !q(r=0)
  q1=1.60 !q(r=1)
  Lz=100.
  Prim(:)=0.
  Prim(1)=1.-0.9*r2(1) 
  Prim(4)=IniDisturbance*r2(1)*(1-r2(1))*SIN(2*PP_Pi*x(3)/Lz)
  Prim(5)=smu_0*2.*PP_Pi*PP_Pi*(q0+q1)/(0.98*Lz*Lz*q0*q1*q1)
  Prim(8)=1.  !Bz=1.
  Bphi= 2*PP_Pi/(Lz*((q1-q0)*r2(1)+q0)) !*r !L_0=100 !
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  dp  = 2*PP_Pi*PP_Pi*r2(1)*((q1-q0)*r2(1)+2*q0)/(Lz*Lz*q0*(((q1-q0)*r2(1)+2*q0)*(q1-q0)*r2(1)+q0**2))

  Prim(5)=Prim(5)-smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
              
  CALL PrimToCons(Prim,Resu)
CASE(910) ! like 91, but with constant density:
         ! q(r)=1/100*(q1*(1-r^2)+q0*r^2) => Bphi(r)=2*pi*r*Bz/(L0*q(r))
  r2(1)=SUM(x(1:2)**2)
  q0=0.73  !q(r=0)
  q1=1.60 !q(r=1)
  Lz=100.
  Prim(:)=0.
  Prim(1)=1.
  Prim(4)=IniDisturbance*r2(1)*(1-r2(1))*SIN(2*PP_Pi*x(3)/Lz)
  Prim(5)=smu_0*2.*PP_Pi*PP_Pi*(q0+q1)/(0.98*Lz*Lz*q0*q1*q1)
  Prim(8)=1.  !Bz=1.
  Bphi= 2*PP_Pi/(Lz*((q1-q0)*r2(1)+q0)) !*r !L_0=100 !
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  dp  = 2*PP_Pi*PP_Pi*r2(1)*((q1-q0)*r2(1)+2*q0)/(Lz*Lz*q0*(((q1-q0)*r2(1)+2*q0)*(q1-q0)*r2(1)+q0**2))

  Prim(5)=Prim(5)-smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
              
  CALL PrimToCons(Prim,Resu)
CASE(911) ! like 91, but with a different q profile 0.82 <= q <= 1.4:
  r2(1)=SUM(x(1:2)**2)
  q0=0.82  !q(r=0)
  q1=1.40 !q(r=1)
  Lz=100.
  Prim(:)=0.
  Prim(1)=1.-0.9*r2(1) 
  Prim(4)=IniDisturbance*r2(1)*(1-r2(1))*SIN(2*PP_Pi*x(3)/Lz)
  Prim(5)=smu_0*2.*PP_Pi*PP_Pi*(q0+q1)/(0.98*Lz*Lz*q0*q1*q1)
  Prim(8)=1.  !Bz=1.
  Bphi= 2*PP_Pi/(Lz*((q1-q0)*r2(1)+q0)) !*r !L_0=100 !
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  dp  = 2*PP_Pi*PP_Pi*r2(1)*((q1-q0)*r2(1)+2*q0)/(Lz*Lz*q0*(((q1-q0)*r2(1)+2*q0)*(q1-q0)*r2(1)+q0**2))

  Prim(5)=Prim(5)-smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
              
  CALL PrimToCons(Prim,Resu)
CASE(100) !tearing mode in a torus around z axis and center (0,0,0), small radius a=1 and large radius R=10
          ! desinty constant, pressure (PsiN)
  R0=10.
  R  = SQRT(x(1)**2+x(2)**2)
  a  = SQRT((R-R0)**2 +x(3)**2)
  psiN= (-0.062829+(1.84896-0.786131*a)*a)*a  !psi=0 for a=0, psi=1 for a=1
  prim=0.
  prim(1)=1.
  Prim(5)=0.0024142343*(0.0497 - 0.0393 * psiN) !rho*T=rho_n*T_n*rho0/(kB*n0)=0.24e-02*rho_n*T_n
  !psi=(psiN-1)*0.20266, dpsi/da
  psi_a= 0.20266*(-0.062829+(2*1.84896-3*0.786131*a)*a) 
  !B_R=-1/R*(dPsi/dZ) = -1/R*(dPsi/da)*(da/dZ) = -Z/(R*a)*(dPsi/da)
  !B_Z=1/R*(dPsi/dR) = 1/R*(dPsi/da)*(da/dR) = (R-R0)/(R*a)*(dPsi/da)
  IF(a.LT.1.0E-12)THEN 
    B_R=-psi_a/R
    Prim(8)=psi_a/R
  ELSE
    B_R=-x(3)/(R*a)*psi_a
    Prim(8)=(R-R0)/(R*a)*psi_a
  END IF
  ! F0=10. B_tor=F0/R
  B_tor=10./R !10/R
  Prim(6)= ( x(2)*B_tor + x(1)*B_R)/R
  Prim(7)= (-x(1)*B_tor + x(2)*B_R)/R
  
  CALL PrimToCons(Prim,Resu)
CASE(101) !internal kink from Jorek in a torus around z axis and center (0,0,0), small radius a=1 and large radius R=10
          ! not using R=R0 as magnetic axis, but the one provided by jorek!! 
  R0=10.
  R  = SQRT(x(1)**2+x(2)**2)
  r2(1)= (R-R0)**2 +x(3)**2
  e  = 0.0310976053 ! dm/(1-dm^2), dm= distance to magnetic axis in xdirection (magnetic axis from Jorek: R0+dm=R0 + 0.03106759)
  !e  = 0.030769957662 ! dm= 0.03074088
  a  = SQRT(((R-R0)-e*(1-r2(1)))**2 +x(3)**2) !correction to idistance from magnetic axis
  !psiN= (1.27024127385*a**2+0.0327486740462*a**3-0.531773141036*a**4+0.128783193*a**5)*a**2/0.9 (from density fit, using a)
  psiN= (1.27024127385+(0.0327486740462+(-0.531773141036+0.128783193*a)*a)*a)*a*a/0.9
  prim=0.
  prim(1)=1.-0.9*psiN
  Prim(5)=2.0e-03*(1-1.7*psiN+0.72*PsiN**2) !rho*T=rho_n*T_n*rho0/(kB*n0)=0.24e-02*rho_n*T_n
  !psi=(psiN-1)*0.4792, dpsi/da, non-normalized psi_min=-0.4792
  psi_a= 0.4792*(2*1.27024127385*a+3*0.0327486740462*a**2-4*0.531773141036*a**3+5*0.128783193*a**4)/0.9
  !psi_a= 0.4792*(2*1.27024127385+(3*0.0327486740462+(-4*0.531773141036+5*0.128783193*a)*a)*a)*a/0.9
  !B_R=1/R*(dPsi/dZ) = 1/R*(dPsi/da)*(da/dZ) 
  !B_Z=-1/R*(dPsi/dR) = -1/R*(dPsi/da)*(da/dR) 
  IF(a.LT.1.0E-12)THEN 
    B_R=psi_a/R
    Prim(8)=-psi_a/R
  ELSE
    B_R=psi_a/(R*a)*  x(3)*(1+2*e*((R-R0)-e*(1-r2(1))))
    Prim(8)=-psi_a/(R*a)* (1+2*e*(R-R0))*((R-R0)-e*(1-r2(1)))
  END IF
  ! F0=10. B_tor=F/R, F=(F0^2-4*0.4792*(Psi_n-Psi_n^2))^0.5
  B_tor=SQRT(100.-4*0.4792*PsiN*(1-0.5*PsiN))/R 
  Prim(6)= ( x(2)*B_tor + x(1)*B_R)/R
  Prim(7)= (-x(1)*B_tor + x(2)*B_R)/R
  
  CALL PrimToCons(Prim,Resu)
CASE(311) ! Orzsag-Tang vortex
  prim    = 0.
  prim(1) = 1.
  prim(2) = -SIN(2.*PP_Pi*x(2))
  prim(3) =  SIN(2.*PP_Pi*x(1))
  prim(5) =  1./kappa
  prim(6) = -SIN(2.*PP_Pi*x(2))/kappa
  prim(7) =  SIN(4.*PP_Pi*x(1))/kappa
  CALL PrimToCons(Prim,Resu)
CASE(24601) ! Alternative insulating Taylor-Green vortex (A) from Brachet et al. Derivation of the pressure initial condition
            ! is found in the appendix of Bohm et al. Constant chosen such that initial Mach number is 0.1
  prim    =  0.
!  r       =  1. so we don't need it, but include it here to see the similarities between these test cases
!
  prim(1) =  1.
  prim(2) =  SIN(x(1))*COS(x(2))*COS(x(3))
  prim(3) = -COS(x(1))*SIN(x(2))*COS(x(3))
  prim(5) =  100./kappa + 0.0625*(COS(2.*x(1))+COS(2.*x(2)))*(2.+COS(2.*x(3))) + &
             0.0625*(COS(4.*x(1))+COS(4.*x(2)))*(2.-COS(4.*x(3)))
  prim(6) =  COS(2.*x(1))*SIN(2.*x(2))*SIN(2.*x(3))
  prim(7) = -SIN(2.*x(1))*COS(2.*x(2))*SIN(2.*x(3))
  CALL PrimToCons(Prim,Resu)
CASE(24602) ! Original insulating Taylor-Green vortex (I) from Brachet et al. Constant chosen such that initial Mach number is 0.1
  prim    =  0.
  r       =  1./SQRT(3.)
!
  prim(1) =  1.
  prim(2) =  SIN(x(1))*COS(x(2))*COS(x(3))
  prim(3) = -COS(x(1))*SIN(x(2))*COS(x(3))
  prim(5) =  100./kappa + 0.0625*(COS(2.*x(1))+COS(2.*x(2)))*(2.+COS(2.*x(3))) + &
             r*0.0625*(COS(2.*x(2))+COS(2.*x(3)))*(2.-COS(2.*x(1))) + &
             r*0.0625*(COS(2.*x(1))+COS(2.*x(3)))*(2.-COS(2.*x(2))) + &
             r*0.25*COS(2.*x(3)) - r*0.125*COS(2.*x(1))*COS(2.*x(2))
  prim(6) =  r*COS(x(1))*SIN(x(2))*SIN(x(3))
  prim(7) =  r*SIN(x(1))*COS(x(2))*SIN(x(3))
  prim(8) = -r*2.*SIN(x(1))*SIN(x(2))*COS(x(3))
  CALL PrimToCons(Prim,Resu)
CASE(24603) ! Conductive Taylor-Green vortex (C) from Brachet et al. Constant chosen such that initial Mach number is 0.1
  prim    =  0.
  r       =  1./SQRT(3.) ! Brachet et al. proposed SQRT(2./3.) but that doesn't make sense to me as the initial kinetic
                         ! and magnetic energies don't match
!
  prim(1) =  1.
  prim(2) =  SIN(x(1))*COS(x(2))*COS(x(3))
  prim(3) = -COS(x(1))*SIN(x(2))*COS(x(3))
  prim(5) =  100./kappa + 0.0625*(COS(2.*x(1))+COS(2.*x(2)))*(2.+COS(2.*x(3))) - &
             r*0.0625*(COS(4.*x(2))+COS(4.*x(3)))*(2.+COS(4.*x(1))) - &
             r*0.0625*(COS(4.*x(1))+COS(4.*x(3)))*(2.+COS(4.*x(2))) - &
             r*0.25*COS(4.*x(3)) - r*0.125*COS(4.*x(1))*COS(4.*x(2))
  prim(6) =  r*SIN(2.*x(1))*COS(2.*x(2))*COS(2.*x(3))
  prim(7) =  r*COS(2.*x(1))*SIN(2.*x(2))*COS(2.*x(3))
  prim(8) = -r*2.*COS(2.*x(1))*COS(2.*x(2))*SIN(2.*x(3))
CALL PrimToCons(Prim,Resu)
END SELECT ! ExactFunction

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
  SELECT CASE(CurrentStage)
  CASE(1)
    ! resu = g(t)
  CASE(2)
    ! resu = g(t) + dt/3*g'(t)
    Resu=Resu + dt*RKc(2)*Resu_t
  CASE(3)
    ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
    Resu=Resu + RKc(3)*dt*Resu_t + RKc(2)*RKb(2)*dt*dt*Resu_tt
  CASE DEFAULT
    ! Stop, works only for 3 Stage O3 LS RK
    CALL abort(__STAMP__,&
               'Exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT !CurrentStage
END IF !fullBoundaryOrder

END SUBROUTINE ExactFunc


!==================================================================================================================================
!> Compute the source terms, called from dg timedrivtaive
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,tIn) 
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc !PP_N
USE MOD_Equation_Vars, ONLY: IniExactFunc,IniFrequency,IniAmplitude
USE MOD_Equation_Vars,ONLY:RefStatePrim,IniRefState
USE MOD_Equation_Vars, ONLY:Kappa,KappaM1
USE MOD_Equation_Vars, ONLY:doCalcSource
USE MOD_Mesh_Vars,     ONLY:Elem_xGP,nElems
#if PARABOLIC
USE MOD_Equation_Vars, ONLY:mu,Pr,eta
#endif
#ifdef PP_GLM
USE MOD_Equation_Vars, ONLY:GLM_scr,GLM_ch,smu_0
USE MOD_Equation_Vars, ONLY:DivBSource
USE MOD_DG_Vars,       ONLY:U
#endif /*PP_GLM*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: tIn        !< evaluation time
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !<time derivative, where source is added 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Omega,a
REAL                            :: Ut_src(PP_nVar)
INTEGER                         :: i,j,k,iElem
REAL                            :: sinXGP,sinXGP2,cosXGP,at
REAL                            :: tmp(6)
REAL                            :: rho,rho_x,rho_xx
#ifdef PP_GLM
REAL                            :: v(3),divB,sGLM_ch
#endif /*PP_GLM*/
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(4) ! navierstokes exact function
  Omega=PP_Pi*IniFrequency
  a=RefStatePrim(IniRefState,2)*2.*PP_Pi
  tmp(1)=-a+3*Omega
  tmp(2)=-a+0.5*Omega*(1.+kappa*5.)
  tmp(3)=IniAmplitude*Omega*KappaM1  
  tmp(4)=0.5*((9.+Kappa*15.)*Omega-8.*a) 
  tmp(5)=IniAmplitude*(3.*Omega*Kappa-a)
#if PARABOLIC
  tmp(6)=3.*mu*Kappa*Omega*Omega/Pr
#else
  tmp(6)=0.
#endif
  tmp=tmp*IniAmplitude
  at=a*tIn
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      cosXGP=COS(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP=SIN(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP2=2.*sinXGP*cosXGP !=SIN(2.*(omega*SUM(Elem_xGP(:,i,j,k,iElem))-a*tIn))
      Ut_src=0.
      Ut_src(1  )       =  tmp(1)*cosXGP
      Ut_src(2:4)       =  tmp(2)*cosXGP + tmp(3)*sinXGP2
      Ut_src(5  )       =  tmp(4)*cosXGP + tmp(5)*sinXGP2 + tmp(6)*sinXGP
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src
    END DO; END DO; END DO ! i,j,k
    !Ut(:,:,:,:,iElem) = Ut(:,:,:,:,iElem)+resu*Amplitude !Original
  END DO ! iElem
CASE(5) ! mhd exact function, KAPPA==2!!!
  Omega=PP_Pi*IniFrequency
#if PARABOLIC
  tmp(1)=6.*mu/Pr
#endif
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      rho    = IniAmplitude*SIN(Omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      rho_x  = IniAmplitude*Omega*COS(Omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      rho    = rho + 2.
      Ut_src = 0.
!
      Ut_src(1  )  =  rho_x
      Ut_src(2:3)  =  rho_x +  4.*rho*rho_x
      Ut_src(4  )  =  4.*rho*rho_x
      Ut_src(5  )  =  rho_x + 12.*rho*rho_x
      Ut_src(6  )  =  rho_x
      Ut_src(7  )  = -rho_x
#if PARABOLIC
      rho_xx       = -Omega*Omega*(rho - 2.)
      Ut_src(5  )  =  Ut_src(5) - tmp(1)*rho_xx - 6.*eta*(rho_x*rho_x+rho*rho_xx)
      Ut_src(6  )  =  Ut_src(6) - 3.*eta*rho_xx
      Ut_src(7  )  =  Ut_src(7) + 3.*eta*rho_xx
#endif /*PARABOLIC*/
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE(6) ! case 5 rotated
  Omega=PP_Pi*IniFrequency
#if PARABOLIC
  tmp(1)=6*mu/Pr
#endif
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      rho    = IniAmplitude*SIN(omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      rho_x  = IniAmplitude*omega*COS(omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      rho_xx =-omega*omega*rho
      rho    = rho+2.
      Ut_src=0.
      Ut_src(1)   =  rho_x 
      Ut_src(2)   =  rho_x +  4.*rho*rho_x 
      Ut_src(3)   =        +  4.*rho*rho_x
      Ut_src(4)   =  rho_x +  4.*rho*rho_x 
      Ut_src(5)   =  rho_x + 12.*rho*rho_x
      Ut_src(6)   =  rho_x 
      Ut_src(8)   =  -rho_x
#if PARABOLIC
      Ut_src(5)   = Ut_src(5) - tmp(1)*rho_xx-6.*eta*(rho_x*rho_x+rho*rho_xx)
      Ut_src(6)   = Ut_src(6) - 3*eta*rho_xx
      Ut_src(8)   = Ut_src(8) + 3*eta*rho_xx
#endif /*PARABOLIC*/
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE DEFAULT
  ! No source -> do nothing
  doCalcSource=.FALSE. 
END SELECT ! ExactFunction

!#ifndef PP_GLM
!CASE DEFAULT
!  ! No source -> do nothing
!  doCalcSource=.FALSE. 
!#endif /*PP_GLM*/
!END SELECT ! ExactFunction
!
!#ifdef PP_GLM
!IF(DivBSource)THEN
!  sGLM_ch=1./GLM_ch
!  DO iElem=1,nElems
!    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
!      v=U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem)
!      ! Ut(9)/ch = -divB
!      divB=-(sGLM_ch*Ut(9,i,j,k,iElem))
!      ! to (rhov)_t : -1/mu_0*B*divB  (from JxB = -1/mu_0* [ div(1/2|B|^2-BB) + B div(B) ] )
!      ! to B_t      : -v*divB
!      ! to energy v.(rhov)_t +1/mu_0(B.B_t) = -1/mu_0* divB ((v*B) + (B*v))
!      ! or, regarding entropy conservation, u*(source_rhou) != source_totE
!      ! so that  energy source is  -1/mu_0* divB (v*B)  
!      Ut_src(2:4) = -smu_0*divB*U(6:8,i,j,k,iElem)  !=-1/mu_0 *divB * B
!      Ut_src(5)   =  SUM(Ut_src(2:4)*v(:))          !=-1/mu_0 *divB * (B.v) 
!      Ut_src(6:8) = -divB*v(:)                      !=-divB * v
!      Ut(2:8,i,j,k,iElem) = Ut(2:8,i,j,k,iElem) +Ut_src(2:8)
!    END DO; END DO; END DO ! i,j,k
!  END DO ! iElem
!END IF !divBsource
!Ut(9,:,:,:,:)=Ut(9,:,:,:,:)-GLM_scr*U(9,:,:,:,:)
!#endif /*PP_GLM*/
END SUBROUTINE CalcSource

!==================================================================================================================================
!> Finalize /deallocate all module variables 
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone,RefStatePrim,RefStateCons
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(RefStatePrim)
SDEALLOCATE(RefStateCons)
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation


!==================================================================================================================================
!> Check Riemann Solver  for consistency, two-point flux for consistency and symmetry
!==================================================================================================================================
SUBROUTINE CheckFluxes()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,ONLY:nAuxVar,PrimToCons 
USE MOD_Flux,         ONLY: EvalAdvectionFlux1D
USE MOD_Flux_Average
USE MOD_Riemann
#if (PP_DiscType==2)
USE MOD_Flux_Average , ONLY: standardDGFluxVec
USE MOD_DG_Vars,       ONLY: DGinitIsDone,nTotal_vol,U
USE MOD_Mesh_Vars,     ONLY: Metrics_fTilde
#endif /*PP_DiscType==2*/
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY:GLM_ch
#endif /*PP_GLM*/
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: PL(PP_nVar),UL(PP_nVar),FrefL(PP_nVar),FrefR(PP_nVar),Fcheck(PP_nVar)
REAL  :: PR(PP_nVar),UR(PP_nVar),Frefsym(PP_nVar)
REAL  :: check,absdiff
INTEGER :: icase,i
PROCEDURE(),POINTER :: fluxProc 
CHARACTER(LEN=255)  :: fluxName
#if PP_DiscType==2
REAL  :: UauxElem(   nAuxVar,0:PP_N,0:PP_N,0:PP_N)
REAL  :: ftildeElem( PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL  :: gtildeElem( PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL  :: htildeElem( PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL  :: metricL(3),metricR(3),mtmp(3),Utmp(PP_nVar)
REAL  :: ULaux(1:nAuxVar),URaux(1:nAuxVar)
LOGICAL :: failed_vol
#endif /*PP_DiscType==2*/
LOGICAL :: failed
!==================================================================================================================================
WRITE(*,*)'    CHECK ALL FLUXES...'
!test consistency, random state:
PL(1:8)=(/1.12794,0.103391,-0.04153,0.097639,74.3605,0.15142,-3.1415,0.5673/)
PR(1:8)=(/0.94325,-0.21058,-0.14351,-0.20958,52.3465,0.32217,-2.0958,-0.243/)
#ifdef PP_GLM
PL(9)= 0.31999469
PR(9)= 0.
GLM_ch = 0.5 !was not set yet
#endif
CALL PrimToCons(PL,UL)
CALL PrimToCons(PR,UR)
!use EvalAdvectionFlux1D as reference Flux
CALL EvalAdvectionFlux1D(UL,FrefL)
CALL EvalAdvectionFlux1D(UR,FrefR)
failed=.FALSE.
DO icase=0,6
  NULLIFY(fluxProc)
  SELECT CASE(icase)
  CASE(0)
    fluxProc => StandardDGFlux
    fluxName = "StandardDGFlux"
  CASE(1)
    fluxProc => RiemannSolverByHLL
    fluxName = "RiemannSolverByHLL"
  CASE(2)
    fluxProc => RiemannSolverByRoe
    fluxName = "RiemannSolverByRoe"
  CASE(3)
    fluxProc => RiemannSolverByHLLC
    fluxName = "RiemannSolverByHLLC"
  CASE(4)
    fluxProc => RiemannSolverByHLLD
    fluxName = "RiemannSolverByHLLD"
!  CASE(5)
!    fluxProc => EntropyAndKinEnergyConservingFlux
!    fluxName = "EntropyAndKinEnergyConservingFlux"
!  CASE(6)
!    fluxProc => EntropyStableFlux
!    fluxName = "EntropyStableFlux"
  CASE DEFAULT
    CYCLE
  END SELECT
  !CONSISTENCY
  CALL fluxProc(UL,UL,Fcheck)
  check=1.0e-12
  DO i=1,PP_nVar
    absdiff=ABS(FrefL(i)-Fcheck(i))
    IF(absdiff.GT.1.0e-12)THEN
      WRITE(*,*)'FrefL /=Fcheck:',i,FrefL(i),Fcheck(i)
      check=max(check,absdiff)
    END IF
  END DO
  CALL fluxProc(UR,UR,Fcheck)
  DO i=1,PP_nVar
    absdiff=ABS(FrefR(i)-Fcheck(i))
    IF(absdiff.GT.1.0e-12)THEN
      WRITE(*,*)'FrefR /=Fcheck:',i,FrefR(i),Fcheck(i)
      check=max(check,absdiff)
    END IF
  END DO
  IF(check.GT.1.0e-12)THEN
    WRITE(*,*) "consistency check for solver "//TRIM(fluxName)//" failed",icase,check
    failed=.TRUE.
  END IF
END DO !icase
#if PP_DiscType==2
#ifdef CARTESIANFLUX
metricL=(/1.5320,0.,0./)
metricR=(/1.5320,0.,0./)
#else
metricL=(/1.5320,-0.05,4.895/)
metricR=(/0.8715,0.594,2.531/)
#endif

!use EvalEulerFluxTilde3D at point (0,0,0) as reference Flux
IF(DGinitIsDone)THEN
  Utmp(:)=U(:,0,0,0,1) !save U
ELSE
  ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1)) !DGinit not yet called!
  nTotal_vol=(PP_N+1)**3
  DO i=1,PP_nVar
    U(i,:,:,:,1)=UL(i)
  END DO
END IF
mtmp(:)=Metrics_ftilde(:,0,0,0,1) !save metric
!overwrite

U(:,0,0,0,1)=UL
Metrics_ftilde(:,0,0,0,1)=metricL
CALL EvalEulerFluxTilde3D(1,ftildeElem,gtildeElem,htildeElem,UauxElem)
ULaux=UauxElem(:,0,0,0)
FrefL = fTildeElem(:,0,0,0)

U(:,0,0,0,1)=UR
Metrics_ftilde(:,0,0,0,1)=metricR
CALL EvalEulerFluxTilde3D(1,ftildeElem,gtildeElem,htildeElem,UauxElem)
URaux=UauxElem(:,0,0,0)
FrefR = fTildeElem(:,0,0,0)

Metrics_fTilde(:,0,0,0,1)=mtmp(:) !put metric back
IF(DGinitIsDone)THEN
  U(:,0,0,0,1)=Utmp(:) !put back U
ELSE
  DEALLOCATE(U)
END IF
failed_vol=.FALSE.
DO icase=0,2
  NULLIFY(fluxProc)
  SELECT CASE(icase)
  CASE(0)
    fluxProc => StandardDGFluxVec
    fluxName = "StandardDGFluxVec"
  CASE(1)
    fluxProc => StandardDGFluxDealiasedMetricVec
    fluxName = "StandardDGFluxDealiasedMetricVec"
!  CASE(2)
!    fluxProc => EntropyandKinEnergyConservingFluxVec
!    fluxName = "EntropyandKinEnergyConservingFluxVec"
  CASE DEFAULT
    CYCLE
  END SELECT
  !CONSISTENCY
  CALL fluxProc(   UL,UL,ULaux,ULaux,metricL  & 
#ifndef CARTESIANFLUX
                                    ,metricL  &
#endif          
                                    ,Fcheck)
  check=1.0e-12
  DO i=1,PP_nVar
    absdiff=ABS(FrefL(i)-Fcheck(i))
    IF(absdiff.GT.1.0e-12)THEN
      WRITE(*,*)'FrefL /=Fcheck:',i,FrefL(i),Fcheck(i)
      check=max(check,absdiff)
    END IF
  END DO
  CALL fluxProc(   UR,UR,URaux,URaux,metricR  & 
#ifndef CARTESIANFLUX
                                    ,metricR  &
#endif          
                                    ,Fcheck)
  DO i=1,PP_nVar
    absdiff=ABS(FrefR(i)-Fcheck(i))
    IF(absdiff.GT.1.0e-12)THEN
      WRITE(*,*)'FrefR /=Fcheck:',i,FrefR(i),Fcheck(i)
      check=max(check,absdiff)
    END IF
  END DO
  IF(check.GT.1.0e-12)THEN
    WRITE(*,*)"consistency check for volume flux "//TRIM(fluxName)//" failed",icase,check
    failed_vol=.TRUE.
  END IF
  !SYMMETRY
  CALL fluxProc(   UL,UR,ULaux,URaux,metricL  & 
#ifndef CARTESIANFLUX
                                    ,metricR  &
#endif
                                    ,Frefsym)
  CALL fluxProc(   UR,UL,URaux,ULaux,metricR  & 
#ifndef CARTESIANFLUX
                                    ,metricL  &
#endif
                                    ,Fcheck)
  check=1.0e-12
  DO i=1,PP_nVar
    absdiff=ABS(Frefsym(i)-Fcheck(i))
    IF(absdiff.GT.1.0e-12)THEN
      WRITE(*,*)'Frefsym /=Fcheck:',i,Frefsym(i),Fcheck(i)
      check=max(check,absdiff)
    END IF
  END DO
  IF(check.GT.1.0e-12)THEN
    WRITE(*,*) "symmetry check for solver "//TRIM(fluxName)//" failed",icase,check
    failed_vol=.TRUE.
  END IF
END DO !iCase
IF(failed_vol) THEN
  CALL ABORT(__STAMP__, &
     "consistency/symmetry check for average volume flux failed")
END IF !failed
#endif /*PP_DiscType==2*/
IF(failed) THEN
  CALL ABORT(__STAMP__, &
     "consistency check for riemann solver failed")
END IF !failed
WRITE(*,*)'    ...SUCESSFULL.'

END SUBROUTINE CheckFluxes

END MODULE MOD_Equation
