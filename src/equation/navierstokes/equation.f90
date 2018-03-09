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
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
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

PUBLIC:: DefineParametersEquation
PUBLIC:: InitEquation
PUBLIC:: FillIni
PUBLIC:: ExactFunc
PUBLIC:: CalcSource
PUBLIC:: FinalizeEquation
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
CALL prms%CreateIntOption(     'IniExactFunc' , " Specifies exactfunc to be used for initialization ")
CALL prms%CreateIntOption(      'IniRefState' , "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState'    , "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE.)
CALL prms%CreateRealArrayOption('AdvVel'      , "for exact function, const velocity.")
CALL prms%CreateRealArrayOption('MachShock'   , "for exact function, Mach shock.")
CALL prms%CreateRealArrayOption('PreShockDens', "for exact function, pre shock density.")
CALL prms%CreateRealArrayOption('IniCenter'   , "for exactfunc, center point.","0.,0.,0.")
CALL prms%CreateRealArrayOption('IniAxis'     , "for exactfunc, center axis.","0.,0.,1.")
CALL prms%CreateRealOption(     'IniFrequency', "for exactfunc, frequency.","1.")
CALL prms%CreateRealOption(     'IniAmplitude', "for exactfunc, Amplitude.","0.1")
CALL prms%CreateRealOption(     'IniHalfwidth', "for exactfunc, Halfwidth.","0.2")
CALL prms%CreateRealOption(     'Kappa'       , "ratio of specific heats.","1.4")
CALL prms%CreateRealOption(     'R'           , "gas constant.","287.058")
CALL prms%CreateRealOption(     'Pr'          , "Prandtl number.","0.72")
#if PP_VISC == 0
CALL prms%CreateRealOption(     'mu0'         , "constant viscosity.","0.")
#elif PP_VISC == 1
! mu-Sutherland
CALL prms%CreateRealOption(     'mu0'         , "sutherland viscosity, prefactor.","1.735e-5")
CALL prms%CreateRealOption(     'Ts'          , "sutherland viscosity, temperature coeff.","110.4")
CALL prms%CreateRealOption(     'Tref'        , "sutherland viscosity, reference temperature.","280.")
CALL prms%CreateRealOption(     'ExpoSuth'    , "sutherland viscosity, exponent.","1.5")
#elif PP_VISC == 2
! mu power-law
CALL prms%CreateRealOption(     'mu0'         , "power-law viscosity, prefactor.","0.")
CALL prms%CreateRealOption(     'Tref'        , "power-law viscosity, reference temperature.","280.")
CALL prms%CreateRealOption(     'ExpoPow'     , "power-law viscosity, exponent.","1.5")
#endif /*PP_VISC==2*/

CALL prms%CreateIntOption(     "Riemann",  " Specifies the riemann flux to be used:"//&
                                           " 1: Lax-Friedrichs"//&
                                           " 2: HLLC"//&
                                           " 3: Roe"//&
                                           " 4: Entropy Stable"//&
                                           " 5: Entropy Conserving"//&
                                           " 6: Kennedy & Gruber"//&
                                           " 7: Ducros"//&
                                           " 8: Morinishi"//&
                                           " 9: EC-KEP"//&
                                           "10: approx. EC-KEP"//&
                                           "11: EC-KEP + press. aver."//&
                                           "12: Kennedy & Gruber (Pirozolli version)"//&
                                           "13: Entropy conservative Ismali and Roe with LLF diss"//&
                                           "14: Kennedy Gruber with LLF diss"//&
                                           "15: Decros Flux with LLF diss"//&
                                           "16: EC+KEP with LLF diss"//&
                                           "17: EC+KEP - pressure with LLF diss"//&
                                           "18: Kennedy Gruber (Pirozilli version)  with LLF diss"//&
                                           "19: Morinishi + LLF diss"//&
                                           "20: Gassner, Winters, Walch flux"//&
                                           "21: Gassner, Winters, Walch flux + LLF diss" &
                                          ,"1")
CALL prms%CreateIntOption(     "VolumeFlux",  " Specifies the two-point flux to be used in split-form flux or Riemann:"//&
                                              "DG volume integral "//&
                                              " 0: Standard DG"//&
                                              " 1: Standard DG with metirc dealiasing"//&
                                              " 2: Kennedy-Gruber"//&
                                              " 3: Ducros"//&
                                              " 4: Morinishi"//&
                                              " 5: EC-KEP"//&
                                              " 6: approx. EC-KEP"//&
                                              " 7: approx. EC-KEP + press. aver."//&
                                              " 8: Kenndy & Gruber (pirozolli version)"//&
                                              " 9: Gassner Winter Walch"//&
                                              "10: Two-Point EC" &
                                             ,"0")
END SUBROUTINE DefineParametersEquation



!==================================================================================================================================
!> initialize equation specific parameters 
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
#if PP_VISC == 2
REAL    :: Tref
#endif /*PP_VISC==2*/
INTEGER :: i,iSide
INTEGER :: MaxBCState,locType,locState,nRefState
!==================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.EquationInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "InitEquation not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT NAVIER-STOKES...'
doCalcSource=.TRUE.

s23=2./3.

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
IniRefState  = 0
SELECT CASE (IniExactFunc)
CASE(1,11,12)
  IniRefState  =GETINT('IniRefState')
CASE(2,21,8) ! synthetic test cases
  AdvVel = GETREALARRAY('AdvVel',3)
CASE(6) ! shock
  MachShock    = GETREAL('MachShock','1.5')
  PreShockDens = GETREAL('PreShockDens','1.0')
END SELECT ! IniExactFunc
IniCenter    = GETREALARRAY('IniCenter',3,'0.,0.,0.')
IniAxis      = GETREALARRAY('IniAxis',3,'0.,0.,1.')
IniAxis      = IniAxis/SQRT(SUM(IniAxis*IniAxis)) !Normalize
IniFrequency = GETREAL('IniFrequency','1.0')
IniAmplitude = GETREAL('IniAmplitude','0.1')
IniHalfwidth = GETREAL('IniHalfwidth','0.2')



! Gas constants
Kappa    =GETREAL('kappa','1.4')
sKappa    =1./Kappa
KappaM1  =Kappa-1.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)
R=GETREAL('R','287.058')
#if PARABOLIC
Pr       =GETREAL('Pr','0.72')
KappasPr =Kappa/Pr

! Viscosity
#if PP_VISC == 0
  mu0=GETREAL('mu0','0.')
#endif /*PP_VISC==0*/
! mu-Sutherland
#if PP_VISC == 1
mu0     =GETREAL('mu0','1.735E-5')
Ts      =GETREAL('Ts','110.4')
sTref   =1./GETREAL('Tref','280.')
ExpoSuth=GETREAL('ExpoSuth','1.5')
Ts      =Ts*sTref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#endif /*PP_VISC==1*/
! mu power-law
#if PP_VISC == 2
mu0=GETREAL('mu0','0.')
Tref    =GETREAL('Tref')
ExpoPow =GETREAL('ExpoPow')
mu0     =mu0/(Tref**ExpoPow)
#endif /*PP_VISC==2*/

#endif /*PARABOLIC*/

! Read Boundary information / RefStates / perform sanity check
nRefState=COUNTOPTION('RefState')

! determine max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  IF((locType.NE.22).AND.(locType.NE.220)) MaxBCState = MAX(MaxBCState,locState) !BCType=22: special BC with different exactfuncs
  IF((locType.EQ.4) .AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No temperature (refstate) defined for BC_TYPE',locType)
  IF((locType.EQ.10).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No pressure (refstate) defined for BC_TYPE',locType)
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
  ALLOCATE(RefStatePrim(nRefState,5))
  ALLOCATE(RefStateCons(nRefState,5))
  DO i=1,nRefState
    RefStatePrim(i,:)  = GETREALARRAY('RefState',5)
    RefStateCons(i,1)  = RefStatePrim(i,1) !cant use primtocons yet
    RefStateCons(i,2:4)= RefStatePrim(i,2:4)*RefStatePrim(i,1)
    RefStateCons(i,5)  = sKappaM1*RefStatePrim(i,5)+0.5*SUM(RefStateCons(i,2:4)*RefStatePrim(i,2:4))
  END DO
END IF

WhichRiemannSolver = GETINT('Riemann','1')
SELECT CASE(WhichRiemannSolver)
  CASE(1)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Lax-Friedrichs'
    SolveRiemannProblem     => RiemannSolverByRusanov
  CASE(2)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLC'
    SolveRiemannProblem     => RiemannSolverByHLLC
  CASE(3)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Roe'
    SolveRiemannProblem     => RiemannSolverByRoe
  CASE(4)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Entropy Stable'
    SolveRiemannProblem     => RiemannSolver_EntropyStable
  CASE(5)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Entropy Conserving'
    VolumeFluxAverage    => TwoPointEntropyConservingFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(6)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Kennedy & Gruber'
    VolumeFluxAverage    => KennedyAndGruberFlux1
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(7)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Ducros'
    VolumeFluxAverage    => DucrosFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(8)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Morinishi'
    VolumeFluxAverage    => MorinishiFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(9)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: EC-KEP'
    VolumeFluxAverage    => EntropyAndEnergyConservingFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(10)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: approx. EC-KEP'
    VolumeFluxAverage    => EntropyAndEnergyConservingFlux2
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(11)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: EC-KEP + press. aver.'
    VolumeFluxAverage    => ggFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(12)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Kennedy & Gruber (Pirozolli version)'
    VolumeFluxAverage    => KennedyAndGruberFlux2
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(13)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Entropy conservative Ismali and Roe with LLF diss'
    VolumeFluxAverage    => TwoPointEntropyConservingFlux
    SolveRiemannProblem     => RiemannSolver_EC_LLF
  CASE(14)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Kennedy Gruber with LLF diss'
    VolumeFluxAverage    => KennedyAndGruberFlux1
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE(15)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Decros Flux with LLF diss'
    VolumeFluxAverage    => DucrosFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE(16)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: EC+KEP with LLF diss'
    VolumeFluxAverage    => EntropyAndEnergyConservingFlux
    SolveRiemannProblem     => RiemannSolver_ECKEP_LLF
  CASE(17)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: EC+KEP - pressure with LLF diss'
    VolumeFluxAverage    => ggFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE(18)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Kennedy Gruber (Pirozilli version)  with LLF diss'
    VolumeFluxAverage    => KennedyAndGruberFlux2
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE(19)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Morinishi + LLF diss'
    VolumeFluxAverage    => MorinishiFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE(20)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Gassner, Winters, Walch flux'
    VolumeFluxAverage    => GassnerWintersWalchFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage
  CASE(21)
    SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Gassner, Winters, Walch flux + LLF diss'
    VolumeFluxAverage    => GassnerWintersWalchFlux
    SolveRiemannProblem     => RiemannSolver_VolumeFluxAverage_LLF
  CASE DEFAULT
    CALL ABORT(__STAMP__,&
         "Riemann solver not implemented")
END SELECT

WhichVolumeFlux = GETINT('VolumeFlux','0')
SELECT CASE(WhichVolumeFlux)
CASE(0)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Standard DG'
  VolumeFluxAverageVec => StandardDGFluxVec
CASE(1)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Standard DG with metirc dealiasing'
  VolumeFluxAverageVec => StandardDGFluxDealiasedMetricVec
CASE(2)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Kennedy-Gruber'
  VolumeFluxAverageVec => KennedyAndGruberFluxVec1
CASE(3)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Ducros'
  VolumeFluxAverageVec => DucrosFluxVec
CASE(4)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Morinishi'
  VolumeFluxAverageVec => MorinishiFluxVec
CASE(5)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: EC-KEP'
  VolumeFluxAverageVec => EntropyAndEnergyConservingFluxVec
CASE(6)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: approx. EC-KEP'
  VolumeFluxAverageVec => EntropyAndEnergyConservingFluxVec2
CASE(7)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: approx. EC-KEP + press. aver.'
  VolumeFluxAverageVec => ggFluxVec
CASE(8)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Kenndy & Gruber (pirozolli version)'
  VolumeFluxAverageVec => KennedyAndGruberFluxVec2
CASE(9)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Gassner Winter Walch'
  VolumeFluxAverageVec => GassnerWintersWalchFluxVec
CASE(10)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Two-Point EC'
  VolumeFluxAverageVec => TwoPointEntropyConservingFluxVec
CASE DEFAULT
  CALL ABORT(__STAMP__,&
         "volume flux not implemented")
END SELECT
EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT NAVIER-STOKES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation


!==================================================================================================================================
!> fill the initial DG solution with a given exactfunction
!==================================================================================================================================
SUBROUTINE FillIni(IniExactFunc_in,U_in)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:Elem_xGP,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: IniExactFunc_in                                !< Exactfunction to be used
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_in(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< Input state
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
!> Collection of analytical function, can represent exact solutions. input is x and t and the exactfunction integer.
!> The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu) 
! MODULES
USE MOD_Preproc
USE MOD_Globals,ONLY:Abort,CROSS
USE MOD_Equation_Vars,ONLY:Kappa,sKappaM1,KappaM1,KappaP1,MachShock,PreShockDens,AdvVel,RefStateCons,RefStatePrim,IniRefState
USE MOD_Equation_Vars,ONLY:IniCenter,IniFrequency,IniHalfwidth,IniAmplitude,IniAxis
USE MOD_Equation_Vars,ONLY:PrimToCons
USE MOD_TimeDisc_Vars,ONLY:dt,CurrentStage,FullBoundaryOrder,RKc,RKb,t
USE MOD_TestCase_ExactFunc,ONLY: TestcaseExactFunc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)                 :: tIn              !< evaluation time 
REAL,INTENT(IN)                 :: x(3)             !< x,y,z position
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: tEval
REAL                            :: Resu_t(5),Resu_tt(5)          ! state in conservative variables
REAL                            :: Frequency,Amplitude
REAL                            :: Omega
REAL                            :: Vel(3),Cent(3),a
REAL                            :: Prim(5)
REAL                            :: r_len,e,nx,ny,sqr,Va,Phi_alv
REAL                            :: Ms,xs
REAL                            :: Resul(5),Resur(5)
REAL                            :: random
REAL                            :: du, dTemp, RT, r2    ! aux var for SHU VORTEX,isentropic vortex case 17,20 
  !                                                  ! vel. perturbation, adiabatic, gasconstant
! needed for the SHU VORTEX 2D case 20
! needed for PAO -spectrum ONLY:
REAL                            :: factors(256)
REAL                            :: sines(256)
INTEGER                         :: f
! needed for cylinder potential flow
REAL                            :: phi
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration
resu_t=0.
resu_tt=0.
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
CASE(2) ! sinus
  Frequency=0.01
  Amplitude=0.3
  Omega=PP_Pi*Frequency
  ! base flow
  prim(1)   = 1.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)*(1.+Amplitude*sin(Omega*SUM(cent(1:3))))
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e 
  ! g'(t)
  Resu_t(1)=-Amplitude*cos(Omega*SUM(cent(1:3)))*Omega*SUM(vel)
  Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
  Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
  ! g''(t)
  Resu_tt(1)=-Amplitude*sin(Omega*SUM(cent(1:3)))*Omega*SUM(vel)*Omega*SUM(vel)
  Resu_tt(2:4)=Resu_tt(1)*prim(2:4) 
  Resu_tt(5)=0.5*Resu_tt(1)*SUM(prim(2:4)*prim(2:4))
CASE(21) ! linear in rho
  ! base flow
  prim(1)   = 100.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)+SUM(cent)
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e
  ! g'(t)
  Resu_t(1)=-SUM(Vel)
  Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
  Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
  ! g''(t)
  Resu_tt=0.
CASE(3) !Alfen Wave without Magnetic Field
  Omega=PP_Pi*IniFrequency
  ! r_len: lenght-variable = lenght of computational domain
  r_len=2.
  ! e: epsilon = 0.2
  e=0.2
  nx  = 1./SQRT(r_len**2+1.)
  ny  = r_len/SQRT(r_len**2+1.)
  sqr = 1.
  Va  = 1./sqr
  Phi_alv = Omega/ny*(nx*(x(1)-0.5*r_len) + ny*(x(2)-0.5*r_len) - Va*tEval)
  ! g(t)
  Resu(1) = 1.
  Resu(2) = 0.*nx - e*ny*COS(Phi_alv)
  Resu(3) = 0.*ny + e*nx*COS(Phi_alv)
  Resu(4) = e*SIN(Phi_alv)
  Prim(5) = 1.
  ! note, for rho=1...resu(2:4) = prim(2:4) !!!
  Resu(5)=prim(5)*sKappaM1+0.5*Resu(1)*SUM(Resu(2:4)*Resu(2:4))
  ! g'(t) 
  Resu_t(1) = 0.
  Resu_t(2) = -e*ny*SIN(Phi_alv)*Omega/ny*Va
  Resu_t(3) =  e*ny*SIN(Phi_alv)*Omega/ny*Va
  Resu_t(4) = -e*COS(Phi_alv)*Omega/ny*Va
  Resu_t(5) = 0.5*Resu(1)*2.*SUM(Resu(2:4)*Resu_t(2:4))
  ! g''(t)
  Resu_tt(1) = 0.
  Resu_tt(2) =  e*ny*COS(Phi_alv)*Omega/ny*Va*Omega/ny*Va
  Resu_tt(3) = -e*ny*COS(Phi_alv)*Omega/ny*Va*Omega/ny*Va
  Resu_tt(4) = -e*SIN(Phi_alv)*Omega/ny*Va*Va*Omega/ny*Va
  Resu_tt(5) = 0.5*Resu(1)*2.*SUM(Resu(2:4)*Resu_tt(2:4) + Resu_t(2:4)*Resu_t(2:4))
CASE(4) ! exact function
  Omega=PP_Pi*IniFrequency

  ! g(t)
  Resu(1:4)=2.+ IniAmplitude*sin(Omega*(SUM(x) - tEval))
  Resu(5)=Resu(1)*Resu(1)
  ! g'(t)
  Resu_t(1:4)=(-omega)*IniAmplitude*cos(Omega*(SUM(x) - tEval))
  Resu_t(5)=2.*Resu(1)*Resu_t(1)
  ! g''(t)
  Resu_tt(1:4)=-omega*omega*IniAmplitude*sin(Omega*(SUM(x) - tEval))
  Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
CASE(5)
  Resu(2:4) = 1.+x(1)+2.*x(2)-7.*x(3) 
  !!Resu(2:4) = 1.+x(1)**2+2.*x(2)**2*x(1)*x(3)-5.*x(3)**3*x(2)**2*x(1)
  Resu(1)=10.
  Resu(5)=10000.
CASE(6) ! shock
  prim=0.

  ! pre-shock
  prim(1) = PreShockDens
  Ms      = MachShock

  prim(5)=prim(1)/Kappa
  CALL PrimToCons(prim,Resur)

  ! post-shock
  prim(3)=prim(1) ! temporal storage of pre-shock density
  prim(1)=prim(1)*((KappaP1)*Ms*Ms)/(KappaM1*Ms*Ms+2.)
  prim(5)=prim(5)*(2.*Kappa*Ms*Ms-KappaM1)/(KappaP1)
  IF (prim(2) .EQ. 0.0) THEN
    prim(2)=Ms*(1.-prim(3)/prim(1))
  ELSE
    prim(2)=prim(2)*prim(3)/prim(1)
  END IF
  prim(3)=0. ! reset temporal storage
  CALL PrimToCons(prim,Resul)
  xs=5.+Ms*tEval ! 5. bei 10x10x10 Rechengebiet
  ! Tanh boundary
  Resu=-0.5*(Resul-Resur)*TANH(5.0*(x(1)-xs))+Resur+0.5*(Resul-Resur)
CASE(7) !TAYLOR GREEN VORTEX
  A=1. ! magnitude of speed
  Ms=0.1  ! maximum Mach number
  prim(1)=1.
  prim(2)= A*SIN(x(1))*COS(x(2))*COS(x(3))
  prim(3)=-A*COS(x(1))*SIN(x(2))*COS(x(3))
  prim(4)=0.
  prim(5)=(A/Ms*A/Ms/Kappa*prim(1))  ! scaling to get Ms
  prim(5)=prim(5)+1./16.*A*A*prim(1)*(COS(2*x(1))*COS(2.*x(3)) + 2.*COS(2.*x(2)) +2.*COS(2.*x(1)) +COS(2*x(2))*COS(2.*x(3)))
  CALL PrimToCons(prim,Resu)
CASE(8) !1D Sinus
  Omega=PP_Pi
  cent(1)=x(1)-AdvVel(1)*tEval
  resu(1)=1.+0.25*SIN(Omega*cent(1))
  resu(2)= resu(1)*AdvVel(1)
  resu(3:4)=0.
  resu(5)=sKappaM1*1.+0.5*Advvel(1)**2*resu(1)
  resu_t(1)=-0.25*Omega*Advvel(1)*COS(Omega*cent(1))
  resu_t(2)=AdvVel(1)*resu_t(1)
  resu_t(3:4)=0.
  resu_t(5)=0.5*AdvVel(1)**2*resu_t(1)
  resu_tt(1)=-0.25*(Omega*Advvel(1))**2*SIN(Omega*cent(1))
  resu_tt(2)=AdvVel(1)*resu_tt(1)
  resu_tt(3:4)=0.
  resu_tt(5)=0.5*AdvVel(1)**2*resu_tt(1)
CASE(9)
  Omega=2.*PP_Pi
  DO f=1,256
    factors(f)=f**(-5./3.)*exp(-1.5*(f*0.02)**(4./3.))
    factors(f)=sqrt(2.*factors(f))
    sines(f) = SIN(f*Omega*x(2))
  END DO
  prim(3) = SUM(factors(1:256)*sines(1:256))
  CALL PrimToCons(prim,Resu)
CASE(10) !Roundjet Bogey Bailly 2002, Re=65000, x-axis is jet axis
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  ! Jet inflow (from x=0, diameter 2.0)
  ! Initial jet radius: rj=1.
  ! Momentum thickness: delta_theta0=0.05=1/20
  ! Re=65000
  ! Uco=0.
  ! Uj=0.9
  r_len=SQRT((x(2)*x(2)+x(3)*x(3)))
  prim(2)=0.9*0.5*(1.+TANH((1.-r_len)*10.))
  CALL RANDOM_NUMBER(random)
  ! Random disturbance +-5%
  random=0.05*2.*(random-0.5)
  prim(2)=prim(2)+random*prim(2)
  prim(3)=x(2)/r_len*0.5*random*prim(2)
  prim(4)=x(3)/r_len*0.5*random*prim(2)
  CALL PrimToCons(prim,ResuL)
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  CALL PrimToCons(prim,ResuR)   
!   after x=10 blend to ResuR
  Resu=ResuL+(ResuR-ResuL)*0.5*(1.+tanh(x(1)-10.))
CASE(11)  ! Cylinder flow
  IF(tEval .EQ. 0.)THEN   ! Initialize potential flow
    prim(1)=RefStatePrim(IniRefState,1)  ! Density
    prim(4)=0.           ! VelocityZ=0. (2D flow)
    ! Calculate cylinder coordinates (0<phi<Pi/2)
    phi=ATAN2(x(2),x(1))
    ! Calculate radius**2
    r2=x(1)*x(1)+x(2)*x(2)
    ! Calculate velocities, radius of cylinder=0.5
    prim(2)=RefStatePrim(IniRefState,2)*(COS(phi)**2*(1.-0.25/r2)+SIN(phi)**2*(1.+0.25/r2))
    prim(3)=RefStatePrim(IniRefState,2)*(-2.)*SIN(phi)*COS(phi)*0.25/r2
    ! Calculate pressure, RefState(2)=u_infinity
    prim(5)=RefStatePrim(IniRefState,5) + &
            0.5*prim(1)*(RefStatePrim(IniRefState,2)*RefStatePrim(IniRefState,2)-prim(2)*prim(2)-prim(3)*prim(3))
  ELSE  ! Use RefState as BC
    prim=RefStatePrim(IniRefState,:)
  END IF  ! t=0
  CALL PrimToCons(prim,resu)
CASE(12) ! SHU VORTEX,isentropic vortex (adapted from HALO)
  ! base flow
  prim=RefStatePrim(IniRefState,:)  ! Density
  ! ini-Parameter of the Example
  vel=prim(2:4)
  RT=prim(PP_nVar)/prim(1) !ideal gas
  cent=(iniCenter+vel*tEval)!centerpoint time dependant
  cent=x-cent ! distance to centerpoint
  cent=cent-IniAxis*SUM(IniAxis*cent)
  cent=cent/iniHalfWidth !Halfwidth is dimension 1
  r2=SUM(cent*cent) !
  du = iniAmplitude/(2.*PP_Pi)*exp(0.5*(1.-r2)) ! vel. perturbation
  dTemp = -kappaM1/(2.*kappa*RT)*du**2 ! adiabatic
  prim(1)=prim(1)*(1.+dTemp)**(1.*skappaM1) !rho
  prim(2:4)=prim(2:4)+du*CROSS(IniAxis,cent) !v
  prim(PP_nVar)=prim(PP_nVar)*(1.+dTemp)**(kappa/kappaM1) !p
  CALL PrimToCons(prim,resu)
CASE(13) ! Sedov-Taylor Circular Blast Wave
  prim(1)   = 1.     ! ambient density
  prim(2:4) = 0.     ! gas at rest initially
  prim(5)   = 0.00001! ambient pressure
  r2 = SQRT(SUM(x*x))! the radius
  IF ((r2.LE.0.1).AND.(r2.NE.0.)) THEN
    du      = 4.*PP_Pi*r2*r2*r2/3. ! the volume of the small sphere
    prim(5) = kappaM1/du ! inject energy into small radius sphere, p = (gamma-1)*E/V
  END IF
  CALL PrimToCons(prim,resu)
END SELECT ! ExactFunction

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
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,tIn)
! MODULES
USE MOD_PreProc
USE MOD_Globals,ONLY:Abort
USE MOD_Equation_Vars, ONLY: IniExactFunc,IniFrequency,IniAmplitude
USE MOD_Equation_Vars, ONLY: Kappa,KappaM1,AdvVel,doCalcSource
USE MOD_Mesh_Vars,     ONLY:Elem_xGP,nElems
#if PARABOLIC
USE MOD_Equation_Vars,ONLY:mu0,Pr
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: tIn              !< evaluation time 
REAL,INTENT(INOUT)              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Omega
INTEGER                         :: i,j,k,iElem
REAL                            :: Ut_src(5)
REAL                            :: sinXGP,sinXGP2,cosXGP
REAL                            :: tmp(6)
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(4) ! exact function
  Omega=PP_Pi*IniFrequency
  tmp(1)=-Omega+3*Omega
  tmp(2)=-Omega+0.5*Omega*(1.+kappa*5.)
  tmp(3)=IniAmplitude*Omega*KappaM1
  tmp(4)=0.5*((9.+Kappa*15.)*Omega-8.*Omega)
  tmp(5)=IniAmplitude*(3.*Omega*Kappa-Omega)
#if PARABOLIC
  tmp(6)=3.*mu0*Kappa*Omega*Omega/Pr
#else
  tmp(6)=0.
#endif
  tmp=tmp*IniAmplitude
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      cosXGP=COS(omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      sinXGP=SIN(omega*(SUM(Elem_xGP(:,i,j,k,iElem))-tIn))
      sinXGP2=2.*sinXGP*cosXGP !=SIN(2.*(omega*SUM(Elem_xGP(:,i,j,k,iElem))-omega*t))
      Ut_src(1)   = tmp(1)*cosXGP
      Ut_src(2:4) = tmp(2)*cosXGP + tmp(3)*sinXGP2
      Ut_src(5)   = tmp(4)*cosXGP + tmp(5)*sinXGP2 + tmp(6)*sinXGP
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE DEFAULT
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource


!==================================================================================================================================
!> Deallocate Equation Variables 
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

END MODULE MOD_Equation
