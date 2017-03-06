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

!===================================================================================================================================
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source 
!===================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
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
CALL prms%CreateRealArrayOption('AdvVel',       "Advection velocity for advection part of LinAdv-Diff.")
CALL prms%CreateRealOption(     'DiffC',        "Diffusion constant for diffusion part of LinAdv-Diff.")
CALL prms%CreateRealArrayOption('IniWaveNumber'," Wave numbers used for exactfunction in linadv.")
CALL prms%CreateIntOption(     'IniExactFunc',  " Specifies exactfunc to be used for initialization ")
#if (PP_DiscType==2)
CALL prms%CreateIntOption(     'VolumeFlux',  " Specifies the two-point flux to be used in the flux of the split-form DG volume integral ")
#endif /*PP_DiscType==2*/
END SUBROUTINE DefineParametersEquation


!==================================================================================================================================
!> initialize equation specific parameters 
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,ONLY:GETREALARRAY,GETREAL,GETINT
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
USE MOD_Equation_Vars
USE MOD_Flux_Average,ONLY: standardDGFluxVec
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
   SWRITE(*,*) "InitLinearScalarAdvection not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SCALAR LINADV...'

! Read the velocity vector from ini file
AdvVel             = GETREALARRAY('AdvVel',3)
! Read the diffusion constant from ini file
DiffC             = GETREAL('DiffC','0.')
IniWavenumber     = GETREALARRAY('IniWavenumber',3,'1.,1.,1.')

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')

#if (PP_DiscType==2)
WhichVolumeFlux = GETINT('VolumeFlux','0')
SELECT CASE(WhichVolumeFlux)
CASE(0)
  SWRITE(UNIT_stdOut,'(A)') 'Flux Average Volume: Standard DG'
  VolumeFluxAverageVec => StandardDGFluxVec
CASE DEFAULT
  CALL ABORT(__STAMP__,&
         "volume flux not implemented")
END SELECT
#endif /*PP_DiscType==2*/

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LINADV DONE!'
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
SUBROUTINE ExactFunc(ExactFunction,tin,x,resu) 
! MODULES
USE MOD_Globals,ONLY:Abort,MPIRoot
USE MOD_Preproc,ONLY:PP_Pi
USE MOD_Equation_Vars,ONLY:AdvVel,DiffC
USE MOD_Equation_Vars,ONLY: IniWavenumber 
USE MOD_TimeDisc_vars,ONLY:dt,CurrentStage,FullBoundaryOrder,t,RKc,RKb
USE MOD_TestCase_ExactFunc,ONLY: TestcaseExactFunc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: tin              !< evaluation time 
REAL,INTENT(IN)                 :: x(3)             !< x,y,z position
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! state in conservative variables
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3)
REAL                            :: x0(3)
REAL                            :: r,theta
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Cent=x-AdvVel*t
SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL TestcaseExactFunc(ExactFunction,tEval,x,Resu,Resu_t,Resu_tt)
CASE(0)        
  CALL TestcaseExactFunc(ExactFunction,tEval,x,Resu,Resu_t,Resu_tt)
CASE(1) !linear
    Resu=0.+SUM(Cent)
    Resu_t=-SUM(Advvel)
    Resu_tt=0.
CASE(21) !linear
    Resu=10.+SUM(Cent)
    Resu_t=-SUM(Advvel)
    Resu_tt=0.
CASE(101) !constant
    Resu=3.5
    Resu_t=0.
    Resu_tt=0.
CASE(102) !linear steadystate in xy
    Resu=SUM(x(1:2))
    Resu_t=0.
    Resu_tt=0.
CASE(103) !radial steadystate in xy
    Resu=SUM(x(1:2)**2)**2
    Resu_t=0.
    Resu_tt=0.
CASE(104) !radial steadystate in xy
    r=SQRT(SUM(x(1:2)**2))
    theta=ATAN2(x(2),x(1))/PP_Pi
    Resu=r**6*(0.1*(3*theta-theta**3))**2
    Resu_t=0.
    Resu_tt=0.
CASE(105) !radial steadystate in xy
    r=SQRT(SUM(x(1:2)**2))
    theta=ATAN2(x(2),x(1))/PP_Pi
    Resu=(1-r**2)*(1.+(r*(3*theta-theta**3))**2)
    Resu_t=0.
    Resu_tt=0.
CASE(106) !sinus steadystate in xy
    r=SQRT(SUM(x(1:2)**2))
    theta=ATAN2(x(2),x(1))
    Resu=(SIN(2*PP_Pi*r)-0.5*SIN(4*PP_Pi*r))*(COS(2*theta)-SIN(4*(theta-0.4)))
    Resu_t=0.
    Resu_tt=0.


CASE(2) !sinus [-1,1] 
  Frequency=0.5
  Amplitude=1.
  Omega=2.*PP_Pi*Frequency
  Resu=Amplitude*SIN(Omega*SUM(Cent))
  Resu_t=-Amplitude*COS(Omega*SUM(Cent))*Omega*SUM(AdvVel)
  Resu_tt=-Amplitude*SIN(Omega*SUM(Cent))*Omega*SUM(AdvVel)*Omega*SUM(AdvVel)
CASE(3) 
  ! not used at the moment
CASE(4) ! quadratic
  ! g(t)
  Resu=5.*(x(1)-AdvVel(1)*teval)**2
  ! g'(t)
  Resu_t=-10.*AdvVel(1)*(x(1)-AdvVel(1)*teval)
  ! g''(t)
  Resu_tt=10.*AdvVel(1)*AdvVel(1)
! This is needed for Taylor Time Discretization
!  SELECT CASE(tDeriv)
!  CASE(0)
!    resu=5.*(x(1)-AdvVel(1)*t)**2
!  CASE(1)
!    Resu=-10.*AdvVel(1)*(x(1)-AdvVel(1)*t)
!  CASE(2)
!    Resu=10.*AdvVel(1)*AdvVel(1)
!  CASE DEFAULT
!    Resu=0.
!  END SELECT
CASE(5) ! Kopriva page 200, advection-diffusion, but for 3D with 1/( (4t+1)^(3/2) )
  x0 = (/-0.5,-0.5,-0.5/)
  Resu=1./((4.*teval+1.)**(1.5))*EXP(-(SUM((x(:)-AdvVel(:)*t-x0(:))**2))/(DiffC*(4.*teval+1.)))
  Resu_t=Resu   *( -6./(4.*teval+1.) &
                   +2./(DiffC*(4.*teval+1.))*SUM(AdvVel(:)*(x(:)-AdvVel(:)*teval-x0(:))) &
                   +4./(DiffC*(4.*teval+1.)**2)*SUM((x(:)-AdvVel(:)*teval-x0(:))**2) )
  Resu_tt=Resu_t*( -6./(4.*t+1.) &
                   +2./(DiffC*(4.*teval+1.))*SUM(AdvVel(:)*(x(:)-AdvVel(:)*teval-x0(:))) &
                   +4./(DiffC*(4.*teval+1.)**2)*SUM((x(:)-AdvVel(:)*teval-x0(:))**2) ) &
          + Resu*( 24./(4.*t+1.)**2 &
                   -8./(DiffC*(4.*teval+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                   -2./(DiffC*(4.*teval+1.))*SUM(AdvVel(:)*AdvVel(:)) &
                   -32./(DiffC*(4.*teval+1.)**3)*SUM((x(:)-AdvVel(:)*teval-x0(:))**2) & 
                   -8./(DiffC*(4.*teval+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*teval-x0(:))) )
CASE(6) !SINUS periodic, angle with IniWavenumber, exp(-|omega|^2*DiffC*t)*SIN(omega(:)*x(:))
         ! note that Omega=IniWavenumber*Pi
  Amplitude=-SUM(IniWavenumber(:)**2)*PP_Pi*PP_Pi*DiffC
  Resu   = EXP(Amplitude*teval)*SIN(PP_Pi*SUM(IniWavenumber(:)*cent(:)))
  Resu_t = Amplitude*Resu
  Resu_tt= Amplitude*Resu_t
CASE(7) ! Kopriva page 200, advection-diffusion, but only 2D in x,y planes  with 1/( (4t+1) )
  x0 = (/-0.,-0.,-0./)
  Resu=1./((4.*teval+1.))*EXP(-(SUM((x(1:2)-AdvVel(1:2)*t-x0(1:2))**2))/(DiffC*(4.*teval+1.)))
  Resu_t=Resu   *( -4./(4.*teval+1.) &
                   +2./(DiffC*(4.*teval+1.))*SUM(AdvVel(1:2)*(x(1:2)-AdvVel(1:2)*teval-x0(1:2))) &
                   +4./(DiffC*(4.*teval+1.)**2)*SUM((x(1:2)-AdvVel(1:2)*teval-x0(1:2))**2) )
  Resu_tt=Resu_t*( -4./(4.*teval+1.) &
                   +2./(DiffC*(4.*teval+1.))*SUM(AdvVel(1:2)*(x(1:2)-AdvVel(1:2)*teval-x0(1:2))) &
                   +4./(DiffC*(4.*teval+1.)**2)*SUM((x(1:2)-AdvVel(1:2)*teval-x0(1:2))**2) ) &
          + Resu*( 16./(4.*t+1.)**2 &
                   -8./(DiffC*(4.*teval+1.)**2)*SUM(AdvVel(1:2)*(x(1:2)-AdvVel(1:2)*teval-x0(1:2))) &
                   -2./(DiffC*(4.*teval+1.))*SUM(AdvVel(1:2)*AdvVel(1:2)) &
                   -32./(DiffC*(4.*teval+1.)**3)*SUM((x(1:2)-AdvVel(1:2)*teval-x0(1:2))**2) & 
                   -8./(DiffC*(4.*teval+1.)**2)*SUM(AdvVel(1:2)*(x(1:2)-AdvVel(1:2)*teval-x0(1:2))) )
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
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
!==================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_PreProc
USE MOD_Mesh_vars,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
REAL,INTENT(INOUT)              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!INTEGER                          :: iElem,i,j,k
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(1) ! constant 
CASE(2) ! sinus
CASE(3)
CASE(4) ! exact function
!  DO iElem=1,nElems
!    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
     ! code
!    END DO; END DO; END DO ! i,j,k
!  END DO
CASE(101)
END SELECT ! ExactFunction
END SUBROUTINE CalcSource



SUBROUTINE FinalizeEquation()
!==================================================================================================================================
! Get the constant advection velocity vector from the ini file
!==================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

