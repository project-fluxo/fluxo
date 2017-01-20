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
!> Computes the DGSEM spatial operator and updates residual Ut

!> Contains the routines to 
!> - initialize and finalize the DG global variables and the DG basis
!> - compute the DG spatial operators/residuals(Ut) using U from the volume, surface and source contribution, incl. 
!> lifting for the gradients and parallelization
!==================================================================================================================================
MODULE MOD_DG
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDG
  MODULE PROCEDURE InitDG
END INTERFACE


INTERFACE DGTimeDerivative
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE

INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE

PUBLIC:: InitDG
PUBLIC:: DGTimeDerivative
PUBLIC:: FinalizeDG
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Allocate all global DG variables like U (solution in volume), U_slave/U_master (solution on faces), Flux, Ut (DG time derivative),
!> also fill the initial solution and call init DG basis. Operator building are also initialized by calling InitDGBasis.  
!==================================================================================================================================
SUBROUTINE InitDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars
USE MOD_Restart_Vars,       ONLY: DoRestart,RestartInitIsDone
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nElems,nSides,MeshInitIsDone
USE MOD_Equation_Vars,      ONLY:IniExactFunc
USE MOD_Equation_Vars,      ONLY:EquationInitIsDone
USE MOD_Equation,           ONLY:FillIni
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.EquationInitIsDone) &
   .OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone)THEN
   CALL abort(__STAMP__, &
   'InitDG not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

CALL initDGbasis(PP_N,xGP,wGP,wBary)
! the local DG solution
ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
! the time derivative computed with the DG scheme
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
U=0.
Ut=0.

nDOFElem=(PP_N+1)**3
nTotalU=PP_nVar*nDOFElem*nElems
nTotal_face=(PP_N+1)*(PP_N+1)
nTotal_vol=nTotal_face*(PP_N+1)
nTotal_IP=nTotal_vol*nElems
nTotalU=PP_nVar*nTotal_vol*nElems




! Allocate the 2D solution vectors on the sides, one array for the data belonging to the proc (the master) 
! and one for the sides which belong to another proc (slaves): side-based
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,1:nSides))
U_master=0.
U_slave=0.

! unique flux per side
ALLOCATE(Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides))
Flux=0.

#if (PP_DiscType==1)
  SWRITE(UNIT_StdOut,'(A)') ' => VolInt in weak form!'
#elif (PP_DiscType==2)
  SWRITE(UNIT_StdOut,'(A)') ' => VolInt in fluxdifferencing form!'
#endif /*PP_DiscType*/
#ifdef CARTESIANFLUX
  SWRITE(UNIT_StdOut,'(A)') ' => CARTESIANFLUX=T: FLUXES ONLY FOR CARTESIAN MESHES!!!'
#endif

IF(.NOT.DoRestart)THEN
  ! U is filled with the ini solution
  CALL FillIni(IniExactFunc,U)
END IF

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDG


!==================================================================================================================================
!> Allocate and initialize the building blocks for the DG operator: Differentiation matrices and prolongation operators
!==================================================================================================================================
SUBROUTINE InitDGbasis(N_in,xGP,wGP,wBary)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Basis,ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,ONLY:PolynomialDerivativeMatrix,LagrangeInterpolationPolys
USE MOD_DG_Vars,ONLY:D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatMinus0,L_HatPlus
USE MOD_DG_Vars,ONLY:DvolSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                             :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)              :: xGP                    !< Gauss/Gauss-Lobatto Nodes 
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wBary                  !< Barycentric weights to evaluate the Gauss/Gauss-Lobatto lagrange basis
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv     
REAL,DIMENSION(0:N_in)                     :: L_Minus,L_Plus
INTEGER                                    :: iMass         
!===================================================================================================================================
ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(Dvolsurf(0:N_in,0:N_in))
ALLOCATE(D(    0:N_in,0:N_in), D_T(    0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

! Build D_Hat matrix. (D^ = M^(-1) * D^T * M
M(:,:)=0.
Minv(:,:)=0.
DO iMass=0,N_in
  M(iMass,iMass)=wGP(iMass)
  Minv(iMass,iMass)=1./wGP(iMass)
END DO
D_Hat(:,:) = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T= TRANSPOSE(D_hat)

!modified D matrix for fluxdifference volint
Dvolsurf=2.*D
Dvolsurf(0,0)=2*D(0,0)+1./wGP(0)
Dvolsurf(N_in,N_in)=2*D(N_in,N_in)-1./wGP(N_in)

! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
L_HatPlus(:) = MATMUL(Minv,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)
L_HatMinus(:) = MATMUL(Minv,L_Minus)
L_HatMinus0 = L_HatMinus(0)
END SUBROUTINE InitDGbasis

!==================================================================================================================================
!> summary: Computes the residual Ut = \f$ \frac {d\vec{U}} {dt} \f$ from the current solution U employing the DG method.
!> Computes the weak DGSEM space operator from surface, volume and source contributions. To do this we need to:
!>
!> - Prolong the solution from the volume to the interface
!> - Invoke the lifting operator to calculate the gradients
!> - Perform the volume integral
!> - Perform the surface integral
!> - If needed, add source terms to the residual
!==================================================================================================================================
SUBROUTINE DGTimeDerivative_weakForm(t)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars             ,ONLY: Ut,U,U_slave,U_master,Flux,L_HatPlus,L_HatMinus
USE MOD_DG_Vars             ,ONLY: D_Hat_T,nDOFElem
USE MOD_DG_Vars             ,ONLY: nTotalU
USE MOD_VolInt              ,ONLY: VolInt
USE MOD_SurfInt             ,ONLY: SurfInt
USE MOD_ProlongToFace       ,ONLY: ProlongToFace
USE MOD_FillFlux            ,ONLY: FillFlux
USE MOD_ApplyJacobian       ,ONLY: ApplyJacobian
USE MOD_Interpolation_Vars  ,ONLY: L_Minus,L_Plus
USE MOD_Testcase            ,ONLY: TestcaseSource
USE MOD_Testcase_Vars       ,ONLY: doTCSource
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Mesh_Vars           ,ONLY: Metrics_fTilde ,Metrics_gTilde ,Metrics_hTilde
USE MOD_Exactfunc           ,ONLY: CalcSource
USE MOD_Equation_Vars       ,ONLY: doCalcSource
USE MOD_FillMortar          ,ONLY: U_Mortar,Flux_Mortar
#if PARABOLIC
USE MOD_Lifting             ,ONLY: Lifting
USE MOD_Lifting_Vars
#endif /*PARABOLIC*/
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI                 ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,           ONLY: nSides
#endif /*MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! -----------------------------------------------------------------------------
! MAIN STEPS       
! -----------------------------------------------------------------------------
! 0.  Convert volume solution to primitive 
! 1.  Pronlong to face (fill U_master/slave)
! 3.  ConsToPrim of face data (U_master/slave)
! 5.  Lifting
! 6.  Viscous volume integral (DG only)
! 7.  Fill flux (Riemann solver) + surface integral
! 8. Ut = -Ut
! 9. Sponge and source terms
! 10. apply Jacobian
! -----------------------------------------------------------------------------

! Nullify arrays
! TODO fix!
! NOTE: UT and U are nullified in DGInit, and Ut is set directly in the volume integral, so in this implementation, 
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
CALL VNullify(nTotalU,Ut)

! 1. Prolong the solution to the face integration points for flux computation (and do overlapping communication)
! -----------------------------------------------------------------------------------------------------------
! General idea: The slave sends its surface data to the master, where the flux is computed and sent back to the slaves.
! Steps:
! * (these steps are done for all slave MPI sides first and then for all remaining sides): 
! 1.1)  Prolong solution to faces and store in U_master/slave. Use them to build mortar data (split into 2/4 smaller sides).
!       Then U_slave can be communicated from the slave to master MPI side.
! 1.2)  Finish all started MPI communications (after step 2. due to latency hiding) 

#if MPI
! Step 1 for all slave MPI sides
! 1.1)
CALL StartReceiveMPIData(U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / U_slave: slave -> master
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / U_slave: slave -> master
#endif /* MPI */

! Step 1 for all remaining sides
! 1.1)
CALL ProlongToFace(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)


#if MPI
! 1.4) complete send / receive of side data from step 1.
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! U_slave: slave -> master 
#endif

#if PARABOLIC
! 5. Lifting 
! Compute the gradients using Lifting (BR1 scheme,BR2 scheme ...)
! The communication of the gradients is initialized within the lifting routines
CALL Lifting(U,U_master,U_slave,t)
#endif /*PARABOLIC*/

! 6. Compute volume integral contribution and add to Ut
CALL VolInt(Ut)


#if PARABOLIC && MPI
! Complete send / receive for gradUx, gradUy, gradUz, started in the lifting routines
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) ! gradUx,y,z: slave -> master
#endif /*PARABOLIC && MPI*/


! 7. Fill flux and Surface integral
! General idea: U_master/slave and gradUx,y,z_master/slave are filled and can be used to compute the Riemann solver
!               and viscous flux at the faces. This is done for the MPI master sides first, to start communication early
!               and then for all other sides.
!               After communication from master to slave the flux can be integrated over the faces.
! Steps:
! * (step 8.2 is done for all MPI master sides first and then for all remaining sides)
! * (step 8.3 and 8.4 are done for all other sides first and then for the MPI master sides) 
! 8.3)  Fill flux (Riemann solver + viscous flux)
! 8.4)  Combine fluxes from the 2/4 small mortar sides to the flux on the big mortar side (when communication finished)
! 8.5)  Compute surface integral 

#if MPI
! 8.3)
CALL StartReceiveMPIData(Flux, DataSizeSide, 1,nSides,MPIRequest_Flux( :,SEND),SendID=1)
                                                                              ! Receive YOUR / Flux_slave: master -> slave
CALL FillFlux(Flux,doMPISides=.TRUE.)
CALL StartSendMPIData(Flux, DataSizeSide, 1,nSides,MPIRequest_Flux( :,RECV),SendID=1)
                                                                              ! Send MINE  /   Flux_slave: master -> slave
#endif /* MPI*/

! 8.3)
CALL FillFlux(Flux,doMPISides=.FALSE.)
! 8.4)
CALL Flux_Mortar(Flux,doMPISides=.FALSE.,weak=.TRUE.)
! 8.5)
CALL SurfInt(PP_N,Flux_master,Flux_slave,Ut,.FALSE.,L_HatMinus,L_hatPlus)

#if MPI
! 8.4)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux )                                        ! Flux_slave: master -> slave 
CALL Flux_Mortar(Flux,doMPISides=.TRUE.,weak=.TRUE.)
! 8.5)
CALL SurfInt(PP_N,Flux,Ut,.TRUE.,L_HatMinus,L_HatPlus)
#endif /*MPI*/


! 10. Swap to right sign :) 
Ut=-Ut

!  Compute source terms and sponge (in physical space, conversion to reference space inside routines)
IF(doCalcSource) CALL CalcSource(Ut,t)
IF(doTCSource)   CALL TestcaseSource(Ut)

! Apply Jacobian 
CALL ApplyJacobian(Ut,toPhysical=.TRUE.)

END SUBROUTINE DGTimeDerivative_weakForm




!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_DG_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(D)
SDEALLOCATE(D_T)
SDEALLOCATE(D_Hat)
SDEALLOCATE(D_Hat_T)
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(Ut)
SDEALLOCATE(U)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(Flux)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG

END MODULE MOD_DG
