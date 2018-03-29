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
USE MOD_Mesh_Vars,          ONLY: firstSlaveSide,LastSlaveSide 
USE MOD_Equation_Vars,      ONLY: IniExactFunc
USE MOD_Equation_Vars,      ONLY: EquationInitIsDone
USE MOD_Equation,           ONLY: FillIni
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
nTotal_face=(PP_N+1)*(PP_N+1)
nTotal_vol=nTotal_face*(PP_N+1)
nTotal_IP=nTotal_vol*nElems
nTotalU=PP_nVar*nTotal_vol*nElems




! Allocate the 2D solution vectors on the sides, one array for the data belonging to the proc (the master) 
! and one for the sides which belong to another proc (slaves): side-based
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
U_master=0.
U_slave=0.

! unique flux per side
ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Flux_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
Flux_master=0.
Flux_slave=0.

#if (PP_DiscType==1)
  SWRITE(UNIT_StdOut,'(A)') ' => VolInt in weak form!'
#elif (PP_DiscType==2)
  SWRITE(UNIT_StdOut,'(A)') ' => VolInt in fluxdifferencing form!'
#endif /*PP_DiscType*/

IF(.NOT.DoRestart)THEN
  ! U is filled with the ini solution
  CALL FillIni(IniExactFunc,U)
END IF

!#if MPI
!!communicate master surface metrics
!ALLOCATE(geo(13,0:PP_N,0:PP_N,1:nSides))
!geo(  1:3,:,:,:)=NormVec( :,:,:,:) 
!geo(  4:6,:,:,:)=TangVec1(:,:,:,:)
!geo(  7:9,:,:,:)=TangVec2(:,:,:,:)
!geo(   10,:,:,:)=SurfElem(  :,:,:)
!geo(11:13,:,:,:)=Face_xGP(:,:,:,:)
!!!start geo communication
!CALL StartReceiveMPIData(geo, 13*(PP_N+1)**2, 1,nSides,MPIRequest_U( :,SEND),SendID=1) ! Receive YOUR  (sendID=1) 
!CALL StartSendMPIData(   geo, 13*(PP_N+1)**2, 1,nSides,MPIRequest_U( :,RECV),SendID=1) ! Send MINE (SendID=1) 
!!finish geo communication
!CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U )  ! Flux, MPI_MINE -> MPI_YOUR 
!ASSOCIATE(low=>FirstMPISide_YOUR,up=>LastMPISide_YOUR)
!NormVec (:,:,:,low:up)=geo(  1:3,:,:,low:up)
!TangVec1(:,:,:,low:up)=geo(  4:6,:,:,low:up)
!TangVec2(:,:,:,low:up)=geo(  7:9,:,:,low:up)
!SurfElem(  :,:,low:up)=geo(   10,:,:,low:up)
!Face_xGP(:,:,:,low:up)=geo(11:13,:,:,low:up)
!END ASSOCIATE
!DEALLOCATE(geo)
!#endif /*MPI*/
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
USE MOD_Basis,              ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys
USE MOD_DG_Vars,            ONLY: D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatMinus0,L_HatPlus
#if PP_DiscType==2
USE MOD_DG_Vars,            ONLY: DvolSurf
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in      !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)  :: xGP       !< Gauss/Gauss-Lobatto Nodes 
REAL,DIMENSION(0:N_in),INTENT(IN)  :: wGP       !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)  :: wBary     !< Barycentric weights to evaluate the Gauss/Gauss-Lobatto lagrange basis
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(0:N_in,0:N_in)      :: M,Minv     
REAL,DIMENSION(0:N_in)             :: L_Minus,L_Plus
INTEGER                            :: iMass         
!===================================================================================================================================
ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
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

#if PP_DiscType==2
ALLOCATE(Dvolsurf(0:N_in,0:N_in))
!modified D matrix for fluxdifference volint
Dvolsurf=2.0*D
Dvolsurf(0,0)=2.0*D(0,0)+1.0/wGP(0)
Dvolsurf(N_in,N_in)=2.0*D(N_in,N_in)-1.0/wGP(N_in)
#endif /*PP_DiscType==2*/

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
! -----------------------------------------------------------------------------
!> Basic steps in the algorithm       
! -----------------------------------------------------------------------------
!> * Prolong to face (fill U_master/slave)
!> * Fill mortar solution
!> * Lifting
!> * volume integral
!> * fill boundary flux
!> * Fill flux
!> * Fill mortar flux
!> * surface integral
!> * apply Jacobian \( U_t=-\frac{1}{J} U_t \)
!> * add pointwise sources (without Jacobian!) 
! -----------------------------------------------------------------------------
!> MPI adds an additional layer, where prolongtoface and mortars are done first for sides involved in the MPI communication
!> only U_slave is used as a send and receive buffer
!> and once the solution is received the flux is computed on the master sides and sent back.
!> in between starting the send/receive process and its finish, we call buffer routines to hide the communication latency.
!==================================================================================================================================
SUBROUTINE DGTimeDerivative_weakForm(tIn)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector              ,ONLY: VNullify,V2D_M_V1D
USE MOD_DG_Vars             ,ONLY: Ut,U,U_slave,U_master,Flux_master,Flux_slave
USE MOD_FillMortar          ,ONLY: U_Mortar,Flux_Mortar
USE MOD_Mesh_Vars           ,ONLY: sJ
USE MOD_DG_Vars             ,ONLY: nTotalU,nTotal_IP
USE MOD_ProlongToFace       ,ONLY: ProlongToFace
USE MOD_VolInt              ,ONLY: VolInt
USE MOD_GetBoundaryFlux     ,ONLY: GetBoundaryFlux
USE MOD_FillFlux            ,ONLY: FillFlux
USE MOD_SurfInt             ,ONLY: SurfInt
USE MOD_Testcase_Vars       ,ONLY: doTCSource
USE MOD_Testcase_Source     ,ONLY: TestcaseSource
USE MOD_Equation_Vars       ,ONLY: doCalcSource
USE MOD_Equation            ,ONLY: CalcSource
#if PARABOLIC
USE MOD_Lifting             ,ONLY: Lifting
#endif /*PARABOLIC*/
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI                 ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars           ,ONLY: firstSlaveSide,LastSlaveSide 
#endif /*MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================


! Nullify arrays
CALL VNullify(nTotalU,Ut)

#if MPI
! Solution is always communicated on the U_Slave array
! start off with the receive command
CALL StartReceiveMPIData(U_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide, &
                         MPIRequest_U(:,SEND),SendID=2) ! Receive MINE (sendID=2) 
! prolong MPI sides and do the mortar on the MPI sides
CALL ProlongToFace(U,U_master,U_slave,doMPISides=.TRUE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
! start the sending command
CALL StartSendMPIData(U_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide, &
                      MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR (sendID=2) 
#endif /* MPI */

! for all remaining sides (buffer routine for latency hiding!)
!!write(*,*)'u in dgtimederivative', U
!!write(*,*)'u_slave before prolong', u_slave
!!write(*,*)'u_master before prolong', u_master

CALL ProlongToFace(U,U_master,U_slave,doMPISides=.FALSE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)


#if MPI
!complete send / receive of side data (WAIT...)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)  ! U_slave: MPI_YOUR -> MPI_MINE (_slave)
#endif

#if PARABOLIC
! Lifting 
! Compute the gradients using Lifting (BR1 scheme,BR2 scheme ...)
! The communication of the gradients is started within the lifting routines
CALL Lifting(tIn)
#endif /*PARABOLIC*/

! Compute volume integral contribution and add to Ut (should buffer latency of gradient communications)
CALL VolInt(Ut)


#if PARABOLIC && MPI
! Complete send / receive for gradUx, gradUy, gradUz, started in the lifting routines
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_Lifting) ! gradUx,y,z: MPI_YOUR -> MPI_MINE (_slave)
#endif /*PARABOLIC && MPI*/


#if MPI
! start off with the receive command
CALL StartReceiveMPIData(Flux_slave, DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_Flux( :,SEND),SendID=1) ! Receive YOUR  (sendID=1) 
! since mortar solutions are already there, we can directly fill the fluxes for all MPI sides 
CALL FillFlux(Flux_master,Flux_slave,doMPISides=.TRUE.)

! start the sending command
CALL StartSendMPIData(Flux_slave, DataSizeSide, firstSlaveSide,lastSlaveSide,MPIRequest_Flux( :,RECV),SendID=1) ! Send MINE (SendID=1) 
#endif /* MPI*/


! fill physical BC, inner side Flux and inner side Mortars (buffer for latency of flux communication)
CALL GetBoundaryFlux(tIn,Flux_master)
CALL FillFlux(Flux_master,Flux_slave,doMPISides=.FALSE.)

! here, weak=T:-F_slave is used, since small sides can be slave and must be added to big sides, which are always master!
CALL Flux_Mortar(Flux_master,Flux_slave,doMPISides=.FALSE.,weak=.TRUE.)

! add inner and BC side surface contibutions to time derivative 
CALL SurfInt(Flux_master,Flux_slave,Ut,doMPISides=.FALSE.)

#if MPI
! Complete send / receive for  Flux array
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux )  ! Flux, MPI_MINE -> MPI_YOUR 

! finally also collect all small side fluxes of MPI sides to big side fluxes
CALL Flux_Mortar(Flux_master,Flux_slave,doMPISides=.TRUE.,weak=.TRUE.) 

! update time derivative with contribution of MPI sides 
CALL SurfInt(Flux_master,Flux_slave,Ut,doMPIsides=.TRUE.)
#endif /*MPI*/


! We have to take the inverse of the Jacobians into account and directly swap the sign
CALL V2D_M_V1D(PP_nVar,nTotal_IP,Ut,(-sJ)) !Ut(:,i)=-Ut(:,i)*sJ(i)

!  Compute source terms and sponge (in physical space, conversion to reference space inside routines)
IF(doCalcSource) CALL CalcSource(Ut,tIn)
IF(doTCSource)   CALL TestcaseSource(Ut,tIn)

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
SDEALLOCATE(Flux_master)
SDEALLOCATE(Flux_slave)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG

END MODULE MOD_DG
