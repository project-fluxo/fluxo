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
#if PARABOLIC
#include "defines.h"

!==================================================================================================================================
!> Contains the BR1 lifting procedure
!>
!> The BR1 scheme has been found to be unstable for purely elliptic equations. But for advection-diffusion systems, it seems to
!> indifferent, not adding any numerical diffusion.
!> 
!==================================================================================================================================
MODULE MOD_Lifting
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitLifting
  MODULE PROCEDURE InitLifting
END INTERFACE

INTERFACE Lifting
  MODULE PROCEDURE Lifting
END INTERFACE

INTERFACE FinalizeLifting
  MODULE PROCEDURE FinalizeLifting
END INTERFACE

PUBLIC:: DefineParametersLifting
PUBLIC:: InitLifting
PUBLIC:: Lifting
PUBLIC:: FinalizeLifting
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersLifting()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Lifting")
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> Initialize the BR1 lifting: allocate the arrays required for the BR1 lifting procedure.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_DG_Vars,              ONLY: DGInitIsDone
USE MOD_Mesh_Vars,            ONLY: nSides,nElems
USE MOD_Mesh_Vars,            ONLY: FirstSlaveSide,LastSlaveSide 
USE MOD_ReadInTools,          ONLY: GETREAL
#if MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.DGInitIsDone).OR.LiftingInitIsDone)THEN
   SWRITE(*,*) "InitDG not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR1 ...'
#if (PP_Lifting_Var==1)
SWRITE(UNIT_stdOut,'(A)') '    USING CONSERVATIVE VARIABLES FOR GRADIENT!'
#elif (PP_Lifting_Var==2)
SWRITE(UNIT_stdOut,'(A)') '    USING PRIMITIVE VARIABLES FOR GRADIENT!'
#elif (PP_Lifting_Var==3)
SWRITE(UNIT_stdOut,'(A)') '    USING ENTROPY VARIABLES FOR GRADIENT!'
#endif


! We store the interior gradients at the each element face
ALLOCATE(gradPx_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradPy_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradPz_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradPx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradPy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradPz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxX        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxY        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxZ        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
gradPx_slave=0.
gradPy_slave=0.
gradPz_slave=0.
gradPx_master=0.
gradPy_master=0.
gradPz_master=0.
FluxX=0.
FluxY=0.
FluxZ=0.

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradPx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradPy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradPz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
gradPx=0.
gradPy=0.
gradPz=0.

LiftingInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LIFTING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLifting



!==================================================================================================================================
!> Computes the DG gradients using the BR1 scheme in x/y/z direction.
!>
!> To calculate the lifted gradients we need to:
!> * Calculate the volume integral in strong form
!> " add the surface flux to the volume gradient
!> * Prolong the volume gradient to the interface (only MPI)
!> * Send gradients on MPI faces
!> * Prolong inner faces
!> * The surface integral of the fluxes is calculated, here the strong form is used 
!>   In the same step, the surface contribution is added to the prolonged volume contribution
!==================================================================================================================================
SUBROUTINE Lifting(tIn)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Vector,            ONLY: VNullify,V2D_M_V1D
USE MOD_DG_Vars,           ONLY: U,nTotalU,nTotal_IP
USE MOD_FillMortar,        ONLY: U_Mortar,Flux_Mortar
USE MOD_Lifting_SurfInt,   ONLY: Lifting_SurfInt
USE MOD_Lifting_VolInt,    ONLY: Lifting_VolInt
USE MOD_ProlongToFace,     ONLY: ProlongToFace
USE MOD_Lifting_FillFlux,  ONLY: Lifting_FillFlux,Lifting_FillFlux_BC
USE MOD_Mesh_Vars,         ONLY: sJ 
USE MOD_Equation_Vars,     ONLY: ConvertToGradPrimVec
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI,               ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,         ONLY: nSides,FirstSlaveSide,LastSlaveSide 
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: tIn        !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! fill the global surface flux list
!fluxX=0. !don't nullify fluxes if not really needed
!fluxY=0. !don't nullify fluxes if not really needed
!fluxZ=0. !don't nullify fluxes if not really needed
CALL VNullify(nTotalU,gradPx)
CALL VNullify(nTotalU,gradPy)
CALL VNullify(nTotalU,gradPz)

#if MPI
! Receive YOUR
CALL StartReceiveMPIData(FluxX,DataSizeSide,1,nSides,MPIRequest_Lifting(:,1,RECV),SendID=1)
CALL StartReceiveMPIData(FluxY,DataSizeSide,1,nSides,MPIRequest_Lifting(:,2,RECV),SendID=1)
CALL StartReceiveMPIData(FluxZ,DataSizeSide,1,nSides,MPIRequest_Lifting(:,3,RECV),SendID=1)
! Compute lifting MPI fluxes
CALL Lifting_FillFlux(FluxX,FluxY,FluxZ,doMPISides=.TRUE.)
! Start Send MINE
CALL StartSendMPIData(   FluxX,DataSizeSide,1,nSides,MPIRequest_Lifting(:,1,SEND),SendID=1)
CALL StartSendMPIData(   FluxY,DataSizeSide,1,nSides,MPIRequest_Lifting(:,2,SEND),SendID=1)
CALL StartSendMPIData(   FluxZ,DataSizeSide,1,nSides,MPIRequest_Lifting(:,3,SEND),SendID=1)
#endif /*MPI*/
! compute volume integral contribution and add to gradP, Jacobian not yet included
! this is onyl the local gradient!
CALL Lifting_VolInt(U,gradPx,gradPy,gradPz)

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(tIn,FluxX, FluxY, FluxZ)
CALL Lifting_FillFlux(FluxX,FluxY,FluxZ,doMPISides=.FALSE.)

!Start now with gradPx
!TODO CALL Flux_Mortar(FluxX,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxX,gradPx,doMPISides=.FALSE.)

#if MPI
! Complete send / receive FluxX
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,1,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
!TODO CALL Flux_Mortar(FluxX,doMPISides=.TRUE.,weak=.FALSE.)
CALL Lifting_SurfInt(FluxX,gradPx,doMPISides=.TRUE.)
#endif /*MPI*/


!apply the Jacobian for gradPx
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradPx,sJ) !gradPx(:,i)=gradPx(:,i)*sJ(i)

!TODO entropy: "gradWx" -> gradPx back, berofre P2F (since GL only copies point data to surfaces)

! now, grad P is the gradient of primivite variables (before, container for cons/prim/entropy vars) 
CALL ConvertToGradPrimVec(nTotal_IP,U,gradPx) !overwrites gradPx!

#if MPI
!prolongtoface and start send/receive gradPx_slave
CALL StartReceiveMPIData(gradPx_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,1,RECV),SendID=2)
CALL ProlongToFace(gradPx,gradPx_master,gradPx_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradPx_master,gradPx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradPx_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,1,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradPx,gradPx_master,gradPx_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradPx_master,gradPx_slave,doMPISides=.FALSE.)

!gradPx sides finished (will be received in DG)
! now gradPy sides
!TODO CALL Flux_Mortar(FluxY,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxY,gradPy,doMPISides=.FALSE.)

#if MPI
! Complete send / receive FluxY
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,2,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
! TODO CALL Flux_Mortar(FluxY,doMPISides=.TRUE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxY,gradPy,doMPISides=.TRUE.)
#endif /*MPI*/

!apply the Jacobian for gradPy
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradPy,sJ) !gradPy(:,i)=gradPy(:,i)*sJ(i)

! now, grad P is the gradient of primivite variables (before, container for cons/prim/entropy vars) 
CALL ConvertToGradPrimVec(nTotal_IP,U,gradPy) !overwrites gradPy!

#if MPI
!ProlongToFace and start send/receive gradPy_slave
CALL StartReceiveMPIData(gradPy_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,2,RECV),SendID=2)
CALL ProlongToFace(gradPy,gradPy_master,gradPy_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradPy_master,gradPy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradPy_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,2,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradPy,gradPy_master,gradPy_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradPy_master,gradPy_slave,doMPISides=.FALSE.)

!gradPy sides finished (will be received in DG)
! now gradPz sides

!TODO CALL Flux_Mortar(FluxZ,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxZ,gradPz,doMPISides=.FALSE.)
#if MPI
! Complete send / receive FluxZ
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,3,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
!TODO CALL Flux_Mortar(FluxZ,doMPISides=.TRUE.,weak=.FALSE.)
CALL Lifting_SurfInt(FluxZ,gradPz,doMPISides=.TRUE.)
#endif /*MPI*/

!apply the Jacobian for gradPz
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradPz,sJ) !gradPz(:,i)=gradPz(:,i)*sJ(i)

! now, grad P is the gradient of primivite variables (before, container for cons/prim/entropy vars) 
CALL ConvertToGradPrimVec(nTotal_IP,U,gradPz) !overwrites gradPz!

#if MPI
!ProlongToFace and start send/receive gradPz_slave
CALL StartReceiveMPIData(gradPz_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,3,RECV),SendID=2)
CALL ProlongToFace(gradPz,gradPz_master,gradPz_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradPz_master,gradPz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradPz_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,3,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradPz,gradPz_master,gradPz_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradPz_master,gradPz_slave,doMPISides=.FALSE.)

!gradPz sides finished (will be received in DG)

END SUBROUTINE Lifting



!==================================================================================================================================
!> Deallocate BR1 arrays (volume and surface gradients and gradient fluxes)
!==================================================================================================================================
SUBROUTINE FinalizeLifting()
! MODULES
USE MOD_Lifting_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(gradPx_slave)
SDEALLOCATE(gradPy_slave)
SDEALLOCATE(gradPz_slave)
SDEALLOCATE(gradPx_master)
SDEALLOCATE(gradPy_master)
SDEALLOCATE(gradPz_master)
SDEALLOCATE(gradPx)
SDEALLOCATE(gradPy)
SDEALLOCATE(gradPz)
SDEALLOCATE(FluxX)
SDEALLOCATE(FluxY)
SDEALLOCATE(FluxZ)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /* PARABOLIC */
