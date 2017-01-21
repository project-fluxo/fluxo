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

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradUy_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradUz_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
ALLOCATE(gradUx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxX        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxY        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxZ        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.
FluxX=0.
FluxY=0.
FluxZ=0.

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
gradUx=0.
gradUy=0.
gradUz=0.

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
SUBROUTINE Lifting(t)
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
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI,               ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,         ONLY: nSides,FirstSlaveSide,LastSlaveSide 
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: t        !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! fill the global surface flux list
!fluxX=0. !don't nullify fluxes if not really needed
!fluxY=0. !don't nullify fluxes if not really needed
!fluxZ=0. !don't nullify fluxes if not really needed
CALL VNullify(nTotalU,gradUx)
CALL VNullify(nTotalU,gradUy)
CALL VNullify(nTotalU,gradUz)

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
! compute volume integral contribution and add to gradU, Jacobian not yet included
! this is onyl the local gradient!
CALL Lifting_VolInt(U,GradUx,GradUy,GradUz)

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,FluxX, FluxY, FluxZ)
CALL Lifting_FillFlux(FluxX,FluxY,FluxZ,doMPISides=.FALSE.)

!Start now with gradUx
CALL Flux_Mortar(FluxX,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxX,gradUx,doMPISides=.FALSE.)

#if MPI
! Complete send / receive FluxX
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,1,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL Flux_Mortar(FluxX,doMPISides=.TRUE.,weak=.FALSE.)
CALL Lifting_SurfInt(FluxX,GradUx,doMPISides=.TRUE.)
#endif /*MPI*/


!apply the Jacobian for gradUx
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradUx,sJ) !gradUx(:,i)=gradUx(:,i)*sJ(i)


#if MPI
!prolongtoface and start send/receive gradUx_slave
CALL StartReceiveMPIData(gradUx_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,1,RECV),SendID=2)
CALL ProlongToFace(gradUx,gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUx_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,1,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradUx,gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.FALSE.)

!gradUx sides finished (will be received in DG)
! now gradUy sides
CALL Flux_Mortar(FluxY,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxY,gradUy,doMPISides=.FALSE.)

#if MPI
! Complete send / receive FluxY
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,2,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL Flux_Mortar(FluxY,doMPISides=.TRUE.,weak=.FALSE.)
CALL Lifting_SurfInt(FluxY,GradUy,doMPISides=.TRUE.)
#endif /*MPI*/

!apply the Jacobian for gradUy
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradUy,sJ) !gradUy(:,i)=gradUy(:,i)*sJ(i)

#if MPI
!ProlongToFace and start send/receive gradUy_slave
CALL StartReceiveMPIData(gradUy_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,2,RECV),SendID=2)
CALL ProlongToFace(gradUy,gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUy_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,2,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradUy,gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.FALSE.)

!gradUy sides finished (will be received in DG)
! now gradUz sides

CALL Flux_Mortar(FluxZ,doMPISides=.FALSE.,weak=.FALSE.)

CALL Lifting_SurfInt(FluxZ,gradUz,doMPISides=.FALSE.)
#if MPI
! Complete send / receive FluxZ
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Lifting(:,3,:))!Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL Flux_Mortar(FluxZ,doMPISides=.TRUE.,weak=.FALSE.)
CALL Lifting_SurfInt(FluxZ,GradUz,doMPISides=.TRUE.)
#endif /*MPI*/

!apply the Jacobian for gradUz
CALL V2D_M_V1D(PP_nVar,nTotal_IP,gradUz,sJ) !gradUz(:,i)=gradUz(:,i)*sJ(i)

#if MPI
!ProlongToFace and start send/receive gradUz_slave
CALL StartReceiveMPIData(gradUz_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,3,RECV),SendID=2)
CALL ProlongToFace(gradUz,gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUz_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide,MPIRequest_Lifting(:,3,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI_MINE sides 
CALL ProlongToFace(gradUz,gradUz_master,gradUz_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.FALSE.)

!gradUz sides finished (will be received in DG)

END SUBROUTINE Lifting



!==================================================================================================================================
!> Deallocate BR1 arrays (volume and surface gradients and gradient fluxes)
!==================================================================================================================================
SUBROUTINE FinalizeLifting()
! MODULES
USE MOD_Lifting_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(gradUx_slave)
SDEALLOCATE(gradUy_slave)
SDEALLOCATE(gradUz_slave)
SDEALLOCATE(gradUx_master)
SDEALLOCATE(gradUy_master)
SDEALLOCATE(gradUz_master)
SDEALLOCATE(gradUx)
SDEALLOCATE(gradUy)
SDEALLOCATE(gradUz)
SDEALLOCATE(FluxX)
SDEALLOCATE(FluxY)
SDEALLOCATE(FluxZ)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /* PARABOLIC */
