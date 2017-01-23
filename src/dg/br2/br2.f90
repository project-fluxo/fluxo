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
!> Contains the BR2 lifting procedure (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi, Rebay et al., "A high-order accurate discontinuous Finite Element method for inviscid an viscous 
!> turbomachinery flows", 1997. The lifted gradients are required for the viscous fluxes.
!>
!> The BR1 scheme has been found to be unstable for purely elliptic equations. Consequently, the BR2 scheme provides a stable method
!> with local lifting operators. In contrast to the BR1 scheme, the BR2 scheme requires a strong form of the lifting operator. The
!> surface gradients are lifted with an additional penalty term \( \eta_{BR2} \). Stability was shown for \( \eta_{BR2} > \)
!> number of element faces.
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

INTERFACE DefineParametersLifting
  MODULE PROCEDURE DefineParametersLifting
END INTERFACE

INTERFACE Lifting
  MODULE PROCEDURE Lifting
END INTERFACE

INTERFACE FinalizeLifting
  MODULE PROCEDURE FinalizeLifting
END INTERFACE

PUBLIC::DefineParametersLifting
PUBLIC::InitLifting
PUBLIC::Lifting
PUBLIC::FinalizeLifting
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
CALL prms%CreateRealOption(   'etaBR2',                "Lifting penalty for BR2. Increase improves stability at the cost of "//&
                                                       "performance and reduces jumps between two cells.", '2.')
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> \brief Initialize the BR2 lifting: get parameters and allocate the arrays required for the BR2 lifting procedure.
!>
!> Important parameters:
!> - etaBR2: Penalty term for the surface contribution of the BR2 lifting. Note, stability is shown only for \( \eta_{BR2} > \)
!>   number of element faces
!>
!> Default is \( eta_{BR2} = 2 \) (1D like)
!>
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_DG_Vars,              ONLY: DGInitIsDone
USE MOD_Mesh_Vars,            ONLY: nSides,nElems
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
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR2 ...'

etaBR2=GETREAL('etaBR2','2.')

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUy_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUz_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides))
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
!> Computes the DG gradients using the BR2 scheme in x/y/z direction.
!>
!> To calculate the lifted gradients we need to:
!> * For reducing the overhead of MPI communication, we calculate the surface contribution, i.e. the numerical fluxes, across the
!>   element faces at MPI boundaries first and communicate them right away. The fluxes are temporarilly stored in the gradient
!>   arrays.
!> * Then, the remainind surface fluxes across the inner element faces are computed
!> * Calculate the volume integral in strong form
!> * Prolong the volume contribution to the interface: Note, this is different from the BR1 scheme. Here we prolong the volume
!>   contribution to the surface before any surface gradient has been applied to the volume terms
!> * The surface integral of the fluxes is calculated, here the strong form is used and the etaBR2 factor is applied.
!>   In the same step, the surface contribution is added to the prolonged volume contribution
!==================================================================================================================================
SUBROUTINE Lifting(t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Vector,            ONLY: VNullify
USE MOD_DG_Vars,           ONLY: nTotalU,U,U_master,U_slave 
USE MOD_FillMortar,        ONLY: U_Mortar,Flux_Mortar
USE MOD_Lifting_SurfInt,   ONLY: Lifting_SurfInt
USE MOD_Lifting_VolInt,    ONLY: Lifting_VolInt
USE MOD_ProlongToFace,     ONLY: ProlongToFace
USE MOD_Lifting_FillFlux,  ONLY: Lifting_FillFlux,Lifting_FillFlux_BC
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI,               ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,         ONLY: nSides
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
CALL StartReceiveMPIData(FluxX,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1)
CALL StartReceiveMPIData(FluxY,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1)
CALL StartReceiveMPIData(FluxZ,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1)
! Compute lifting MPI fluxes
CALL Lifting_FillFlux(U_master,U_slave,FluxX,FluxY,FluxZ,doMPISides=.TRUE.)
! Start Send MINE
CALL StartSendMPIData(   FluxX,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1)
CALL StartSendMPIData(   FluxY,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1)
CALL StartSendMPIData(   FluxZ,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1)
#endif /*MPI*/

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,FluxX, FluxY, FluxZ)
CALL Lifting_FillFlux(U_master,U_slave,FluxX,FluxY,FluxZ,doMPISides=.FALSE.)

CALL Flux_Mortar(FluxX,doMPISides=.FALSE.,weak=.FALSE.)
CALL Flux_Mortar(FluxY,doMPISides=.FALSE.,weak=.FALSE.)
CALL Flux_Mortar(FluxZ,doMPISides=.FALSE.,weak=.FALSE.)

! compute volume integral contribution and add to gradU, Jacobian already included in BR2 volint
! this is onyl the local gradient!
CALL Lifting_VolInt(U,GradUx,GradUy,GradUz)
!

! The local gradient is now interpolated to the face of the grid cells
#if MPI
! Prolong to face for MPI sides - send direction
CALL ProlongToFace(gradUx,gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL ProlongToFace(gradUy,gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL ProlongToFace(gradUz,gradUz_master,gradUz_slave,doMPISides=.TRUE.)
#endif /*MPI*/
! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(gradUx,gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL ProlongToFace(gradUy,gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL ProlongToFace(gradUz,gradUz_master,gradUz_slave,doMPISides=.FALSE.)

#if MPI
! Complete send / receive
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU)
CALL Flux_Mortar(FluxX,doMPISides=.TRUE.,weak=.FALSE.)
CALL Flux_Mortar(FluxY,doMPISides=.TRUE.,weak=.FALSE.)
CALL Flux_Mortar(FluxZ,doMPISides=.TRUE.,weak=.FALSE.)

CALL StartReceiveMPIData(gradUx_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=2)
CALL StartReceiveMPIData(gradUy_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=2)
CALL StartReceiveMPIData(gradUz_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=2)

!update the local volume gradient together with the surface gradients with penalty etaBR2
!big mortar sides need to be finished and interpolated to small mortars before sending
CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUx_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=2)

CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUy_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=2)

CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUz_slave,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=2)
#endif /*MPI*/

! Add the surface lifting flux to the prolonged volume contributions of the gradients and computes the surface integral
CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.FALSE.)

END SUBROUTINE Lifting



!==================================================================================================================================
!> Deallocate BR2 arrays (volume and surface gradients and gradient fluxes)
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
