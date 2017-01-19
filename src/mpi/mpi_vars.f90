!=====================================================================================================================================
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

!==================================================================================================================================
!> Contains the variables that are used to control non-blocking communication
!==================================================================================================================================
MODULE MOD_MPI_Vars
#if MPI
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: MPIRequest_U(:,:)        !< communication handle for the surface solution
INTEGER,ALLOCATABLE :: MPIRequest_Flux(:,:)     !< communication handle for the surface flux
#ifdef EDDYVISCOSITY
INTEGER,ALLOCATABLE :: MPIRequest_DeltaS(:,:)   !< communication handle for the surface flux used for overintegration
#endif
#if PARABOLIC
INTEGER,ALLOCATABLE :: MPIRequest_gradU(:,:,:)  !< communication handle for the surface gradients
#endif /*PARABOLIC*/
INTEGER             :: nSendVal                
INTEGER             :: nRecVal
INTEGER             :: DataSizeSide             !< datasize of one face, =PP_nVar*(PP_N+1)**2

INTEGER             :: SideID_start,SideID_end
INTEGER             :: nNbProcs                 !< number of neighbor procs, is set in ReadMesh
INTEGER,ALLOCATABLE :: NbProc(:)                !< list of neighbor procs; allocated from 1:nNbProcs, is set in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_Proc(:)        !< number of mpisides for all neighbor procs, is set in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_MINE_Proc(:)   !< number of MINE mpisides for all neighbor procs, is set in setLocalSideIDs
                                                !< (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: nMPISides_YOUR_Proc(:)   !< number of YOUR mpisides for all neighbor procs, is set in setLocalSideIDs
                                                !< (prepare_mesh.f90)

INTEGER,ALLOCATABLE :: offsetMPISides_MINE(:)   !< gives position of send/recv block in *_MINE arrays,allocated from 0:nNbProcs, is
                                                !< set in setLocalSideIDs (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: offsetMPISides_YOUR(:)   !< gives position of send/recv block in *_YOUR arrays,allocated from 0:nNbProcs,
                                                !< is set in setLocalSideIDs (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: offsetElemMPI(:)         !< gives offsetposotion of elements of all procs, allocated from 1:nProcessors set
                                                !< in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_send(:,:)      !< number of sides to send, (1:nNbProcs,1:2), last index: 1: SEND MINE, 2: SEND YOUR
INTEGER,ALLOCATABLE :: nMPISides_rec(:,:)       !< number of sides to receive, (1:nNbProcs,1:2), last index: 1: RECEIVE YOUR,
                                                !< 2: RECEIVE MINE
INTEGER,ALLOCATABLE :: OffsetMPISides_send(:,:) !< offset of sides to send,(1:nNbProcs,1:2), last index: 1: SEND MINE, 2: SEND YOUR
INTEGER,ALLOCATABLE :: OffsetMPISides_rec(:,:)  !< offset of sides to receive, (1:nNbProcs,1:2), last index: 1: RECEIVE YOUR,
                                                !< 2: RECEIVE MINE
!==================================================================================================================================
#endif
END MODULE MOD_MPI_Vars
