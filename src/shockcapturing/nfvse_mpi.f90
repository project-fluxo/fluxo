!==================================================================================================================================
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!> Module containing utilities for the NFVSE shock-capturing method
!> Attention: this module also contains the MPI routines for the communication of the blending coefficient in the NFVSE method:
!> -> Done here and not in DGTimeDerivative_weakForm for modularity!
!> -> Here the MPI approach is a bit different than in the main code: Both master and slave sides send and receive, then the
!>    blending coefficient correction is computed on both sides of the MPI interface. This is needed since we want to correct the DG
!>    volume integral BEFORE adding the diffusive contribution.
!===================================================================================================================================
module MOD_NFVSE_MPI
  implicit none
  
  private
  public ProlongBlendingCoeffToFaces, PropagateBlendingCoeff, Get_externalU
contains
!===================================================================================================================================
!> Prolong the blending coefficient to the faces
!===================================================================================================================================
  subroutine ProlongBlendingCoeffToFaces()
    use MOD_NFVSE_Vars, only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars , only: SideToElem, firstMortarInnerSide, nSides
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    !-------------------------------------------------------------------------------------------------------------------------------
    
    do sideID=firstMortarInnerSide, nSides
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      !master sides(ElemID,locSide and flip =-1 if not existing)
      if(ElemID.NE.-1) then ! element belonging to master side is on this processor
        alpha_Master(SideID) = alpha(ElemID) 
      end if
      
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      !slave side (nbElemID,nblocSide and flip =-1 if not existing)
      if(nbElemID.NE.-1) then! element belonging to slave side is on this processor
        alpha_Slave (SideID) = alpha(nbElemID) 
      end if
    end do
    
#if MPI
    call Start_BlendCoeff_MPICommunication()
#endif /*MPI*/
    
    call Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides=.FALSE.)
  end subroutine ProlongBlendingCoeffToFaces
!===================================================================================================================================
!> Propagate the blending coefficient to the neighbor elements
!> -> First the routine makes sure that all the MPI communication is done
!> -> Then, the blending coefficient is set to alpha = max (alpha, alpha_neighbor)
!===================================================================================================================================
  subroutine PropagateBlendingCoeff()
    use MOD_NFVSE_Vars         , only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide, SideToElem
    use MOD_Mesh_Vars          , only: firstMortarInnerSide, lastMortarInnerSide, firstMortarMPISide,lastMortarMPISide
    use MOD_NFVSE_Vars         , only: SpacePropFactor
    use MOD_Mesh_Vars          , only: MortarType,MortarInfo
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_MPI                , only: FinishExchangeMPIData
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha
#endif /*MPI*/
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    real    :: minAlpha (firstSlaveSide:LastSlaveSide)
    integer :: MortarSideID, nMortars, locSide, iMortar
    !-------------------------------------------------------------------------------------------------------------------------------
    
!   Finish the MPI exchange
!   -----------------------
#if MPI
    call FinishExchangeMPIData(4*nNbProcs,MPIRequest_alpha)
#endif /*MPI*/
    
!   Compute the minimum alpha (for each pair master/slave) and enforce it on non-mortar sides
!   -----------------------------------------------------------------------------------------
    do sideID=firstSlaveSide, lastSlaveSide
      
      minAlpha(sideID) = max( alpha_Master(sideID), alpha_Slave(sideID) )
      minAlpha(sideID) = SpacePropFactor * minAlpha(sideID)
      
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side

      !master sides(ElemID,locSide and flip =-1 if not existing)
      if(ElemID.NE.-1) then ! element belonging to master side is on this processor
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      end if
      
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      !slave side (nbElemID,nblocSide and flip =-1 if not existing)
      if(nbElemID.NE.-1) then! element belonging to slave side is on this processor
        if (alpha(nbElemID) < minAlpha(sideID)) alpha(nbElemID) = minAlpha(sideID)
      end if
    end do

!   Enforce minimum alpha on mortar sides
!   -------------------------------------
    
    ! Inner mortar sides
    DO MortarSideID=firstMortarInnerSide, lastMortarInnerSide

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      ElemID =SideToElem(S2E_ELEM_ID,MortarSideID) !element belonging to master side
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      END DO !iMortar
    END DO !MortarSideID
    
    ! MPI mortar sides
    DO MortarSideID=firstMortarMPISide,lastMortarMPISide

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      ElemID =SideToElem(S2E_ELEM_ID,MortarSideID) !element belonging to master side
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      END DO !iMortar
    END DO !MortarSideID
  end subroutine
#if MPI
!===================================================================================================================================
!> Start the MPI communication of the blending coefficient (send and receive)
!===================================================================================================================================
  subroutine Start_BlendCoeff_MPICommunication()
    use MOD_NFVSE_Vars         , only: alpha_Master, alpha_Slave
    use MOD_MPI                , only: StartReceiveMPIData,StartSendMPIData
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    
!   Start off with the receive command
!   **********************************
    ! receive the slave
    call StartReceiveMPIData(alpha_Slave , 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,1), SendID=2) ! Receive MINE (sendID=2)
    
    ! receive the master
    call StartReceiveMPIData(alpha_Master(firstSlaveSide:lastSlaveSide), 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,2), SendID=1) ! Receive YOUR  (sendID=1) 
    
    ! transfer the mortar to the right MPI faces 
    CALL Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides=.TRUE.)

    ! Send the slave
    call StartSendMPIData   (alpha_Slave , 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,3), SendID=2) ! SEND YOUR (sendID=2) 
    
    ! Send the master
    call StartSendMPIData   (alpha_Master(firstSlaveSide:lastSlaveSide), 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,4),SendID=1) 
    
  end subroutine Start_BlendCoeff_MPICommunication
#endif /*MPI*/
!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side
!==================================================================================================================================
  subroutine Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides)
    ! MODULES
    use MOD_Preproc
    use MOD_Mesh_Vars,   only: MortarType,MortarInfo
    use MOD_Mesh_Vars,   only: firstMortarInnerSide,lastMortarInnerSide
    use MOD_Mesh_Vars,   only: firstMortarMPISide,lastMortarMPISide
    use MOD_Mesh_Vars,   only: firstSlaveSide,lastSlaveSide
    use MOD_Mesh_Vars,   only: nSides 
    implicit none
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    real,INTENT(INOUT) :: alpha_Master(firstMortarInnerSide:nSides)  !< (INOUT) 
    real,INTENT(INOUT) :: alpha_Slave (firstSlaveSide:lastSlaveSide) !< (INOUT) 
    LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER      :: iMortar,nMortars
    INTEGER      :: firstMortarSideID,lastMortarSideID
    INTEGER      :: MortarSideID,SideID,locSide,flip
    !==================================================================================================================================
    IF(doMPISides)THEN
      firstMortarSideID = firstMortarMPISide
      lastMortarSideID =  lastMortarMPISide
    ELSE
      firstMortarSideID = firstMortarInnerSide
      lastMortarSideID =  lastMortarInnerSide
    END IF !doMPISides

    DO MortarSideID=firstMortarSideID,lastMortarSideID

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        flip  = MortarInfo(MI_FLIP,iMortar,locSide)
        SELECT CASE(flip)
          CASE(0) ! master side
            alpha_Master(SideID)=Alpha_Master(MortarSideID)
          CASE(1:4) ! slave side
            alpha_Slave (SideID)=Alpha_Master(MortarSideID)
        END SELECT !flip(iMortar)
      END DO !iMortar
    END DO !MortarSideID
  end subroutine Alpha_Mortar
!===================================================================================================================================
!> Gets outer solution for the TVD reconstruction
!> ATTENTION 1: Must be called after FinishExchangeMPIData
!> ATTENTION 2: Mortar faces are not considered (TODO?)
!===================================================================================================================================
  subroutine Get_externalU(nVar,U_ext,U,U_master,U_slave,fluxSign)
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: TanDirs1, TanDirs2
    use MOD_Mesh_Vars , only: nBCSides, SideToElem, S2V, firstInnerSide, lastMPISide_YOUR, nElems, nSides, firstSlaveSide, LastSlaveSide
    implicit none
    !-arguments-----------------------------------------------
    integer, intent(in) :: nVar
    real, intent(out) :: U_ext   (nVar,0:PP_N,0:PP_N,6     ,nElems)
    real, intent(in)  :: U       (nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real, intent(in)  :: U_master(nVar,0:PP_N,0:PP_N,nSides)
    real, intent(in)  :: U_slave (nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide)
    logical, optional :: fluxSign
    !-local-variables-----------------------------------------
    integer :: SideID, ElemID, locSide, p, q, ijk(3), flip, nbElemID, nblocSide, nbFlip
    integer :: ax1, ax2
    integer :: signIdx
    integer, parameter :: signo(2) = [1, -1]
    !---------------------------------------------------------
    
    if (present(fluxSign)) then
      if (fluxSign) then
        signIdx=2
      else
        signIdx=1
      end if
    else
      signIdx=1
    end if
    
    ! First assign the inner value on boundaries (we don't use an "external state")
    do SideID=1, nBCSides
      ElemID  = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip    = 0 ! flip is 0 in master faces!
      
      ax1 = TanDirs1(locSide)
      ax2 = TanDirs2(locSide)
      DO q=0,PP_N; DO p=0,PP_N
        ijk(:)=S2V(:,0,p,q,flip,locSide)
        U_ext(:,ijk(ax1),ijk(ax2),locSide,ElemID)=U(:,ijk(1),ijk(2),ijk(3),ElemID)
      END DO; END DO !p,q=0,PP_N
    end do !SideID=1, nBCSides
    
    ! TODO: Include routines for mortar faces here....
    
    ! Now do the rest
    do SideID=firstInnerSide, lastMPISide_YOUR
      
!     Master side
!     -----------
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      !master sides(ElemID,locSide and flip =-1 if not existing)
      IF(ElemID.NE.-1)THEN ! element belonging to master side is on this processor
        locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
        flip    = 0 ! flip is 0 in master faces!
        ax1 = TanDirs1(locSide)
        ax2 = TanDirs2(locSide)
        DO q=0,PP_N; DO p=0,PP_N
          ijk(:)=S2V(:,0,p,q,flip,locSide)
          U_ext(:,ijk(ax1),ijk(ax2),locSide,ElemID)=U_slave(:,p,q,SideID) ! Slave side is force-aligned with master
        END DO; END DO !p,q=0,PP_N
      end if !(ElemID.NE.-1)
      
!     Slave side
!     ----------
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      IF(nbElemID.NE.-1)THEN! element belonging to slave side is on this processor
        nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
        nbFlip    = SideToElem(S2E_FLIP,SideID)
        
        ax1 = TanDirs1(nblocSide)
        ax2 = TanDirs2(nblocSide)
        DO q=0,PP_N; DO p=0,PP_N
          ijk(:)=S2V(:,0,p,q,nbFlip,nblocSide)
          U_ext(:,ijk(ax1),ijk(ax2),nblocSide,nbElemID)=signo(signIdx)*U_master(:,p,q,SideID) ! Slave side is force-aligned with master
        END DO; END DO !p,q=0,PP_N
      end if !(nbElemID.NE.-1)
      
    end do ! SideID=firstInnerSide, lastMPISide_YOUR
    
  end subroutine Get_externalU
end module MOD_NFVSE_MPI
