#include "defines.h"
!===================================================================================================================================
!> Module containing utilities for the NFVSE shock-capturing method
!> Attention: this module also contains the MPI routines for the communication of the blending coefficient in the NFVSE method:
!> -> Done here and not in DGTimeDerivative_weakForm for modularity!
!> -> Here the MPI approach is a bit different than in the main code: Both master and slave sides send and receive, then the
!>    blending coefficient correction is computed on both sides of the MPI interface. This is needed since we want to correct the DG
!>    volume integral BEFORE adding the diffusive contribution.
!> -> This does not work yet for mortar sides across MPI boundaries (big TODO!)
!===================================================================================================================================
module MOD_NFVSE_MPI
  implicit none
  
  private
  public ProlongBlendingCoeffToFaces, PropagateBlendingCoeff
contains
!===================================================================================================================================
!> Prolong the blending coefficient to the faces
!===================================================================================================================================
  subroutine ProlongBlendingCoeffToFaces()
    use MOD_ShockCapturing_Vars, only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide, SideToElem
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    !-------------------------------------------------------------------------------------------------------------------------------
    
    do sideID=firstSlaveSide, lastSlaveSide
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
    
  end subroutine ProlongBlendingCoeffToFaces
!===================================================================================================================================
!> Propagate the blending coefficient to the neighbor elements
!> -> First the routine makes sure that all the MPI communication is done
!> -> Then, the blending coefficient is set to alpha = max (alpha, alpha_neighbor)
!===================================================================================================================================
  subroutine PropagateBlendingCoeff()
    use MOD_ShockCapturing_Vars, only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide, SideToElem
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha, SpacePropFactor
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_MPI                , only: FinishExchangeMPIData
#endif /*MPI*/
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    real    :: minAlpha
    !-------------------------------------------------------------------------------------------------------------------------------
    
#if MPI
    call FinishExchangeMPIData(4*nNbProcs,MPIRequest_alpha) ! gradUx,y,z: MPI_YOUR -> MPI_MINE (_slave)
#endif /*MPI*/
    
    do sideID=firstSlaveSide, lastSlaveSide
      
      minAlpha = max( alpha_Master(sideID), alpha_Slave(sideID) )
      minAlpha = SpacePropFactor * minAlpha
      
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side

      !master sides(ElemID,locSide and flip =-1 if not existing)
      if(ElemID.NE.-1) then ! element belonging to master side is on this processor
        if (alpha(ElemID) < minAlpha) alpha(ElemID) = minAlpha
      end if
      
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      !slave side (nbElemID,nblocSide and flip =-1 if not existing)
      if(nbElemID.NE.-1) then! element belonging to slave side is on this processor
        if (alpha(nbElemID) < minAlpha) alpha(nbElemID) = minAlpha
      end if
    end do
    
  end subroutine
#if MPI
!===================================================================================================================================
!> Start the MPI communication of the blending coefficient (send and receive)
!===================================================================================================================================
  subroutine Start_BlendCoeff_MPICommunication()
    use MOD_ShockCapturing_Vars, only: alpha_Master, alpha_Slave
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
    call StartReceiveMPIData(alpha_Master, 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,2), SendID=1) ! Receive YOUR  (sendID=1) 
    
    ! TODO: transfer the mortar to the right MPI faces 
    !CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)

    ! Send the slave
    call StartSendMPIData   (alpha_Slave , 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,3), SendID=2) ! SEND YOUR (sendID=2) 
    
    ! Send the master
    call StartSendMPIData   (alpha_Master, 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,4),SendID=1) 
    
  end subroutine Start_BlendCoeff_MPICommunication
#endif /*MPI*/
  
end module MOD_NFVSE_MPI
