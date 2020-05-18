!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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
!> Contains routines to interpolate the interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace
IMPLICIT NONE
PRIVATE

INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFace

CONTAINS

!==================================================================================================================================
!> Interpolates the interior volume element data to surface data 
!> This is one of the main routines depending on the choice of Gauss / Gauss-Lobatto nodes
!> It uses a specific side to volume mapping (S2V) built in mesh/mappings.f90
!> where the side index (p,q) is mapped to the 1D line (i,j,k) in the volume with an additional index l. for l=0, its the point 
!> on the 1D line that is closest to the surface (or actuall the same point for Gauss-Lobatto).
!> Example: 
!> * xi_minus side: (p,q) -> (j,k) depending on the flip and l=0,1,...N -> i=0,1,...N
!> * xi_plus side:  (p,q) -> (j,k) depending on the flip and l=0,1,...N -> i=N,N-1,...0
!> Note for Gauss: one can use the interpolation of the basis functions L_Minus(l), since the node distribution is symmetric
!==================================================================================================================================
SUBROUTINE ProlongToFace(nVars,Uvol,Uface_master,Uface_slave,doMPISides)
! MODULES
USE MOD_Preproc          
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_MINE, nSides
USE MOD_Mesh_Vars,          ONLY: firstSlaveSide,LastSlaveSide
USE MOD_Mesh_Vars,          ONLY: S2V  !magic mapping of side to volume
#if (PP_NodeType==1)
USE MOD_Interpolation_Vars, ONLY: L_Minus
#endif /*PP_NodeType*/ 
#if SHOCK_ARTVISC && mhd
USE MOD_ShockCapturing_Vars,ONLY:nu_Master,nu_Slave,nu
#endif /*SHOCK_ARTVISC && mhd*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
integer,intent(in)    :: nVars
LOGICAL,INTENT(IN)    :: doMPISides  !< either fill MPI sides (=.true.) or local sides (=.false.)
REAL,INTENT(IN)       :: Uvol (nVars,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< input volume data
REAL,INTENT(INOUT)    :: Uface_master(nVars,0:PP_N,0:PP_N,1:nSides)!< output master side data
REAL,INTENT(INOUT)    :: Uface_slave( nVars,0:PP_N,0:PP_N,FirstSlaveSide:LastSlaveSide)!< output slave side data
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_NodeType==1)
INTEGER               :: l
#endif /*PP_NodeType*/
INTEGER               :: ijk(3),p,q,firstSideID,lastSideID
INTEGER               :: ElemID,locSide,SideID,flip
INTEGER               :: nbElemID,nblocSide,nbFlip
!==================================================================================================================================
IF(doMPISides)THEN
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  firstSideID = 1
   lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side

  !master sides(ElemID,locSide and flip =-1 if not existing)
  IF(ElemID.NE.-1)THEN ! element belonging to master side is on this processor
    locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip=0
#if (PP_NodeType==1)
    !gauss nodes
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,flip,locSide) !maps from first point 0,p,q to i,j,k line of the volume (master: flip=0)
      Uface_master(:,p,q,SideID)=Uvol(:,ijk(1),ijk(2),ijk(3),ElemID)*L_Minus(0)
      DO l=1,PP_N
        ijk(:)=S2V(:,l,p,q,flip,locSide) !0: flip=0
        Uface_master(:,p,q,SideID)=Uface_master(:,p,q,SideID) + &
                                   Uvol(:,ijk(1),ijk(2),ijk(3),ElemID)*L_Minus(l)
      END DO !l=1,PP_N
    END DO; END DO !p,q=0,PP_N
#elif (PP_NodeType==2)
    !gauss-lobatto nodes
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,flip,locSide)
      Uface_master(:,p,q,SideID)=Uvol(:,ijk(1),ijk(2),ijk(3),ElemID)
    END DO; END DO !p,q=0,PP_N
#endif /*PP_NodeType*/
#if SHOCK_ARTVISC && mhd
   nu_Master(SideID)=nu(ElemID) 
#endif /*SHOCK_ARTVISC && mhd*/
  END IF !master ElemID > 0

  nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
  !slave side (nbElemID,nblocSide and flip =-1 if not existing)
  IF(nbElemID.NE.-1)THEN! element belonging to slave side is on this processor
    nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    nbFlip      = SideToElem(S2E_FLIP,SideID)
#if (PP_NodeType==1)
    !gauss nodes
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,nbFlip,nblocSide)
      Uface_slave(:,p,q,SideID)=Uvol(:,ijk(1),ijk(2),ijk(3),nbElemID)*L_Minus(0)
      DO l=1,PP_N
        ijk(:)=S2V(:,l,p,q,nbFlip,nblocSide)
        Uface_slave(:,p,q,SideID)=Uface_slave(:,p,q,SideID) + &
                                   Uvol(:,ijk(1),ijk(2),ijk(3),nbElemID)*L_Minus(l)
      END DO !l=1,PP_N
    END DO; END DO !p,q=0,PP_N
#elif (PP_NodeType==2)
    !gauss-lobatto nodes
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,nbFlip,nblocSide)
      Uface_slave(:,p,q,SideID)=Uvol(:,ijk(1),ijk(2),ijk(3),nbElemID)
    END DO; END DO !p,q=0,PP_N
#endif /*PP_NodeType*/
#if SHOCK_ARTVISC && mhd
   nu_Slave(SideID)=nu(nbElemID) 
#endif /*SHOCK_ARTVISC && mhd*/
  END IF !slave nbElemID > 0
END DO !SideID=firstSideID,lastSideID

END SUBROUTINE ProlongToFace


END MODULE MOD_ProlongToFace
