#if PARABOLIC
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
!> Contains the Surface integral for the BR1 lifting scheme
!===================================================================================================================================
MODULE MOD_Lifting_SurfInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_SurfInt
  MODULE PROCEDURE BR1_SurfInt
END INTERFACE

PUBLIC::Lifting_SurfInt
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> specific implementation of the BR1 Lifting: 
!> input volume and surface gradients are only the local gradients, but already transformed to x,y,z gradients, 
!> without the Jacobian and updates only the volume gradient. Note that because of the strong volume integral, slave and master
!> updates have both a positive sign.
!>
!> It uses a specific side to volume mapping (S2V) built in mesh/mappings.f90
!> where the side index (p,q) is mapped to the 1D line (i,j,k) in the volume with an additional index l. for l=0, its the point 
!> on the 1D line that is closest to the surface (or actuall the same point for Gauss-Lobatto).
!> Example: 
!> * xi_minus side: (p,q) -> (j,k) depending on the flip and l=0,1,...N -> i=0,1,...N
!> * xi_plus side:  (p,q) -> (j,k) depending on the flip and l=0,1,...N -> i=N,N-1,...0
!> Note for Gauss: one can use the interpolation of the basis functions L_Minus(l), since the node distribution is symmetric
!===================================================================================================================================
SUBROUTINE BR1_SurfInt(Flux,gradP,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if (PP_NodeType==1)
USE MOD_DG_Vars,            ONLY: L_HatMinus
#elif (PP_NodeType==2)
USE MOD_DG_Vars,            ONLY: L_HatMinus0 
#endif /*PP_NodeType*/
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nElems,nSides
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE 
USE MOD_Mesh_Vars,          ONLY: S2V !Side To volume
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE  
REAL,INTENT(IN)    :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: gradP(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_NodeType==1)
INTEGER                :: l
#endif /*PP_NodeType*/ 
INTEGER                :: ijk(3),p,q,firstSideID,lastSideID
INTEGER                :: ElemID,locSide,SideID,flip
INTEGER                :: nbElemID,nblocSide,nbFlip
!===================================================================================================================================
IF(doMPISides)THEN
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  firstSideID = 1
   lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ! update DG time derivative with corresponding SurfInt contribution
  ! master side, flip=0 =>gradP_master
  ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side

  !master sides(ElemID,locSide and flip =-1 if not existing)
  IF(ElemID.NE.-1)THEN ! element belonging to master side is on this processor
    locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip=0
#if (PP_NodeType==1)
    DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
      ijk(:)=S2V(:,l,p,q,flip,locSide) !flip=0
      gradP(:,ijk(1),ijk(2),ijk(3),ElemID) = gradP(:,ijk(1),ijk(2),ijk(3),ElemID) +Flux(:,p,q,SideID)*L_hatMinus(l)
    END DO; END DO; END DO ! l,p,q
#elif (PP_NodeType==2)
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,flip,locSide) !flip=0, l=0
      gradP(:,ijk(1),ijk(2),ijk(3),ElemID) = gradP(:,ijk(1),ijk(2),ijk(3),ElemID) +Flux(:,p,q,SideID)*L_hatMinus0
    END DO; END DO ! p,q
#endif /*PP_NodeType*/
  END IF !master ElemID .NE. 1

  nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
  !slave side (nbElemID,nblocSide and flip =-1 if not existing)
  IF(nbElemID.NE.-1)THEN! element belonging to slave side is on this processor
    nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    nbFlip    = SideToElem(S2E_FLIP,SideID)
#if (PP_NodeType==1)
    DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
      ijk(:)=S2V(:,l,p,q,nbFlip,nblocSide) 
      gradP(:,ijk(1),ijk(2),ijk(3),nbElemID) = gradP(:,ijk(1),ijk(2),ijk(3),nbElemID) +Flux(:,p,q,SideID)*L_hatMinus(l)
    END DO; END DO; END DO ! l,p,q
#elif (PP_NodeType==2)
    DO q=0,PP_N; DO p=0,PP_N
      ijk(:)=S2V(:,0,p,q,nbFlip,nblocSide) !l=0 
      gradP(:,ijk(1),ijk(2),ijk(3),nbElemID) = gradP(:,ijk(1),ijk(2),ijk(3),nbElemID) +Flux(:,p,q,SideID)*L_hatMinus0
    END DO; END DO ! p,q
#endif /*PP_NodeType*/
  END IF !nbElemID.NE.-1
END DO ! SideID=1,nSides
END SUBROUTINE BR1_SurfInt

END MODULE MOD_Lifting_SurfInt
#endif /*PARABOLIC*/
