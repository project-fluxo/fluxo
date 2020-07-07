!==================================================================================================================================
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
!> \brief Routines providing support for geometric features of non-conforming meshes (generally a preprocessing step)
!==================================================================================================================================
MODULE MOD_Mortar_Metrics
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Mortar_CalcSurfMetrics
  MODULE PROCEDURE Mortar_CalcSurfMetrics
END INTERFACE

PUBLIC::Mortar_CalcSurfMetrics

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Interpolates surface metrics Ja and xGP on face from master to slave (small) mortar sides.
!> 1D interpolation operators M_0_1,M_0_2 are built locally per polynomial degree.
!>
!> Already existing surface metrics are overwritten, metrics for small sides are built from
!> big (master) side, i.e. all small sides belonging to a mortar interface are slave sides 
!> (with inward pointing normal vector). NOTE THAT THIS IS NOT THE CASE FOR MPI_YOUR MORTAR SIDES!
!> In an MPI setting if the big sides are not present on a CPU and this CPU has small master sides
!> they are not rebuilt and fluxes need to be rotated at the big mortar.
!>
!>~~~~~~~~~~~~~~~~~~~~
!>       Type 1               Type 2              Type3
!>        eta                  eta                 eta
!>         ^                    ^                   ^
!>         |                    |                   |
!>     +---+---+            +---+---+           +---+---+
!>     | 3 | 4 |            |   2   |           |   |   |
!>     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>     | 1 | 2 |            |   1   |           |   |   |
!>     +---+---+            +---+---+           +---+---+
!>~~~~~~~~~~~~~~~~~~~~
!>
!>
!>==================================================================================================================================
SUBROUTINE Mortar_CalcSurfMetrics(SideID,Nloc,Face_Ja,Face_xGP,&
                                  Mortar_Ja,Mortar_xGP,nbSideID)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_FillMortar,  ONLY: InterpolateBigToSmall
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                         !< SideID of mortar master side
INTEGER,INTENT(IN) :: Nloc                           !< polynomial degree
REAL,INTENT(IN)    :: Face_Ja(  3,3,0:Nloc,0:Nloc)   !< surface metrics of side
REAL,INTENT(IN)    :: Face_xGP(   3,0:Nloc,0:Nloc)   !< face xGP
REAL,INTENT(OUT)   :: Mortar_Ja(3,3,0:Nloc,0:Nloc,4) !< mortarized surface metrics of side
REAL,INTENT(OUT)   :: Mortar_xGP( 3,0:Nloc,0:Nloc,4) !< mortarized face xGP
INTEGER,INTENT(OUT):: nbSideID(4)                    !< index of neighbour sideIDs
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iNb,jNb,ind,SideIDMortar
REAL     :: Ja_small( 1:3,1:3,0:Nloc,0:Nloc,-2:4)
REAL     :: xGP_small(1:3,0:Nloc,0:Nloc,-2:4)
!==================================================================================================================================

nbSideID=-1

! Surface metrics derived from big sides are only built for inner sides and MPI_MINE sides!
SideIDMortar=MortarType(2,SideID)

CALL InterpolateBigToSmall(3*3,MortarType(1,SideID),Face_Ja ,Ja_small )
CALL InterpolateBigToSmall(  3,MortarType(1,SideID),Face_xGP,xGP_small)

SELECT CASE(MortarType(1,SideID))
CASE(1) !1->4
  DO jNb=1,2
    DO iNb=1,2
      ind=iNb+2*(jNb-1)
      IF(.NOT.(MortarInfo(E2S_FLIP,ind,SideIDMortar).GT.0)) & !no slave sides (MPI)
        nbSideID(ind)=MortarInfo(E2S_SIDE_ID,ind,SideIDMortar)
      Mortar_Ja(:,:,:,:,ind)=0.25*Ja_small(:,:,:,:,ind)
      Mortar_xGP( :,:,:,ind)=     xGP_small(:,:,:,ind)
    END DO !iNb=1,2
  END DO !jNb=1,2

CASE(2) !1->2 in eta
  DO jNb=1,2
    IF(.NOT.(MortarInfo(E2S_FLIP,jNb,SideIDMortar).GT.0)) & !no slave sides (MPI)
      nbSideID(jNb)=MortarInfo(E2S_SIDE_ID,jNb,SideIDMortar)
    Mortar_Ja(:,:,:,:,jNb)=0.5*Ja_small(:,:,:,:,jNb)
    Mortar_xGP( :,:,:,jNb)=   xGP_small(  :,:,:,jNb)
  END DO !jNb

CASE(3) !1->2 in xi
  DO iNb=1,2
    IF(MortarInfo(E2S_FLIP,iNb,SideIDMortar).GT.0) CYCLE !no slave sides (MPI)
    IF(.NOT.(MortarInfo(E2S_FLIP,iNb,SideIDMortar).GT.0)) & !no slave sides (MPI)
      nbSideID(iNb)=MortarInfo(E2S_SIDE_ID,iNb,SideIDMortar)
    Mortar_Ja(:,:,:,:,iNb)=0.5*Ja_small(:,:,:,:,iNb)
    Mortar_xGP( :,:,:,iNb)=   xGP_small(  :,:,:,iNb)
  END DO !iNb

END SELECT !MortarType
END SUBROUTINE Mortar_CalcSurfMetrics

END MODULE MOD_Mortar_Metrics
