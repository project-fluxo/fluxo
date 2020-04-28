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
!> Fills the inner and MPI side fluxes
!===================================================================================================================================
MODULE MOD_Lifting_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_FillFlux
  MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC
  MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

PUBLIC::Lifting_FillFlux
PUBLIC::Lifting_FillFlux_BC
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Computes the BR2 Surface Fluxes at inner faces for all x,y,z directions, since the flux is in strong form, 
!>  ( 1/2(U_outer+U_inner) - U_inner) *sHat * nvec_x/y/z .
!> Note that for both slave and master updates, the sign of this flux does not change!
!===================================================================================================================================
SUBROUTINE Lifting_FillFlux(FluxX,FluxY,FluxZ,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Flux,            ONLY: EvalLiftingSurfFlux
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: FluxX(1:PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: FluxY(1:PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: FluxZ(1:PP_nVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
REAL               :: F_loc(1:PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
  lastSideID  = lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
  lastSideID  = lastInnerSide
END IF

DO SideID = firstSideID,lastSideID
  CALL EvalLiftingSurfFlux(SideID,F_loc) !strong flux (1/2(U_s+U_m) -U_m)*surfElem
  DO q=0,PP_N; DO p=0,PP_N
    FluxX(:,p,q,SideID)=F_loc(:,p,q)*NormVec(1,p,q,SideID)
    FluxY(:,p,q,SideID)=F_loc(:,p,q)*NormVec(2,p,q,SideID)
    FluxZ(:,p,q,SideID)=F_loc(:,p,q)*NormVec(3,p,q,SideID)
  END DO; END DO !p,q
END DO ! SideID

END SUBROUTINE Lifting_FillFlux



!===================================================================================================================================
!> Computes the BR2 Surface Fluxes at boundary faces for all x,y,z directions
!> boundary flux must always be computed in strong form
!> surfelem contribution is considered as well
!===================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(t,FluxX,FluxY,FluxZ)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec
USE MOD_Mesh_Vars,       ONLY: nBCSides
USE MOD_GetBoundaryFlux, ONLY: Lifting_GetBoundaryFlux
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: FluxX(PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(OUT)   :: FluxY(PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(OUT)   :: FluxZ(PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q
!===================================================================================================================================
! fill flux for boundary sides (surface element is already added in getboundaryflux
CALL Lifting_GetBoundaryFlux(t,FluxZ) ! uses the current U_master at the boundary , includes surfelem

DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    FluxX(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(1,p,q,SideID)
    FluxY(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(2,p,q,SideID)
    FluxZ(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(3,p,q,SideID)
  END DO; END DO !p,q
END DO !SideID
END SUBROUTINE Lifting_FillFlux_BC

END MODULE MOD_Lifting_FillFlux
#endif /* PARABOLIC */
