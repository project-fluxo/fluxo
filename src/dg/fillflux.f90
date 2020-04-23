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

MODULE MOD_FillFlux
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!>Fills the inner and the mpi sides fluxes, BC fluxes are separately called in dg
!===================================================================================================================================
SUBROUTINE FillFlux(Flux_master,Flux_slave,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY: U_Master,U_Slave
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: firstSlaveSide,LastSlaveSide
USE MOD_Riemann,         ONLY: Riemann
#if NONCONS
USE MOD_Riemann,         ONLY: AddNonConsFlux
#endif /*NONCONS*/
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY: gradPx_Master,gradPy_Master,gradPz_Master
USE MOD_Lifting_Vars,    ONLY: gradPx_Slave ,gradPy_Slave ,gradPz_Slave
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: Flux_slave(1:PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
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

DO SideID=firstSideID,lastSideID
  CALL Riemann(Flux_master(:,:,:,SideID),     U_Master(:,:,:,SideID),     U_Slave(:,:,:,SideID), &
#if PARABOLIC
                                  gradPx_Master(:,:,:,SideID),gradPx_Slave(:,:,:,SideID), &
                                  gradPy_Master(:,:,:,SideID),gradPy_Slave(:,:,:,SideID), &
                                  gradPz_Master(:,:,:,SideID),gradPz_Slave(:,:,:,SideID), &
#endif /*PARABOLIC*/
                                  NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))

  !conservative flux:
  Flux_slave(:,:,:,SideID)=Flux_master(:,:,:,SideID)
#if NONCONS
  !add nonconservative fluxes
  CALL AddNonConsFlux(Flux_master(:,:,:,SideID), &
                      U_Master(:,:,:,SideID),     U_Slave(:,:,:,SideID),&
                      NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
  CALL AddNonConsFlux(Flux_slave(:,:,:,SideID), &
                      U_Slave(:,:,:,SideID) ,U_Master(:,:,:,SideID),&
                      NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#endif /*NONCONS*/
  DO q=0,PP_N
    DO p=0,PP_N
      Flux_master(:,p,q,SideID)=Flux_master(:,p,q,SideID)*SurfElem(p,q,SideID)
      Flux_slave( :,p,q,SideID)=Flux_slave( :,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO ! SideID


END SUBROUTINE FillFlux


END MODULE MOD_FillFlux
