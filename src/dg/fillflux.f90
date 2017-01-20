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
!>Fills the inner, periodic and bc fluxes
!===================================================================================================================================
SUBROUTINE FillFlux(Flux,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY: U_Master,U_Slave
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Riemann,         ONLY: Riemann
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY: gradUx_Slave,gradUx_Master,gradUy_Slave,gradUy_Master,gradUz_Slave,gradUz_Master
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  firstSideID = nBCSides+nInnerSides+1
  lastSideID  = firstSideID-1+nMPISides_MINE 
ELSE
  ! fill only InnerSides
  firstSideID = nBCSides+1
  lastSideID  = firstSideID-1+nInnerSides 
END IF
!firstSideID=nBCSides+1
!lastSideID  =nBCSides+nInnerSides+nMPISides_MINE
DO SideID=firstSideID,lastSideID
  CALL Riemann(Flux(:,:,:,SideID),     U_Master(:,:,:,SideID),     U_Slave(:,:,:,SideID), &
#if PARABOLIC
                                  gradUx_Master(:,:,:,SideID),gradUx_Slave(:,:,:,SideID), &
                                  gradUy_Master(:,:,:,SideID),gradUy_Slave(:,:,:,SideID), &
                                  gradUz_Master(:,:,:,SideID),gradUz_Slave(:,:,:,SideID), &
#endif /*PARABOLIC*/
                                  NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO ! SideID

END SUBROUTINE FillFlux


END MODULE MOD_FillFlux
