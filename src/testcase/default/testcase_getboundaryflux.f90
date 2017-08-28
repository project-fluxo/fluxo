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

!==================================================================================================================================
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!==================================================================================================================================
MODULE MOD_Testcase_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE TestcaseGetBoundaryFlux
  MODULE PROCEDURE TestcaseGetBoundaryFlux
END INTERFACE

#if PARABOLIC
INTERFACE TestcaseLiftingGetBoundaryFlux
  MODULE PROCEDURE TestcaseLiftingGetBoundaryFlux
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::TestcaseGetBoundaryFlux
#if PARABOLIC
PUBLIC::TestcaseLiftingGetBoundaryFlux
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes the boundary values for the main equation part. HERE ONLY FLUX WITHOUT SURFACE METRIC!!!
!==================================================================================================================================
SUBROUTINE TestcaseGetBoundaryFlux(iBC,tIn,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nSides,BoundaryType
USE MOD_Equation_Vars,ONLY: nBCByType
!USE MOD_Equation_Vars,ONLY: BCdata,BCSideID
!USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2
!USE MOD_DG_Vars      ,ONLY: U_Master
!#if PARABOLIC
!USE MOD_Lifting_Vars ,ONLY: gradUx_master,gradUy_master,gradUz_master
!#endif /*PARABOLIC*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iBC       !< current BC
REAL,INTENT(IN)                      :: tIn       !<current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nSides)  !<Navierstokes boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: BCType,BCState,nBCLoc
!INTEGER                              :: iSide,p,q,SideID
!==================================================================================================================================
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
nBCLoc =nBCByType(iBC)
SELECT CASE(BCType)
!  CASE(10001)
!    DO iSide=1,nBCLoc
!      SideID=BCSideID(iBC,iSide)
!      DO q=0,PP_N
!        DO p=0,PP_N
!          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
!        END DO ! p
!      END DO ! q
!      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
!#if PARABOLIC
!                   gradUx_Master(:,:,:,SideID),gradUx_Master(:,:,:,SideID), &
!                   gradUy_Master(:,:,:,SideID),gradUy_Master(:,:,:,SideID), &
!                   gradUz_Master(:,:,:,SideID),gradUz_Master(:,:,:,SideID), &
!#endif /*PARABOLIC*/
!                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
!    END DO !iSide=1,nBCloc
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in default testcase_getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE TestcaseGetBoundaryFlux


#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary fluxes for the lifting procedure for all sides. HERE ONLY FLUX WITHOUT SURFACE METRICS!
!==================================================================================================================================
SUBROUTINE TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
! 
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,BoundaryType
USE MOD_Equation_Vars,ONLY: nBCByType
!USE MOD_Equation_Vars,ONLY: BCdata,BCSideID
!USE MOD_DG_Vars      ,ONLY: U_Master
!USE MOD_Mesh_Vars    ,ONLY: NormVec
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iBC       !< current BC
REAL,INTENT(IN)                      :: tIn       !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nBCSides) !<lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: BCType,BCState,nBCLoc
!INTEGER                              :: iSide,p,q,SideID
!==================================================================================================================================
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
nBCLoc =nBCByType(iBC)
SELECT CASE(BCType)
!CASE(10001)
!  DO iSide=1,nBCLoc
!    SideID=BCSideID(iBC,iSide)
!    DO q=0,PP_N
!      DO p=0,PP_N
!        CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
!      END DO ! p
!    END DO ! q
!  END DO !iSide=1,nBCloc
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in default /lifting testcase_getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE TestcaseLiftingGetBoundaryFlux
#endif /*PARABOLIC*/


END MODULE MOD_Testcase_GetBoundaryFlux
