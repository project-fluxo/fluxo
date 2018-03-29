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
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

#if PARABOLIC
INTERFACE Lifting_GetBoundaryFlux
  MODULE PROCEDURE Lifting_GetBoundaryFlux
END INTERFACE
#endif /*PARABOLIC*/

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE


PUBLIC::InitBC
PUBLIC::GetBoundaryFlux
#if PARABOLIC
PUBLIC::Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/
PUBLIC::FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iBC,SideID
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,&
     "InitBC not ready to be called or already called.")
END IF

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO SideID=1,nBCSides
  DO iBC=1,nBCs
    IF(BC(SideID).EQ.iBC) nBCByType(iBC)=nBCByType(iBC)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO SideID=1,nBCSides
  DO iBC=1,nBCs
    IF(BC(SideID).EQ.iBC)THEN
      nBCByType(iBC)=nBCByType(iBC)+1
      BCSideID(iBC,nBCByType(iBC))=SideID
    END IF
  END DO
END DO

END SUBROUTINE InitBC



!==================================================================================================================================
!> Computes the boundary values for all sides
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCsideID
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_master
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradPx_master,gradPy_master,gradPz_master
#endif /*PARABOLIC*/
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseGetBoundaryFlux
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)         :: tIn       !< current time (provided by time integration scheme)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)        :: Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides) !< boundary flux 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iBC,BCType,BCstate
INTEGER                 :: p,q,iSide,SideID,nBCloc
REAL                    :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
#if PARABOLIC
REAL                    :: gradPx_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL                    :: gradPy_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL                    :: gradPz_Face_loc(PP_nVar,0:PP_N,0:PP_N)
#endif /*PARABOLIC*/
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) ! exact BC = Dirichlet BC !!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                   gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                   gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(22) ! exact BC = Dirichlet BC !!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tin,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                   gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                   gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(3) !Neumann Flux(U_inside)
#if PARABOLIC
    gradPx_Face_loc(:,:,:) = 0. 
    gradPy_Face_loc(:,:,:) = 0. 
    gradPz_Face_loc(:,:,:) = 0. 
#endif /*PARABOLIC*/
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_master(:,:,:,SideID),         &
#if PARABOLIC
                   gradPx_Face_loc(:,:,:),gradPx_Face_loc(:,:,:), &
                   gradPy_Face_loc(:,:,:),gradPy_Face_loc(:,:,:), &
                   gradPz_Face_loc(:,:,:),gradPz_Face_loc(:,:,:), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide
  CASE(4) !Neumann Flux=0.
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      Flux(:,:,:,SideID)=0.
    END DO !iSide
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO !iBC=1,nBCs
! Integrate over the surface
DO SideID=1,nBCSides
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO! SideID
END SUBROUTINE GetBoundaryFlux


#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary fluxes for the lifting procedure for all sides.
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCsideID
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_master
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseLiftingGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: tIn                                  !< current time (provided by time integration scheme)
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides) !< lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iBC,BCType,BCstate
INTEGER                        :: p,q,iSide,SideID,nBCloc
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) ! exact BC = Dirichlet BC !!
    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,tIn,face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCstate,tIn,face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(3,4) !Neumann
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      Flux(:,:,:,SideID)=U_master(:,:,:,SideID)
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO ! iBC
!for BR1 and BR2, lifting is in strong form: Flux=Flux-U_master...
DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=(Flux(:,p,q,SideID)-U_master(:,p,q,SideID))*SurfElem(p,q,SideID)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/


!==================================================================================================================================
!> Finalize arrays used for boundary conditions.
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC

END MODULE MOD_GetBoundaryFlux
