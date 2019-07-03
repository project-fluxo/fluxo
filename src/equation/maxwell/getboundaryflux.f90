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
!> Contains FillBoundary (which depends on the considered equation)
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

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

! Public Part ---------------------------------------------------------------------------------------------------------------------

PUBLIC::InitBC
PUBLIC::GetBoundaryFlux
PUBLIC::FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID,BCdata
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,nBCs,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iBC,SideID,BCType
LOGICAL :: fillBC
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL abort(__STAMP__,&
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

fillBC=.FALSE.
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  SELECT CASE(BCType)
  CASE(20,220)
    fillBC=.TRUE.
  END SELECT
END DO
IF(fillBC)THEN
  ! Allocate buffer array to store temp data for all BC sides
  ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
  ! Fill  BC data for steady BCs (BCtype = 20, 22)
  CALL FillBCdata(0.,BCdata)
END IF !fillBC

END SUBROUTINE InitBC


!==================================================================================================================================
!> Compute steadyState BC data for BCtype 20,220 (Dirichlet BC like 2,22)
!==================================================================================================================================
SUBROUTINE FillBCdata(t,BCdata)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t                                      !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: BCdata(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< boundary data
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iBC,iSide,p,q,SideID
INTEGER          :: BCType,BCState,nBCLoc
!==================================================================================================================================
! some BCdata is not changing over time, store into BCdata
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(20)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(220) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,t,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  END SELECT ! BCType
END DO !iBC=1,nBCs
END SUBROUTINE FillBCdata


!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!> Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Mesh_Vars    ,ONLY: nSides,nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: tIn                                !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Flux(PP_nVar,0:PP_N,0:PP_N,nSides) !< boundary flux array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iBC,iSide,p,q,SideID
INTEGER          :: BCType,BCState,nBCLoc
REAL             :: n_loc(3),resul(PP_nVar)
REAL             :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(20,220) !steadyStateBCs
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),BCdata(:,:,:,SideID), &
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(2) ! exact BC = Dirichlet BC !!
    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
  
    END DO !iSide=1,nBCloc
  CASE(3) ! 1st order absorbing BC
    U_Face_loc=0.
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(4) ! perfectly conducting wall (MunzOmnesSchneider 2000, pp. 97-98)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      ! Determine the exact BC state
      DO q=0,PP_N
        DO p=0,PP_N
          resul=U_Master(:,p,q,SideID)
          n_loc=NormVec(:,p,q,SideID)
          U_Face_loc(1:3,p,q) = -resul(1:3) + 2.*SUM(resul(1:3)*n_loc)*n_loc
          U_Face_loc(4:6,p,q) =  resul(4:6) - 2.*SUM(resul(4:6)*n_loc)*n_loc
          U_Face_loc(  7,p,q) =  resul(  7)
          U_Face_loc(  8,p,q) = -resul(  8)
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
  
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO !iBC=1,nBCs
! Integrate over the surface
DO SideID=1,nBCSides
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO ! p
  END DO! q
END DO! SideID
END SUBROUTINE GetBoundaryFlux


!==================================================================================================================================
!> Finalize boundary conditions
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC


END MODULE MOD_GetBoundaryFlux
