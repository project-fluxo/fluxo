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
!> Fills the flux of all boundary faces
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
SUBROUTINE FillBCdata(tIn,BCdata)
! MODULES
USE MOD_PreProc
USE MOD_Globals           ,ONLY: Abort
USE MOD_Mesh_Vars         ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars         ,ONLY: Face_xGP
USE MOD_Equation          ,ONLY: ExactFunc
USE MOD_Equation_Vars     ,ONLY: RefStateCons
USE MOD_Equation_Vars     ,ONLY: IniExactFunc
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: tIn       !< evaluation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: BCdata(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< data on BC face
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iBC,iSide,p,q,SideID
INTEGER                   :: BCType,BCState,nBCLoc
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
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /= 0
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            BCdata(:,p,q,SideID) = RefStateCons(BCState,:)
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    END IF !BCState=0
  CASE(220) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  END SELECT ! BCType
END DO !iBC=1,nBCs
END SUBROUTINE FillBCdata


!==================================================================================================================================
!> Computes the flux for a all boundary faces
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_DG_Vars      ,ONLY: U_master
USE MOD_Mesh_Vars    ,ONLY: nSides,nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: s2mu_0
USE MOD_Equation_Vars,ONLY: ConsToPrim
USE MOD_Equation_Vars,ONLY: FastestWave1D
USE MOD_Equation_Vars,ONLY: FastestWave1D_Roe
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY: GLM_ch 
#endif /*PP_GLM*/
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradUx_master,gradUy_master,gradUz_master
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
#endif /*PARABOLIC*/
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: tIn       !< evaluation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nSides) !< boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL,DIMENSION(1:PP_nVar)            :: PrimL
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
      CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),BCdata(:,:,:,SideID), &
#if PARABOLIC
                   gradUx_master(:,:,:,SideID),gradUx_master(:,:,:,SideID), &
                   gradUy_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), &
                   gradUz_master(:,:,:,SideID),gradUz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(2)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
          END DO ! p
        END DO ! q
        CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradUx_master(:,:,:,SideID),gradUx_master(:,:,:,SideID), &
                     gradUy_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), &
                     gradUz_master(:,:,:,SideID),gradUz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
      END DO !iSide=1,nBCloc
    ELSE !BCstate /= 0
      DO q=0,PP_N
        DO p=0,PP_N
          U_Face_loc(:,p,q) = RefStateCons(BCState,:)
        END DO ! p
      END DO ! q
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradUx_master(:,:,:,SideID),gradUx_master(:,:,:,SideID), &
                     gradUy_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), &
                     gradUz_master(:,:,:,SideID),gradUz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
      END DO !iSide=1,nBCloc
    END IF !BCState=0
  
  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradUx_master(:,:,:,SideID),gradUx_master(:,:,:,SideID), &
                   gradUy_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), &
                   gradUz_master(:,:,:,SideID),gradUz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(9) ! Euler Wall, slip wall, perfectly conducting, symmetry BC, vn=0 Bn=0
          ! pressure and density from inside,no viscous contributions
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ASSOCIATE(nvec =>NormVec(:,p,q,SideID))
          CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
          PrimL(2:4)=PrimL(2:4) - SUM(PrimL(2:4)*nvec(:))*nvec(:) !only tangential velocities
          PrimL(6:8)=PrimL(6:8) - SUM(PrimL(6:8)*nvec(:))*nvec(:) !only tangential magn. field
          Flux(1,  p,q,SideID) = 0.
          Flux(2:4,p,q,SideID) = (PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2))*nvec(:)
          Flux(5,  p,q,SideID) = 0. !(E+phat)*(v.n)-1/mu_0(B.v)*(B.n)
          Flux(6:8,p,q,SideID) = 0.   ! (v.n)B - (B.n)v
      
#ifdef PP_GLM
          Flux(9,p,q,SideID) = 0.5*GLM_ch*PrimL(9) !outflow for psi
#endif /*PP_GLM*/
          END ASSOCIATE !nvec
        
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCLoc
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
!> Computes the lifting boudnary fluxes
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_master
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: ConsToPrim,PrimToCons
USE MOD_Equation_Vars,ONLY: FastestWave1D
USE MOD_Equation_Vars,ONLY: FastestWave1D_Roe
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseLiftingGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: tIn       !< evaluation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: Flux(PP_nVar,0:PP_N,0:PP_N,1:nBCSides) !< lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(1:PP_nVar)            :: PrimL,PrimR
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(20,220)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          !Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+BCdata(:,p,q,SideID))
          Flux(:,p,q,SideID)=BCdata(:,p,q,SideID)
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(2)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    ! Dirichlet means that we use the gradients from inside the grid cell
    ! BR1 uses arithmetic mean value of states for the Riemann flux
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            !CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:))
            !Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+U_Face_loc(:))
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /=0
      ! use BCState as refstate number
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
              !Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+RefStateCons(BCState,:))
              Flux(:,p,q,SideID)=RefStateCons(BCState,:)
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    END IF !BCstate=0
  
  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    ! Dirichlet means that we use the gradients from inside the grid cell
    ! BR1 uses arithmetic mean value of states for the Riemann flux
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          !CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:))
          !Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+U_Face_loc(:))
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  
  CASE(9) ! Euler Wall, slip wall, perfectly conducting, symmetry BC, vn=0 Bn=0
          ! pressure and density from inside,no viscous contributions
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ASSOCIATE(nvec =>NormVec(:,p,q,SideID))
          CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
          PrimR(  1)=PrimL(1)
          PrimR(2:4)=PrimL(2:4) - SUM(PrimL(2:4)*nvec(:))*nvec(:) !only tangential velocities
          PrimR(  5)=PrimL(5)
          PrimR(6:8)=PrimL(6:8) - SUM(PrimL(6:8)*nvec(:))*nvec(:) !only tangential magn. field
#ifdef PP_GLM
          PrimR(9)=0.   
#endif /* PP_GLM */
          CALL PrimToCons(PrimR(:),Flux(:,p,q,SideID))
          END ASSOCIATE !nvec
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO ! iBC
!for BR1 and BR2, lifting is in strong form: Flux=Flux-U_master

DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=(Flux(:,p,q,SideID)-U_master(:,p,q,SideID))*SurfElem(p,q,SideID)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/


SUBROUTINE FinalizeBC()
!==================================================================================================================================
! Initialize boundary conditions
!==================================================================================================================================
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
