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
#if NONCONS
USE MOD_Riemann,         ONLY: AddNonConsFlux
#endif /*NONCONS*/
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradPx_master,gradPy_master,gradPz_master
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
REAL                                 :: PrimL(1:PP_nVar),lambda_max
#ifdef NONCONS
REAL                                 :: phi_L(2:8) 
#endif
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
                   gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                   gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                   gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#if NONCONS
      CALL AddNonConsFlux(Flux(:,:,:,SideID),U_master(:,:,:,SideID),  BCdata(:,:,:,SideID), &
                       NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#endif 
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
                     gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                     gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                     gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#if NONCONS
        CALL AddNonConsFlux(Flux(:,:,:,SideID),U_master(:,:,:,SideID),     U_Face_loc, &
                         NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#endif 
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
                     gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                     gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                     gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#if NONCONS
        CALL AddNonConsFlux(Flux(:,:,:,SideID),U_master(:,:,:,SideID),     U_Face_loc, &
                         NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#endif 
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
                   gradPx_master(:,:,:,SideID),gradPx_master(:,:,:,SideID), &
                   gradPy_master(:,:,:,SideID),gradPy_master(:,:,:,SideID), &
                   gradPz_master(:,:,:,SideID),gradPz_master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#if NONCONS
        CALL AddNonConsFlux(Flux(:,:,:,SideID),U_master(:,:,:,SideID),     U_Face_loc, &
                         NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#endif 
    END DO !iSide=1,nBCloc
  CASE(9) ! inviscid slip wall, perfectly conducting, symmetry BC, vn=0 Bn=0
          ! WITHOUT VISCOUS CONTRIBUTIONS!!
          ! pressure and density from inside, but flux is exactly the Lax-Friedrichs flux with a mirrored outer state
          ! (rho, rho(v-v_n nvec), E, (B- B_n nvec),psi)
          ! shown to be entropy stable
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ASSOCIATE(nvec =>NormVec(:,p,q,SideID),t1vec=>TangVec1(:,p,q,SideID),t2vec=>TangVec2(:,p,q,SideID))
          CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
          !rotate to normal system
          PrimL(2:4)=(/SUM(PrimL(2:4)*nvec(:)),SUM(PrimL(2:4)*t1vec(:)),SUM(PrimL(2:4)*t2vec(:))/)
          PrimL(6:8)=(/SUM(PrimL(6:8)*nvec(:)),SUM(PrimL(6:8)*t1vec(:)),SUM(PrimL(6:8)*t2vec(:))/)
          CALL FastestWave1D(PrimL(:),lambda_max) !=c_f
          lambda_max=ABS(PrimL(2))+lambda_max
          Flux(1,  p,q,SideID) = 0. !no mass flux
          Flux(5,  p,q,SideID) = 0. ! no energy flux
                                 !p* = p + rho vn( vn+lambda_max) + 1/(2 mu_0) |B|^2 - 1/mu_0 B_n^2 
          Flux(2:4,p,q,SideID) = (PrimL(5)+PrimL(1)*PrimL(2)*(PrimL(2)+lambda_max)  &
                                  + s2mu_0*(PrimL(7)**2+PrimL(8)**2-PrimL(6)**2)    ) *nvec(:)
          Flux(6:8,p,q,SideID) = lambda_max*PrimL(6)*nvec(:)  ! lambda_max*B_n n

          !! EC TEST: 
          !! set numerical flux to physical flux, entropy contribution should be zero!
          !Flux(2:4,p,q,SideID) = (PrimL(5) + s2mu_0*SUM(PrimL(6:8)**2) ) *nvec(:)
          !Flux(6:8,p,q,SideID) = 0. 
      
#ifdef PP_GLM
          Flux(6:8,p,q,SideID) = Flux(6:8,p,q,SideID)+ GLM_ch*PrimL(9)*nvec(:)
          Flux(  9,p,q,SideID) = 0.
#endif /*PP_GLM*/
          
#if NONCONS
          !NONCONS: Powell term, B*nvec=0, should have no contribution here.
          !CAREFUL! the "weak" implementation of the nonconservative volume term already includes -1/2(phi B_n)^inner
          ! but since we want the non-conservative flux to be zero, another -1/2(phi B_n)^inner is missing to get full 
          ! inner surface contribution.
          phi_L(2:4)=U_master(6:8,p,q,SideID)
          phi_L(6:8)=U_master(2:4,p,q,SideID)/U_master(1,p,q,SideID)
          phi_L(5)  =SUM(phi_L(2:4)*phi_L(6:8))
          Flux(2:8,p,q,SideID)=Flux(2:8,p,q,SideID)-0.5*phi_L(2:8)*PrimL(6)   !1/2 phi B_n
#ifdef PP_GLM
          ! galilean invariance for GLM
          Flux((/5,PP_nVar/),p,q,SideID)=Flux((/5,PP_nVar/),p,q,SideID)-0.5*PrimL(2)*(/PrimL(9)**2,PrimL(9)/)  !1/2 v_n (psi^2,psi) 
#endif /*PP_GLM*/
#endif 
          END ASSOCIATE !nvec,t1vec,t2vec
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
USE MOD_DG_Vars      ,ONLY: U_master,nTotal_face
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: ConsToPrim,PrimToCons
USE MOD_Equation_Vars,ONLY: FastestWave1D
USE MOD_Equation_Vars,ONLY: FastestWave1D_Roe
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
#if PP_Lifting_Var==2
USE MOD_Equation_Vars,ONLY: ConsToPrimVec
#elif PP_Lifting_Var==3
USE MOD_Equation_Vars,ONLY: ConsToEntropyVec
#endif /*PP_Lifting_Var**/
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
INTEGER                   :: iBC,iSide,p,q,SideID
INTEGER                   :: BCType,BCState,nBCLoc
REAL,DIMENSION(1:PP_nVar) :: PrimL,PrimR
REAL                      :: P_m(PP_nVar,0:PP_N,0:PP_N) !< conservative / primitive /entropy variable
REAL                      :: Floc(PP_nVar,0:PP_N,0:PP_N) !< conservative / primitive /entropy variable
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
  
  CASE(9) ! inviscid slip wall, perfectly conducting, symmetry BC, vn=0 Bn=0
          ! use pressure and density from inside (only consistent with BC, since no viscous contribution)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ASSOCIATE(nvec =>NormVec(:,p,q,SideID))
          CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
          PrimR(  :)=PrimL(:)
          PrimR(2:4)=PrimL(2:4) - SUM(PrimL(2:4)*nvec(:))*nvec(:) !only tangential velocities
          PrimR(6:8)=PrimL(6:8) - SUM(PrimL(6:8)*nvec(:))*nvec(:) !only tangential magn. field
          CALL PrimToCons(PrimR(:),Flux(:,p,q,SideID))
          END ASSOCIATE !nvec
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO ! iBC

DO SideID=1,nBCSides
  !for BR1 and BR2, lifting is in strong form: Flux=Flux-U_master
#if PP_Lifting_Var==1
  P_m(:,:,:)=U_master(:,:,:,SideID)
  Floc(:,:,:)=Flux(:,:,:,SideID)
#elif PP_Lifting_Var==2
  !convert BC lifting flux to primitive variables
  CALL ConsToPrimVec(nTotal_face,P_m,U_master(:,:,:,SideID)) !prim_var
  CALL ConsToPrimVec(nTotal_face,Floc,Flux(:,:,:,SideID)) !prim_var
#elif PP_Lifting_Var==3
  !convert BC lifting flux to entropy variables, 
  CALL ConsToEntropyVec(nTotal_face,P_m,U_master(:,:,:,SideID)) !entropy_var
  CALL ConsToEntropyVec(nTotal_face,Floc,Flux(:,:,:,SideID)) !entropy_var
#endif /*PP_Lifting_Var**/
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=(Floc(:,p,q)-P_m(:,p,q))*SurfElem(p,q,SideID)
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
