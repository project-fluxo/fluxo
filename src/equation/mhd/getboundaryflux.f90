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
  !CASE(21) like case(20), but for testcase mhd_equilibrium, will be initialized there
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
REAL,DIMENSION(1:PP_nVar)            :: U_L,U_R,PrimL,PrimR
REAL                                 :: pstar,Bystar2,Bzstar2
REAL                                 :: SL,cf,cf_Roe,RoeVelx
REAL                                 :: N_loc(1:3,1:3)
#if PARABOLIC
REAL                                 :: BCGradMat(1:3,1:3)
REAL                                 :: Fd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Gd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Hd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUx_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUy_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUz_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
#endif /*PARABOLIC*/
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(20,21,220) !steadyStateBCs
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
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_L(:) = U_master(:,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_L(2)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,1))
          U_L(3)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,2))
          U_L(4)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,3))
          U_L(6)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,1))
          U_L(7)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,2))
          U_L(8)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,3))
          ! use approximate HLLC solver and deduct wall pressure 
          ! symmetry state
          CALL ConsToPrim(PrimL(:),U_L(:))
          U_R(1:8)=(/U_L(1),-U_L(2),U_L(3:5),-U_L(6),U_L(7:8)/)
          PrimR(1:8)=(/PrimL(1),-PrimL(2),PrimL(3:5),-PrimL(6),PrimL(7:8)/)
#ifdef PP_GLM
          PrimR(9)=0.
          U_R(9)=0.
#endif /*PP_GLM*/
          CALL FastestWave1D(PrimL,cf) !waves left = waves right
          CALL FastestWave1D_Roe(U_L,U_R,PrimL,PrimR,RoeVelx,cf_Roe)
          SL = MIN(PrimL(2)-cf,-cf_Roe) !(roeVelx=0 for rho_l=rho_r and u_l=-u_r)
          ! SR = -SL => HLL state (SR*Ur-SL*UL+FL-FR)/(SR-SL) becomes 0.5*(UR+UL+ (FR-FL)/SL)
          ! SM=0. (since U_HLL(2)=0)
          ! MomemtumX Flux_HLLC(2) = rho_l*u_l*(u_l-SL)+p_l+s2mu_0*(-Bx_l^2 + By_l^2+Bz_l^2)
          ! and we want the MomentumX Flux = pstar+s2mu_0(Bystar^2+Bzstar^2) 
          ! if we use Bstar=BHLL 
          Bystar2=(PrimL(7)- (PrimL(2)*PrimL(7)-PrimL(3)*PrimL(6))/SL)**2
          Bzstar2=(PrimL(8)- (PrimL(2)*PrimL(8)-PrimL(4)*PrimL(6))/SL)**2
          ! so pstar=p_l+ rho_l*u_l*(u_l-SL) - s2mu0*(Bx_l^2+(Bystar^2-By_l^2)+(Bzstar^2-Bz_l^2)
          ! in addition, due to the discontinuity in Bx, we get 
          ! Flux_HLLC(6)=-SL*BxL
          pstar = PrimL(5)+U_L(2)*(PrimL(2)-SL)-s2mu_0*(PrimL(6)**2+(Bystar2-PrimL(7)**2)+(Bzstar2-PrimL(8)**2))
          ! Now we compute the 1D Euler flux, but use the info that the normal component u=0 and Bn=0
          ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0,n1*Fbx,n2*Fbx,n3*Fbx)^T
          Flux(1,  p,q,SideID) = 0.
          Flux(2:4,p,q,SideID) = (pstar+s2mu_0*(Bystar2+Bzstar2))*N_loc(:,1)
          Flux(5,  p,q,SideID) = 0.
          Flux(6:8,p,q,SideID) = (-SL*PrimL(6))*N_loc(:,1)
#ifdef PP_GLM
          Flux(9,p,q,SideID) = GLM_ch*PrimL(6)
#endif /*PP_GLM*/
        
#if PARABOLIC
          !! Diffusion
          ! We prepare the gradients and set the normal derivative to zero (symmetry condition!)
          ! BCGradMat = I - n * n^T = (gradient -normal component of gradient)
          BCGradMat(1,1) = 1. - N_loc(1,1)*N_loc(1,1)
          BCGradMat(2,2) = 1. - N_loc(2,1)*N_loc(2,1)
          BCGradMat(3,3) = 1. - N_loc(3,1)*N_loc(3,1)
          BCGradMat(1,2) = -N_loc(1,1)*N_loc(2,1)
          BCGradMat(1,3) = -N_loc(1,1)*N_loc(3,1)
          BCGradMat(3,2) = -N_loc(3,1)*N_loc(2,1)
          BCGradMat(2,1) = BCGradMat(1,2)
          BCGradMat(3,1) = BCGradMat(1,3)
          BCGradMat(2,3) = BCGradMat(3,2)
          gradUx_Face_loc(1:8,p,q) = BCGradMat(1,1) * gradUx_master(1:8,p,q,SideID) &
                                   + BCGradMat(1,2) * gradUy_master(1:8,p,q,SideID) &
                                   + BCGradMat(1,3) * gradUz_master(1:8,p,q,SideID)
          gradUy_Face_loc(1:8,p,q) = BCGradMat(2,1) * gradUx_master(1:8,p,q,SideID) &
                                   + BCGradMat(2,2) * gradUy_master(1:8,p,q,SideID) &
                                   + BCGradMat(2,3) * gradUz_master(1:8,p,q,SideID)
          gradUz_Face_loc(1:8,p,q) = BCGradMat(3,1) * gradUx_master(1:8,p,q,SideID) &
                                   + BCGradMat(3,2) * gradUy_master(1:8,p,q,SideID) &
                                   + BCGradMat(3,3) * gradUz_master(1:8,p,q,SideID)
#ifdef PP_GLM
          !gradient from inside for divcorr (dirichlet = 0 for state)
          gradUx_Face_loc(9,p,q) = gradUx_master(9,p,q,SideID)
          gradUy_Face_loc(9,p,q) = gradUx_master(9,p,q,SideID)
          gradUz_Face_loc(9,p,q) = gradUx_master(9,p,q,SideID)
#endif /*PP_GLM*/
#endif /*PARABOLIC*/
        END DO ! p
      END DO ! q
#if PARABOLIC
      ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
      CALL EvalDiffFlux3D(Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_master(:,:,:,SideID), &
                          gradUx_Face_loc,gradUy_Face_loc,gradUz_Face_loc)
      ! Sum up Euler and Diffusion Flux
      DO q=0,PP_N
        DO p=0,PP_N
        Flux(:,p,q,SideID) = Flux(:,p,q,SideID)                          + & 
                                NormVec(1,p,q,SideID)*Fd_Face_loc(:,p,q) + &
                                NormVec(2,p,q,SideID)*Gd_Face_loc(:,p,q) + &
                                NormVec(3,p,q,SideID)*Hd_Face_loc(:,p,q)
        END DO ! p
      END DO ! q
#endif /*PARABOLIC*/
    END DO !iSide=1,nBCLoc
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in getBoundaryFlux in mhd/getboundaryflux.f90!',BCType,999.)
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
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_master
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: s2mu_0
USE MOD_Equation_Vars,ONLY: ConsToPrim,PrimToCons
USE MOD_Equation_Vars,ONLY: FastestWave1D
USE MOD_Equation_Vars,ONLY: FastestWave1D_Roe
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
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
REAL                                 :: U_Face_loc(PP_nVar)
REAL,DIMENSION(1:PP_nVar)            :: U_L,U_R,PrimL,PrimR
REAL                                 :: pstar,Bystar,Bzstar
REAL                                 :: SL,cf,cf_Roe,RoeVelx
REAL                                 :: N_loc(1:3,1:3)
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(20,21,220)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+BCdata(:,p,q,SideID))
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
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:))
            Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+U_Face_loc(:))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /=0
      ! use BCState as refstate number
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
              Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+RefStateCons(BCState,:))
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
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:))
          Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+U_Face_loc(:))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  
  CASE(9) ! Euler Wall, slip wall, perfectly conducting, symmetry BC, vn=0 Bn=0
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_L(:) = U_master(:,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_L(2)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,1))
          U_L(3)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,2))
          U_L(4)   = SUM(U_master(2:4,p,q,SideID)*N_loc(:,3))
          U_L(6)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,1))
          U_L(7)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,2))
          U_L(8)   = SUM(U_master(6:8,p,q,SideID)*N_loc(:,3))
          ! use approximate HLLC solver and deduct wall pressure 
          ! symmetry state
          U_R(:)=(/U_L(1),-U_L(2),U_L(3:5),-U_L(6),U_L(7:PP_nVar)/)
          CALL ConsToPrim(PrimL(:),U_L(:))
          PrimR(:)=(/PrimL(1),-PrimL(2),PrimL(3:5),-PrimL(6),PrimL(7:PP_nVar)/)
          CALL FastestWave1D(PrimL,cf) !waves left = waves right
          CALL FastestWave1D_Roe(U_L,U_R,PrimL,PrimR,RoeVelx,cf_Roe)
          SL = MIN(PrimL(2)-cf,-cf_Roe) !(roeVelx=0 for rho_l=rho_r and u_l=-u_r)
          ! SR = -SL => HLL state (SR*Ur-SL*UL+FL-FR)/(SR-SL) becomes 0.5*(UR+UL+ (FR-FL)/SL)
          ! SM=0. (since U_HLL(2)=0)
          ! MomemtumX Flux_HLLC(2) = rho_l*u_l*(u_l-SL)+p_l+s2mu_0*(-Bx_l^2 + By_l^2+Bz_l^2)
          ! and we want the MomentumX Flux = pstar+s2mu_0(Bystar^2+Bzstar^2) 
          ! if we use Bstar=BHLL 
          Bystar=(PrimL(7)- (PrimL(2)*PrimL(7)-PrimL(3)*PrimL(6))/SL)
          Bzstar=(PrimL(8)- (PrimL(2)*PrimL(8)-PrimL(4)*PrimL(6))/SL)
          ! so pstar=p_l+ rho_l*u_l*(u_l-SL) - s2mu0*(Bx_l^2+(Bystar^2-By_l^2)+(Bzstar^2-Bz_l^2)
          ! in addition, due to the discontinuity in Bx, we get 
          ! Flux_HLLC(6)=-SL*BxL
          pstar = PrimL(5)+U_L(2)*(PrimL(2)-SL)-s2mu_0*(PrimL(6)**2+(Bystar**2-PrimL(7)**2)+(Bzstar**2-PrimL(8)**2))
          !overwrite PrimL and UL
          PrimL(2)=0.
          PrimL(5)=pstar
          PrimL(6)=0.
          PrimL(7)=Bystar
          PrimL(8)=Bzstar
          CALL PrimToCons(PrimL,U_L)
          ! Compute Flux
          Flux(1,p,q,SideID)   = U_L(1)
          Flux(5,p,q,SideID)   = 0.5*(U_L(5)+U_master(5,p,q,SideID))
          ! Rotate momentum and Bfield back in x,y,z coords
          Flux(2:4,p,q,SideID) = U_L(2)*N_loc(:,1)+U_L(3)*N_loc(:,2)+U_L(4)*N_loc(:,3)
          Flux(2:4,p,q,SideID) = 0.5 * (Flux(2:4,p,q,SideID)+U_master(2:4,p,q,SideID))
          Flux(6:8,p,q,SideID) = U_L(6)*N_loc(:,1)+U_L(7)*N_loc(:,2)+U_L(8)*N_loc(:,3)
          Flux(6:8,p,q,SideID) = 0.5 * (Flux(6:8,p,q,SideID)+U_master(6:8,p,q,SideID))
#ifdef PP_GLM
          Flux(9,p,q,SideID)   = 0.
#endif
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in Lifting_getBoundaryFlux in mhd/getboundaryflux.f90!',BCType,999.)
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
