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
  BCData=0.
  
  ! Fill  BC data for steady BCs (BCtype = 20, 22)
  CALL FillBCdata(0.,BCdata)
END IF !fillBC

END SUBROUTINE InitBC


!==================================================================================================================================
!> precompute and store steadyState BC data for BCtype 20,220 (Dirichlet BC like 2,22) 
!==================================================================================================================================
SUBROUTINE FillBCdata(tIn,BCdata)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: tIn                                  !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: BCdata(PP_nVar,0:PP_N,0:PP_N,nBCSides) !<boundary date filled
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
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
!> Computes the boundary values for the navierstokes part
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Mesh_Vars    ,ONLY: nSides,nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: ConstoPrim_aux,ConstoPrim,PrimtoCons
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Flux         ,ONLY: EvalEulerFlux1D
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradUx_Master,gradUy_Master,gradUz_Master
USE MOD_Flux         ,ONLY: EvalDiffFlux3D,EvalDiffFlux1D_Outflow
#endif /*PARABOLIC*/
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: tIn                                 !<current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nSides)  !<Navierstokes boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,iVar,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N),N_loc(1:3,1:3)
REAL                                 :: Prim(1:8),ar,br
REAL                                 :: U_loc(PP_nVar)
#if PARABOLIC
REAL                                 :: BCGradMat(1:3,1:3)
REAL                                 :: Fd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Gd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Hd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUx_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUy_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUz_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradUn_loc(PP_nVar),gradUt1_loc(PP_nVar),gradUt2_loc(PP_nVar)
REAL                                 :: gradRho(3,0:PP_N,0:PP_N),gradVel(3,0:PP_N,0:PP_N)
#endif /*PARABOLIC*/
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
#if PARABOLIC
                   gradUx_Master(:,:,:,SideID),gradUx_Master(:,:,:,SideID), &
                   gradUy_Master(:,:,:,SideID),gradUy_Master(:,:,:,SideID), &
                   gradUz_Master(:,:,:,SideID),gradUz_Master(:,:,:,SideID), &
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
        CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradUx_Master(:,:,:,SideID),gradUx_Master(:,:,:,SideID), &
                     gradUy_Master(:,:,:,SideID),gradUy_Master(:,:,:,SideID), &
                     gradUz_Master(:,:,:,SideID),gradUz_Master(:,:,:,SideID), &
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
        CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradUx_Master(:,:,:,SideID),gradUx_Master(:,:,:,SideID), &
                     gradUy_Master(:,:,:,SideID),gradUy_Master(:,:,:,SideID), &
                     gradUz_Master(:,:,:,SideID),gradUz_Master(:,:,:,SideID), &
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
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradUx_Master(:,:,:,SideID),gradUx_Master(:,:,:,SideID), &
                   gradUy_Master(:,:,:,SideID),gradUy_Master(:,:,:,SideID), &
                   gradUz_Master(:,:,:,SideID),gradUz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
    END DO !iSide=1,nBCloc
  CASE(4) ! Isothermal Wall, Euler = Slip Wall (see CASE(9)), Diffusion: density=inside, velocity=0, rhoE=c_v * rho_L *T_wall
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          !! Advection
          ! Compute the Euler state: tangential component of v=0, density from inside and pressure with a 1D riemann problem
          U_loc(1) = U_Master(1,p,q,SideID)
          U_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives 
          CALL ConsToPrim_aux(Prim,U_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          ! Now we compute the 1D Euler flux, but use the info that the normal component u=0
          ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0)^T
          Flux(  1,p,q,SideID) = 0.
          Flux(2:4,p,q,SideID) = Prim(5)*N_loc(:,1)
          Flux(  5,p,q,SideID) = 0.
    
          !! Diffusion
          ! For isothermal wall, all gradients are from interior
          ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall =  rho_L*C_v*Twall
          !U_Face_loc(1,p,q)   = U_Face(1,p,q) !old version, problem since density gradient already from inside
          U_Face_loc(1,p,q)   = prim(5)/RefStatePrim(BCState,5)*RefStatePrim(BCState,1) !p=(R*Twall) !pressure from outside 
          U_Face_loc(2:4,p,q) = 0.
          !U_Face_loc(5,p,q)   = U_Face_loc(1,p,q)*R*sKappaM1*Twall !old version, problem since density gradient already from insid
          U_Face_loc(5,p,q)   = prim(5)*sKappaM1 !pressure from outside 
        END DO ! p
      END DO ! q
#if PARABOLIC
      ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
      CALL EvalDiffFlux3D(Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Face_loc(:,:,:),&
                          gradUx_Master(:,:,:,SideID),                           &
                          gradUy_Master(:,:,:,SideID),                           &
                          gradUz_Master(:,:,:,SideID))                           
      ! Sum up Euler and Diffusion Flux
      DO iVar=2,PP_nVar
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID)               + & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:)
      END DO ! ivar
#endif /*PARABOLIC*/
    END DO !iSide=1,nBCLoc
  
  CASE(9) ! Euler Wall, slip wall, symmetry BC, v=0 strategy a la HALO (is very perfect)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_loc(1) = U_Master(1,p,q,SideID)
          U_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives
          CALL ConsToPrim_aux(Prim,U_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          ! Now we compute the 1D Euler flux, but use the info that the normal component u=0
          ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0)^T
          Flux(  1,p,q,SideID) = 0.
          Flux(2:4,p,q,SideID) = Prim(5)*N_loc(:,1)
          Flux(  5,p,q,SideID) = 0.
  
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
          gradUx_Face_loc(:,p,q) = BCGradMat(1,1) * gradUx_Master(:,p,q,SideID) &
                                 + BCGradMat(1,2) * gradUy_Master(:,p,q,SideID) &
                                 + BCGradMat(1,3) * gradUz_Master(:,p,q,SideID)
          gradUy_Face_loc(:,p,q) = BCGradMat(2,1) * gradUx_Master(:,p,q,SideID) &
                                 + BCGradMat(2,2) * gradUy_Master(:,p,q,SideID) &
                                 + BCGradMat(2,3) * gradUz_Master(:,p,q,SideID)
          gradUz_Face_loc(:,p,q) = BCGradMat(3,1) * gradUx_Master(:,p,q,SideID) &
                                 + BCGradMat(3,2) * gradUy_Master(:,p,q,SideID) &
                                 + BCGradMat(3,3) * gradUz_Master(:,p,q,SideID)
#endif /*PARABOLIC*/
        END DO ! p
      END DO ! q
#if PARABOLIC
      ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
      CALL EvalDiffFlux3D(Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Master(:,:,:,SideID), &
                          gradUx_Face_loc,gradUy_Face_loc,gradUz_Face_loc)
      ! Sum up Euler and Diffusion Flux
      DO iVar=2,PP_nVar
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID)                       + & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:)
      END DO ! ivar
#endif /*PARABOLIC*/
    END DO !iSide=1,nBCLoc
  CASE(10) ! outflow with gradient zero and solution intern
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! get pressure from refstate
          U_loc = U_Master(:,p,q,SideID)
          CALL ConsToPrim(Prim(1:5),U_loc)
          Prim(5) = RefStatePrim(BCState,5)
          ! U_loc contains now the state with pressure from outside (refstate)
          CALL PrimToCons(Prim(1:5),U_loc(:))
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! transform state into normal system
          U_Face_loc(1,p,q)= U_loc(1)
          U_Face_loc(5,p,q)= U_loc(5)
          U_Face_loc(2,p,q)= SUM(U_loc(2:4)*N_loc(:,1))
          U_Face_loc(3,p,q)= SUM(U_loc(2:4)*N_loc(:,2))
          U_Face_loc(4,p,q)= SUM(U_loc(2:4)*N_loc(:,3))
#if PARABOLIC
          ! for diffusion, we use the rotational invariance of the diffusion fluxes
          !   for this, we need to transform the gradients into the normal system (see GG Diss for details)
          !     First step: tranform the derivatives (dx,dy,dz) into normal system 
          gradUn_loc  = N_loc(1,1)*gradUx_Master(:,p,q,SideID) &
                      + N_loc(2,1)*gradUy_Master(:,p,q,SideID) &
                      + N_loc(3,1)*gradUz_Master(:,p,q,SideID)
          gradUt1_loc = N_loc(1,2)*gradUx_Master(:,p,q,SideID) &
                      + N_loc(2,2)*gradUy_Master(:,p,q,SideID) &
                      + N_loc(3,2)*gradUz_Master(:,p,q,SideID)
          gradUt2_loc = N_loc(1,3)*gradUx_Master(:,p,q,SideID) &
                      + N_loc(2,3)*gradUy_Master(:,p,q,SideID) &
                      + N_loc(3,3)*gradUz_Master(:,p,q,SideID)
          !     Second step: tranform the state of the gradients into normal system
   
          !TODO: actually only gradUx,y,z(1) and gradUx..(2,), gradUy..(3,), gradUz..(4) are needed
          ! maybe save some time
          !gradUx_Face_loc(1,p,q)=gradUn_loc(1)
          !gradUx_Face_loc(5,p,q)=gradUn_loc(5)
          !gradUx_Face_loc(2,p,q)=SUM(gradUn_loc(2:4)*N_loc(:,1))
          !gradUx_Face_loc(3,p,q)=SUM(gradUn_loc(2:4)*N_loc(:,2))
          !gradUx_Face_loc(4,p,q)=SUM(gradUn_loc(2:4)*N_loc(:,3))
   
          !gradUy_Face_loc(1,p,q)=gradUt1_loc(1)
          !gradUy_Face_loc(5,p,q)=gradUt1_loc(5)
          !gradUy_Face_loc(2,p,q)=SUM(gradUt1_loc(2:4)*N_loc(:,1))
          !gradUy_Face_loc(3,p,q)=SUM(gradUt1_loc(2:4)*N_loc(:,2))
          !gradUy_Face_loc(4,p,q)=SUM(gradUt1_loc(2:4)*N_loc(:,3))
   
          !gradUz_Face_loc(1,p,q)=gradUt2_loc(1)
          !gradUz_Face_loc(5,p,q)=gradUt2_loc(5)
          !gradUz_Face_loc(2,p,q)=SUM(gradUt2_loc(2:4)*N_loc(:,1))
          !gradUz_Face_loc(3,p,q)=SUM(gradUt2_loc(2:4)*N_loc(:,2))
          !gradUz_Face_loc(4,p,q)=SUM(gradUt2_loc(2:4)*N_loc(:,3))
   
          gradRho(:,p,q) = (/gradUn_loc(1),gradUt1_loc(1),gradUt2_loc(1)/)
          gradVel(:,p,q) = (/SUM(gradUn_loc(2:4) *N_loc(:,1)),&
                             SUM(gradUt1_loc(2:4)*N_loc(:,2)),&
                             SUM(gradUt2_loc(2:4)*N_loc(:,3)) /)
#endif /*PARABOLIC*/
        END DO ! p
      END DO ! q
      ! Compute 1D Euler Flux Fe_Face_loc
      CALL EvalEulerFlux1D(U_Face_loc,Flux(:,:,:,SideID))
#if PARABOLIC
      ! Compute 1D diffusion Flux Fd_Face_loc
      !   Use: tau_12 = 0, tau_13 = 0, q1 = 0 (Paper POINSOT and LELE, JCP 1992, page 113, Table IV)
      !   We use special evalflux routine
      CALL EvalDiffFlux1D_Outflow(Fd_Face_loc,U_Face_loc,gradRho,gradVel)
      ! Sum up Euler and Diffusion Flux and tranform back into Cartesian system
      Flux(:,:,:,SideID) = Flux(:,:,:,SideID) + Fd_Face_loc
#endif /*PARABOLIC*/
      DO q=0,PP_N
        DO p=0,PP_N
          Flux(2:4,p,q,SideID)=NormVec( :,p,q,SideID)*Flux(2,p,q,SideID) &
                              +TangVec1(:,p,q,SideID)*Flux(3,p,q,SideID) &
                              +TangVec2(:,p,q,SideID)*Flux(4,p,q,SideID)
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
!> Computes the boundary fluxes for the lifting procedure for all sides.
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: ConsToPrim_aux,ConsToPrim,PrimToCons
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseLiftingGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: tIn       !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nBCSides) !<lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: ar,br,N_loc(1:3,1:3),Prim(8)
REAL                                 :: U_loc(PP_nVar)
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
            CALL ExactFunc(IniExactFunc,tIn,face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /=0
      ! use BCState as refstate number
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
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
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  
  CASE(4) ! Isothermal Wall, Density from interior, velocity=0, energy = c_v * R* sKappaM1 * T_wall
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          u_loc(1) = U_Master(1,p,q,SideID)
          u_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives
          CALL ConsToPrim_aux(Prim,u_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          !  Get Wall Temperature:  Twall = RefState(BCState,5)/(R*RefState(BCState,1))
          ! Use Density from interior
          !F_Face(1,p,q)=U_Face(1,p,q) ! OLD version, problem, since density gradient already from inside
          Flux(1,p,q,SideID)=prim(5)/RefStatePrim(BCState,5)*RefStatePrim(BCState,1)!pressure from outside: rho=p/(Twall*R/rho_wall
          ! Velocity is zero, but is forced via the Riemann flux of the lifting Operator (u_riemann = 0.5(u_l + u_r)
          Flux(2:4,p,q,SideID) = 0.
          ! Energy is computed with density from inside and Twall and zero(!) velocity: rhoE_wall
          !F_Face(5,p,q) = U_Face(1,p,q)*R*sKappaM1*Twall! OLD version, problem, since density gradient already from inside 
          Flux(5,p,q,SideID) = prim(5)*sKappaM1  !pressure from outside
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(9) ! Euler Wall, slip wall, symmetry BC, v=0 strategy a la HALO (is very perfect)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          u_loc(1) = U_Master(1,p,q,SideID)
          u_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives
          CALL ConsToPrim_aux(Prim,u_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          Prim(2)=0.
          CALL PrimToCons(prim(1:5),u_loc)
          ! Compute Flux
          Flux(1,p,q,SideID)=u_loc(1)
          Flux(5,p,q,SideID)=0.5*(u_loc(5)+U_Master(5,p,q,SideID))
          ! Rotate momentum back in x,y,z coords
          Flux(2:4,p,q,SideID) = 0.5*(u_loc(2)*N_loc(:,1)+u_loc(3)*N_loc(:,2)+u_loc(4)*N_loc(:,3) &
                                 +U_Master(2:4,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(10) ! outflow with gradient zero and solution intern
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_Loc = U_Master(:,p,q,SideID)
          CALL ConsToPrim(Prim(1:5),U_Loc)
          Prim(5) = RefStatePrim(BCState,5)
          CALL PrimToCons(Prim(1:5),U_Loc)
          Flux(:,p,q,SideID) = U_Loc
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO ! iBC
!for BR1 and BR2, lifting is in strong form, flux=U-Uminus...
DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=(Flux(:,p,q,SideID)-U_Master(:,p,q,SideID))*SurfElem(p,q,SideID)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/


!==================================================================================================================================
!> Finalize arrays used for boundary conditions.
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
