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
USE MOD_Equation_Vars,ONLY: BCdata,BCSideID
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_DG_Vars      ,ONLY: U_Master
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradUx_master,gradUy_master,gradUz_master
#endif /*PARABOLIC*/
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_Equation_Vars,ONLY: ConsToPrim,PrimToCons
USE MOD_Equation_Vars,ONLY: smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY: GLM_ch 
#endif /*PP_GLM*/
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
INTEGER                              :: iSide,p,q,SideID
REAL,DIMENSION(1:PP_nVar)            :: PrimL,PrimR
REAL                                 :: BR_n
!==================================================================================================================================
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
nBCLoc =nBCByType(iBC)
SELECT CASE(BCType)
CASE(21) !steadyStateBCs
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
CASE(29) !steadyState BC using State from mhd_equilibriumTestcase, 
         ! but with slip v*n=0, rho=rho_inner, p=p_equi, B=B_equi (B.n) /=0 !!
         !and no viscous contribution 
  DO iSide=1,nBCLoc
    SideID=BCSideID(iBC,iSide)
    DO q=0,PP_N
      DO p=0,PP_N
        ASSOCIATE(nvec =>NormVec(:,p,q,SideID))
        CALL ConsToPrim(PrimR(:),BCdata(:,p,q,SideID))
        CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
        PrimL(2:4)=PrimL(2:4) - SUM(PrimL(2:4)*nvec(:))*nvec(:) !only tangential velocities
        BR_n=SUM(PrimR(6:8)*nvec(:))
        Flux(1,  p,q,SideID) = 0.
        Flux(2:4,p,q,SideID) = (PrimR(5)+s2mu_0*SUM(PrimR(6:8)**2))*nvec(:)-smu_0*BR_n*PrimR(6:8)
        Flux(5,  p,q,SideID) =  -smu_0*BR_n*SUM(PrimR(6:8)*PrimL(2:4)) !-1/mu_0(B.v)*(B.n)
        Flux(6:8,p,q,SideID) =        -BR_n*PrimL(2:4) ! (v.n)B - (B.n)v
#ifdef PP_GLM
        Flux(9,p,q,SideID) =  GLM_ch*BR_n !0.5*GLM_ch*PrimL(9) !outflow for GLM (psi_outer=0, LF jump)
#endif /* PP_GLM */
        END ASSOCIATE !nvec
      END DO !p=0,PP_N
    END DO !q=0,PP_N
  END DO !iSide=1,nBCLoc
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
USE MOD_Equation_Vars,ONLY: BCdata,BCSideID
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Mesh_Vars    ,ONLY: NormVec
USE MOD_Equation_Vars,ONLY: ConsToPrim,PrimToCons
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
INTEGER                              :: iSide,p,q,SideID
REAL,DIMENSION(1:PP_nVar)            :: PrimL,PrimR
!==================================================================================================================================
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
nBCLoc =nBCByType(iBC)
SELECT CASE(BCType)
CASE(21)
  DO iSide=1,nBCLoc
    SideID=BCSideID(iBC,iSide)
    DO q=0,PP_N
      DO p=0,PP_N
        Flux(:,p,q,SideID)=BCdata(:,p,q,SideID)
      END DO ! p
    END DO ! q
  END DO !iSide=1,nBCloc
CASE(29) !steadyState BC using State from mhd_equilibriumTestcase, 
         ! but with slip v*n=0, rho=rho_inner, p=p_equi, B=B_equi (B.n) /=0 !!
         !and no viscous contribution 
  DO iSide=1,nBCLoc
    SideID=BCSideID(iBC,iSide)
    DO q=0,PP_N
      DO p=0,PP_N
        ASSOCIATE(nvec =>NormVec(:,p,q,SideID))
        CALL ConsToPrim(PrimR(:),BCdata(:,p,q,SideID))
        CALL ConsToPrim(PrimL(:),U_master(:,p,q,SideID))
        PrimR(1  )=PrimL(1)
        PrimR(2:4)=PrimL(2:4) - SUM(PrimL(2:4)*nvec(:))*nvec(:) !only tangential velocities
        !PrimR(5:8) = PrimR(5:8)
#ifdef PP_GLM
        PrimR(9)=0. 
#endif /* PP_GLM */
        CALL PrimToCons(PrimR(:),Flux(:,p,q,SideID))

        END ASSOCIATE !nvec
      END DO !p=0,PP_N
    END DO !q=0,PP_N
  END DO !iSide=1,nBCLoc
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in default /lifting testcase_getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE TestcaseLiftingGetBoundaryFlux
#endif /*PARABOLIC*/


END MODULE MOD_Testcase_GetBoundaryFlux
