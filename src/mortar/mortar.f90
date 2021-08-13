!==================================================================================================================================
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
! Copyright (c) 2020 - 2021 Florian Hindenlang
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
!> \brief Build mortar interpolation/projection operators for (2->1) and (1->2) non-conforming interfaces.
!> Contains the routines to initialize and finalize the mortar operator matrices M_0_1, M_0_2, etc (see module Mortar_Vars)
!==================================================================================================================================
MODULE MOD_Mortar
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersMortar 
  MODULE PROCEDURE DefineParametersMortar
END INTERFACE

INTERFACE InitMortar
  MODULE PROCEDURE InitMortar
END INTERFACE

INTERFACE InitMortarArrays
  MODULE PROCEDURE InitMortarArrays
END INTERFACE

INTERFACE FinalizeMortar
  MODULE PROCEDURE FinalizeMortar
END INTERFACE

PUBLIC::DefineParametersMortar,InitMortarBase,InitMortar,InitMortarArrays,FinalizeMortar,FinalizeMortarArrays 

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersMortar()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mortar")
CALL prms%CreateIntOption(     'whichMortar',           "0: projection mortar, 1: collocation mortar. ")
END SUBROUTINE DefineParametersMortar

!==================================================================================================================================
!> Basic Mortar initialization.
!==================================================================================================================================
SUBROUTINE InitMortarBase()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Interpolation     ,ONLY: getNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone,NodeType
USE MOD_Basis,             ONLY: buildLegendreVdm 
USE MOD_Mortar_Vars
USE MOD_ReadInTools, ONLY: GETINT
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
REAL                          :: error
REAL,DIMENSION(0:PP_N,0:PP_N) :: Vdm_Leg,sVdm_Leg
REAL,DIMENSION(0:PP_N)        :: test1,test2,xi_GP,w_GP
INTEGER                       :: whichMortar
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MORTAR...'
IF(MortarInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,&
     'InitMortar not ready to be called or already called.')
END IF
!index 1:2/1:4 interpolation to small sides, index -2:-1 intermediate interpolation, index 0: big side
! DG interfaces
ALLOCATE(MInt(0:PP_N,0:PP_N,2))
ALLOCATE(MInt_h(0:PP_N,0:PP_N,2))
ALLOCATE(MProj(0:PP_N,0:PP_N,2))
ALLOCATE(MProj_h(0:PP_N,0:PP_N,2))
CALL MortarBasis_BigToSmall(PP_N,NodeType,   Mint(:,:,1),   Mint(:,:,2))
MInt_h=0.5*MInt
#ifdef JESSE_MORTAR
whichMortar = 1  
SWRITE(UNIT_StdOut,'(A)')"Compiled with Jesse's mortar,  Mortar set to collocation!"
#else
whichMortar = GETINT('whichMortar','0')
#endif
SELECT CASE (whichMortar)
CASE(0)
  CALL MortarBasis_SmallToBig_Projection(PP_N,NodeType,   Mproj(:,:,1),   Mproj(:,:,2))
  SWRITE(UNIT_StdOut,'(A)')'Projection Mortar chosen.'
CASE(1)
  CALL MortarBasis_SmallToBig_Collocation(PP_N,NodeType,  Mproj(:,:,1),   Mproj(:,:,2))
  SWRITE(UNIT_StdOut,'(A)')'Collocation Mortar chosen.'
CASE DEFAULT
  CALL abort(__STAMP__,&
    'which mortar either 0 (projection) or 1 (collocation)')
END SELECT
Mproj_h=0.5*Mproj
  

ASSOCIATE(M_0_1=>MInt(:,:,1),M_0_2=>Mint(:,:,2),M_1_0=>MProj(:,:,1),M_2_0=>Mproj(:,:,2))

!> TODO: Make a unit test out of this one
!Test mean value property 0.5*(0.5+1.5)=1.  !ONLY GAUSS
test1=0.5
test2=1.5
CALL GetNodesAndWeights(PP_N,NodeType,xi_GP,w_GP) !Gauss nodes and integration weights
error=ABS(0.25*SUM((MATMUL(TRANSPOSE(M_1_0),test1)+MATMUL(TRANSPOSE(M_2_0),test2))*w_GP)-1.)

IF(error.GT. 100.*PP_RealTolerance) THEN
  CALL abort(__STAMP__,&
    'problems in building Mortar 1',999,error)
END IF

! freestream test
test2 = 1.33d0
test1 = MATMUL(TRANSPOSE(M_0_1),test2) !interpolate to mortar1
error = SUM(ABS(test1-1.33d0))/REAL(PP_N+1)
SWRITE(UNIT_StdOut,*)'Error of interpolate constant to mortar 1:',error
test2 = MATMUL(TRANSPOSE(M_0_2),test2) !interpolate to mortar2
error = SUM(ABS(test2-1.33d0))/REAL(PP_N+1)
SWRITE(UNIT_StdOut,*)'Error of interpolate constant to mortar 2:',error
test2 = MATMUL(TRANSPOSE(M_1_0),test1)+MATMUL(TRANSPOSE(M_2_0),test2)
error = SUM(ABS(0.5d0*test2-1.33d0))/REAL(PP_N+1)
SWRITE(UNIT_StdOut,*)'Error of project constant back to big side:',error

IF(error.GT. 100.*PP_RealTolerance) THEN
  CALL abort(__STAMP__,&
    'problems in building Mortar 2',999,error)
END IF

CALL buildLegendreVdm(PP_N,xi_GP,Vdm_Leg,sVdm_Leg)
test2 = 1.0d0
test2(0) = -10.33d0
error  = SUM(test2*w_GP)  !save mean value of big side
!SWRITE(*,*)'big side',test2
!SWRITE(*,*)'big side modes',MATMUL(sVdm_Leg,test2)
!SWRITE(*,*)'MV',SUM(test2*w_GP)

test1 = MATMUL(TRANSPOSE(M_0_1),test2) !interpolate to mortar1
test2 = MATMUL(TRANSPOSE(M_0_2),test2) !interpolate to mortar2
test2 = 0.5*(MATMUL(TRANSPOSE(M_1_0),test1)+MATMUL(TRANSPOSE(M_2_0),test2)) !project  back
!SWRITE(*,*)'big side',test2
!SWRITE(*,*)'big side modes',MATMUL(sVdm_Leg,test2)
!SWRITE(*,*)'MV',SUM(test2*w_GP)

error  = error - SUM(test2*w_GP)  !difference to initial mean value
SWRITE(UNIT_StdOut,*)'Error of mean value of polynomial, projected back to big side:',error

END ASSOCIATE !M_0_1,M_0_2,M_1_0,M_2_0

IF(error.GT. 100.*PP_RealTolerance) THEN
  CALL abort(__STAMP__,&
    'problems in building Mortar 3',999,error)
ELSE
  SWRITE(UNIT_StdOut,'(A)')'Mortar operators build successfully.'
END IF
SWRITE(UNIT_stdOut,'(A)')' INIT MORTAR DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitMortarBase

!==================================================================================================================================
!> Basic Mortar initialization.
!==================================================================================================================================
SUBROUTINE InitMortar()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY: MeshInitIsDone
USE MOD_Mortar_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
IF(MortarInitIsDone.OR.(.NOT.MeshInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,&
     'InitMortar not ready to be called or already called.')
END IF
CALL InitMortarArrays()

MortarInitIsDone=.TRUE.
END SUBROUTINE InitMortar

!==================================================================================================================================
!> Initialize arrays needed for mortars that depend on the mesh.
!==================================================================================================================================
SUBROUTINE InitMortarArrays()
! MODULES
USE MOD_Preproc
USE MOD_Globals
#ifdef JESSE_MORTAR
USE MOD_Mesh_Vars,  ONLY: nMortarSides
USE MOD_Mesh_Vars,  ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,  ONLY: NormVec,SurfElem
USE MOD_Mortar_Vars
#endif /*JESSE_MORTAR*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef JESSE_MORTAR
INTEGER      :: MortarSideID,iSide
REAL         :: Ns_loc(1:3,0:PP_N,0:PP_N)
#endif /*JESSE_MORTAR*/
!==================================================================================================================================

#ifdef JESSE_MORTAR
!index 1:2/1:4 interpolation to small sides, index -2:-1 intermediate interpolation, index 0: big side
ALLOCATE(U_small(PP_nVar,0:PP_N,0:PP_N,-2:4,nMortarSides)) 
U_small=-HUGE(1.)

ALLOCATE(delta_flux_jesse(PP_nVar,0:PP_N,0:PP_N,nMortarSides)) 
delta_flux_jesse=-HUGE(1.)

ALLOCATE(Ns_small(1:3   ,0:PP_N,0:PP_N,-2:4,nMortarSides))
Ns_small=-HUGE(1.)
DO iSide=1,nMortarSides
    MortarSideID=MortarInfo(MI_SIDEID,0,iSide)
    Ns_loc(1,:,:)=NormVec(1,:,:,MortarSideID)*SurfElem(:,:,MortarSideID)
    Ns_loc(2,:,:)=NormVec(2,:,:,MortarSideID)*SurfElem(:,:,MortarSideID)
    Ns_loc(3,:,:)=NormVec(3,:,:,MortarSideID)*SurfElem(:,:,MortarSideID)
    CALL InterpolateBigToSmall_ALL(3,MortarType(1,MortarSideID),Ns_loc,Ns_small(:,:,:,:,iSide))
END DO !nMortarSides
#endif /*JESSE_MORTAR*/
END SUBROUTINE InitMortarArrays

!==================================================================================================================================
!> Deallocate mortar interpolation matrices.
!==================================================================================================================================
SUBROUTINE FinalizeMortar()
! MODULES
USE MOD_Mortar_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(Mint)
SDEALLOCATE(Mint_h)
SDEALLOCATE(Mproj)
SDEALLOCATE(Mproj_h)
CALL FinalizeMortarArrays()

END SUBROUTINE FinalizeMortar

!==================================================================================================================================
!> Deallocate mortar arrays that depend on the mesh.
!==================================================================================================================================
SUBROUTINE FinalizeMortarArrays()
! MODULES
USE MOD_Mortar_Vars
IMPLICIT NONE
!==================================================================================================================================
#ifdef JESSE_MORTAR
SDEALLOCATE(U_small)
SDEALLOCATE(Ns_small)
SDEALLOCATE(delta_flux_jesse)
#endif /*JESSE_MORTAR*/

END SUBROUTINE FinalizeMortarArrays

END MODULE MOD_Mortar
