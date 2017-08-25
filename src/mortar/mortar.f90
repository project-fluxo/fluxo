!==================================================================================================================================
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
!> \brief Build mortar interpolation/projection operators for (2->1) and (1->2) non-conforming interfaces.
!> Contains the routines to initialize and finalize the mortar operator matrices M_0_1, M_0_2, etc (see module Mortar_Vars)
!==================================================================================================================================
MODULE MOD_Mortar
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitMortar
  MODULE PROCEDURE InitMortar
END INTERFACE

INTERFACE MortarBasis_BigToSmall
  MODULE PROCEDURE MortarBasis_BigToSmall
END INTERFACE

INTERFACE MortarBasis_SmallToBig
  MODULE PROCEDURE MortarBasis_SmallToBig_Collocation
END INTERFACE

INTERFACE FinalizeMortar
  MODULE PROCEDURE FinalizeMortar
END INTERFACE

PUBLIC::InitMortar,FinalizeMortar,MortarBasis_BigToSmall,MortarBasis_SmallToBig

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Basic Mortar initialization.
!==================================================================================================================================
SUBROUTINE InitMortar()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Interpolation     ,ONLY: getNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone,NodeType
USE MOD_Mortar_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
REAL                          :: error
REAL,DIMENSION(0:PP_N)        :: test1,test2,xi_GP,w_GP
!==================================================================================================================================
IF(MortarInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,&
     'InitMortar not ready to be called or already called.')
END IF

! DG interfaces
ALLOCATE(M_0_1(0:PP_N,0:PP_N))
ALLOCATE(M_0_2(0:PP_N,0:PP_N))
ALLOCATE(M_1_0(0:PP_N,0:PP_N))
ALLOCATE(M_2_0(0:PP_N,0:PP_N))
CALL MortarBasis_BigToSmall(PP_N,NodeType,   M_0_1,   M_0_2)
CALL MortarBasis_SmallToBig(PP_N,NodeType,   M_1_0,   M_2_0)

!> TODO: Make a unit test out of this one
!Test mean value property 0.5*(0.5+1.5)=1.  !ONLY GAUSS
test1=0.5
test2=1.5
CALL GetNodesAndWeights(PP_N,NodeType,xi_GP,w_GP) !Gauss nodes and integration weights
error=ABS(0.25*SUM((MATMUL(TRANSPOSE(M_1_0),test1)+MATMUL(TRANSPOSE(M_2_0),test2))*w_GP)-1.)

IF(error.GT. 100.*PP_RealTolerance) THEN
  CALL abort(__STAMP__,&
    'problems in building Mortar',999,error)
ELSE
  SWRITE(UNIT_StdOut,'(A)')'Mortar operators build successfully.'
END IF

MortarInitIsDone=.TRUE.
END SUBROUTINE InitMortar


!==================================================================================================================================
!> Build 1D operators for non-conforming interfaces:
!>    M_0_1(:,:)  interpolation from full  interval 0: [-1,1] to left  interval 1: [-1,0]
!>    M_0_2(:,:)  interpolation from full  interval 0: [-1,1] to right interval 2: [0, 1]
!> see doc/mortar for details...
!==================================================================================================================================
SUBROUTINE MortarBasis_BigToSmall(N_In,NodeType_In,M_0_1,M_0_2)
! MODULES
USE MOD_Basis,             ONLY: InitializeVandermonde
USE MOD_Interpolation     ,ONLY: getNodesAndWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                      :: N_In        !< polynomial degree
CHARACTER(LEN=255),INTENT(IN)                           :: NodeType_In !< nodetype
!> precomputed mortar matrices: big to small
REAL,DIMENSION(0:N_In,0:N_in),INTENT(OUT)               :: M_0_1,M_0_2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in)        :: xi_In,wBary_In
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_In,xi_In,wIPBary=wBary_In)

!build interpolation operators M 0->1,M 0->2
CALL InitializeVandermonde(N_In,N_In,wBary_In,xi_In,0.5*(xi_In-1.),M_0_1)
CALL InitializeVandermonde(N_In,N_In,wBary_In,xi_In,0.5*(xi_In+1.),M_0_2)

! later the transposed version is mostly used
! ATTENTION: MortarBasis_BigToSmall computes the transposed matrices, which is useful when they are used
!            in hand-written matrix multiplications. For the use with the intrinsic MATMUL, they must be transposed.
M_0_1=TRANSPOSE(M_0_1)
M_0_2=TRANSPOSE(M_0_2)
END SUBROUTINE MortarBasis_BigToSmall



!==================================================================================================================================
!> Build 1D operators for non-conforming interfaces:
!>    M_1_0(:,:)  discrete projection    from left  interval 1: [-1,0] to full  interval 0: [-1,1]
!>    M_2_0(:,:)  discrete projection    from right interval 1: [0, 1] to full  interval 0: [-1,1]
!> integrals are approximated with the same collocation points (gauss/gauss-Lob) as the scheme. for gauss, collocation is exact.
!> for GL, collocation makes scheme FSP and angular momentum conservation on cartesian mortars (exact projection doesnt!!!). 
!==================================================================================================================================
SUBROUTINE MortarBasis_SmallToBig_Collocation(N_In,NodeType_In,M_1_0,M_2_0)
! MODULES
USE MOD_Interpolation     ,ONLY: getNodesAndWeights
USE MOD_Basis,             ONLY: InitializeVandermonde 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                      :: N_In  !< polynomial degree
CHARACTER(LEN=255),INTENT(IN)                           :: NodeType_In !< nodetype
!> precomputed mortar matrices: small to big
REAL,DIMENSION(0:N_In,0:N_in),INTENT(OUT)  :: M_1_0,M_2_0
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,j
REAL,DIMENSION(0:N_in)        :: xi_In,w_in,wBary_In
REAL,DIMENSION(0:N_in,0:N_in) :: Vdm_0_1,Vdm_0_2
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_In,xi_In,wIP=w_in,wIPBary=wBary_In)

! Vdm^{01}_{ij} = l^0_j (xi^1_i)
! Vdm^{02}_{ij} = l^0_j (xi^2_i)
CALL InitializeVandermonde(N_In,N_In,wBary_In,xi_In,0.5*(xi_In-1.),Vdm_0_1)
CALL InitializeVandermonde(N_In,N_In,wBary_In,xi_In,0.5*(xi_In+1.),Vdm_0_2)

!M^{10}_{ij} = w_j/w_i l^0_i(xi^1_j) !*0.5 accounted in shat
!M^{20}_{ij} = w_j/w_i l^0_i(xi^2_j) !*0.5 accounted in shat

DO i=0,N_in
  DO j=0,N_in
    M_1_0(i,j)=w_in(j)/w_in(i)*Vdm_0_1(j,i)
    M_2_0(i,j)=w_in(j)/w_in(i)*Vdm_0_2(j,i)
  END DO
END DO
! later the transposed version is mostly used
! ATTENTION: MortarBasis_SmallToBig computes the transposed matrices, which is useful when they are used
!            in hand-written matrix multiplications. For the use with the intrinsic MATMUL, they must be transposed.
M_1_0=TRANSPOSE(M_1_0)
M_2_0=TRANSPOSE(M_2_0)
END SUBROUTINE MortarBasis_SmallToBig_Collocation

!==================================================================================================================================
!> Build 1D operators for non-conforming interfaces:
!>    M_1_0(:,:)  projection    from left  interval 1: [-1,0] to full  interval 0: [-1,1]
!>    M_2_0(:,:)  projection    from right interval 1: [0, 1] to full  interval 0: [-1,1]
!> see doc/mortar for details...
!==================================================================================================================================
SUBROUTINE MortarBasis_SmallToBig_projection(N_In,NodeType_In,M_1_0,M_2_0)
! MODULES
USE MOD_Basis,             ONLY: LegendrePolynomialAndDerivative
USE MOD_Interpolation     ,ONLY: getNodesAndWeights,GetVandermonde
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                      :: N_In  !< polynomial degree
CHARACTER(LEN=255),INTENT(IN)                           :: NodeType_In !< nodetype
!> precomputed mortar matrices: small to big
REAL,DIMENSION(0:N_In,0:N_in),INTENT(OUT)  :: M_1_0,M_2_0
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: dummy
INTEGER                       :: i,j
REAL,DIMENSION(0:N_in,0:N_in) :: VGP,W,Vphi1,Vphi2,Vdm_Leg
REAL,DIMENSION(0:N_in)        :: xi_In,xi_Gauss,w_Gauss  ! Gauss Nodes
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_In,xi_In)
CALL GetNodesAndWeights(N_in,'GAUSS',xi_Gauss,w_Gauss) !Gauss nodes and integration weights

!build projection operators M 1->0,M 2->0
!evaluate nodal basis (depends on NodeType, for Gauss: unity matrix)
CALL GetVandermonde(N_In,NodeType_In,N_In,'GAUSS',VGP)

!multiply with W_jj=wGP_j
W=0.
DO i=0,N_In
  W(i,i)=w_Gauss(i)
END DO
VGP=MATMUL(W,VGP)
!compute the Vandermonde on xGP (Depends on NodeType)
DO i=0,N_In
  DO j=0,N_In
    CALL LegendrePolynomialAndDerivative(j,xi_In(i),Vdm_Leg(i,j),dummy)
    CALL LegendrePolynomialAndDerivative(j,0.5*(xi_Gauss(i)-1.),Vphi1(i,j),dummy) ! evaluate Legendre in [-1,0]
    CALL LegendrePolynomialAndDerivative(j,0.5*(xi_Gauss(i)+1.),Vphi2(i,j),dummy) ! evaluate Legendre in [ 0,1]
  END DO !i
END DO !j
! final Mortar: Vphi1
M_1_0=MATMUL(Vdm_Leg,MATMUL(TRANSPOSE(Vphi1),VGP))
M_2_0=MATMUL(Vdm_Leg,MATMUL(TRANSPOSE(Vphi2),VGP))

! later the transposed version is mostly used
! ATTENTION: MortarBasis_SmallToBig computes the transposed matrices, which is useful when they are used
!            in hand-written matrix multiplications. For the use with the intrinsic MATMUL, they must be transposed.
M_1_0=TRANSPOSE(M_1_0)
M_2_0=TRANSPOSE(M_2_0)
END SUBROUTINE MortarBasis_SmallToBig_projection

!==================================================================================================================================
!> Deallocate mortar interpolation matrices.
!==================================================================================================================================
SUBROUTINE FinalizeMortar()
! MODULES
USE MOD_Mortar_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(M_0_1)
SDEALLOCATE(M_0_2)
SDEALLOCATE(M_1_0)
SDEALLOCATE(M_2_0)
END SUBROUTINE FinalizeMortar

END MODULE MOD_Mortar
