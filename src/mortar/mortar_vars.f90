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
!> Variables used for mortars: mortar interpolation and projection matrices
!==================================================================================================================================
MODULE MOD_Mortar_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

REAL,ALLOCATABLE :: Mint(:,:,:)          !< 1D-Mortar Operator: interpolation full interval 0: [-1,1] to left interval 1: [-1,0] 
                                         !< and right intervall 2: [0,1], size is (0:N,0:N,1:2), 
                                        
REAL,ALLOCATABLE :: Mint_h(:,:,:)        !< Mint*0.5
REAL,ALLOCATABLE :: Mproj(:,:,:)         !< 1D-Mortar Operator: projection left interval 1: [-1,0] 
                                         !< and right intervall 2: [0,1] to full intervall 0: [-1,1],size is (0:N,0:N,1:2)
REAL,ALLOCATABLE :: Mproj_h(:,:,:)       !< Mproj*0.5
#ifdef JESSE_MORTAR
REAL,ALLOCATABLE :: U_small(:,:,:,:,:)   !< interpolated solution for the big mortar side to small.
                                         !< two sides for 2-1 mortars, and 4+2 sides for  4-1 mortars
REAL,ALLOCATABLE :: delta_flux_jesse(:,:,:,:)  !< contribution of the two-point flux correction to the big mortar side 
REAL,ALLOCATABLE :: Ns_small(:,:,:,:,:)  !< Ja normal of big face interpolated to small faces (scaled normal vector nvec*surfelem
                                         !< with factor 4 (4-1) or factor 2 (2-1).), also with intermediate solution for 4-1 mortar
#endif /*JESSE_MORTAR*/
LOGICAL          :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!==================================================================================================================================

INTERFACE MortarBasis_BigToSmall
  MODULE PROCEDURE MortarBasis_BigToSmall
END INTERFACE

INTERFACE MortarBasis_SmallToBig_collocation
  MODULE PROCEDURE MortarBasis_SmallToBig_collocation
END INTERFACE

INTERFACE MortarBasis_SmallToBig_projection
  MODULE PROCEDURE MortarBasis_SmallToBig_projection
END INTERFACE


!INTERFACE InterpolateBigToSmall
!  MODULE PROCEDURE InterpolateBigToSmall
!END INTERFACE

CONTAINS

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
!> interpolates the data from the big  mortar  side to the small (and intermediate for 4-1) sides, stored in "small" 
!>
!>  Type 1,1st step    Type 1 ,2nd step        Type 2              Type3
!>
!>       eta                eta                  eta                 eta
!>        ^                  ^                    ^                   ^
!>        |                  |                    |                   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>    |       |          | 3 | 4 |            |   2   |           |   |   |
!>    +---+---+ --->     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>    |       |          | 1 | 2 |            |   1   |           |   |   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE InterpolateBigToSmall(ndim1,whichMortarType,Big,Small)
! MODULES
USE MOD_Preproc
USE MOD_Globals, ONLY: abort
!USE MOD_Mortar_Vars, ONLY: Mint
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)  :: ndim1  !< size of first dimension of array
INTEGER, INTENT(IN)  :: whichMortarType   !< either 1,2,3
REAL,INTENT(IN)      :: Big(  1:ndim1,0:PP_N,0:PP_N) !< solution on the big side 
REAL,INTENT(INOUT)   :: small(1:ndim1,0:PP_N,0:PP_N,1:4)
                                                    !< 4-1 mortar: sol. on small (1:4)
                                                    !< 2-1 mortar: sol. on small (1:2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l,iNb,jNb
REAL         :: tmp_jNb(1:ndim1,0:PP_N,0:PP_N)
!==================================================================================================================================

SELECT CASE(WhichMortarType)
CASE(1) !1->4
  !first  split 1 side into two, in eta direction
  DO jNb=1,2
    tmp_jNb=0.
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          tmp_jNb(:,p,q)=tmp_jNb(:,p,q) +MInt(l,q,jNb)*Big(:,p,l)
        END DO
      END DO
    END DO
    ! then split each side again into two, now in xi direction
    DO iNb=1,2
      small(:,:,:,iNb+2*(jNb-1))=0.
      DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
        DO p=0,PP_N
          DO l=0,PP_N
            small(:,p,q,iNb+2*(jNb-1))=small(:,p,q,iNb+2*(jNb-1)) +MInt(l,p,iNb)*tmp_jNb(:,l,q)
          END DO !l=1,PP_N
        END DO
      END DO 
    END DO !iNb=1,2
  END DO !jNb=1,2

CASE(2) !1->2 in eta
  DO jNb=1,2
    small(:,:,:,jNb)=0.
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          small(:,p,q,jNb)=small(:,p,q,jNb) +MInt(l,q,jNb)*Big(:,p,l)
        END DO
      END DO
    END DO
  END DO !jNb=1,2

CASE(3) !1->2 in xi
  DO iNb=1,2
    small(:,:,:,iNb)=0.
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      DO p=0,PP_N
        DO l=0,PP_N
          small(:,p,q,iNb)=small(:,p,q,iNb) +MInt(l,p,iNb)*Big(:,l,q)
        END DO
      END DO
    END DO
  END DO !iNb=1,2
CASE DEFAULT
  CALL abort(__STAMP__,&
               'Wrong mortar type')
END SELECT ! mortarType(SideID)

END SUBROUTINE InterpolateBigToSmall

!==================================================================================================================================
!> interpolates the data from the big  mortar  side to the small (and intermediate for 4-1) sides, stored in "small" 
!> ALL means that the big side and the intermediate sides are stored as well in 'small' array
!>
!>  Type 1,1st step    Type 1 ,2nd step        Type 2              Type3
!>
!>       eta                eta                  eta                 eta
!>        ^                  ^                    ^                   ^
!>        |                  |                    |                   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>    |  -1   |          | 3 | 4 |            |   2   |           |   |   |
!>    +---+---+ --->     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>    |  -2   |          | 1 | 2 |            |   1   |           |   |   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE InterpolateBigToSmall_ALL(ndim1,whichMortarType,Big,Small)
! MODULES
USE MOD_Globals, ONLY: abort
USE MOD_Preproc
!USE MOD_Mortar_Vars, ONLY: Mint
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)  :: ndim1  !< size of first dimension of array
INTEGER, INTENT(IN)  :: whichMortarType   !< either 1,2,3
REAL,INTENT(IN)      :: Big(  1:ndim1,0:PP_N,0:PP_N) !< solution on the big side 
REAL,INTENT(INOUT)   :: small(1:ndim1,0:PP_N,0:PP_N,-2:4)
                                                    !< 4-1 mortar: sol. on intermediate level (-2:1) big (0) and small (1:4)
                                                    !< 2-1 mortar: sol. on big (0) and small (1:2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l,iNb,jNb
!==================================================================================================================================
small(:,:,:,0)=Big(:,:,:) !save big mortar solution, too

SELECT CASE(WhichMortarType)
CASE(1) !1->4
  !first  split 1 side into two, in eta direction
  DO jNb=1,2
    small(:,:,:,jNb-3)=0.
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          small(:,p,q,jNb-3)=small(:,p,q,jNb-3) +MInt(l,q,jNb)*Big(:,p,l)
        END DO
      END DO
    END DO
    ! then split each side again into two, now in xi direction
    DO iNb=1,2
      small(:,:,:,iNb+2*(jNb-1))=0.
      DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
        DO p=0,PP_N
          DO l=0,PP_N
            small(:,p,q,iNb+2*(jNb-1))=small(:,p,q,iNb+2*(jNb-1)) +MInt(l,p,iNb)*small(:,l,q,jNb-3)
          END DO !l=1,PP_N
        END DO
      END DO 
    END DO !iNb=1,2
  END DO !jNb=1,2

CASE(2) !1->2 in eta
  DO jNb=1,2
    small(:,:,:,jNb)=0.
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          small(:,p,q,jNb)=small(:,p,q,jNb) +MInt(l,q,jNb)*Big(:,p,l)
        END DO
      END DO
    END DO
  END DO !jNb=1,2

CASE(3) !1->2 in xi
  DO iNb=1,2
    small(:,:,:,iNb)=0.
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      DO p=0,PP_N
        DO l=0,PP_N
          small(:,p,q,iNb)=small(:,p,q,iNb) +MInt(l,p,iNb)*Big(:,l,q)
        END DO
      END DO
    END DO
  END DO !iNb=1,2
CASE DEFAULT
  CALL abort(__STAMP__,&
               'Wrong mortar type')
END SELECT ! mortarType(SideID)

END SUBROUTINE InterpolateBigToSmall_ALL

END MODULE MOD_Mortar_Vars
