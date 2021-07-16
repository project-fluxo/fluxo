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


!INTERFACE InterpolateBigToSmall
!  MODULE PROCEDURE InterpolateBigToSmall
!END INTERFACE

CONTAINS

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
END SELECT ! mortarType(SideID)

END SUBROUTINE InterpolateBigToSmall_ALL

END MODULE MOD_Mortar_Vars
