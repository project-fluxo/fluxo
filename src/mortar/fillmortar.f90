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
!> Routines that perform the projection operation between nonconforming interfaces using the operators set up in module
!> mortar
!> 
!> Contains the routines to
!> - interpolate the solution at the large sides to the small ones, which are used for flux computation
!> - project the flux from the small sides back to the large ones
!==================================================================================================================================
MODULE MOD_FillMortar
IMPLICIT NONE
PRIVATE

INTERFACE U_Mortar
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE U_Mortar_Eqn
  MODULE PROCEDURE U_Mortar_Eqn
END INTERFACE

INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

INTERFACE Flux_Cons_Mortar
  MODULE PROCEDURE Flux_Cons_Mortar
END INTERFACE

PUBLIC::U_Mortar,U_Mortar_Eqn,Flux_Mortar,Flux_Cons_mortar

#ifdef JESSE_MORTAR
PUBLIC::Fill_delta_flux_jesse
#endif

CONTAINS

!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side, using 1D interpolation operators MInt(:,:,1:2).
!> This is used to obtain the face solution for flux computation.
!>
!> NOTE: that input arrays can be both normal solution or gradient data.
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
SUBROUTINE U_Mortar(Uface_master,Uface_slave,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,   ONLY: FS2M,nSides 
USE MOD_Mortar_Vars, ONLY: InterpolateBigToSmall
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Uface_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: Uface_slave( 1:PP_nVar,0:PP_N,0:PP_N,FirstSlaveSide:LastSlaveSide) !< (INOUT) can be U or Grad_Ux/y/z_master
LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
REAL         :: U_small_loc(1:PP_nVar,0:PP_N,0:PP_N,1:4)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
  lastMortarSideID =  lastMortarMPISide
ELSE
  firstMortarSideID = firstMortarInnerSide
  lastMortarSideID =  lastMortarInnerSide
END IF !doMPISides


DO MortarSideID=firstMortarSideID,lastMortarSideID
  iSide=MortarType(2,MortarSideID) !ID in list 1:nMortarSides

  CALL InterpolateBigToSmall(PP_nVar,MortarType(1,MortarSideID),Uface_master(:,:,:,MortarSideID),U_small_loc(:,:,:,:))
  !Now save the small sides into master/slave arrays
  IF(MortarType(1,MortarSideID).EQ.1)THEN
    nMortars=4
  ELSE
    nMortars=2
  END IF !MortarType
  !iSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,iSide)
    flip  = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
      CASE(0) ! small master side
        Uface_master(:,:,:,SideID)=U_small_loc(:,:,:,iMortar)
      CASE(1:4) ! small slave side
        DO q=0,PP_N; DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=U_small_loc(:,FS2M(1,p,q,flip), &
                                              FS2M(2,p,q,flip),iMortar)
        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar

!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side, using 1D interpolation operators MInt(:,:,1:2).
!> This is used to obtain the face solution for flux computation.
!>
!> NOTE: input arrays can be only conservative face solution, 
!> U_small array is filled here containing only big mortar face interpolations, needed for JESSE_MORTAR
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
SUBROUTINE U_Mortar_Eqn(Uface_master,Uface_slave,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,   ONLY: FS2M,nSides 
#if defined(JESSE_MORTAR)
USE MOD_Mortar_Vars, ONLY: U_small
#else
USE MOD_Mortar_Vars, ONLY: InterpolateBigToSmall
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Uface_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: Uface_slave( 1:PP_nVar,0:PP_N,0:PP_N,FirstSlaveSide:LastSlaveSide) !< (INOUT) can be U or Grad_Ux/y/z_master
LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
REAL         :: U_small_loc(1:PP_nVar,0:PP_N,0:PP_N,1:4)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
  lastMortarSideID =  lastMortarMPISide
ELSE
  firstMortarSideID = firstMortarInnerSide
  lastMortarSideID =  lastMortarInnerSide
END IF !doMPISides


DO MortarSideID=firstMortarSideID,lastMortarSideID
  iSide=MortarType(2,MortarSideID) !ID in list 1:nMortarSides

#if defined(JESSE_MORTAR)
  CALL InterpolateBigToSmall_ALL_Eqn(MortarType(1,MortarSideID),Uface_master(:,:,:,MortarSideID),U_small(:,:,:,:,iSide))
  U_small_loc=U_small(:,:,:,1:4,iSide)
#else
  CALL InterpolateBigToSmall(PP_nVar,MortarType(1,MortarSideID),Uface_master(:,:,:,MortarSideID),U_small_loc(:,:,:,:))
#endif
  !Now save the small sides into master/slave arrays
  IF(MortarType(1,MortarSideID).EQ.1)THEN
    nMortars=4
  ELSE
    nMortars=2
  END IF !MortarType
  !iSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,iSide)
    flip  = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
      CASE(0) ! small master side
        Uface_master(:,:,:,SideID)=U_small_loc(:,:,:,iMortar)
      CASE(1:4) ! small slave side
        DO q=0,PP_N; DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=U_small_loc(:,FS2M(1,p,q,flip), &
                                                  FS2M(2,p,q,flip),iMortar)
        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar_Eqn

!==================================================================================================================================
!>  Fills master side from small non-conforming sides, using 1D projection operators Mproj(:,:,1:2)
!>
!> This routine is used to project the numerical flux at the small sides of the nonconforming interface to the corresponding large
!>  ones.
!>
!>  Type 1,1st step    Type 1 ,2nd step        Type 2              Type3
!>
!>       eta                eta                  eta                 eta
!>        ^                  ^                    ^                   ^
!>        |                  |                    |                   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>    | 3 | 4 |          |  -1   |            |   2   |           |   |   |
!>    +---+---+ --->     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>    | 1 | 2 |          |  -2   |            |   1   |           |   |   |
!>    +---+---+          +---+---+            +---+---+           +---+---+
!>
!> flag weak changes the sign of incoming flux. 
!> when used for lifting, weak=.false, since volint of lifting is in strong form
!> and hence Flux_L=1/2*(u_R-u_L)*outwardnormal_L  =  Flux_R = 1/2*(u_L-u_R)*outwardnormal_R
!==================================================================================================================================
SUBROUTINE Flux_Mortar(Flux_master,Flux_slave,doMPISides,weak)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: MProj
#ifdef JESSE_MORTAR
USE MOD_Mortar_Vars, ONLY: delta_flux_jesse
#endif
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide,FS2M
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,LastSlaveSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Flux_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< on input: has flux from small mortar sides 
                                                                      !< on output: flux on big mortar sides filled
REAL,INTENT(IN   )   :: Flux_slave(1:PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide) !<has flux from small mortar sides,
                                                                      !< set -F_slave in call if surfint is weak (dg.f90)
                                                                      !< set +F_slave in call if surfint is strong (lifting)
                                                            
LOGICAL,INTENT(IN) :: doMPISides                                    !< flag whether MPI sides are processed
LOGICAL,INTENT(IN) :: weak                                          !< flag whether strong(lifting) or weak(equation) form is used 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l,iNb,jNb
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
REAL         :: Flux_small(PP_nVar,0:PP_N,0:PP_N,1:4)
REAL         :: Flux_tmp(PP_nVar,0:PP_N,0:PP_N)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
  lastMortarSideID =  lastMortarMPISide
ELSE
  firstMortarSideID = firstMortarInnerSide
  lastMortarSideID =  lastMortarInnerSide
END IF !doMPISides

DO MortarSideID=firstMortarSideID,lastMortarSideID

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! small master side
      Flux_small(:,:,:,iMortar)=Flux_master(:,:,:,SideID)
    CASE(1:4) ! slave sides (should only occur for MPI)
      IF(weak)THEN
        DO q=0,PP_N; DO p=0,PP_N
          Flux_small(:,FS2M(1,p,q,flip),FS2M(2,p,q,flip),iMortar)=-Flux_slave(:,p,q,SideID)
        END DO; END DO !p,q
      ELSE    ! do not change sign if strong form when used for lifting! 
        DO q=0,PP_N; DO p=0,PP_N
          Flux_small(:,FS2M(1,p,q,flip),FS2M(2,p,q,flip),iMortar)= Flux_slave(:,p,q,SideID)
        END DO; END DO !p,q
      END IF !weak
    END SELECT !slave sides
  END DO

  Flux_master(:,:,:,MortarSideID)=0.
#ifdef JESSE_MORTAR
  IF(weak)THEN !use 'weak' flag to identify that its called from eqn, not for lifting flux!!!
    !add previously computed two-point correction flux (Fill_delta_flux_jesse must be called before!!) 
    Flux_master(:,:,:,MortarSideID)=delta_Flux_Jesse(:,:,:,iSide)
  END IF
#endif

  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    ! first in xi
    DO jNb=1,2
      Flux_tmp(:,:,:)=0.
      DO iNb=1,2
        DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
          DO p=0,PP_N
            DO l=0,PP_N
              Flux_tmp(:,p,q)=Flux_tmp(:,p,q) + Mproj(l,p,iNb)*Flux_small(:,l,q,iNb+2*(jNb-1))
            END DO !l=0,PP_N
          END DO !p=0,PP_N
        END DO !q=0,PP_N
      END DO !iNb=1,2

      !then in eta (l,q)
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO q=0,PP_N
!          Flux_master(:,p,q,MortarSideID)= 0. 
          DO l=0,PP_N
            Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) &
                                                  + Mproj(l,q,jNb)*Flux_tmp(:,p,l) 
          END DO !l=1,PP_N
        END DO !q=0,PP_N
      END DO !p=0,PP_N
    END DO !jNb=1,2

  CASE(2) !1->2 in eta
    DO jNb=1,2
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO q=0,PP_N !index big side 
          DO l=0,PP_N !index small side
            Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) &
                                            + Mproj(l,q,jNb)*Flux_small(:,p,l,jNb)
          END DO !l=0,PP_N
        END DO !q=0,PP_N
      END DO !p=0,PP_N
    END DO !jNb=1,2

  CASE(3) !1->2 in xi
    DO iNb=1,2
      DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
        DO p=0,PP_N !index big side
!          Flux_master(:,p,q,MortarSideID)= 0.
          DO l=0,PP_N !index small side
            Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID)  &
                                            + Mproj(l,p,iNb)*Flux_small(:,l,q,iNb)
          END DO !l=0,PP_N
        END DO !p=0,PP_N
      END DO !q=0,PP_N
    END DO !iNb=1,2

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID

END SUBROUTINE Flux_Mortar



!==================================================================================================================================
!> Calls flux_mortar if flux is conservative, slave=master, with array 
!>   Flux_master(:,:,:,1:nSides)=Flux(:,:,:,1:nSides)
!>   Flux_slave(:,:,:,FirstSlaveSide:LastSlaveSide)=Flux(:,:,:,FirstSlaveSide:LastSlaveSide)
!==================================================================================================================================
SUBROUTINE Flux_Cons_Mortar(Flux,doMPISides,weak)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: nSides
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,LastSlaveSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< on input: has flux from small mortar sides 
LOGICAL,INTENT(IN) :: doMPISides                                    !< flag whether MPI sides are processed
LOGICAL,INTENT(IN) :: weak                                          !< flag whether strong or weak form is used
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL Flux_Mortar(Flux(:,:,:,1:nSides),Flux(:,:,:,FirstSlaveSide:LastSlaveSide),doMPIsides,weak)

END SUBROUTINE Flux_Cons_mortar


#ifdef JESSE_MORTAR
!==================================================================================================================================
!>  if useEntropyMortar=.TRUE. , Data is transformed to entropy variables before interpolation 
!>  and transformed back to conservative after interpolation 
!>  else conservative variables are used
!>
!==================================================================================================================================
SUBROUTINE InterpolateBigToSmall_ALL_Eqn(whichMortarType,BigCons,SmallCons)
! MODULES
USE MOD_Preproc
#if defined(navierstokes) || defined(mhd)
USE MOD_Equation_Vars, ONLY: useEntropyMortar,ConsToEntropyVec,EntropyToConsVec
#endif
USE MOD_Mortar_vars, ONLY: InterpolateBigToSmall_ALL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)  :: whichMortarType   !< either 1,2,3
REAL,INTENT(IN)      :: BigCons(  1:PP_nVar,0:PP_N,0:PP_N) !< solution on the big side 
REAL,INTENT(INOUT)   :: smallCons(1:PP_nVar,0:PP_N,0:PP_N,-2:4)
                                                    !< 4-1 mortar: sol. on intermediate level (-2:1) big (0) and small (1:4)
                                                    !< 2-1 mortar: sol. on big (0) and small (1:2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(navierstokes) || defined(mhd)
REAL         :: BigEnt(  1:PP_nVar,0:PP_N,0:PP_N) !< solution on the big side 
REAL         :: smallEnt(1:PP_nVar,0:PP_N,0:PP_N,-2:4)
!==================================================================================================================================
IF(useEntropyMortar)THEN
  smallCons(:,:,:,0)=BigCons(:,:,:) !save big mortar Cons solution, too
  CALL ConsToEntropyVec((PP_N+1)*(PP_N+1),BigEnt,BigCons)
  CALL InterpolateBigToSmall_ALL(PP_nVar,whichMortarType,BigEnt,SmallEnt)
  SELECT CASE(WhichMortarType)
  CASE(1) !1->4
    CALL EntropyToConsVec((PP_N+1)*(PP_N+1)*2,smallEnt(:,:,:,-2:1), SmallCons(:,:,:,-2:1))
    CALL EntropyToConsVec((PP_N+1)*(PP_N+1)*4,smallEnt(:,:,:, 1:4), SmallCons(:,:,:, 1:4))
  CASE(2,3) !1->2 in eta,1->2 in xi
    CALL EntropyToConsVec((PP_N+1)*(PP_N+1)*2,smallEnt(:,:,:, 1:2), SmallCons(:,:,:, 1:2))
  END SELECT ! mortarType(SideID)
ELSE
  CALL InterpolateBigToSmall_ALL(PP_nVar,whichMortarType,BigCons,SmallCons)
END IF
#else 
CALL InterpolateBigToSmall_ALL(PP_nVar,whichMortarType,BigCons,SmallCons)
#endif /* (defined(navierstokes) || defined(mhd)) */

END SUBROUTINE InterpolateBigToSmall_ALL_Eqn

!==================================================================================================================================
!> Computes the two-point flux correction (Jesse's paper: https://arxiv.org/pdf/2005.03237.pdf) and stores it in "delta_flux_jesse" 
!> small side surface metric is already scaled by factor of 2 for 2-1 mortar or 4 for 4-1 mortar, to be the same polynomial
!> and then 0.5 must be put in the projection matrix Mproj=>Mproj_h
!==================================================================================================================================
SUBROUTINE Fill_delta_Flux_Jesse()
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: Mint,Mproj_h,U_small,Ns_small
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nMortarSides
USE MOD_Equation_Vars,ONLY: MortarFluxAverageVec
USE MOD_Mortar_Vars, ONLY: delta_flux_jesse !<<<<== is filled
#if defined(navierstokes) || defined(mhd)
USE MOD_Equation_Vars,ONLY: nAuxVar
USE MOD_Flux_Average, ONLY: EvalUaux
#endif /*navierstokes || mhd */
#if defined(NONCONS)
USE MOD_Flux_Average, ONLY: AddNonConsFluxVec
#endif /*NONCONS*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l,r,iNb,jNb
INTEGER      :: MortarSideID,iSide
REAL         :: Flux_tmp(PP_nVar,0:PP_N,0:PP_N)
REAL         :: Flux_tp_l(PP_nVar,0:PP_N)
REAL         :: Flux_corr_l(PP_nVar)
#if defined(navierstokes) || defined(mhd)
REAL         :: Uaux_small(1:nAuxVar,0:PP_N,0:PP_N,-2:4)
#endif /*navierstokes || mhd */
#ifdef NONCONS
REAL         :: Flux_tp(PP_nVar,0:PP_N,0:PP_N)
REAL         :: Flux_tp_T(PP_nVar,0:PP_N,0:PP_N)
#endif
!==================================================================================================================================

delta_Flux_jesse(:,:,:,:)= 0. 
DO iSide=1,nMortarSides

  MortarSideID=MortarInfo(MI_SIDEID,0,iSide)
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
#if defined(navierstokes) || defined(mhd)
    CALL EvalUaux((PP_N+1)*(PP_N+1)*7,U_small(:,:,:,-2:4,iSide),Uaux_small(:,:,:,-2:4))
#endif /*navierstokes || mhd */

    ! first in xi
    DO jNb=1,2
      Flux_tmp(:,:,:)= 0.
      DO iNb=1,2
        DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
#ifdef NONCONS
          DO l=0,PP_N !index small side                       
            DO p=0,PP_N ! index big side                                     !<this is the intermediate side!>
              CALL MortarFluxAverageVec(   U_small(:,l,q,iNb+2*(jNb-1),iSide),   U_small(:,p,q,jNb-3,iSide), &
                                        Uaux_small(:,l,q,iNb+2*(jNb-1))      ,Uaux_small(:,p,q,jNb-3)      , &
                                          Ns_small(:,l,q,iNb+2*(jNb-1),iSide),  Ns_small(:,p,q,jNb-3,iSide),Flux_tp(:,l,p))
              Flux_tp_T(:,p,l)=Flux_tp(:,l,p)
              CALL AddNonConsFluxVec(      U_small(:,l,q,iNb+2*(jNb-1),iSide),   U_small(:,p,q,jNb-3,iSide), &
                                        Uaux_small(:,l,q,iNb+2*(jNb-1))      ,Uaux_small(:,p,q,jNb-3)      , &
                                          Ns_small(:,l,q,iNb+2*(jNb-1),iSide),  Ns_small(:,p,q,jNb-3,iSide), Flux_tp_T(:,p,l))
              CALL AddNonConsFluxVec(      U_small(:,p,q,jNb-3,iSide),   U_small(:,l,q,iNb+2*(jNb-1),iSide), &
                                        Uaux_small(:,p,q,jNb-3)      ,Uaux_small(:,l,q,iNb+2*(jNb-1))      , &
                                          Ns_small(:,p,q,jNb-3,iSide),  Ns_small(:,l,q,iNb+2*(jNb-1),iSide), Flux_tp(:,l,p))
            END DO !p
          END DO !l
          DO l=0,PP_N !index small side                       
            Flux_corr_l(:)=0.
            DO r=0,PP_N  !index big side
              Flux_corr_l(:)=Flux_corr_l(:) + Flux_tp_T(:,r,l)*Mint(r,l,iNb) !is the reduction with *lagbaseBig_r(ximortar_l)
            END DO !r=0,PP_N
            DO p=0,PP_N !index big side
              Flux_tmp(:,p,q)=Flux_tmp(:,p,q)  &
                                  + Mproj_h(l,p,iNb)*(Flux_tp(:,l,p)-Flux_corr_l(:))
            END DO !p=0,PP_N
          END DO !l=0,PP_N
#else 
          DO l=0,PP_N !index small side                       
            DO p=0,PP_N ! index big side                                     !<this is the intermediate side!>
              CALL MortarFluxAverageVec( U_small(:,l,q,iNb+2*(jNb-1),iSide),   U_small(:,p,q,jNb-3,iSide), &
#if defined(navierstokes) || defined(mhd)
                                      Uaux_small(:,l,q,iNb+2*(jNb-1))      ,Uaux_small(:,p,q,jNb-3)      , &
#endif /*navierstokes || mhd */
                                        Ns_small(:,l,q,iNb+2*(jNb-1),iSide),  Ns_small(:,p,q,jNb-3,iSide),Flux_tp_l(:,p))
            END DO
            Flux_corr_l(:)=0.
            DO r=0,PP_N  !index big side
              Flux_corr_l(:)=Flux_corr_l(:) + Flux_tp_l(:,r)*Mint(r,l,iNb) !is the reduction with *lagbaseBig_r(ximortar_l)
            END DO !r=0,PP_N
            DO p=0,PP_N !index big side
              Flux_tmp(:,p,q)=Flux_tmp(:,p,q)  &
                                  + Mproj_h(l,p,iNb)*(Flux_tp_l(:,p)-Flux_corr_l(:))
            END DO !p=0,PP_N
          END DO !l=0,PP_N
#endif /*NONCONS*/
        END DO !q=0,PP_N
      END DO !iNb=1,2
      !then in eta (l,q)
     
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        !mortar flux in eta (l,q), for fixed p 
#ifdef NONCONS
        DO l=0,PP_N !index small side                      
          DO q=0,PP_N  !index big side                              !<this is the big side>
            CALL MortarFluxAverageVec(   U_small(:,p,l,jNb-3,iSide),   U_small(:,p,q,    0,iSide), &
                                      Uaux_small(:,p,l,jNb-3)      ,Uaux_small(:,p,q,    0)      , &
                                        Ns_small(:,p,l,jNb-3,iSide),  Ns_small(:,p,q,    0,iSide),Flux_tp(:,l,q))
            Flux_tp_T(:,q,l)=Flux_tp(:,l,q) !symmetric flux
            CALL AddNonConsFluxVec(      U_small(:,p,l,jNb-3,iSide),   U_small(:,p,q,    0,iSide), &
                                      Uaux_small(:,p,l,jNb-3)      ,Uaux_small(:,p,q,    0)      , &
                                        Ns_small(:,p,l,jNb-3,iSide),  Ns_small(:,p,q,    0,iSide), Flux_tp_T(:,q,l))
            CALL AddNonConsFluxVec(      U_small(:,p,q,    0,iSide),   U_small(:,p,l,jNb-3,iSide), &
                                      Uaux_small(:,p,q,    0)      ,Uaux_small(:,p,l,jNb-3)      , &
                                        Ns_small(:,p,q,    0,iSide),  Ns_small(:,p,l,jNb-3,iSide), Flux_tp(:,l,q))
          END DO !q
        END DO !l
        DO l=0,PP_N !index small side                      
          Flux_corr_l(:)=0.
          DO r=0,PP_N !index big side
            Flux_corr_l(:)= Flux_corr_l(:) + Flux_tp_T(:,r,l)*Mint(r,l,jNb) !is the reduction with *lagbaseBig_r(etamortar_l)
          END DO !r
          DO q=0,PP_N !index big side 
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                            + Mproj_h(l,q,jNb)*(Flux_tmp(:,p,l)+ Flux_tp(:,l,q) - Flux_corr_l(:))
          END DO !q=0,PP_N
        END DO !l=0,PP_N
#else
        DO l=0,PP_N !index small side                      
          DO q=0,PP_N  !index big side                             !<this is the big side>
            CALL MortarFluxAverageVec( U_small(:,p,l,jNb-3,iSide),   U_small(:,p,q,0,iSide), &
#if defined(navierstokes) || defined(mhd)
                                    Uaux_small(:,p,l,jNb-3)      ,Uaux_small(:,p,q,0)      , &
#endif /*navierstokes*/
                                      Ns_small(:,p,l,jNb-3,iSide),  Ns_small(:,p,q,0,iSide),Flux_tp_l(:,q))
          END DO
          Flux_corr_l(:)=0.
          DO r=0,PP_N !index big side
            Flux_corr_l(:)= Flux_corr_l(:) + Flux_tp_l(:,r)*Mint(r,l,jNb) !is the reduction with *lagbaseBig_r(etamortar_l)
          END DO !r
          DO q=0,PP_N !index big side 
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                            + Mproj_h(l,q,jNb)*(Flux_tmp(:,p,l)+ Flux_tp_l(:,q) - Flux_corr_l(:))
          END DO !q=0,PP_N
        END DO !l=0,PP_N
#endif /*NONCONS*/
      END DO !p=0,PP_N
    END DO !jNb=1,2

  CASE(2) !1->2 in eta
#if defined(navierstokes) || defined(mhd)
    CALL EvalUaux((PP_N+1)*(PP_N+1)*3,U_small(:,:,:,0:2,iSide),Uaux_small(:,:,:,0:2))
#endif /*navierstokes || mhd */
    DO jNb=1,2
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        !mortar flux in eta (l,q), for fixed p 
#ifdef NONCONS
        DO l=0,PP_N !index small side
          DO q=0,PP_N ! index big side                             !<this is the big side!>
            CALL MortarFluxAverageVec(   U_small(:,p,l,jNb,iSide),   U_small(:,p,q,  0,iSide), &
                                      Uaux_small(:,p,l,jNb)      ,Uaux_small(:,p,q,  0)      , &
                                        Ns_small(:,p,l,jNb,iSide),  Ns_small(:,p,q,  0,iSide),Flux_tp(:,l,q))
            Flux_tp_T(:,q,l)=Flux_tp(:,l,q) !symmetric flux
            CALL AddNonConsFluxVec(      U_small(:,p,l,jNb,iSide),   U_small(:,p,q,  0,iSide), &
                                      Uaux_small(:,p,l,jNb)      ,Uaux_small(:,p,q,  0)      , &
                                        Ns_small(:,p,l,jNb,iSide),  Ns_small(:,p,q,  0,iSide), Flux_tp_T(:,q,l))
            CALL AddNonConsFluxVec(      U_small(:,p,q,  0,iSide),   U_small(:,p,l,jNb,iSide), &
                                      Uaux_small(:,p,q,  0)      ,Uaux_small(:,p,l,jNb)      , &
                                        Ns_small(:,p,q,  0,iSide),  Ns_small(:,p,l,jNb,iSide), Flux_tp(:,l,q))
          END DO !q
        END DO !l
        DO l=0,PP_N !index small side
          Flux_corr_l(:)=0.
          DO r=0,PP_N !index big side
            Flux_corr_l(:)= Flux_corr_l(:) + Flux_tp_T(:,r,l)*Mint(r,l,jNb) !is the reduction with *lagbaseBig_r(etamortar_l)
          END DO !r
          DO q=0,PP_N !index big side 
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                            + Mproj_h(l,q,jNb)*(Flux_tp(:,l,q) - Flux_corr_l(:))
          END DO !q=0,PP_N
        END DO !l=0,PP_N
#else
        DO l=0,PP_N !index small side                      
          DO q=0,PP_N  !index big side                         !<this is the big side>
            CALL MortarFluxAverageVec( U_small(:,p,l,jNb,iSide),   U_small(:,p,q,0,iSide), &
#if defined(navierstokes) || defined(mhd)
                                    Uaux_small(:,p,l,jNb)      ,Uaux_small(:,p,q,0)      , &
#endif /*navierstokes || mhd */
                                      Ns_small(:,p,l,jNb,iSide),  Ns_small(:,p,q,0,iSide),Flux_tp_l(:,q))
          END DO
          Flux_corr_l(:)=0.
          DO r=0,PP_N !index big side
            Flux_corr_l(:)= Flux_corr_l(:) + Flux_tp_l(:,r)*Mint(r,l,jNb) !is the reduction with *lagbaseBig_r(etamortar_l)
          END DO !r
          DO q=0,PP_N !index big side 
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                            + Mproj_h(l,q,jNb)*(Flux_tp_l(:,q) - Flux_corr_l(:))
          END DO !q=0,PP_N
        END DO !l=0,PP_N
#endif /*NONCONS*/
      END DO !p=0,PP_N
    END DO !jNb=1,2

  CASE(3) !1->2 in xi
#if defined(navierstokes) || defined(mhd)
    CALL EvalUaux((PP_N+1)*(PP_N+1)*3,U_small(:,:,:,0:2,iSide),Uaux_small(:,:,:,0:2))
#endif /*navierstokes || mhd */
    DO iNb=1,2
      DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
#ifdef NONCONS
        DO l=0,PP_N !index small side
          DO p=0,PP_N ! index big side                             !<this is the big side!> 
            CALL MortarFluxAverageVec(   U_small(:,l,q,iNb,iSide),   U_small(:,p,q,  0,iSide), &
                                      Uaux_small(:,l,q,iNb)      ,Uaux_small(:,p,q,  0)      , &
                                        Ns_small(:,l,q,iNb,iSide),  Ns_small(:,p,q,  0,iSide),Flux_tp(:,l,p))
            Flux_tp_T(:,p,l)=Flux_tp(:,l,p) !symmetric flux
            CALL AddNonConsFluxVec(      U_small(:,l,q,iNb,iSide),   U_small(:,p,q,  0,iSide), &
                                      Uaux_small(:,l,q,iNb)      ,Uaux_small(:,p,q,  0)      , &
                                        Ns_small(:,l,q,iNb,iSide),  Ns_small(:,p,q,  0,iSide), Flux_tp_T(:,p,l))
            CALL AddNonConsFluxVec(      U_small(:,p,q,  0,iSide),   U_small(:,l,q,iNb,iSide), &
                                      Uaux_small(:,p,q,  0)      ,Uaux_small(:,l,q,iNb)      , &
                                        Ns_small(:,p,q,  0,iSide),  Ns_small(:,l,q,iNb,iSide), Flux_tp(:,l,p))
          END DO !p
        END DO !l
        DO l=0,PP_N !index small side
          Flux_corr_l(:)=0.
          DO r=0,PP_N  !index big side
            Flux_corr_l(:)=Flux_corr_l(:) + Flux_tp_T(:,r,l)*Mint(r,l,iNb) !is the reduction with *lagbaseBig_r(ximortar_l)
          END DO !r=0,PP_N
          DO p=0,PP_N !index big side
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                              + Mproj_h(l,p,iNb)*(Flux_tp(:,l,p)-Flux_corr_l(:))
          END DO !p=0,PP_N
        END DO !l=0,PP_N
#else
        DO l=0,PP_N !index small side
          DO p=0,PP_N ! index big side                              !<this is the big side!>
            CALL MortarFluxAverageVec(   U_small(:,l,q,iNb,iSide),   U_small(:,p,q,0,iSide), &
#if defined(navierstokes) || defined(mhd)
                                      Uaux_small(:,l,q,iNb)      ,Uaux_small(:,p,q,0)      , &
#endif /*navierstokes || mhd */
                                        Ns_small(:,l,q,iNb,iSide),  Ns_small(:,p,q,0,iSide),Flux_tp_l(:,p))
          END DO
          Flux_corr_l(:)=0.
          DO r=0,PP_N  !index big side
            Flux_corr_l(:)=Flux_corr_l(:) + Flux_tp_l(:,r)*Mint(r,l,iNb) !is the reduction with *lagbaseBig_r(ximortar_l)
          END DO !r=0,PP_N
          DO p=0,PP_N !index big side
            delta_flux_jesse(:,p,q,iSide)=delta_flux_jesse(:,p,q,iSide) &
                                              + Mproj_h(l,p,iNb)*(Flux_tp_l(:,p)-Flux_corr_l(:))
          END DO !p=0,PP_N
        END DO !l=0,PP_N
#endif /* NONCONS*/
      END DO !q=0,PP_N
    END DO !iNb=1,2

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID

END SUBROUTINE Fill_delta_Flux_Jesse
#endif /*JESSE_MORTAR*/


END MODULE MOD_FillMortar
