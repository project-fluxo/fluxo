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

INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_Mortar,Flux_Mortar

CONTAINS

!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side, using 1D interpolation operators M_0_1,M_0_2.
!> This is used to obtain the face solution for flux computation.
!>
!> NOTE: that input arrays can be both normal solution or gradient data.
!> NOTE2: fillmortar is only built for PP_N as even in case of overint mortarized data
!>
!>       Type 1               Type 2              Type3
!>
!>        eta                  eta                 eta
!>         ^                    ^                   ^
!>         |                    |                   |
!>     +---+---+            +---+---+           +---+---+
!>     | 3 | 4 |            |   2   |           |   |   |
!>     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>     | 1 | 2 |            |   1   |           |   |   |
!>     +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE U_Mortar(Uface_master,Uface_slave,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2,U_small
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,   ONLY: FS2M,nSides 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Uface_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: Uface_slave( 1:PP_nVar,0:PP_N,0:PP_N,FirstSlaveSide:LastSlaveSide) !< (INOUT) can be U or Grad_Ux/y/z_master
LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
!REAL         :: U_tmp( PP_nVar,0:PP_N,0:PP_N,1:4)
!REAL         :: U_tmp2(PP_nVar,0:PP_N,0:PP_N,1:2)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
  lastMortarSideID =  lastMortarMPISide
ELSE
  firstMortarSideID = firstMortarInnerSide
  lastMortarSideID =  lastMortarInnerSide
END IF !doMPISides

U_small=0.

DO MortarSideID=firstMortarSideID,lastMortarSideID
  iSide=MortarType(2,MortarSideID)
  U_small(:,:,:,0,iSide)=Uface_master(:,:,:,MortarSideID) !save solution on big mortar
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    !first  split 1 side into two, in eta direction
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    U_tmp2(iVar,p,:,1)  =  M1 * Uface_master(iVar,p,:,MortarSideID)
    !    U_tmp2(iVar,p,:,2)  =  M2 * Uface_master(iVar,p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          U_small(:,p,q,1-3,iSide)=U_small(:,p,q,1-3,iSide)+M_0_1(l,q)*Uface_master(:,p,l,MortarSideID)
          U_small(:,p,q,2-3,iSide)=U_small(:,p,q,2-3,iSide)+M_0_2(l,q)*Uface_master(:,p,l,MortarSideID)
        END DO
      END DO
    END DO
    ! then split each side again into two, now in xi direction
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    U_tmp(iVar,:,q,1)  =  M1 * U_tmp2(iVar,:,q,1)
      !    U_tmp(iVar,:,q,2)  =  M2 * U_tmp2(iVar,:,q,1)
      !    U_tmp(iVar,:,q,3)  =  M1 * U_tmp2(iVar,:,q,2)
      !    U_tmp(iVar,:,q,4)  =  M2 * U_tmp2(iVar,:,q,2)
      DO p=0,PP_N
        DO l=0,PP_N
          U_small(:,p,q,1,iSide)=U_small(:,p,q,1,iSide)+M_0_1(l,p)*U_small(:,l,q,1-3,iSide)
          U_small(:,p,q,2,iSide)=U_small(:,p,q,2,iSide)+M_0_2(l,p)*U_small(:,l,q,1-3,iSide)
          U_small(:,p,q,3,iSide)=U_small(:,p,q,3,iSide)+M_0_1(l,p)*U_small(:,l,q,2-3,iSide)
          U_small(:,p,q,4,iSide)=U_small(:,p,q,4,iSide)+M_0_2(l,p)*U_small(:,l,q,2-3,iSide)
        END DO !l=1,PP_N
      END DO
    END DO 

  CASE(2) !1->2 in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    U_tmp(iVar,p,:,1)  =  M1 * Uface_master(iVar,p,:,MortarSideID)
    !    U_tmp(iVar,p,:,2)  =  M2 * Uface_master(iVar,p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          U_small(:,p,q,1,iSide)=U_small(:,p,q,1,iSide)+M_0_1(l,q)*Uface_master(:,p,l,MortarSideID)
          U_small(:,p,q,2,iSide)=U_small(:,p,q,2,iSide)+M_0_2(l,q)*Uface_master(:,p,l,MortarSideID)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    U_tmp(iVar,:,q,1)  =  M1 * Uface_master(iVar,:,q,MortarSideID)
      !    U_tmp(iVar,:,q,2)  =  M2 * Uface_master(iVar,:,q,MortarSideID)
      DO p=0,PP_N
        DO l=0,PP_N
          U_small(:,p,q,1,iSide)=U_small(:,p,q,1,iSide)+M_0_1(l,p)*Uface_master(:,l,q,MortarSideID)
          U_small(:,p,q,2,iSide)=U_small(:,p,q,2,iSide)+M_0_2(l,p)*Uface_master(:,l,q,MortarSideID)
        END DO
      END DO
    END DO
  END SELECT ! mortarType(SideID)
 
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
        Uface_master(:,:,:,SideID)=U_small(:,:,:,iMortar,iSide)
      CASE(1:4) ! small slave side
        DO q=0,PP_N; DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=U_small(:,FS2M(1,p,q,flip), &
                                              FS2M(2,p,q,flip),iMortar,iSide)
        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar


!==================================================================================================================================
!>  Fills master side from small non-conforming sides, using 1D projection operators M_1_0,M_2_0
!>
!> This routine is used to project the numerical flux at the small sides of the nonconforming interface to the corresponding large
!>  ones.
!>
!>        Type 1               Type 2              Type3
!>         eta                  eta                 eta
!>          ^                    ^                   ^
!>          |                    |                   |
!>      +---+---+            +---+---+           +---+---+
!>      | 3 | 4 |            |   2   |           |   |   |
!>      +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>      | 1 | 2 |            |   1   |           |   |   |
!>      +---+---+            +---+---+           +---+---+
!>
!> flag weak changes the sign of incoming flux. 
!> when used for lifting, weak=.false, since volint of lifting is in strong form
!> and hence Flux_L=1/2*(u_R-u_L)*outwardnormal_L  =  Flux_R = 1/2*(u_L-u_R)*outwardnormal_R
!==================================================================================================================================
SUBROUTINE Flux_Mortar(Flux_master,Flux_slave,doMPISides,weak)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_1_0,M_2_0
#ifdef JESSE_MORTAR
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2,U_small,Ns_small
USE MOD_Equation_Vars,ONLY: MortarFluxAverageVec
#endif /*JESSE_MORTAR*/
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
LOGICAL,INTENT(IN) :: weak                                          !< flag whether strong or weak form is used
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
REAL         :: Flux_small(PP_nVar,0:PP_N,0:PP_N,1:4)
REAL         :: Flux_tmp(PP_nVar,0:PP_N,0:PP_N,1:2)
#ifdef JESSE_MORTAR
REAL         :: Flux_tp1(PP_nVar,0:PP_N,0:PP_N,1:2)
REAL         :: Flux_corr1(PP_nVar,0:PP_N,1:2)
#endif /*JESSE_MORTAR*/
REAL         :: M1(0:PP_N,0:PP_N),M2(0:PP_N,0:PP_N)
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
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    M1=0.5*M_1_0; M2=0.5*M_2_0
    Flux_small = 4.*Flux_small
    ! first in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_tmp(iVar,:,q,1)  =  M1 * Flux_small(iVar,:,q,1) + M2 * Flux_small(iVar,:,q,2)
      !    Flux_tmp(iVar,:,q,2)  =  M1 * Flux_small(iVar,:,q,1) + M2 * Flux_small(iVar,:,q,2)
      DO p=0,PP_N
        Flux_tmp(:,p,q,1)=                      M1(0,p)*Flux_small(:,0,q,1)+M2(0,p)*Flux_small(:,0,q,2)
        Flux_tmp(:,p,q,2)=                      M1(0,p)*Flux_small(:,0,q,3)+M2(0,p)*Flux_small(:,0,q,4)
        DO l=1,PP_N
          Flux_tmp(:,p,q,1)=Flux_tmp(:,p,q,1) + M1(l,p)*Flux_small(:,l,q,1)+M2(l,p)*Flux_small(:,l,q,2)
          Flux_tmp(:,p,q,2)=Flux_tmp(:,p,q,2) + M1(l,p)*Flux_small(:,l,q,3)+M2(l,p)*Flux_small(:,l,q,4)
        END DO !l=1,PP_N
      END DO !p=0,PP_N
    END DO !q=0,PP_N
    !then in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    Flux_m(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp(iVar,p,:,1) + M2 * Flux_tmp(iVar,p,:,2)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        Flux_master(:,p,q,MortarSideID)=    M1(0,q)*Flux_tmp(:,p,0,1)+M2(0,q)*Flux_tmp(:,p,0,2)
        DO l=1,PP_N
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) &
                                   + M1(l,q)*Flux_tmp(:,p,l,1)+M2(l,q)*Flux_tmp(:,p,l,2)
        END DO !l=1,PP_N
      END DO !p=0,PP_N
    END DO !q=0,PP_N

  CASE(2) !1->2 in eta
    
    !small side surface metric must be scaled by factor of 2, to be the same polynomial, same for the Flux_small,
    ! and then 0.5 must be put in the projection matrix M_1_0, and M_2_0
    M1=0.5*M_1_0; M2=0.5*M_2_0
    Flux_small(:,:,:,1:2) = 2.*Flux_small(:,:,:,1:2)
     
    DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
#ifdef JESSE_MORTAR
      !mortar flux in eta (l,q), for fixed p 
      DO q=0,PP_N  !index big side
        DO l=0,PP_N !index small side                       !<this is the big side>
          CALL MortarFluxAverageVec( U_small(:,p,l,1,iSide), U_small(:,p,q,0,iSide), &
                                    Ns_small(:,p,l,1,iSide),Ns_small(:,p,q,0,iSide),Flux_tp1(:,l,q,1))
          CALL MortarFluxAverageVec( U_small(:,p,l,2,iSide), U_small(:,p,q,0,iSide), &
                                    Ns_small(:,p,l,2,iSide),Ns_small(:,p,q,0,iSide),Flux_tp1(:,l,q,2))
        END DO
      END DO
      DO l=0,PP_N !index small side
        Flux_corr1(:,l,1)=0.
        Flux_corr1(:,l,2)=0.
        DO q=0,PP_N !index big side
          Flux_corr1(:,l,1)= Flux_corr1(:,l,1) + Flux_tp1(:,l,q,1)*M_0_1(q,l) !is the reduction with *lagbaseBig_q(etamortar_l)
          Flux_corr1(:,l,2)= Flux_corr1(:,l,2) + Flux_tp1(:,l,q,2)*M_0_2(q,l) !is the reduction with *lagbaseBig_q(etamortar_l)
        END DO
      END DO
#endif /*JESSE_MORTAR*/

      DO q=0,PP_N !index big side 
        Flux_master(:,p,q,MortarSideID)=0.
        DO l=0,PP_N !index small side
#ifdef JESSE_MORTAR
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) &
                                   + M1(l,q)*(Flux_small(:,p,l,1) + Flux_tp1(:,l,q,1) - Flux_corr1(:,l,1)) &
                                   + M2(l,q)*(Flux_small(:,p,l,2) + Flux_tp1(:,l,q,2) - Flux_corr1(:,l,2))
#else
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) &
                                   + M1(l,q)*(Flux_small(:,p,l,1)) &
                                   + M2(l,q)*(Flux_small(:,p,l,2))
#endif /*JESSE_MORTAR*/
        END DO !l=0,PP_N
      END DO !q=0,PP_N
    END DO !p=0,PP_N

  CASE(3) !1->2 in xi
    !small side surface metric must be scaled by factor of 2, to be the same polynomial, same for the Flux_small,
    ! and then 0.5 must be put in the projection matrix M_1_0, and M_2_0
    M1=0.5*M_1_0; M2=0.5*M_2_0
    Flux_small(:,:,:,1:2) = 2.*Flux_small(:,:,:,1:2)

    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
#ifdef JESSE_MORTAR
      ASSOCIATE(U_big=>U_small(:,:,:,0,iSide),Ns_big=>Ns_small(:,:,:,0,iSide))
      !mortar flux in xi (p,l), for fixed q
      DO p=0,PP_N ! index big side 
        DO l=0,PP_N !index small side                       !<this is the big side!>
          CALL MortarFluxAverageVec( U_small(:,l,q,1,iSide), U_small(:,p,q,0,iSide), &
                                    Ns_small(:,l,q,1,iSide),Ns_small(:,p,q,0,iSide),Flux_tp1(:,l,p,1))
          CALL MortarFluxAverageVec( U_small(:,l,q,2,iSide), U_small(:,p,q,0,iSide), &
                                    Ns_small(:,l,q,2,iSide),Ns_small(:,p,q,0,iSide),Flux_tp1(:,l,p,2))
        END DO
      END DO
      END ASSOCIATE
      DO l=0,PP_N !index small side
        Flux_corr1(:,l,1)=0.
        Flux_corr1(:,l,2)=0.
        DO p=0,PP_N  !index big side
          Flux_corr1(:,l,1)=Flux_corr1(:,l,1) + Flux_tp1(:,l,p,1)*M_0_1(p,l) !is the reduction with *lagbaseBig_p(ximortar_l)
          Flux_corr1(:,l,2)=Flux_corr1(:,l,2) + Flux_tp1(:,l,p,2)*M_0_2(p,l) !is the reduction with *lagbaseBig_p(ximortar_l)
        END DO
      END DO
#endif /*JESSE_MORTAR*/

      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_m(iVar,:,q,MortarSideID)  =   M1 * Flux_m(iVar,:,q,1) + M2 * Flux_m(iVar,:,q,2)
      DO p=0,PP_N !index big side
        Flux_master(:,p,q,MortarSideID)= 0.
        DO l=0,PP_N !index small side
#ifdef JESSE_MORTAR
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID)  &
                                   + M1(l,p)*(Flux_small(:,l,q,1) +Flux_tp1(:,l,p,1)-Flux_corr1(:,l,1)) &
                                   + M2(l,p)*(Flux_small(:,l,q,2) +Flux_tp1(:,l,p,2)-Flux_corr1(:,l,2))
#else
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID)  &
                                   + M1(l,p)*(Flux_small(:,l,q,1)) &
                                   + M2(l,p)*(Flux_small(:,l,q,2))
#endif /*JESSE_MORTAR*/
        END DO !l=0,PP_N
      END DO !p=0,PP_N
    END DO !q=0,PP_N

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID
END SUBROUTINE Flux_Mortar


END MODULE MOD_FillMortar
