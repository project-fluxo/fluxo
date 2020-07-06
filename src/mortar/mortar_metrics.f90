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
!> \brief Routines providing support for geometric features of non-conforming meshes (generally a preprocessing step)
!==================================================================================================================================
MODULE MOD_Mortar_Metrics
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Mortar_CalcSurfMetrics
  MODULE PROCEDURE Mortar_CalcSurfMetrics
END INTERFACE

PUBLIC::Mortar_CalcSurfMetrics

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Interpolates surface metrics Ja and xGP on face from master to slave (small) mortar sides.
!> 1D interpolation operators M_0_1,M_0_2 are built locally per polynomial degree.
!>
!> Already existing surface metrics are overwritten, metrics for small sides are built from
!> big (master) side, i.e. all small sides belonging to a mortar interface are slave sides 
!> (with inward pointing normal vector). NOTE THAT THIS IS NOT THE CASE FOR MPI_YOUR MORTAR SIDES!
!> In an MPI setting if the big sides are not present on a CPU and this CPU has small master sides
!> they are not rebuilt and fluxes need to be rotated at the big mortar.
!>
!>~~~~~~~~~~~~~~~~~~~~
!>       Type 1               Type 2              Type3
!>        eta                  eta                 eta
!>         ^                    ^                   ^
!>         |                    |                   |
!>     +---+---+            +---+---+           +---+---+
!>     | 3 | 4 |            |   2   |           |   |   |
!>     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>     | 1 | 2 |            |   1   |           |   |   |
!>     +---+---+            +---+---+           +---+---+
!>~~~~~~~~~~~~~~~~~~~~
!>
!>
!>==================================================================================================================================
SUBROUTINE Mortar_CalcSurfMetrics(SideID,Nloc,Face_Ja,Face_xGP,&
                                  Mortar_Ja,Mortar_xGP,nbSideID)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar,      ONLY: MortarBasis_BigToSmall
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
#ifdef JESSE_MORTAR
USE MOD_Mortar_Vars, ONLY: Ns_small !,M_0_1,M_0_2
USE MOD_Mesh_Vars,   ONLY: SideToElem,NormalDirs,NormalSigns
#endif /*JESSE_MORTAR*/
USE MOD_Interpolation_Vars,ONLY: NodeType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                         !< SideID of mortar master side
INTEGER,INTENT(IN) :: Nloc                           !< polynomial degree
REAL,INTENT(IN)    :: Face_Ja(  3,3,0:Nloc,0:Nloc)   !< surface metrics of side
REAL,INTENT(IN)    :: Face_xGP(   3,0:Nloc,0:Nloc)   !< face xGP
REAL,INTENT(OUT)   :: Mortar_Ja(3,3,0:Nloc,0:Nloc,4) !< mortarized surface metrics of side
REAL,INTENT(OUT)   :: Mortar_xGP( 3,0:Nloc,0:Nloc,4) !< mortarized face xGP
INTEGER,INTENT(OUT):: nbSideID(4)                    !< index of neighbour sideIDs
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,dir1,dir2,iNb,jNb,ind,SideIDMortar
REAL     :: M_0_12(0:Nloc,0:Nloc,2),M_0_12_h(0:Nloc,0:Nloc,2)
REAL     :: Mortar_Ja2(1:3,1:3,0:Nloc,0:Nloc)
REAL     :: Mortar_xGP2 (  1:3,0:Nloc,0:Nloc)
#ifdef JESSE_MORTAR
INTEGER  :: iLocSide,NormalDir !,l
REAL     :: NormalSign 
#endif /*JESSE_MORTAR*/
!==================================================================================================================================
CALL MortarBasis_BigToSmall(Nloc,NodeType,M_0_12(:,:,1),M_0_12(:,:,2))
! ATTENTION: MortarBasis_BigToSmall computes the transposed matrices, which is useful when they are used
!            in hand-written matrix multiplications. For the use with the intrinsic MATMUL, they must be transposed.
M_0_12_h(:,:,1)=0.5*TRANSPOSE(M_0_12(:,:,1))
M_0_12_h(:,:,2)=0.5*TRANSPOSE(M_0_12(:,:,2))

nbSideID=-1

! Surface metrics derived from big sides are only built for inner sides and MPI_MINE sides!
SideIDMortar=MortarType(2,SideID)
#ifdef JESSE_MORTAR
IF(Nloc.NE.PP_N) CALL abort(__STAMP__, &
     'Nloc ne PP_N in mortar metrics')
iLocSide   = SideToElem(S2E_LOC_SIDE_ID,SideID)
NormalDir  = NormalDirs(iLocSide)
NormalSign = NormalSigns(iLocSide)
#endif /*JESSE_MORTAR*/

#ifdef JESSE_MORTAR
Ns_small(:,:,:,0,SideIDMortar) = NormalSign*Face_Ja(NormalDir,:,:,:) !big mortar metric
#endif /*JESSE_MORTAR*/
SELECT CASE(MortarType(1,SideID))
CASE(1) !1->4
  !inb=1,jNb=1 > Nb=1
  !inb=2,jNb=1 > Nb=2
  !inb=1,jNb=2 > Nb=3
  !inb=2,jNb=2 > Nb=4
  !first in eta
  DO iNb=1,2
    DO p=0,Nloc
      DO dir1=1,3
        DO dir2=1,3
          Mortar_Ja2(dir1,dir2,p,:)=MATMUL(M_0_12_h(:,:,iNb),Face_Ja(dir1,dir2,p,:))
        END DO !dir2=1,3
        Mortar_xGP2(dir1,p,:)      =MATMUL(TRANSPOSE(M_0_12(  :,:,iNb)),Face_xGP(dir1,p,:))
      END DO !dir1=1,3
    END DO !p=0,Nloc
#ifdef JESSE_MORTAR
   Ns_small(:,:,:,iNb-3,SideIDMortar) = 2.*NormalSign*Mortar_Ja2(NormalDir,:,:,:) !revert 0.5 from M_0_12_h
#endif /*JESSE_MORTAR*/
    !now in xi
    DO jNb=1,2
      ind=iNb+2*(jNb-1)
!      IF(MortarInfo(E2S_FLIP,ind,SideIDMortar).GT.0) CYCLE !no slave sides (MPI)
      IF(.NOT.(MortarInfo(E2S_FLIP,ind,SideIDMortar).GT.0)) & !no slave sides (MPI)
        nbSideID(ind)=MortarInfo(E2S_SIDE_ID,ind,SideIDMortar)

      DO q=0,Nloc
        DO dir1=1,3
          DO dir2=1,3
            Mortar_Ja(dir1,dir2,:,q,ind)=MATMUL(M_0_12_h(:,:,jNb),Mortar_Ja2(dir1,dir2,:,q))
          END DO !dir2=1,3
          Mortar_xGP(dir1,:,q,ind)      =MATMUL(TRANSPOSE(M_0_12(  :,:,jNb)),Mortar_xGP2(dir1,:,q))
        END DO !dir1=1,3
      END DO !q=0,Nloc
#ifdef JESSE_MORTAR
      Ns_small(:,:,:,ind,SideIDMortar) = 4.*NormalSign*Mortar_Ja(NormalDir,:,:,:,ind) !revert 0.5*0.5 from twice M_0_12_h
#endif /*JESSE_MORTAR*/
    END DO !jNb
  END DO !iNb
!#ifdef JESSE_MORTAR
!    Ns_small(:,:,:,-2:-1,SideIDMortar)=0.
!    Ns_small(:,:,:,1:4,SideIDMortar)=0.
!    !first  split 1 side into two, in eta direction
!    DO q=0,PP_N
!      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
!        DO l=0,PP_N
!          Ns_small(:,p,q,1-3,SideIDMortar)=Ns_small(:,p,q,1-3,SideIDMortar)+M_0_1(l,q)*Ns_small(:,p,l,0,SideIDMortar)
!          Ns_small(:,p,q,2-3,SideIDMortar)=Ns_small(:,p,q,2-3,SideIDMortar)+M_0_2(l,q)*Ns_small(:,p,l,0,SideIDMortar)
!        END DO
!      END DO
!    END DO
!    ! then split each side again into two, now in xi direction
!    DO iNb=1,2
!      DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
!        DO p=0,PP_N
!          DO l=0,PP_N
!            Ns_small(:,p,q,1+2*(iNb-1),SideIDMortar)=Ns_small(:,p,q,1+2*(iNb-1),SideIDMortar)+M_0_1(l,p)*Ns_small(:,l,q,iNb-3,SideIDMortar)
!            Ns_small(:,p,q,2+2*(iNb-1),SideIDMortar)=Ns_small(:,p,q,2+2*(iNb-1),SideIDMortar)+M_0_2(l,p)*Ns_small(:,l,q,iNb-3,SideIDMortar)
!          END DO !l=1,PP_N
!        END DO
!      END DO 
!    END DO !iNb=1,2
!#endif /*JESSE_MORTAR*/

CASE(2) !1->2 in eta
  DO jNb=1,2
!    IF(MortarInfo(E2S_FLIP,jNb,SideIDMortar).GT.0) CYCLE !no slave sides (MPI)
    IF(.NOT.(MortarInfo(E2S_FLIP,jNb,SideIDMortar).GT.0)) & !no slave sides (MPI)
      nbSideID(jNb)=MortarInfo(E2S_SIDE_ID,jNb,SideIDMortar)

    DO p=0,Nloc
      DO dir1=1,3
        DO dir2=1,3
          Mortar_Ja(dir1,dir2,p,:,jNb)=MATMUL(M_0_12_h(:,:,jNb),Face_Ja(dir1,dir2,p,:))
        END DO !dir2=1,3
        Mortar_xGP(dir1,p,:,jNb)      =MATMUL(TRANSPOSE(M_0_12(  :,:,jNb)),Face_xGP(dir1,p,:))
      END DO !dir1=1,3
    END DO !p=0,Nloc
#ifdef JESSE_MORTAR
  Ns_small(:,:,:,jNb,SideIDMortar) = 2.*NormalSign*Mortar_Ja(NormalDir,:,:,:,jNb) !revert 0.5 from M_0_12_h
#endif /*JESSE_MORTAR*/
  END DO !jNb
!#ifdef JESSE_MORTAR
!  Ns_small(:,:,:,1:2,SideIDMortar)=0.
!  DO q=0,PP_N
!    DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
!      DO l=0,PP_N
!        Ns_small(:,p,q,1,SideIDMortar)=Ns_small(:,p,q,1,SideIDMortar)+M_0_1(l,q)*Ns_small(:,p,l,0,SideIDMortar)
!        Ns_small(:,p,q,2,SideIDMortar)=Ns_small(:,p,q,2,SideIDMortar)+M_0_2(l,q)*Ns_small(:,p,l,0,SideIDMortar)
!      END DO
!    END DO
!  END DO
!#endif /*JESSE_MORTAR*/

CASE(3) !1->2 in xi
  DO iNb=1,2
!    IF(MortarInfo(E2S_FLIP,iNb,SideIDMortar).GT.0) CYCLE !no slave sides (MPI)
    IF(.NOT.(MortarInfo(E2S_FLIP,iNb,SideIDMortar).GT.0)) & !no slave sides (MPI)
      nbSideID(iNb)=MortarInfo(E2S_SIDE_ID,iNb,SideIDMortar)

    DO q=0,Nloc
      DO dir1=1,3
        DO dir2=1,3
          Mortar_Ja(dir1,dir2,:,q,iNb)=MATMUL(M_0_12_h(:,:,iNb),Face_Ja(dir1,dir2,:,q))
        END DO !dir2=1,3
        Mortar_xGP(dir1,:,q,iNb)      =MATMUL(TRANSPOSE(M_0_12(  :,:,iNb)),Face_xGP(dir1,:,q))
      END DO !dir1=1,3
    END DO !q=0,Nloc
#ifdef JESSE_MORTAR
    Ns_small(:,:,:,iNb,SideIDMortar) = 2.*NormalSign*Mortar_Ja(NormalDir,:,:,:,iNb) !revert 0.5 from M_0_12_h
#endif /*JESSE_MORTAR*/
  END DO !iNb
!#ifdef JESSE_MORTAR
!  Ns_small(:,:,:,1:2,SideIDMortar)=0.
!  DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
!    DO p=0,PP_N
!      DO l=0,PP_N
!        Ns_small(:,p,q,1,SideIDMortar)=Ns_small(:,p,q,1,SideIDMortar)+M_0_1(l,p)*Ns_small(:,l,q,0,SideIDMortar)
!        Ns_small(:,p,q,2,SideIDMortar)=Ns_small(:,p,q,2,SideIDMortar)+M_0_2(l,p)*Ns_small(:,l,q,0,SideIDMortar)
!      END DO
!    END DO
!  END DO
!#endif /*JESSE_MORTAR*/

END SELECT !MortarType
END SUBROUTINE Mortar_CalcSurfMetrics

END MODULE MOD_Mortar_Metrics
