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
!> Routines for the constraint transport 

!> Contains the routines to 
!> - initialize and finalize the CT global variables
!> - compute the curl of a projected At
!> 
!==================================================================================================================================
MODULE MOD_CT
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitCT
  MODULE PROCEDURE InitCT
END INTERFACE

INTERFACE CTTimeDerivative
  MODULE PROCEDURE CTTimeDerivative
END INTERFACE

INTERFACE swapB
  MODULE PROCEDURE swapB
END INTERFACE

INTERFACE FinalizeCT
  MODULE PROCEDURE FinalizeCT
END INTERFACE

PUBLIC:: InitCT
PUBLIC:: CTTimeDerivative
PUBLIC:: swapB
PUBLIC:: FinalizeCT
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Allocate all global CT variables
!> 
!==================================================================================================================================
SUBROUTINE InitCT()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CT_Vars
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: MeshInitIsDone,nElems
USE MOD_Restart_Vars,       ONLY: RestartInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone))THEN
   CALL abort(__STAMP__, &
   'InitCT not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT CT...'

ALLOCATE(curlA( 3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(curlAt(3,0:PP_N,0:PP_N,0:PP_N,nElems))
curlA=0.
curlAt=0.

nTotalcurlA=3*(PP_N+1)**3*nElems


IF(.NOT.DoRestart)THEN
  CALL FillCurlA(IniExactFunc,curlA)
ELSE
!  CALL ReadCurlAFromRestart()
END IF

CTInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT CT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitCT

!==================================================================================================================================
!> 
!==================================================================================================================================
SUBROUTINE CTTimeDerivative()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: U       !input!
USE MOD_CT_Vars    ,ONLY: curlAt  !output!
USE MOD_Mesh_Vars  ,ONLY: nElems
USE MOD_Mesh_Vars, ONLY: dXGL_N
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: At(3),AtCov(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_N;DO j=0,PP_N; DO i=0,PP_N
     At =  CROSS(U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem),U(6:8,i,j,k,iElem))  !=v x B = - B x v
     AtCov(:,i,j,k,iElem)= MATMUL(dXGL_N(:,:,i,j,k,iElem),At(:)) !covariant components
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems

CALL ProjectToNedelec(AtCov)
 
!compute curl in logical space, then transform pointwise to cartesian components
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    dAdxi=0.
    dAdeta=0.
    dAdzeta=0.
    DO l=0,PP_N
      dAdxi(:)   = dAdxi(:)   +D(i,l)*Atcov(:,l,j,k,iElem)
      dAdeta(:)  = dAdeta(:)  +D(j,l)*Atcov(:,i,l,k,iElem)
      dAdzeta(:) = dAdzeta(:) +D(k,l)*Atcov(:,i,j,l,iElem)
    END DO !l
    curlAt_ref(1)=(dAdeta( 3)-dAdzeta(2))
    curlAt_ref(2)=(dAdzeta(1)-dAdxi(  3))
    curlAt_ref(3)=(dAdxi(  2)-dAdeta( 1))
    !cartesian components (J Bt)^i = (Ja^i.Bcart) scaled contravariant components
    InvJacMat(1,:)=Metrics_ftilde(:,i,j,k,iElem)
    InvJacMat(2,:)=Metrics_gtilde(:,i,j,k,iElem)
    InvJacMat(3,:)=Metrics_htilde(:,i,j,k,iElem)
    CALL INV33(InvJacMat,JacMat) !JacMat=(InvJacMat)^-1
    curlAt(:,i,j,k,iElem)= MATMUL(JacMat,curlAt_ref(:))
END DO !iElem=1,nElems

END SUBROUTINE CTTimeDerivative


!==================================================================================================================================
!> project discontinuous covariant vector represented by LGL polynomials to Nedelec Finite Element LGL polynomial basis
!> 1 PROC & CARTESIAN ONLY & PERIODIC !!!
!==================================================================================================================================
SUBROUTINE projectToNedelec(vec)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: vec(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL       :: vec_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!==================================================================================================================================
vec_in=vec
DO iElem=1,nElems
  !A3: C0 in x,y
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_MINUS,iElem))       !    y- 
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,nb(0,-1)))    ! x-,y-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,nb(0,-1)))    ! x+,y-
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,iElem))       ! x-
  nb( 0, 0) =iElem  ! 
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,iElem))       ! x+
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_PLUS ,iElem))       !    y+ 
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,nb(0, 1)))    ! x-,y+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,nb(0, 1)))    ! x+,y+

  DO k=0,PP_N
    !corner x-,y-
    vec(3,   0,   0,k,iElem)= 0.25*( vec_in(3,   0,   0,k,nb( 0, 0)) &
                                    +vec_in(3,PP_N,   0,k,nb(-1, 0)) &
                                    +vec_in(3,   0,PP_N,k,nb( 0,-1)) &
                                    +vec_in(3,PP_N,PP_N,k,nb(-1,-1)) ) 
    !corner x-,y+
    vec(3,   0,PP_N,k,iElem)= 0.25*( vec_in(3,   0,PP_N,k,nb( 0, 0)) &
                                    +vec_in(3,PP_N,PP_N,k,nb(-1, 0)) &
                                    +vec_in(3,   0,   0,k,nb( 0, 1)) &
                                    +vec_in(3,PP_N,   0,k,nb(-1, 1)) ) 
    !corner x+,y-
    vec(3,PP_N,   0,k,iElem)= 0.25*( vec_in(3,PP_N,   0,k,nb( 0, 0)) &
                                    +vec_in(3,   0,   0,k,nb( 1, 0)) &
                                    +vec_in(3,PP_N,PP_N,k,nb( 0,-1)) &
                                    +vec_in(3,   0,PP_N,k,nb( 1,-1)) ) 
    !corner x+,y+
    vec(3,PP_N,PP_N,k,iElem)= 0.25*( vec_in(3,PP_N,PP_N,k,nb( 0, 0)) &
                                    +vec_in(3,   0,PP_N,k,nb( 1, 0)) &
                                    +vec_in(3,PP_N,   0,k,nb( 0, 1)) &
                                    +vec_in(3,   0,   0,k,nb( 1, 1)) ) 
    !inner surface points
    DO j=1,PP_N-1 
      vec(3,   0, j,k,iElem)=0.5*( vec_in(3,   0,   j,k,nb( 0, 0)) &
                                  +vec_in(3,PP_N,   j,k,nb(-1, 0)) ) !x-
      vec(3,PP_N, j,k,iElem)=0.5*( vec_in(3,PP_N,   j,k,nb( 0, 0)) &
                                  +vec_in(3,   0,   j,k,nb( 1, 0)) ) !x+
    END DO
    DO i=1,PP_N-1 
      vec(3, i,   0,k,iElem)=0.5*( vec_in(3,   i,   0,k,nb( 0, 0)) &
                                  +vec_in(3,   i,PP_N,k,nb( 0,-1)) ) !y-
      vec(3, i,PP_N,k,iElem)=0.5*( vec_in(3,   i,PP_N,k,nb( 0, 0)) &
                                  +vec_in(3,   i,   0,k,nb( 0, 1)) ) !y+
    END DO
  END DO!k
  !A2: C0 in x,z
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,ZETA_MINUS,iElem))       !    z- 
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,nb(0,-1)))    ! x-,z-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,nb(0,-1)))    ! x+,z-
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,iElem))       ! x-
  nb( 0, 0) =iElem  ! 
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,iElem))       ! x+
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,ZETA_PLUS ,iElem))       !    z+ 
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_MINUS,nb(0, 1)))    ! x-,z+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,  XI_PLUS ,nb(0, 1)))    ! x+,z+

  DO j=0,PP_N
    !corner x-,z-
    vec(2,   0,j,   0,iElem)= 0.25*( vec_in(2,   0,j,   0,nb( 0, 0)) &
                                    +vec_in(2,PP_N,j,   0,nb(-1, 0)) &
                                    +vec_in(2,   0,j,PP_N,nb( 0,-1)) &
                                    +vec_in(2,PP_N,j,PP_N,nb(-1,-1)) ) 
    !corner x-,z+                                   
    vec(2,   0,j,PP_N,iElem)= 0.25*( vec_in(2,   0,j,PP_N,nb( 0, 0)) &
                                    +vec_in(2,PP_N,j,PP_N,nb(-1, 0)) &
                                    +vec_in(2,   0,j,   0,nb( 0, 1)) &
                                    +vec_in(2,PP_N,j,   0,nb(-1, 1)) ) 
    !corner x+,z-                                   
    vec(2,PP_N,j,   0,iElem)= 0.25*( vec_in(2,PP_N,j,   0,nb( 0, 0)) &
                                    +vec_in(2,   0,j,   0,nb( 1, 0)) &
                                    +vec_in(2,PP_N,j,PP_N,nb( 0,-1)) &
                                    +vec_in(2,   0,j,PP_N,nb( 1,-1)) ) 
    !corner x+,z+                                   
    vec(2,PP_N,j,PP_N,iElem)= 0.25*( vec_in(2,PP_N,j,PP_N,nb( 0, 0)) &
                                    +vec_in(2,   0,j,PP_N,nb( 1, 0)) &
                                    +vec_in(2,PP_N,j,   0,nb( 0, 1)) &
                                    +vec_in(2,   0,j,   0,nb( 1, 1)) ) 
    !inner surface points
    DO k=1,PP_N-1 
      vec(2,   0, j,k,iElem)=0.5*( vec_in(2,   0,   j,k,nb( 0, 0)) &
                                  +vec_in(2,PP_N,   j,k,nb(-1, 0)) ) !x-
      vec(2,PP_N, j,k,iElem)=0.5*( vec_in(2,PP_N,   j,k,nb( 0, 0)) &
                                  +vec_in(2,   0,   j,k,nb( 1, 0)) ) !x+
    END DO
    DO i=1,PP_N-1 
      vec(2, i,j,   0,iElem)=0.5*( vec_in(2,   i,j,   0,nb( 0, 0)) &
                                  +vec_in(2,   i,j,PP_N,nb( 0,-1)) ) !z-
      vec(2, i,j,PP_N,iElem)=0.5*( vec_in(2,   i,j,PP_N,nb( 0, 0)) &
                                  +vec_in(2,   i,j,   0,nb( 0, 1)) ) !z+
    END DO
  END DO!j
  !A1: C0 in y,z
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,ZETA_MINUS,iElem))       !    z- 
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_MINUS,nb(0,-1)))    ! y-,z-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_PLUS ,nb(0,-1)))    ! y+,z-
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_MINUS,iElem))       ! y-
  nb( 0, 0) =iElem  ! 
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_PLUS ,iElem))       ! y+
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID,ZETA_PLUS ,iElem))       !    z+ 
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_MINUS,nb(0, 1)))    ! y-,z+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SideID, ETA_PLUS ,nb(0, 1)))    ! y+,z+

  DO i=0,PP_N
    !corner y-,z-
    vec(1,i,   0,   0,iElem)= 0.25*( vec_in(1,i,   0,   0,nb( 0, 0)) &
                                    +vec_in(1,i,PP_N,   0,nb(-1, 0)) &
                                    +vec_in(1,i,   0,PP_N,nb( 0,-1)) &
                                    +vec_in(1,i,PP_N,PP_N,nb(-1,-1)) ) 
    !corner y-,z+                                   
    vec(1,i,   0,PP_N,iElem)= 0.25*( vec_in(1,i,   0,PP_N,nb( 0, 0)) &
                                    +vec_in(1,i,PP_N,PP_N,nb(-1, 0)) &
                                    +vec_in(1,i,   0,   0,nb( 0, 1)) &
                                    +vec_in(1,i,PP_N,   0,nb(-1, 1)) ) 
    !corner y+,z-                                   
    vec(1,i,PP_N,   0,iElem)= 0.25*( vec_in(1,i,PP_N,   0,nb( 0, 0)) &
                                    +vec_in(1,i,   0,   0,nb( 1, 0)) &
                                    +vec_in(1,i,PP_N,PP_N,nb( 0,-1)) &
                                    +vec_in(1,i,   0,PP_N,nb( 1,-1)) ) 
    !corner y+,z+                                   
    vec(1,i,PP_N,PP_N,iElem)= 0.25*( vec_in(1,i,PP_N,PP_N,nb( 0, 0)) &
                                    +vec_in(1,i,   0,PP_N,nb( 1, 0)) &
                                    +vec_in(1,i,PP_N,   0,nb( 0, 1)) &
                                    +vec_in(1,i,   0,   0,nb( 1, 1)) ) 
    !inner surface points
    DO k=1,PP_N-1 
      vec(1, i,   0,k,iElem)=0.5*( vec_in(1,   i,   0,k,nb( 0, 0)) &
                                  +vec_in(1,   i,PP_N,k,nb(-1, 0)) ) !y-
      vec(1, i,PP_N,k,iElem)=0.5*( vec_in(1,   i,PP_N,k,nb( 0, 0)) &
                                  +vec_in(1,   i,   0,k,nb( 1, 0)) ) !y+
    END DO
    DO j=1,PP_N-1 
      vec(1, i,j,   0,iElem)=0.5*( vec_in(1,   i,j,   0,nb( 0, 0)) &
                                  +vec_in(1,   i,j,PP_N,nb( 0,-1)) ) !z-
      vec(1, i,j,PP_N,iElem)=0.5*( vec_in(1,   i,j,PP_N,nb( 0, 0)) &
                                  +vec_in(1,   i,j,   0,nb( 0, 1)) ) !z+
    END DO
  END DO!i
  
END DO !iElem=1,nElems

END SUBROUTINE projecTToNedelec


!==================================================================================================================================
!> Replaces B of current U by curlA, and changing total energy and keeping the same pressure
!==================================================================================================================================
SUBROUTINE swapB()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: U      !input/output!
USE MOD_CT_Vars    ,ONLY: curlA  !input!
USE MOD_Mesh_Vars  ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: EmagU,EmagCurlA
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_N;DO j=0,PP_N; DO i=0,PP_N
     EmagU=s2mu0*SUM(U(6:8,i,j,k,iElem)**2)
     EmagcurlA=s2mu0*SUM(curlA(:,i,j,k,iElem)**2)
     U(5,i,j,k,iElem)=U(5,i,j,k,iElem)+EmagcurlA-EmagU
     U(6:8,i,j,k,iElem)=curlA(:,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems

END SUBROUTINE swapB



!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeCT()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_CT_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(CurlAt)
SDEALLOCATE(CurlA)
CTInitIsDone = .FALSE.
END SUBROUTINE FinalizeCT

END MODULE MOD_CT
