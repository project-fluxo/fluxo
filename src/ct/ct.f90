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

INTERFACE CT_TimeDerivative
  MODULE PROCEDURE CT_TimeDerivative
END INTERFACE

!INTERFACE swapB
!  MODULE PROCEDURE swapB
!END INTERFACE

INTERFACE FinalizeCT
  MODULE PROCEDURE FinalizeCT
END INTERFACE

PUBLIC:: InitCT
PUBLIC:: CT_TimeDerivative
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
USE MOD_Mesh_Vars,          ONLY: Metrics_ftilde,Metrics_gtilde,Metrics_htilde,Elem_xGP
USE MOD_Restart_Vars,       ONLY: RestartInitIsDone,doRestart
USE MOD_Basis,              ONLY: Inv33
USE MOD_Equation_Vars,      ONLY: IniExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i,j,k,iElem
REAL    :: InvJacMat(3,3)
REAL    :: Acart(3,0:PP_N,0:PP_N,0:PP_N,nElems)
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

ALLOCATE(JacMat(3,3,0:PP_N,0:PP_N,0:PP_N,nElems))
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    InvJacMat(1,:)=Metrics_ftilde(:,i,j,k,iElem)
    InvJacMat(2,:)=Metrics_gtilde(:,i,j,k,iElem)
    InvJacMat(3,:)=Metrics_htilde(:,i,j,k,iElem)
    CALL INV33(InvJacMat,JacMat(:,:,i,j,k,iElem)) !JacMat=(InvJacMat)^-1
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems

#if PP_NodeType != 2
STOP 'constraint transport not implemented for Gauss points!'
#endif

IF(.NOT.DoRestart)THEN
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      CALL EvalMagneticVectorPotential(IniExactFunc,Elem_xGP(:,i,j,k,iElem),Acart(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k=0,PP_N
  END DO !iElem=1,nElems
  !since computed at GL points, already continuous
  CALL ComputeCartCurl(Acart,curlA,vec_in_Base=0)
ELSE
!  CALL ReadCurlAFromRestart()
END IF

CTInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT CT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitCT

!==================================================================================================================================
!> compute ODE time derivative of A, A_t=-(B x v [+eta J] )
!> project to Necelec base (in covariant components!), compute curl back to cartesian components. 
!==================================================================================================================================
SUBROUTINE CT_TimeDerivative()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: U       !! input  !!
USE MOD_CT_Vars    ,ONLY: curlAt  !! output !!
USE MOD_Mesh_Vars  ,ONLY: nElems,dXGL_N
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,l,iElem
REAL    :: AtCov(3,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL    :: At(3)
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_N;DO j=0,PP_N; DO i=0,PP_N
     At =  CROSS(U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem),U(6:8,i,j,k,iElem))  !=v x B = - B x v
     AtCov(:,i,j,k,iElem)= MATMUL(dXGL_N(:,:,i,j,k,iElem),At(:)) !covariant components
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems

CALL ProjectToNedelec(AtCov)
 
CALL ComputeCartCurl(AtCov,curlAt,vec_in_Base=1)

END SUBROUTINE CT_TimeDerivative

!==================================================================================================================================
!> compute the cartesian components of the curl, only locally per element,
!> such that the transformed divergence in the volume is discretely zero (1/J sum_i d/dxi^i Ja^i . v), using the inverse of Ja^i .
!> input can be cartesian or covariant 
!==================================================================================================================================
SUBROUTINE ComputeCartCurl(vec_in,curl_cart,vec_in_Base)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: D
USE MOD_CT_Vars    ,ONLY: JacMat
USE MOD_Mesh_Vars  ,ONLY: nElems,dXGL_N
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: vec_in_Base                           !< base of input vector: 0: cartesian components
                                                              !< 1: covariant components v_i = v . acov_i
REAL,INTENT(IN)      :: vec_in(3,0:PP_N,0:PP_N,0:PP_N,nElems) !< vector in covariant components 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)     :: curl_cart(3,0:PP_N,0:PP_N,0:PP_N,nElems) !< curl of vector in cartesian components 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,l,iElem
REAL    :: dvdxi(3),dvdeta(3),dvdzeta(3)
REAL    :: vec_cov(3,0:PP_N,0:PP_N,0:PP_N) 
!==================================================================================================================================
DO iElem=1,nElems
  SELECT CASE(vec_in_Base)
  CASE(0) !cartesian components are input
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      vec_cov(:,i,j,k)=MATMUL(dXGL_N(:,:,i,j,k,iElem),vec_in(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k=0,PP_N                      
  CASE(1) !covariant components are input
    vec_cov(:,:,:,:)=vec_in(:,:,:,:,iElem)
  END SELECT !vec_in_Base
  
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    dvdxi(:)   = D(i,0)*vec_cov(:,0,j,k)
    dvdeta(:)  = D(j,0)*vec_cov(:,i,0,k)
    dvdzeta(:) = D(k,0)*vec_cov(:,i,j,0)
    DO l=1,PP_N
      dvdxi(:)   = dvdxi(:)   +D(i,l)*vec_cov(:,l,j,k)
      dvdeta(:)  = dvdeta(:)  +D(j,l)*vec_cov(:,i,l,k)
      dvdzeta(:) = dvdzeta(:) +D(k,l)*vec_cov(:,i,j,l)
    END DO !l
    !cartesian components (J Bt)^i = (Ja^i.Bcart) scaled contravariant components
    curl_cart(:,i,j,k,iElem)= MATMUL(JacMat(:,:,i,j,k,iElem),(/(dvdeta( 3)-dvdzeta(2)), &
                                                               (dvdzeta(1)-dvdxi(  3)), &
                                                               (dvdxi(  2)-dvdeta( 1)) /))
  END DO; END DO; END DO !i,j,k=0,PP_N                      
END DO !iElem=1,nElems

END SUBROUTINE ComputeCartCurl


!==================================================================================================================================
!> project discontinuous covariant vector represented by LGL polynomials to Nedelec Finite Element LGL polynomial basis
!> !!! WORKS ONLY FOR 1 PROC &  PERIODIC STRUCTURED BLOCK MESH !!!
!==================================================================================================================================
SUBROUTINE projectToNedelec(vec)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY: nElems,SideToElem,ElemToSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: vec(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL       :: vec_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
INTEGER    :: i,j,k,iElem,nb(-1:1,-1:1),xx,yy,zz
!==================================================================================================================================
vec_in=vec
DO iElem=1,nElems
  !A1: C0 in y,z
  nb( 0, 0) =iElem  ! 
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))       !    z- 
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,ZETA_PLUS ,iElem))       !    z+ 
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_MINUS,iElem))       ! y-
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_PLUS ,iElem))       ! y+
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_MINUS,nb(0,-1)))    ! y-,z-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_PLUS ,nb(0,-1)))    ! y+,z-
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_MINUS,nb(0, 1)))    ! y-,z+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_PLUS ,nb(0, 1)))    ! y+,z+
  !corners
  DO yy=0,1; DO zz=0,1
    vec(1,:,yy*PP_N,zz*PP_N,iElem)= 0.25*( vec_in(1,:,   0,   0,nb(yy  ,zz  )) &
                                          +vec_in(1,:,   0,PP_N,nb(yy  ,zz-1)) &
                                          +vec_in(1,:,PP_N,   0,nb(yy-1,zz  )) &
                                          +vec_in(1,:,PP_N,PP_N,nb(yy-1,zz-1)) )
  END DO; END DO !yy,zz
  !inner surface points
  DO k=1,PP_N-1 
    vec(1, :,   0,k,iElem)=0.5*( vec_in(1,   :,   0,k,nb( 0, 0)) &
                                +vec_in(1,   :,PP_N,k,nb(-1, 0)) ) !y-
    vec(1, :,PP_N,k,iElem)=0.5*( vec_in(1,   :,PP_N,k,nb( 0, 0)) &
                                +vec_in(1,   :,   0,k,nb( 1, 0)) ) !y+
  END DO
  DO j=1,PP_N-1 
    vec(1, :,j,   0,iElem)=0.5*( vec_in(1,   :,j,   0,nb( 0, 0)) &
                                +vec_in(1,   :,j,PP_N,nb( 0,-1)) ) !z-
    vec(1, :,j,PP_N,iElem)=0.5*( vec_in(1,   :,j,PP_N,nb( 0, 0)) &
                                +vec_in(1,   :,j,   0,nb( 0, 1)) ) !z+
  END DO
  !A2: C0 in x,z
  nb( 0, 0) =iElem  ! 
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))       !    z- 
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,ZETA_PLUS ,iElem))       !    z+ 
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,iElem))       ! x-
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,iElem))       ! x+
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,nb(0,-1)))    ! x-,z-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,nb(0,-1)))    ! x+,z-
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,nb(0, 1)))    ! x-,z+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,nb(0, 1)))    ! x+,z+
  !corners
  DO xx=0,1; DO zz=0,1
    vec(2,xx*PP_N,:,zz*PP_N,iElem)= 0.25*( vec_in(2,   0,:,   0,nb(xx  ,zz  )) &
                                          +vec_in(2,   0,:,PP_N,nb(xx  ,zz-1)) &
                                          +vec_in(2,PP_N,:,   0,nb(xx-1,zz  )) &
                                          +vec_in(2,PP_N,:,PP_N,nb(xx-1,zz-1)) )
  END DO; END DO !xx,zz
  !inner surface points
  DO k=1,PP_N-1 
    vec(2,   0, :,k,iElem)=0.5*( vec_in(2,   0,   :,k,nb( 0, 0)) &
                                +vec_in(2,PP_N,   :,k,nb(-1, 0)) ) !x-
    vec(2,PP_N, :,k,iElem)=0.5*( vec_in(2,PP_N,   :,k,nb( 0, 0)) &
                                +vec_in(2,   0,   :,k,nb( 1, 0)) ) !x+
  END DO
  DO i=1,PP_N-1 
    vec(2, i,:,   0,iElem)=0.5*( vec_in(2,   i,:,   0,nb( 0, 0)) &
                                +vec_in(2,   i,:,PP_N,nb( 0,-1)) ) !z-
    vec(2, i,:,PP_N,iElem)=0.5*( vec_in(2,   i,:,PP_N,nb( 0, 0)) &
                                +vec_in(2,   i,:,   0,nb( 0, 1)) ) !z+
  END DO
  !A3: C0 in x,y
  nb( 0, 0) =iElem  ! 
  nb( 0,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_MINUS,iElem))       !    y- 
  nb( 0, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID, ETA_PLUS ,iElem))       !    y+ 
  nb(-1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,iElem))       ! x-
  nb( 1, 0) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,iElem))       ! x+
  nb(-1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,nb(0,-1)))    ! x-,y-
  nb( 1,-1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,nb(0,-1)))    ! x+,y-
  nb(-1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_MINUS,nb(0, 1)))    ! x-,y+
  nb( 1, 1) =SideToElem(S2E_NB_ELEM_ID,ElemToSide(E2S_SIDE_ID,  XI_PLUS ,nb(0, 1)))    ! x+,y+
  !corners
  DO xx=0,1; DO yy=0,1
    vec(3,xx*PP_N,yy*PP_N,:,iElem)= 0.25*( vec_in(3,   0,   0,:,nb(xx  ,yy  )) &
                                          +vec_in(3,   0,PP_N,:,nb(xx  ,yy-1)) &
                                          +vec_in(3,PP_N,   0,:,nb(xx-1,yy  )) &
                                          +vec_in(3,PP_N,PP_N,:,nb(xx-1,yy-1)) )
  END DO; END DO !xx,yy
  !inner surface points
  DO j=1,PP_N-1 
    vec(3,   0, j,:,iElem)=0.5*( vec_in(3,   0,   j,:,nb( 0, 0)) &
                                +vec_in(3,PP_N,   j,:,nb(-1, 0)) ) !x-
    vec(3,PP_N, j,:,iElem)=0.5*( vec_in(3,PP_N,   j,:,nb( 0, 0)) &
                                +vec_in(3,   0,   j,:,nb( 1, 0)) ) !x+
  END DO
  DO i=1,PP_N-1 
    vec(3, i,   0,:,iElem)=0.5*( vec_in(3,   i,   0,:,nb( 0, 0)) &
                                +vec_in(3,   i,PP_N,:,nb( 0,-1)) ) !y-
    vec(3, i,PP_N,:,iElem)=0.5*( vec_in(3,   i,PP_N,:,nb( 0, 0)) &
                                +vec_in(3,   i,   0,:,nb( 0, 1)) ) !y+
  END DO
  
END DO !iElem=1,nElems

END SUBROUTINE projectToNedelec


!==================================================================================================================================
!> evaluate the magnetic vector potential A, belonging to the magnetic field
!> of a given IniExactfunction 
!> B_1= d(A_3)/dy - d(A_2)/dz
!> B_2= d(A_1)/dz - d(A_3)/dx
!> B_3= d(A_2)/dx - d(A_1)/dy
!==================================================================================================================================
SUBROUTINE EvalMagneticVectorPotential(ExactFunc,x,Acart)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0,IniHalfwidth
USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: ExactFunc
REAL   ,INTENT(IN)  :: x(3) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Acart(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: B(3)
!==================================================================================================================================
SELECT CASE(ExactFunc)
CASE(73)
  B(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
  !B2=0.
  B(3)=1.0
  Acart(1) = 0.
  Acart(2) = x(1)*B(3) - x(3)*B(1)
  Acart(3) = 0. 
CASE(74)
  B(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
  !B2=0.
  B(3)=SQRT(1.0-B(1)*B(1))
  Acart(1) = 0.
  Acart(2) = x(1)*B(3)-x(3)*B(1)
  Acart(3) = 0. 

CASE DEFAULT
  STOP 'this exactfunc is not implemented in EvalMagneticVectorPotential'
END SELECT !ExactFunc

END SUBROUTINE EvalMagneticVectorPotential


!==================================================================================================================================
!> Replaces B of current U by curlA, and changing total energy and keeping the same pressure
!==================================================================================================================================
SUBROUTINE swapB(U_inout,curlA_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0
USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: curlA_in(3,ntotal_IP)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_inout(PP_nVar,ntotal_IP)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: EmagU,EmagCurlA
!==================================================================================================================================
DO i=1,nTotal_IP
     EmagU=SUM(U_inout(6:8,i)**2)
     EmagcurlA=SUM(curlA_in(:,i)**2)
     U_inout(5,i)=U_inout(5,i)+s2mu_0*(EmagcurlA-EmagU) !add new magn. energy, substract old energy
     U_inout(6:8,i)=curlA_in(:,i)
END DO !i=1,nTotal_IP

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
SDEALLOCATE(JacMat)
CTInitIsDone = .FALSE.
END SUBROUTINE FinalizeCT

END MODULE MOD_CT
