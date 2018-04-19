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

INTERFACE InitB_CT
  MODULE PROCEDURE InitB_CT
END INTERFACE

INTERFACE CT_TimeDerivative
  MODULE PROCEDURE CT_TimeDerivative
END INTERFACE

!INTERFACE swapBt
!  MODULE PROCEDURE swapBt
!END INTERFACE

INTERFACE FinalizeCT
  MODULE PROCEDURE FinalizeCT
END INTERFACE

PUBLIC:: InitCT
PUBLIC:: InitB_CT
PUBLIC:: CT_TimeDerivative
PUBLIC:: swapBt
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
USE MOD_Mesh_Vars,          ONLY: Metrics_ftilde,Metrics_gtilde,Metrics_htilde
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

doEvalSource_A=.TRUE. !default

CTInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT CT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitCT


!==================================================================================================================================
!> Initialize divergence free magnetic field
!> 
!==================================================================================================================================
SUBROUTINE InitB_CT(ExactFunc,U_inout)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CT_Vars
USE MOD_Mesh_Vars,          ONLY: nElems,Elem_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: exactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: U_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i,j,k,iElem
REAL    :: Acart(3,0:PP_N,0:PP_N,0:PP_N,nElems),constB(3)
LOGICAL :: discreteDivB
!===================================================================================================================================
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      CALL EvalExactFunc_A(ExactFunc,Elem_xGP(:,i,j,k,iElem),Acart(:,i,j,k,iElem), &
                                                  curlA(:,i,j,k,iElem),constB,discreteDivB)
    END DO; END DO; END DO !i,j,k=0,PP_N
  END DO !iElem=1,nElems
  IF(.NOT.discreteDivB)THEN
    !TEST
      curlA=Acart
      CALL ProjectToNedelec(Acart)
      WRITE(*,*)'TESTTEST: (Acart - Acart_nedelec)',MINVAL(curlA-Acart),MAXVAL(curlA-Acart)
    !TEST
    !since computed at GL points, already continuous
    CALL ComputeCartCurl(Acart,curlA,vec_in_Base=0)
    curlA(1,:,:,:,:)=curlA(1,:,:,:,:)+constB(1)
    curlA(2,:,:,:,:)=curlA(2,:,:,:,:)+constB(2)
    curlA(3,:,:,:,:)=curlA(3,:,:,:,:)+constB(3)
  END IF
  WRITE(*,*)'CHECK B=curlA_init...'
  CALL CheckB(curlA)
  !replace B in U_in
  CALL swapB(U_inout,curlA)
   
END SUBROUTINE InitB_CT

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
USE MOD_CT_Vars    ,ONLY: doEvalSource_A 
USE MOD_Equation_Vars    ,ONLY: IniExactFunc
USE MOD_Mesh_Vars  ,ONLY: nElems,dXGL_N,Elem_xGP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,l,iElem
REAL    :: AtCov(3,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL    :: At(3),At_source(3)
!==================================================================================================================================
At_source=0.
DO iElem=1,nElems
  DO k=0,PP_N;DO j=0,PP_N; DO i=0,PP_N
     IF(doEvalSource_A) CALL EvalSource_A(IniExactFunc,Elem_xGP(:,i,j,k,iElem),At_source)
     At = CROSS(U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem),U(6:8,i,j,k,iElem))+At_source  !=v x B = - B x v
     AtCov(:,i,j,k,iElem)= MATMUL(dXGL_N(:,:,i,j,k,iElem),At(:)) !covariant components
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems

CALL ProjectToNedelec(AtCov)
 
CALL ComputeCartCurl(AtCov,curlAt,vec_in_Base=1)

!WRITE(*,*)'check B=curlA_t  !!!!!!!!!'
!CALL checkB(curlAt)
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
USE MOD_Mesh_Vars,ONLY: nElems,SideToElem,ElemToSide,Elem_xGP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: vec(3,0:PP_N,0:PP_N,0:PP_N,nElems) !< covariant components of a generally discontinuous vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL       :: vec_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
INTEGER    :: i,j,k,iElem,nb(-1:1,-1:1),xx,yy,zz
!==================================================================================================================================
vec_in=vec
DO iElem=1,nElems
  !A1: C0 in y,z
  nb( 0, 0) =iElem  ! 
  nb( 0,-1) =GetnbElemID(ZETA_MINUS,iElem)      !    z- 
  nb( 0, 1) =GetnbElemID(ZETA_PLUS ,iElem)      !    z+ 
  nb(-1, 0) =GetnbElemID( ETA_MINUS,iElem)      ! y-
  nb( 1, 0) =GetnbElemID( ETA_PLUS ,iElem)      ! y+
  nb(-1,-1) =GetnbElemID( ETA_MINUS,nb(0,-1))   ! y-,z-
  nb( 1,-1) =GetnbElemID( ETA_PLUS ,nb(0,-1))   ! y+,z-
  nb(-1, 1) =GetnbElemID( ETA_MINUS,nb(0, 1))   ! y-,z+
  nb( 1, 1) =GetnbElemID( ETA_PLUS ,nb(0, 1))   ! y+,z+
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
  nb( 0,-1) =GetnbElemID(ZETA_MINUS,iElem)    !    z- 
  nb( 0, 1) =GetnbElemID(ZETA_PLUS ,iElem)    !    z+ 
  nb(-1, 0) =GetnbElemID(  XI_MINUS,iElem)    ! x-
  nb( 1, 0) =GetnbElemID(  XI_PLUS ,iElem)    ! x+
  nb(-1,-1) =GetnbElemID(  XI_MINUS,nb(0,-1)) ! x-,z-
  nb( 1,-1) =GetnbElemID(  XI_PLUS ,nb(0,-1)) ! x+,z-
  nb(-1, 1) =GetnbElemID(  XI_MINUS,nb(0, 1)) ! x-,z+
  nb( 1, 1) =GetnbElemID(  XI_PLUS ,nb(0, 1)) ! x+,z+
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
  nb( 0,-1) =GetnbElemID( ETA_MINUS,iElem)       !    y- 
  nb( 0, 1) =GetnbElemID( ETA_PLUS ,iElem)       !    y+ 
  nb(-1, 0) =GetnbElemID(  XI_MINUS,iElem)       ! x-
  nb( 1, 0) =GetnbElemID(  XI_PLUS ,iElem)       ! x+
  nb(-1,-1) =GetnbElemID(  XI_MINUS,nb(0,-1))    ! x-,y-
  nb( 1,-1) =GetnbElemID(  XI_PLUS ,nb(0,-1))    ! x+,y-
  nb(-1, 1) =GetnbElemID(  XI_MINUS,nb(0, 1))    ! x-,y+
  nb( 1, 1) =GetnbElemID(  XI_PLUS ,nb(0, 1))    ! x+,y+
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


CONTAINS

  FUNCTION GETnbElemID(locSide,elemID)
    !----------------------------------
    INTEGER,INTENT(IN) :: locSide
    INTEGER,INTENT(IN) :: elemID
    INTEGER            :: GETnbElemID
    !----------------------------------
    INTEGER            :: SideID,Flip
    !----------------------------------
    SideID=ElemToSide(E2S_SIDE_ID,locSide,ElemID)
    Flip  =ElemToSide(E2S_FLIP   ,locSide,ElemID)
    IF(flip.EQ.0)THEN
      GETnbElemID=SideToElem(S2E_NB_ELEM_ID,SideID)
    ELSE
      GETnbElemID=SideToElem(S2E_ELEM_ID,SideID)
    END IF

  END FUNCTION GETnbElemID

END SUBROUTINE projectToNedelec


!==================================================================================================================================
!> evaluate the magnetic vector potential A, belonging to the magnetic field
!> of a given IniExactfunction 
!> B_1= d(A_3)/dy - d(A_2)/dz
!> B_2= d(A_1)/dz - d(A_3)/dx
!> B_3= d(A_2)/dx - d(A_1)/dy
!==================================================================================================================================
SUBROUTINE EvalExactFunc_A(ExactFunc,x,Acart,Bcart,constB,discreteDivB)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0,IniHalfwidth,IniFrequency,IniAmplitude
USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: ExactFunc
REAL   ,INTENT(IN)  :: x(3) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Acart(3)
REAL,INTENT(OUT)    :: Bcart(3)
REAL,INTENT(OUT)    :: constB(3)
LOGICAL,INTENT(OUT) :: discreteDivB   !< set to true if divergence of Bcart is already DISCRETELY zero!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: Omega,e,nx,ny,phi_alv
!==================================================================================================================================
constB=0.
SELECT CASE(ExactFunc)
CASE(3) !alfven wave domain [-1,1]^3, A is periodic
  Omega=2.*PP_Pi*IniFrequency
  !r=2,sqr=1,teval=0.
  e   = 0.2
  nx  = 1./SQRT(2**2+1.)
  ny  = 2./SQRT(2**2+1.)
  phi_alv = omega*(nx/ny*(x(1)-1.) + (x(2)-1.))


  Acart(1) = 0.
  Acart(2) = e*ny/(nx*Omega)*COS(phi_alv)
  Acart(3) = e*ny/Omega*SIN(phi_alv) 

  constB(1)=nx
  constB(2)=ny
  constB(3)=0.

  discreteDivB=.FALSE.

  Bcart(1) =nx+ e*ny*COS(phi_alv)
  Bcart(2) =ny -e*nx*COS(phi_alv)
  Bcart(3) =   -e*SIN(phi_alv) 

CASE(3001) !alfven wave domain [-1,1]^3, A is periodic, rotated
  Omega=2.*PP_Pi*IniFrequency
  !r=2,sqr=1,teval=0.
  e   = 0.2
  nx  = 1./SQRT(2**2+1.)
  ny  = 2./SQRT(2**2+1.)
  phi_alv = omega*(nx/ny*(x(1)-1.) + (x(3)-1.))


  Acart(1) = 0.
  Acart(2) =-e*ny/Omega*SIN(phi_alv) 
  Acart(3) =-e*ny/(nx*Omega)*COS(phi_alv)

  constB(1)=nx
  constB(2)=0
  constB(3)=ny

  discreteDivB=.FALSE.

  Bcart(1) =nx+ e*ny*COS(phi_alv)
  Bcart(2) =   -e*SIN(phi_alv) 
  Bcart(3) =ny -e*nx*COS(phi_alv)

CASE(7)
  Omega=PP_Pi*IniFrequency
  Acart(1)= 2.*IniAmplitude*SIN(Omega*SUM(x(:)))
  Acart(2)=    IniAmplitude*SIN(Omega*SUM(x(:)))
  Acart(3)=  - IniAmplitude*SIN(Omega*SUM(x(:)))
  discreteDivB=.FALSE.

  !B1=dA3/dy-dA2/dz = -2*IniAmplitude*Omega*COS(Omega*SUM(x))
  !B2=dA1/dz-dA3/dx =  3*IniAmplitude*Omega*COS(Omega*SUM(x))
  !B3=dA2/dx-dA1/dy = -1*IniAmplitude*Omega*COS(Omega*SUM(x))

  Bcart(1) = -2*IniAmplitude*Omega*COS(Omega*SUM(x))
  Bcart(2) =  3*IniAmplitude*Omega*COS(Omega*SUM(x))
  Bcart(3) = -  IniAmplitude*Omega*COS(Omega*SUM(x))
CASE(73) !IS NOT PERIODIC!!!
  Bcart(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
  Bcart(2)=0.
  Bcart(3)=1.
  discreteDivB=.TRUE. !A,constB are not used...
  Acart(1) = 0.
  Acart(2) = - x(3)*Bcart(1) !NOT PERIODIC!!!!
  Acart(3) = 0. 
  !constant Bfield,to be added
  constB(3)=1.0
CASE(74) !IS NOT PERIODIC !!!
  Bcart(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
  Bcart(2)=0.
  Bcart(3)=SQRT(1.0-Bcart(1)**2)
  discreteDivB=.TRUE.  ! A,constB are not used...
  Acart(1) = 0.
  Acart(2) = x(1)*Bcart(3)-x(3)*Bcart(1) !NOT PERIODIC
  Acart(3) = 0. 

CASE DEFAULT
  STOP 'this exactfunc is not implemented in EvalExactFunc_A'
END SELECT !ExactFunc

END SUBROUTINE EvalExactFunc_A

!==================================================================================================================================
!> evaluate the magnetic vector potential A, belonging to the magnetic field
!> of a given IniExactfunction 
!> B_1= d(A_3)/dy - d(A_2)/dz
!> B_2= d(A_1)/dz - d(A_3)/dx
!> B_3= d(A_2)/dx - d(A_1)/dy
!==================================================================================================================================
SUBROUTINE EvalSource_A(ExactFunc,x,Atcart)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_CT_Vars,        ONLY: doEvalSource_A
USE MOD_Equation_Vars  ,ONLY: s2mu_0,IniHalfwidth,IniFrequency,IniAmplitude
USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: ExactFunc
REAL   ,INTENT(IN)  :: x(3) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Atcart(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: Omega
!==================================================================================================================================
SELECT CASE (ExactFunc)
CASE(4,5,6) ! navierstokes exact function
  STOP 'you need to provide a source either for the magnetic vector potential At'
CASE DEFAULT
  ! No source -> do nothing
  doEvalSource_A=.FALSE. 
END SELECT ! ExactFunction

END SUBROUTINE EvalSource_A

!==================================================================================================================================
!> Replaces B of current U by curlA, and changing total energy and keeping the same pressure
!==================================================================================================================================
SUBROUTINE swapB(U_inout,curlA_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0
USE MOD_Mesh_Vars      ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: curlA_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: EmagU,EmagCurlA
REAL    :: diffB(3),mindiffB(3),maxdiffB(3) 
!==================================================================================================================================
maxdiffB=-1.0e20
mindiffB=1.0e20
DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
  diffB=U_inout(6:8,i,j,k,iElem)-curlA_in(:,i,j,k,iElem)
  mindiffB=MIN(mindiffB,diffB)
  maxdiffB=MAX(maxdiffB,diffB)
END DO; END DO; END DO; END DO ! i,j,k,iElem

!WRITE(*,*)'SWAPB: min(B-curlA)',mindiffB
!WRITE(*,*)'       max(B-curlA)',maxdiffB


!WRITE(*,*)'CHECK B(U)...'
!CALL CheckB(U_inout(6:8,:,:,:,:))
!WRITE(*,*)'CHECK B=curlA...'
!CALL CheckB(curlA_in)

DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
  EmagU=SUM(U_inout(6:8,i,j,k,iElem)**2)
  EmagcurlA=SUM(curlA_in(:,i,j,k,iElem)**2)
  U_inout(5,i,j,k,iElem)=U_inout(5,i,j,k,iElem)+s2mu_0*(EmagcurlA-EmagU) !add new magn. energy, substract old energy
  U_inout(6:8,i,j,k,iElem)=curlA_in(:,i,j,k,iElem)
END DO; END DO; END DO; END DO ! i,j,k,iElem

END SUBROUTINE swapB


!==================================================================================================================================
!> Replaces Bt of current Ut by curlAt, and change total energy Et to remain consistent with the entropy
!> inequality:
!> S_t=w^T Ut = ((gamma-s)/(gamma-1)-(rho/p)1/2|v|^2,(rho/p)vvec,-(rho/p),(rho/p)Bvec) ^T 
!>               (rho_t, (rho*vvec)_t,E_t,Bvec_t)
!>            =  ... (rho/p) (-E_t + Bvec*Bvec_t  + [ Bvec*(curlA_t-curlA_t) ] ) 
!>
!> in new Ut, the B_t is replaced with curlA_t, yielding
!>        =>  ... (rho/p) (-E*_t + Bvec*curlA_t)
!> if we choose:
!> E*_t =E_t + Bvec*(curlA_t-B_t)
!> we recover the original form, such that (w^T *U_t) remains the same!

!==================================================================================================================================
SUBROUTINE swapBt(U_in,Ut_inout)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: smu_0
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_CT_Vars      ,ONLY: curlAt
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Ut_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: EmagUt,EmagCurlAt
!==================================================================================================================================

DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
  EmagUt=SUM(U_in(6:8,i,j,k,iElem)*Ut_inout(6:8,i,j,k,iElem))
  EmagcurlAt=SUM(U_in(6:8,i,j,k,iElem)*curlAt(:,i,j,k,iElem))
  Ut_inout(5,i,j,k,iElem)=Ut_inout(5,i,j,k,iElem)+smu_0*(EmagcurlAt-EmagUt) !to keep entropy inequality
  Ut_inout(6:8,i,j,k,iElem)=curlAt(:,i,j,k,iElem)
END DO; END DO; END DO; END DO ! i,j,k,iElem

END SUBROUTINE swapBt

!==================================================================================================================================
!> checks divergence and surface jumps of cartesian vector field 
!> ONLY 1 PROC & STRUCTURED MESH 
!==================================================================================================================================
SUBROUTINE checkB(B_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_Mesh_Vars      ,ONLY: sJ,metrics_ftilde,metrics_gtilde,metrics_htilde
USE MOD_Mesh_Vars      ,ONLY: ElemToSide,nSides,SideToElem
USE MOD_DG_Vars        ,ONLY: D 
!USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: B_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: i
REAL    :: Bn_plus(0:PP_N,0:PP_N,nSides),Bn_minus(0:PP_N,0:PP_N,nSides)
INTEGER :: i,j,k,l,iElem
REAL    :: mindivB,maxDivB,divB_loc,Btilde(3,0:PP_N,0:PP_N,0:PP_N)
!==================================================================================================================================

mindivB=1.0e20
maxdivB=-1.0e20
Bn_plus=9999.
Bn_minus=9999.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Btilde(1,i,j,k)=SUM(Metrics_ftilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
    Btilde(2,i,j,k)=SUM(Metrics_gtilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
    Btilde(3,i,j,k)=SUM(Metrics_htilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
  Bn_plus( :,:,ElemToSide(E2S_SIDE_ID,   XI_PLUS,iElem))=Btilde(1,PP_N,:,:)
  Bn_plus( :,:,ElemToSide(E2S_SIDE_ID,  ETA_PLUS,iElem))=Btilde(2,:,PP_N,:)
  Bn_plus( :,:,ElemToSide(E2S_SIDE_ID, ZETA_PLUS,iElem))=Btilde(3,:,:,PP_N)
  Bn_minus(:,:,ElemToSide(E2S_SIDE_ID,  XI_MINUS,iElem))=Btilde(1,0,:,:)
  Bn_minus(:,:,ElemToSide(E2S_SIDE_ID, ETA_MINUS,iElem))=Btilde(2,:,0,:)
  Bn_minus(:,:,ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))=Btilde(3,:,:,0)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    divB_loc=0.
    DO l=0,PP_N
      divB_loc=divB_loc + D(i,l)*Btilde(1,l,j,k) + &
                          D(j,l)*Btilde(2,i,l,k) + &
                          D(k,l)*Btilde(3,i,j,l)
    END DO ! l
    divB_loc = sJ(i,j,k,iElem) * divB_loc 
    mindivB=MIN(mindivB,divB_loc)
    maxdivB=MAX(maxdivB,divB_loc)

  END DO; END DO; END DO ! i,j,k
END DO ! iElem

WRITE(*,*)'..check ..min/max divB    : ',mindivB,maxdivB
WRITE(*,*)'        ..min/max [[B*n]] : ',MINVAL((Bn_plus-Bn_minus)),MAXVAL((Bn_plus-Bn_minus))
!IF(MAXVAL(ABS(Bn_plus-Bn_minus)).GT.1.0e-10)THEN
!  DO i=1,nSides
!    WRITE(*,*)'SideID',i,'locSide master',SideToElem(S2E_LOC_SIDE_ID,i) &
!                        ,'locSide slave ',SideToElem(S2E_NB_LOC_SIDE_ID,i)               
!                                                 
!    WRITE(*,*)i,'     ..min/max [[B*n]] : ',MINVAL((Bn_plus(:,:,i)-Bn_minus(:,:,i))),MAXVAL((Bn_plus(:,:,i)-Bn_minus(:,:,i))),MAXVAL(Bn_plus(:,:,i))
!  END DO
!  READ(*,*)
!END IF

END SUBROUTINE checkB


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
