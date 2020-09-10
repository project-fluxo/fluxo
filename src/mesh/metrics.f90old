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
!> \brief This module contains routines for computing the geometries volume and surface metric terms.
!>
!> Compute the volume and surface metric terms:
!>     Metrics_fTilde(n=1:3,i,j,k,iElem)=Ja_n^1
!>     Metrics_gTilde(n=1:3,i,j,k,iElem)=Ja_n^2
!>     Metrics_hTilde(n=1:3,i,j,k,iElem)=Ja_n^3
!> 
!>   Per Element we do:
!>   1.) a.) Preparation: the geometry (equidistant nodal basis, NGeo+1 points/dir) is interpolated to a high precision
!>           mapping X_n(xi_i) using a Gauss-Lobatto basis and stored in XGL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!>       b.) Computing the gradients: compute the derivative of the mapping XGL_NGeo in \f$ (xi_1,xi_2,xi_3) \f$ direction,
!>           using a polynomial derivative Matrix at degree NGeo.
!>       c.) Computing the Jacobian: compute Jacobian JRef at a degree of NGeoRef=3*NGeo (exact). 
!>                                   For this gradients have to be interpolated to NGeoRef first.
!>                                   Then project JRef down to degree N. Finally check for negative Jacobians.
!>       d.) For computing Ja the gradients at degree N are required: if N>=NGeo directly interpolate dXGL_NGeo to dXGL_N,
!>                                                                    else compute dXGL_N from XGL_N directly.
!>
!>   2.) for each direction n
!>       a.) compute the nth vector and for each Gauss-Lobatto point (:,i,j,k)
!>          \f$(dXGL_n^1,dXGL_n^2,dXGL_n^3)^T=(X_l grad_xi (X_m) )\f$ for n=1,2,3 and (n,m,l) cyclic
!>       b.) interpolate the dXGL_n vector defined primarily on (NGeo+1)x(NGeo+1)x(Ngeo+1) Gauss-Lobatto points to
!>             (N+1)x(N+1)x(N+1) Gauss-Lobatto points and write to Ja_n(1:3,i,j,k) i,j,k=[0:N]
!>       c.) compute the curl of vector Ja_n(1:3,i,j,k) using the derivative Matrix DGL_N [NxN]
!>       d.) interpolate from (N+1)x(N+1)x(N+1) Gauss-Lobatto points to  Gauss-Points (N+1)x(N+1)x(N+1) (exact!)
!>       e.) store Ja_n in the Metrics arrays
!>
!>   3.) Compute the surface metrics (normal/tangential vectors, surface area) from volume metrics for each side.
!> 
!==================================================================================================================================
MODULE MOD_Metrics
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CalcMetrics
  MODULE PROCEDURE CalcMetrics
END INTERFACE

INTERFACE CalcSurfMetrics
  MODULE PROCEDURE CalcSurfMetrics
END INTERFACE

INTERFACE SurfMetricsFromJa
  MODULE PROCEDURE SurfMetricsFromJa
END INTERFACE

PUBLIC::CalcMetrics
PUBLIC::CalcSurfMetrics
PUBLIC::SurfMetricsFromJa
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine computes the geometries volume metric terms.
!==================================================================================================================================
SUBROUTINE CalcMetrics()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,     ONLY:NGeo,NgeoRef,nElems,offsetElem,crossProductMetrics,NodeTypeMesh
USE MOD_Mesh_Vars,     ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars,     ONLY:sJ,detJac_Ref
USE MOD_Mesh_Vars,     ONLY:Vdm_GLN_N,dXGL_N
USE MOD_Mesh_Vars,     ONLY:NodeCoords,Elem_xGP
USE MOD_Interpolation_Vars
USE MOD_Interpolation, ONLY:GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_ChangeBasis,   ONLY:changeBasis3D
USE MOD_Basis,         ONLY:LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,q,iElem
INTEGER :: ll
! Jacobian on GL N and NGeoRef
REAL    :: DetJac_N( 1,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on GL N
REAL    :: XGL_N(      3,  0:PP_N,0:PP_N,0:PP_N)          ! mapping X(xi) P\in N
REAL    :: XGL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: dXGL_Ngeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on GL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

REAL    :: R_GL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! buffer for metric terms, uses XGL_N,dXGL_N
REAL    :: JaGL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! metric terms P\in N
REAL    :: scaledJac(2)

! Polynomial derivativion matrices
REAL    :: DGL_NGeo(0:Ngeo,0:Ngeo)
REAL    :: DGL_N(   0:PP_N,0:PP_N)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNgeo_GLNgeo( 0:Ngeo   ,0:Ngeo)
REAL    :: Vdm_GLNGeo_NgeoRef(0:NgeoRef,0:Ngeo)
REAL    :: Vdm_NgeoRef_N(     0:PP_N   ,0:NgeoRef)
REAL    :: Vdm_GLNGeo_GLN(    0:PP_N   ,0:Ngeo)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)
REAL    :: xiGL_N(0:PP_N)   ,wBaryGL_N(0:PP_N)
!==================================================================================================================================
! Prerequisites
Metrics_fTilde=0.
Metrics_gTilde=0.
Metrics_hTilde=0.

! Initialize Vandermonde and D matrices
! Only use modal Vandermonde for terms that need to be conserved as Jacobian if N_out>N_in
! Always use interpolation for the rest!

! 1.a) NodeCoords: EQUI Ngeo to GLNgeo and GLN
CALL GetVandermonde(    Ngeo   , NodeTypeMesh, Ngeo    , NodeTypeGL, Vdm_EQNgeo_GLNgeo , modal=.FALSE.)

! 1.b) dXGL_Ngeo:
CALL GetDerivativeMatrix(Ngeo  , NodeTypeGL  , DGL_Ngeo)

! 1.c) Jacobian: GLNgeo to NgeoRef, GLNgeoRef to N
CALL GetVandermonde(    Ngeo   , NodeTypeGL  , NgeoRef , NodeType  , Vdm_GLNgeo_NgeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NgeoRef_N     , modal=.TRUE.)
CALL GetNodesAndWeights(NgeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)

! 1.d) derivatives (dXGL) by projection or by direct derivation (D_GL):
CALL GetVandermonde(    Ngeo   , NodeTypeGL  , PP_N    , NodeTypeGL, Vdm_GLNgeo_GLN    , modal=.FALSE.)
CALL GetDerivativeMatrix(PP_N  , NodeTypeGL  , DGL_N)

! 2.d) derivatives (dXGL) by projection or by direct derivation (D_GL):
CALL GetNodesAndWeights(PP_N   , NodeTypeGL  , xiGL_N  , wIPBary=wBaryGL_N)

! Outer loop over all elements
DO iElem=1,nElems
  !1.a) Transform from EQUI_Ngeo to GL points on Ngeo and N
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_GLNGeo,NodeCoords(:,:,:,:,iElem) ,XGL_Ngeo               )
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_GLNGeo_GLN   ,XGL_Ngeo                  ,XGL_N                  )
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GLN_N        ,XGL_N                     ,Elem_xGP(:,:,:,:,iElem))

  !1.b) Jacobi Matrix of d/dxi_dd(X_nn): dXGL_NGeo(dd,nn,i,j,k))
  dXGL_NGeo=0.
  DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
    ! Matrix-vector multiplication
    DO ll=0,Ngeo
      dXGL_Ngeo(1,:,i,j,k)=dXGL_Ngeo(1,:,i,j,k) + DGL_Ngeo(i,ll)*XGL_Ngeo(:,ll,j,k)
      dXGL_Ngeo(2,:,i,j,k)=dXGL_Ngeo(2,:,i,j,k) + DGL_Ngeo(j,ll)*XGL_Ngeo(:,i,ll,k)
      dXGL_Ngeo(3,:,i,j,k)=dXGL_Ngeo(3,:,i,j,k) + DGL_Ngeo(k,ll)*XGL_Ngeo(:,i,j,ll)
    END DO !l=0,N
  END DO; END DO; END DO !i,j,k=0,Ngeo

  ! 1.c)Jacobians! grad(X_1) (grad(X_2) x grad(X_3))
  ! Compute Jacobian on NGeo and then interpolate:
  ! required to guarantee conservativity when restarting with N<NGeo
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_GLNGeo_NgeoRef,dXGL_NGeo(:,1,:,:,:),dX_NgeoRef(:,1,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_GLNGeo_NgeoRef,dXGL_NGeo(:,2,:,:,:),dX_NgeoRef(:,2,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_GLNGeo_NgeoRef,dXGL_NGeo(:,3,:,:,:),dX_NgeoRef(:,3,:,:,:))
  !detJac_Ref(:,:,:,:,iElem)=0. !set above
  detJac_Ref(1,:,:,:,iElem)=0.
  DO k=0,NgeoRef; DO j=0,NgeoRef; DO i=0,NgeoRef
    ASSOCIATE(dX_Ngeo_R => dX_NgeoRef(:,:,i,j,k))
    detJac_Ref(1,i,j,k,iElem)=detJac_Ref(1,i,j,k,iElem) & 
                              + dX_Ngeo_R(1,1)*(dX_Ngeo_R(2,2)*dX_Ngeo_R(3,3) - dX_Ngeo_R(3,2)*dX_Ngeo_R(2,3))  &
                              + dX_Ngeo_R(2,1)*(dX_Ngeo_R(3,2)*dX_Ngeo_R(1,3) - dX_Ngeo_R(1,2)*dX_Ngeo_R(3,3))  &
                              + dX_Ngeo_R(3,1)*(dX_Ngeo_R(1,2)*dX_Ngeo_R(2,3) - dX_Ngeo_R(2,2)*dX_Ngeo_R(1,3))  
    END ASSOCIATE !dX_Ngeo_R => dX_NgeoRef
  END DO; END DO; END DO !i,j,k=0,NgeoRef
!
  ! projection detJacref from NgeoRef to N (mean value=Volume independant of the polynomial degree)
  CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N)

  ! assign to global Variable sJ
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    sJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
  END DO; END DO; END DO !i,j,k=0,PP_N

  ! check for negative Jacobians
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IF(detJac_N(1,i,j,k).LE.0.)&
      WRITE(Unit_StdOut,*) 'Negative Jacobian found on Gauss point. Coords:', Elem_xGP(:,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,N
  ! check scaled Jacobians
  scaledJac(2)=MINVAL(detJac_N(1,:,:,:))/MAXVAL(detJac_N(1,:,:,:))
  IF(scaledJac(2).LT.0.01) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (GL/Gauss):', scaledJac
    CALL abort(__STAMP__,&
      'Scaled Jacobian lower then tolerance in global element:',iElem+offsetElem)
  END IF

  !2.a) Jacobi Matrix of d/dxi_dd(X_nn): dXGL_N(dd,nn,i,j,k))
  ! N>=Ngeo: interpolate from dXGL_Ngeo (default)
  ! N< Ngeo: directly derive XGL_N
  IF(PP_N.GE.NGeo)THEN !compute first derivative on Ngeo and then interpolate
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_GLNGeo_GLN,dXGL_NGeo(:,1,:,:,:),dXGL_N(:,1,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_GLNGeo_GLN,dXGL_NGeo(:,2,:,:,:),dXGL_N(:,2,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_GLNGeo_GLN,dXGL_NGeo(:,3,:,:,:),dXGL_N(:,3,:,:,:,iElem))
  ELSE  !N<Ngeo: first interpolate and then compute derivative (important if curved&periodic)
    dXGL_N(:,:,:,:,:,iElem)=0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! Matrix-vector multiplication
      DO ll=0,PP_N
        dXGL_N(1,:,i,j,k,iElem)=dXGL_N(1,:,i,j,k,iElem) + DGL_N(i,ll)*XGL_N(:,ll,j,k)
        dXGL_N(2,:,i,j,k,iElem)=dXGL_N(2,:,i,j,k,iElem) + DGL_N(j,ll)*XGL_N(:,i,ll,k)
        dXGL_N(3,:,i,j,k,iElem)=dXGL_N(3,:,i,j,k,iElem) + DGL_N(k,ll)*XGL_N(:,i,j,ll)
      END DO !l=0,N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=Ngeo
  !save to covar array

  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXGL => dXGL_N(:,:,i,j,k,iElem))
      ! exact (cross-product) form
      ! Ja(:)^nn = ( d/dxi_(nn+1) XGL_N(:) ) x (d/xi_(nn+2) XGL_N(:))
      !
      ! JaGL_N(dd,nn) = dXGL_N(dd+1,nn+1)*dXGL_N(dd+2,nn+2) -dXGL_N(dd+1,nn+2)*dXGL_N(dd+2,nn+1)
      JaGL_N(1,1,i,j,k)=dXGL(2,2)*dXGL(3,3) - dXGL(2,3)*dXGL(3,2)  
      JaGL_N(2,1,i,j,k)=dXGL(3,2)*dXGL(1,3) - dXGL(3,3)*dXGL(1,2)  
      JaGL_N(3,1,i,j,k)=dXGL(1,2)*dXGL(2,3) - dXGL(1,3)*dXGL(2,2)  
      JaGL_N(1,2,i,j,k)=dXGL(2,3)*dXGL(3,1) - dXGL(2,1)*dXGL(3,3)  
      JaGL_N(2,2,i,j,k)=dXGL(3,3)*dXGL(1,1) - dXGL(3,1)*dXGL(1,3)  
      JaGL_N(3,2,i,j,k)=dXGL(1,3)*dXGL(2,1) - dXGL(1,1)*dXGL(2,3)  
      JaGL_N(1,3,i,j,k)=dXGL(2,1)*dXGL(3,2) - dXGL(2,2)*dXGL(3,1)  
      JaGL_N(2,3,i,j,k)=dXGL(3,1)*dXGL(1,2) - dXGL(3,2)*dXGL(1,1)  
      JaGL_N(3,3,i,j,k)=dXGL(1,1)*dXGL(2,2) - dXGL(1,2)*dXGL(2,1)  
      END ASSOCIATE !dXGL => dXGL_N
    END DO; END DO; END DO !i,j,k=0,N
  ELSE ! curl metrics
    ! invariant curl form, as cross product: R^dd = 1/2( XGL_N(:) x (d/dxi_dd XGL_N(:)))
    !
    !R_GL_N(dd,nn)=1/2*( XGL_N(nn+2)* d/dxi_dd XGL_N(nn+1) - XGL_N(nn+1)* d/dxi_dd XGL_N(nn+2))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXGL => dXGL_N(:,:,i,j,k,iElem))
      R_GL_N(:,1,i,j,k)=0.5*(XGL_N(3,i,j,k)*dXGL(:,2) - XGL_N(2,i,j,k)*dXGL(:,3) )
      R_GL_N(:,2,i,j,k)=0.5*(XGL_N(1,i,j,k)*dXGL(:,3) - XGL_N(3,i,j,k)*dXGL(:,1) )
      R_GL_N(:,3,i,j,k)=0.5*(XGL_N(2,i,j,k)*dXGL(:,1) - XGL_N(1,i,j,k)*dXGL(:,2) ) 
      END ASSOCIATE !dXGL => dXGL_N
    END DO; END DO; END DO !i,j,k=0,N
    ! Metrics are the curl of R:  Ja(:)^nn = -(curl R_GL(:,nn))
    ! JaGL_N(dd,nn)= -[d/dxi_(dd+1) RGL(dd+2,nn) - d/dxi_(dd+2) RGL(dd+1,nn) ]
    !              =   d/dxi_(dd+2) RGL(dd+1,nn) - d/dxi_(dd+1) RGL(dd+2,nn) 
    JaGL_N=0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      DO q=0,PP_N
        JaGL_N(1,:,i,j,k)=JaGL_N(1,:,i,j,k) - DGL_N(j,q)*R_GL_N(3,:,i,q,k)
        JaGL_N(2,:,i,j,k)=JaGL_N(2,:,i,j,k) - DGL_N(k,q)*R_GL_N(1,:,i,j,q)
        JaGL_N(3,:,i,j,k)=JaGL_N(3,:,i,j,k) - DGL_N(i,q)*R_GL_N(2,:,q,j,k)
      END DO!q=0,PP_N
      DO q=0,PP_N
        JaGL_N(1,:,i,j,k)=JaGL_N(1,:,i,j,k) + DGL_N(k,q)*R_GL_N(2,:,i,j,q) 
        JaGL_N(2,:,i,j,k)=JaGL_N(2,:,i,j,k) + DGL_N(i,q)*R_GL_N(3,:,q,j,k) 
        JaGL_N(3,:,i,j,k)=JaGL_N(3,:,i,j,k) + DGL_N(j,q)*R_GL_N(1,:,i,q,k) 
      END DO!q=0,PP_N
! same with only one loop, gives different roundoff ...
!      DO q=0,PP_N
!        JaGL_N(1,:,i,j,k)=JaGL_N(1,:,i,j,k) - DGL_N(j,q)*R_GL_N(3,:,i,q,k) + DGL_N(k,q)*R_GL_N(2,:,i,j,q)
!        JaGL_N(2,:,i,j,k)=JaGL_N(2,:,i,j,k) - DGL_N(k,q)*R_GL_N(1,:,i,j,q) + DGL_N(i,q)*R_GL_N(3,:,q,j,k)
!        JaGL_N(3,:,i,j,k)=JaGL_N(3,:,i,j,k) - DGL_N(i,q)*R_GL_N(2,:,q,j,k) + DGL_N(j,q)*R_GL_N(1,:,i,q,k)
!      END DO!q=0,PP_N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !crossProductMetrics


  ! interpolate Metrics from Gauss-Lobatto N onto GaussPoints N
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GLN_N,JaGL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GLN_N,JaGL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GLN_N,JaGL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem))
  CALL CalcSurfMetrics(PP_N,JaGL_N,XGL_N,Vdm_GLN_N,iElem)
END DO !iElem=1,nElems

END SUBROUTINE CalcMetrics



!==================================================================================================================================
!> Prepares computation of the faces' normal, tangential vectors, surface area and Gauss points from volume metrics.
!> Input is JaGL_N, the 3D element metrics on Cebychev-Lobatto points.
!> For each side the volume metrics are interpolated to the surface and rotated into the side reference frame. 
!==================================================================================================================================
SUBROUTINE CalcSurfMetrics(Nloc,JaGL_N,XGL_N,Vdm_GLN_N,iElem)
! MODULES
USE MOD_Globals,        ONLY:CROSS
USE MOD_Mesh_Vars,      ONLY:ElemToSide,MortarType
USE MOD_Mesh_Vars,      ONLY:NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Mesh_Vars,      ONLY:NormalDirs,TangDirs,NormalSigns
USE MOD_Mappings,       ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis,    ONLY:ChangeBasis2D
USE MOD_Mortar_Metrics, ONLY:Mortar_CalcSurfMetrics
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                !< (IN) polynomial degree
INTEGER,INTENT(IN) :: iElem                               !< (IN) element index
REAL,INTENT(IN)    :: JaGL_N(  3,3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) volume metrics of element
REAL,INTENT(IN)    :: XGL_N(     3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) element geo. interpolation points (GL)
REAL,INTENT(IN)    :: Vdm_GLN_N(   0:Nloc,0:Nloc)         !< (IN) Vandermonde matrix from Gauss-Lob on N to final nodeset on N
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID,SideID2,iMortar,nbSideIDs(4)
INTEGER            :: NormalDir,TangDir
REAL               :: NormalSign
REAL               :: Ja_Face(3,3,0:Nloc,0:Nloc)
REAL               :: Mortar_Ja(3,3,0:Nloc,0:Nloc,4)
REAL               :: Mortar_xGP( 3,0:Nloc,0:Nloc,4)
REAL               :: tmp(        3,0:Nloc,0:Nloc)
REAL               :: tmp2(       3,0:Nloc,0:Nloc)
!==================================================================================================================================

DO iLocSide=1,6
  IF(ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) CYCLE ! only master sides with flip=0
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)

  SELECT CASE(iLocSide)
  CASE(XI_MINUS)
    tmp=XGL_N(1:3,0   ,:   ,:   )
  CASE(XI_PLUS)
    tmp=XGL_N(1:3,Nloc,:   ,:   )
  CASE(ETA_MINUS)
    tmp=XGL_N(1:3,:   ,0   ,:   )
  CASE(ETA_PLUS)
    tmp=XGL_N(1:3,:   ,Nloc,:   )
  CASE(ZETA_MINUS)
    tmp=XGL_N(1:3,:   ,:   ,0   )
  CASE(ZETA_PLUS)
    tmp=XGL_N(1:3,:   ,:   ,Nloc)
  END SELECT
  CALL ChangeBasis2D(3,Nloc,Nloc,Vdm_GLN_N,tmp,tmp2)
  ! turn into right hand system of side
  DO q=0,Nloc; DO p=0,Nloc
    pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
    ! Compute Face_xGP for sides
    Face_xGP(1:3,p,q,sideID)=tmp2(:,pq(1),pq(2))
  END DO; END DO ! p,q

  DO dd=1,3
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=JaGL_N(dd,1:3,0   ,:   ,:   )
    CASE(XI_PLUS)
      tmp=JaGL_N(dd,1:3,Nloc,:   ,:   )
    CASE(ETA_MINUS)
      tmp=JaGL_N(dd,1:3,:   ,0   ,:   )
    CASE(ETA_PLUS)
      tmp=JaGL_N(dd,1:3,:   ,Nloc,:   )
    CASE(ZETA_MINUS)
      tmp=JaGL_N(dd,1:3,:   ,:   ,0   )
    CASE(ZETA_PLUS)
      tmp=JaGL_N(dd,1:3,:   ,:   ,Nloc)
    END SELECT
    CALL ChangeBasis2D(3,Nloc,Nloc,Vdm_GLN_N,tmp,tmp2)
    ! turn into right hand system of side
    DO q=0,Nloc; DO p=0,Nloc
      pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
      Ja_Face(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd


  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)
  CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face,&
                         NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),&
                         TangVec2(:,:,:,SideID),SurfElem(:,:,SideID))

  !compute metrics for mortar faces, interpolate Ja_Face to small sides
  IF(MortarType(1,SideID).GT.0)THEN
    CALL Mortar_CalcSurfMetrics(SideID,Nloc,Ja_Face,Face_xGP(:,:,:,SideID),&
                                            Mortar_Ja,Mortar_xGP,nbSideIDs)
    DO iMortar=1,4
      SideID2=nbSideIDs(iMortar)
      IF(SideID2.LT.1) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
      Face_xGP(:,:,:,SideID2) = Mortar_xGP(:,:,:,iMortar)
      CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Mortar_Ja(:,:,:,:,iMortar),&
                             NormVec(:,:,:,SideID2),TangVec1(:,:,:,SideID2),&
                             TangVec2(:,:,:,SideID2),SurfElem(:,:,SideID2))
    END DO

  END IF
END DO

END SUBROUTINE CalcSurfMetrics

!==================================================================================================================================
!> Computes surface normal and tangential vectors and surface area from surface metrics Ja_Face.
!==================================================================================================================================
SUBROUTINE SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face,NormVec,TangVec1,TangVec2,SurfElem)
! MODULES
USE MOD_Globals,     ONLY: CROSS
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                       !< polynomial degree
INTEGER,INTENT(IN) :: NormalDir                  !< direction of normal vector
INTEGER,INTENT(IN) :: TangDir                    !< direction of 1. tangential vector
REAL,INTENT(IN)    :: NormalSign                 !< sign of normal vector
REAL,INTENT(IN)    :: Ja_Face(3,3,0:Nloc,0:Nloc) !< face metrics
REAL,INTENT(OUT)   ::   NormVec(3,0:Nloc,0:Nloc) !< element face normal vectors
REAL,INTENT(OUT)   ::  TangVec1(3,0:Nloc,0:Nloc) !< element face tangential vectors
REAL,INTENT(OUT)   ::  TangVec2(3,0:Nloc,0:Nloc) !< element face tangential vectors
REAL,INTENT(OUT)   ::  SurfElem(  0:Nloc,0:Nloc) !< element face surface area
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,Nloc; DO p=0,Nloc
  SurfElem(  p,q) = SQRT(SUM(Ja_Face(NormalDir,:,p,q)**2))
  NormVec( :,p,q) = NormalSign*Ja_Face(NormalDir,:,p,q)/SurfElem(p,q)
  TangVec1(:,p,q) = Ja_Face(TangDir,:,p,q) - SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q)) &
                    *NormVec(:,p,q)
  TangVec1(:,p,q) = TangVec1(:,p,q)/SQRT(SUM(TangVec1(:,p,q)**2))
  TangVec2(:,p,q) = CROSS(NormVec(:,p,q),TangVec1(:,p,q))
END DO; END DO ! p,q
END SUBROUTINE SurfMetricsFromJa

END MODULE MOD_Metrics
